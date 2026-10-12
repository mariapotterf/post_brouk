gc()

library(terra)
library(data.table)
library(dplyr)
library(readr)
library(writexl)
library(stringi)

# =============================================================================
# INPUTS
# =============================================================================
spei_file            <- "raw/SPEI_2025/SPEI12_harg1_2025/SPEI12_harg1_2025.tif"
okres_shp            <- terra::vect("raw/CR_administrativa/OKRESY_P.shp")
okres_csv_file       <- "raw/SPEI_2025/R_OKRES.csv"     # File from Roman to link 

out_dir              <- "outSPEIRoman"
out_month            <- file.path(out_dir, "SPEI12_zonal")
out_year             <- file.path(out_dir, "SPEI12_zonal_year")
dir.create(out_month, showWarnings = FALSE, recursive = TRUE)
dir.create(out_year,  showWarnings = FALSE, recursive = TRUE)

# optional manual join key if the automatic detection (section 9) picks wrongly,
# e.g. key_override <- c(shp = "KOD_OKRES", csv = "KOD")
key_override <- NULL

# =============================================================================
# 1. READ RASTER
# =============================================================================
r_stack <- rast(spei_file)
crs(r_stack) <- "EPSG:3035"     # file is tagged "unknown" LAEA; parameters = ETRS89-LAEA

# overlay check (tilted footprint is expected: original grid was lon/lat)
if (crs(okres_shp) == "") crs(okres_shp) <- "EPSG:5514"   # S-JTSK if .prj is missing
okres_3035 <- project(okres_shp, "EPSG:3035")
plot(r_stack[[1]]); lines(okres_3035)

r_stack[is.infinite(r_stack)] <- NA   # ±Inf would break mean/min/max
# optional, your methodological call:
# r_stack[r_stack < -3 | r_stack > 3] <- NA

# band names are ym_YYYYMM -> year / month
meta <- tibble(name = names(r_stack)) |>
  mutate(year  = as.integer(substr(name, 4, 7)),
         month = as.integer(substr(name, 8, 9)),
         date  = sprintf("%04d_%02d", year, month))
names(r_stack) <- meta$date

# =============================================================================
# 2. PREPARE VECTORS
# =============================================================================
okres_shp$okres_id <- seq_len(nrow(okres_shp))
v_ok <- okres_shp
if (!same.crs(r_stack, v_ok)) v_ok <- project(v_ok, crs(r_stack))

# =============================================================================
# 3. RASTERIZE OKRES
# =============================================================================
r_template <- r_stack[[1]]
r_okres <- rasterize(v_ok, r_template, field = "okres_id")
names(r_okres) <- "okres_id"

# =============================================================================
# 4. EXTRACT VALUES
# =============================================================================
dt <- as.data.table(values(c(r_okres, r_stack)))
dt_filt <- dt[!is.na(okres_id)]
dt_long <- melt(dt_filt, id.vars = "okres_id",
                variable.name = "date", value.name = "value")

zonal_stats <- dt_long[!is.na(value),
                       .(mean = mean(value), sd = sd(value), median = median(value),
                         min = min(value), max = max(value), n_cells = .N),
                       by = .(okres_id, date)]

zonal_stats[, c("year", "month") := tstrsplit(date, "_", fixed = TRUE)]
zonal_stats[, `:=`(year = as.integer(year), month = as.integer(month))]

# =============================================================================
# 5. JOIN ATTRIBUTES (shapefile)
# =============================================================================
okres_attrib <- as_tibble(v_ok) |> select(okres_id, NAZEV)
df_final <- left_join(zonal_stats, okres_attrib, by = "okres_id") |>
  relocate(NAZEV, .after = okres_id) |>
  arrange(year, month, okres_id)

# =============================================================================
# 6. COUNTRY-LEVEL ANNUAL SUMMARY (of okres monthly means)
# =============================================================================
df_annual <- as.data.table(df_final)[
  , .(country_mean = mean(mean, na.rm = TRUE),
      country_sd = sd(mean, na.rm = TRUE),
      country_median = median(mean, na.rm = TRUE),
      country_min = min(mean, na.rm = TRUE),
      country_max = max(mean, na.rm = TRUE),
      n_okres = uniqueN(okres_id)),
  by = year][order(year)]

write_csv(df_annual, file.path(out_dir, "country_SPEI12_annual_stats.csv"))
write_xlsx(df_annual, file.path(out_dir, "country_SPEI12_annual_stats.xlsx"))
write_csv(df_final,  file.path(out_dir, "okres_SPEI12_monthly_stats.csv"))
write_xlsx(df_final, file.path(out_dir, "okres_SPEI12_monthly_stats.xlsx"))

# =============================================================================
# 7. HELPER: okres value -> raster
# =============================================================================
okres_to_raster <- function(ids, values) {
  lookup <- rep(NA_real_, max(values(r_okres), na.rm = TRUE))
  lookup[ids] <- values
  classify(r_okres, cbind(seq_along(lookup), lookup))
}

# =============================================================================
# 8. RASTERS: MONTHLY + YEARLY ZONAL MEAN
# =============================================================================
dt_zonal <- as.data.table(df_final)
for (d in unique(dt_zonal$date)) {
  vals  <- dt_zonal[date == d]
  r_out <- okres_to_raster(vals$okres_id, vals$mean)
  names(r_out) <- paste0("spei12_mean_", d)
  writeRaster(r_out, file.path(out_month, sprintf("spei12_mean_%s.tif", d)),
              overwrite = TRUE)
}

df_zonal_year <- df_final |>
  group_by(okres_id, NAZEV, year, NUTS3) |>
  summarise(mean_spei   = mean(mean, na.rm = TRUE),
            median_spei = median(mean, na.rm = TRUE),
            min_spei    = min(mean, na.rm = TRUE),
            max_spei    = max(mean, na.rm = TRUE),
            .groups = "drop")

for (yr in sort(unique(df_zonal_year$year))) {
  vals  <- filter(df_zonal_year, year == yr)
  r_out <- okres_to_raster(vals$okres_id, vals$mean_spei)
  names(r_out) <- paste0("spei12_year_mean_", yr)
  writeRaster(r_out, file.path(out_year, sprintf("spei12_year_mean_%04d.tif", yr)),
              overwrite = TRUE)
}

plot(rast(file.path(out_year, "spei12_year_mean_2025.tif")),
     range = c(-3, 3), col = hcl.colors(100, "RdBu"))

# =============================================================================
# 9. ATTACH 2025 SPEI12 TO R_OKRES.csv  ->  R_OKRES_with_2025_attached.csv
# =============================================================================
okres_csv <- read_delim(okres_csv_file, 
                        delim = ";",
                        locale = locale(encoding = "windows-1250")) #,
                        #col_types = cols(.default = "c"))

# need to be joined by names, 
# as the okres_id and Item, neither the NUTS3 (kraj)/ NUTS4 fits (NUTS4 is missing from csv)
norm <- function(x) stringi::stri_replace_all_regex(
  tolower(stringi::stri_trans_general(as.character(x), "Latin-ASCII")), "[^a-z0-9]", "")

# Praha is named differently in the two tables, so give it a common key
k_csv <- okres_csv |>
  mutate(.k = if_else(NUTS4 == "CZ0100", "praha", norm(REGNAME)))

k_spei <- df_zonal_year |>
  mutate(NAZEV = if_else(NAZEV == "Ostrava-město", "Ostrava", NAZEV),
         .k    = if_else(NUTS3 == "CZ010", "praha", norm(NAZEV))) |>
  select(.k, NUTS3, 
         SPEI12_2025_mean = mean_spei, 
         SPEI12_2025_median = median_spei,
         SPEI12_2025_min = min_spei, 
         SPEI12_2025_max = max_spei,
         df_year = year)

out <- k_csv |>
  left_join(k_spei, by = c(".k", "NUTS3"), relationship = "one-to-one") |> #
  select(-.k, -df_year)



# head(okres_csv)
# head(df_zonal_year)
# out <- okres_csv |>
#   #mutate(.k = norm(REGNAME)) |>
#   left_join(df_zonal_year,by = c("ITEM" = "okres_id", "NUTS3")) 
# 
head(out)
tail(out)
nrow(out)

write.table(out, file.path(out_dir, "R_OKRES_spei_2025_added.csv"),
            sep = ";", dec = ".", row.names = FALSE, quote = FALSE,
            na = "", fileEncoding = "windows-1250")


# =============================================================================
# 10. JPEG MAPS OF ZONAL MEANS (plots, not rasters)
# =============================================================================
out_jpg <- file.path(out_dir, "SPEI12_maps_jpeg")
dir.create(out_jpg, showWarnings = FALSE, recursive = TRUE)

spei_cols <- hcl.colors(100, "RdBu")       # red = dry, blue = wet
spei_rng  <- c(-3, 3)                      # same scale on every map, so maps are comparable

plot_okres_map <- function(r, title, file) {
  jpeg(file, width = 2400, height = 1600, res = 300, quality = 95)
  par(mar = c(1, 1, 3, 1))
  plot(r, range = spei_rng, col = spei_cols, ext = ext(v_ok),
       main = title, axes = FALSE, box = FALSE,
       plg = list(title = "SPEI12", title.cex = 0.8))
  lines(v_ok, col = "grey25", lwd = 0.6)    # district outlines
  dev.off()
}

# --- monthly maps ------------------------------------------------------------
for (d in unique(dt_zonal$date)) {
  vals <- dt_zonal[date == d]
  r_m  <- okres_to_raster(vals$okres_id, vals$mean)
  plot_okres_map(r_m, paste0("SPEI12, mean per okres, ", sub("_", "-", d)),
                 file.path(out_jpg, sprintf("spei12_mean_%s.jpg", d)))
}

# --- yearly map --------------------------------------------------------------
for (yr in sort(unique(df_zonal_year$year))) {
  vals <- filter(df_zonal_year, year == yr)
  r_y  <- okres_to_raster(vals$okres_id, vals$mean_spei)
  plot_okres_map(r_y, paste0("SPEI12, Yearly mean per okres, ", yr),
                 file.path(out_jpg, sprintf("spei12_year_mean_%04d.jpg", yr)))
}


# =============================================================================
# 11. GIF OF MONTHLY MAPS
# =============================================================================
library(magick)

frames <- list.files(out_jpg, pattern = "^spei12_mean_\\d{4}_\\d{2}\\.jpg$",
                     full.names = TRUE)      # sorted by name = chronological
stopifnot(length(frames) > 0)

gif <- image_read(frames) |>
  image_scale("1200") |>                     # smaller width = smaller file
  image_animate(fps = 2, optimize = TRUE)  # 1.5 frames per second

image_write(gif, file.path(out_dir, "SPEI12_mean_okres_monthly.gif"))
