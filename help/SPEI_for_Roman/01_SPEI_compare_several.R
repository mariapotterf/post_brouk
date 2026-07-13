
#
#
# Loop over all SPEIs and make zonal stats and gifs
#
#

# --- PACKAGES ---------------------------------------------------------------
library(terra)
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(magick)
library(rasterVis)
library(RColorBrewer)
library(writexl)
library(fs)
library(purrr)

# --- PARAMETERS & INPUTS ---------------------------------------------------
spei_names <- c("spei1", "spei3", "spei6", "spei12")
method_suffix <- "harg1_split"  # folder suffix
base_spei_folder <- "raw"       # base parent folder for all SPEI subfolders
okres_shp_path <- "raw/CR_administrativa/OKRESY_P.shp"

# Output directories (you may want to create these in advance)
out_zonal_folder <- "outDataSPEI"              # for per‑SPEI zonal raster output subfolders
out_gif_folder <- "outGifsSPEI"                # all GIFs per SPEI
out_table_folder <- "outTableSPEI"             # for CSV / XLSX results

fs::dir_create(out_zonal_folder)
fs::dir_create(out_gif_folder)
fs::dir_create(out_table_folder)

# Read district shapefile
v_ok0 <- terra::vect(okres_shp_path)
ok_attrib0 <- as_tibble(v_ok0) %>% mutate(okres_id = row_number())
v_ok0$okres_id <- ok_attrib0$okres_id

# --- FUNCTIONS --------------------------------------------------------------

# Function to process one SPEI type
process_one_spei <- function(spei_name) {
  cat("Processing:", spei_name, "\n")
  
  spei_folder <- file.path(base_spei_folder, toupper(spei_name), paste0(spei_name, "_", method_suffix))
  
  # list all tif files
  files <- list.files(spei_folder, pattern = "\\.tif$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No files found for ", spei_name, " in folder: ", spei_folder)
  }
  
  # metadata table
  meta <- tibble(file = files) %>%
    mutate(
      name = basename(file),
      year = as.integer(str_extract(name, "(?<=_)(\\d{4})(?=_)")),
      month = as.integer(str_extract(name, "(?<=_\\d{4}_)(\\d{2})(?=_)"))
      # you can also parse day if needed
    ) %>%
    filter(!is.na(year), !is.na(month)) %>%
    arrange(year, month)
  
  # prepare district raster template & transform
  r0 <- rast(meta$file[1])
  v_ok <- v_ok0
  if (!terra::same.crs(r0, v_ok)) {
    v_ok <- terra::project(v_ok, r0)
  }
  r_okres <- terra::rasterize(v_ok, r0, field = "okres_id")
  
  # build raster stack of SPEI rasters
  r_stack <- rast(meta$file)
  names(r_stack) <- paste0(meta$year, "_", meta$month)
  r_all <- c(r_okres, r_stack)
  
  # extract values
  m <- values(r_all)
  dt <- as.data.table(m)
  # name first column
  setnames(dt, 1, "okres_id")
  
  # filter out NA in okres
  dt_f <- dt[!is.na(okres_id)]
  dt_long <- melt(
    dt_f,
    id.vars = "okres_id",
    variable.name = "date",
    value.name = "spei_value"
  )
  # split date
  dt_long[, c("year", "month") := tstrsplit(date, "_", fixed = TRUE)]
  dt_long[, `:=`(
    year = as.integer(year),
    month = as.integer(month)
  )]
  # compute zonal stats
  zonal <- dt_long[!is.na(spei_value),
                   .(
                     mean = mean(spei_value),
                     sd = sd(spei_value),
                     median = median(spei_value),
                     min = min(spei_value),
                     max = max(spei_value)
                   ),
                   by = .(okres_id, year, month)
  ]
  # join district names
  zonal <- zonal %>%
    left_join(ok_attrib0 %>% select(okres_id, NAZEV), by = "okres_id") %>%
    relocate(NAZEV, .after = okres_id) %>%
    arrange(year, month, okres_id)
  
  # also make annual zonal means per okres (optional)
  zonal_annual <- zonal %>%
    group_by(okres_id, NAZEV, year) %>%
    summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")
  
  # Save zonal table
  out_csv <- file.path(out_table_folder, paste0("zonal_stats_", spei_name, ".csv"))
  write.csv(zonal, out_csv, row.names = FALSE)
  write_xlsx(zonal, file.path(out_table_folder, paste0("zonal_stats_", spei_name, ".xlsx")))
  
  # --- Export zonal rasters (monthly) ---
  out_subfold <- file.path(out_zonal_folder, spei_name)
  fs::dir_create(out_subfold)
  
  unique_dates <- unique(zonal[, .(year, month)])
  for (i in seq_len(nrow(unique_dates))) {
    yr <- unique_dates$year[i]
    mo <- unique_dates$month[i]
    date_label <- sprintf("%04d_%02d", yr, mo)
    vals <- zonal[year == yr & month == mo, .(okres_id, mean)]
    
    # prepare lookup vector
    max_id <- max(r_okres[], na.rm = TRUE)
    lookup <- rep(NA_real_, max_id)
    lookup[vals$okres_id] <- vals$mean
    
    mat <- cbind(1:length(lookup), lookup)
    r_out <- classify(r_okres, mat)
    names(r_out) <- date_label
    
    out_path <- file.path(out_subfold, sprintf("%s_mean.tif", date_label))
    writeRaster(r_out, out_path, overwrite = TRUE)
  }
  
  # --- Make monthly GIF ---
  gif_out_path <- file.path(out_gif_folder, paste0("zonal_monthly_", spei_name, ".gif"))
  ras_files <- list.files(out_subfold, pattern = "_mean\\.tif$", full.names = TRUE)
  ras_files <- sort(ras_files)
  date_labels <- str_extract(basename(ras_files), "\\d{4}_\\d{2}")
  
  make_plot <- function(r_path, date_lab) {
    r_ter <- rast(r_path)
    r_r <- raster::raster(r_ter)
    img_tmp <- tempfile(fileext = ".png")
    png(img_tmp, width = 800, height = 800)
    plt <- levelplot(
      r_r,
      margin = FALSE,
      col.regions = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100),
      at = seq(-3, 3, length.out = 101),
      scales = list(draw = FALSE),
      main = paste0(spei_name, " mean ", date_lab),
      colorkey = list(space = "right", labels = list(at = seq(-3, 3, by = 1)))
    )
    print(plt)
    dev.off()
    magick::image_read(img_tmp)
  }
  
  frames <- mapply(make_plot, ras_files, date_labels, SIMPLIFY = FALSE)
  anim <- image_animate(image_join(frames), fps = 5)
  image_write(anim, gif_out_path)
  
  # --- Make yearly GIF ---
  gif_year_path <- file.path(out_gif_folder, paste0("zonal_yearly_", spei_name, ".gif"))
  # build yearly rasters
  yrs <- sort(unique(zonal$year))
  year_rasters <- lapply(yrs, function(yr) {
    dfyr <- zonal %>% filter(year == yr)
    # build vector per cell
    vals <- rep(NA_real_, terra::ncell(r_okres))
    match_id <- match(values(r_okres), dfyr$okres_id)
    valid <- !is.na(match_id)
    vals[valid] <- dfyr$mean[match_id[valid]]
    r_out <- setValues(r_okres, vals)
    names(r_out) <- as.character(yr)
    r_out
  })
  
  make_plot_year <- function(r, yr_lab) {
    r_r <- raster::raster(r)
    img_tmp <- tempfile(fileext = ".png")
    png(img_tmp, width = 800, height = 800)
    plt <- levelplot(
      r_r,
      margin = FALSE,
      col.regions = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100),
      at = seq(-3, 3, length.out = 101),
      scales = list(draw = FALSE),
      main = paste0(spei_name, " mean year ", yr_lab),
      colorkey = list(space = "right", labels = list(at = seq(-3, 3, by = 1)))
    )
    print(plt)
    dev.off()
    magick::image_read(img_tmp)
  }
  frames_y <- mapply(make_plot_year, year_rasters, yrs, SIMPLIFY = FALSE)
  anim_y <- image_animate(image_join(frames_y), fps = 2)
  image_write(anim_y, gif_year_path)
  
  # Return the zonal (monthly) table, with an added column indicating this SPEI name
  zonal <- zonal %>% rename(!!spei_name := mean)
  return(zonal)
}

# --- MAIN LOOP & MERGE -----------------------------------------------------

all_list <- vector("list", length(spei_names))
for (i in seq_along(spei_names)) {
  nm <- spei_names[i]
  all_list[[i]] <- process_one_spei(nm)
}

# Merge all monthly zonal stats by (okres_id, NAZEV, year, month)
df_all <- reduce(all_list, full_join, by = c("okres_id", "NAZEV", "year", "month"))

# Optionally reorder columns: e.g. NAZEV, okres_id, year, month, spei1, spei3, spei6, spei12
df_all <- df_all %>%
  dplyr::select(okres_id, NAZEV, year, month, all_of(spei_names))

# Save merged table
write.csv(df_all, file.path(out_table_folder, "zonal_stats_all_spei.csv"))
write_xlsx(df_all, file.path(out_table_folder, "zonal_stats_all_spei.xlsx"))

cat("All done. Merged zonal stats saved to:", file.path(out_table_folder, "zonal_stats_all_spei.csv"), "\n")

