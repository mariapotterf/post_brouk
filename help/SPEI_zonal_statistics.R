# --- PACKAGES ---------------------------------------------------------------
library(terra)
library(dplyr)
library(stringr)
library(tibble)
library(readr)
library(fs)

# --- INPUTS -----------------------------------------------------------------
spei_folder <- "raw/SPEI12/spei12_harg1_split"  # 
okres_shp   <- "raw/CR_administrativa/OKRESY_P.shp"

# --- HELPERS ----------------------------------------------------------------
# Parse YYYY_MM_DD from filenames like "spei12_2000_01_01.tif"
get_year_month_day <- function(x) {
  m <- str_match(basename(x), "([12][0-9]{3})[-_ ]?(0?[1-9]|1[0-2])[-_ ]?(0?[1-9]|[12][0-9]|3[01])")
  tibble(
    year  = suppressWarnings(as.integer(m[,2])),
    month = suppressWarnings(as.integer(m[,3])),
    day   = suppressWarnings(as.integer(m[,4])),
    file  = x
  )
}

# Single-stat extractor (fast & memory-safe). fun can be mean, sd, median, min, max
zonal_stat <- function(r, v, fun) {
  terra::extract(r, v, fun = fun, na.rm = TRUE, ID = TRUE) |>
    as_tibble() |>
    rename(okres_id = ID, value = 2)
}

# --- LOAD DATA --------------------------------------------------------------
# list rasters (one per month)
files <- list.files(spei_folder, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)

meta <- get_year_month_day(files) |>
  dplyr::filter(!is.na(year), !is.na(month)) |>
  arrange(year, month, day)

# read polygons
v_ok <- vect(okres_shp)

# make sure CRS aligns: project polygons to the first raster's CRS if needed
r0 <- rast(meta$file[1])
if (crs(v_ok) != crs(r0)) {
  v_ok <- project(v_ok, r0)
}

# keep polygon attributes for later join
okres_attrib <- as_tibble(v_ok) |>
  mutate(okres_id = row_number())  # will match extract(ID=TRUE)

# Try to pick a human-readable name column (fallback to none)
name_col <- names(okres_attrib)[which(names(okres_attrib) %in% c("NAZEV", "NAZEV_OKRESU", "NAZ_LAU", "NAZEVNUTS", "NAZ_OKRES"))]
name_col <- if (length(name_col)) name_col[1] else NULL

# --- MONTHLY ZONAL STATS ----------------------------------------------------
monthly_stats <- lapply(seq_len(nrow(meta)), function(i) {
  r <- rast(meta$file[i])
  
  df_mean   <- zonal_stat(r, v_ok, mean)
  df_sd     <- zonal_stat(r, v_ok, sd)
  df_med    <- zonal_stat(r, v_ok, median)
  df_min    <- zonal_stat(r, v_ok, min)
  df_max    <- zonal_stat(r, v_ok, max)
  
  df <- df_mean |>
    rename(mean = value) |>
    left_join(rename(df_sd,   sd = value),   by = "okres_id") |>
    left_join(rename(df_med,  median = value), by = "okres_id") |>
    left_join(rename(df_min,  min = value), by = "okres_id") |>
    left_join(rename(df_max,  max = value), by = "okres_id") |>
    mutate(year = meta$year[i], month = meta$month[i])
  
  df
})

monthly_stats <- bind_rows(monthly_stats) |>
  left_join(okres_attrib, by = "okres_id") |>
  relocate(any_of(name_col), .after = okres_id) |>
  arrange(year, month, okres_id)

# Save monthly table
write_csv(monthly_stats, "outTable/okres_SPEI12_monthly_stats.csv")

# --- YEARLY ZONAL STATS (from monthly mean rasters) -------------------------
# If memory is tight, we compute yearly means on-the-fly
years <- sort(unique(meta$year))

yearly_stats <- lapply(years, function(yr) {
  f_yr <- meta$file[meta$year == yr]
  r    <- rast(f_yr)
  r_m  <- mean(r, na.rm = TRUE)      # yearly mean raster
  
  df_mean   <- zonal_stat(r_m, v_ok, mean)
  df_sd     <- zonal_stat(r_m, v_ok, sd)
  df_med    <- zonal_stat(r_m, v_ok, median)
  df_min    <- zonal_stat(r_m, v_ok, min)
  df_max    <- zonal_stat(r_m, v_ok, max)
  
  df <- df_mean |>
    rename(mean = value) |>
    left_join(rename(df_sd,   sd = value),   by = "okres_id") |>
    left_join(rename(df_med,  median = value), by = "okres_id") |>
    left_join(rename(df_min,  min = value), by = "okres_id") |>
    left_join(rename(df_max,  max = value), by = "okres_id") |>
    mutate(year = yr)
  
  df
})

yearly_stats <- bind_rows(yearly_stats) |>
  left_join(okres_attrib, by = "okres_id") |>
  relocate(any_of(name_col), .after = okres_id) |>
  arrange(year, okres_id)

# Save yearly table
write_csv(yearly_stats, "outTable/okres_SPEI12_yearly_stats.csv")

message("Written:",
        "\n - outTable/okres_SPEI12_monthly_stats.csv",
        "\n - outTable/okres_SPEI12_yearly_stats.csv")
