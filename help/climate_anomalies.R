# ===============================
# Temperature anomaly script
# Baseline: 1991–2010 mean raster
# Annual means: 2000–2024 (one raster per year)
# Output: anomalies for 2018–2024 (annual - baseline): 
# In climate science, anomalies are almost always defined as: anomaly =  observed (year) − reference_mean (result)
# ===============================

# ---- Packages ----
library(terra)          # spatial rasters
library(dplyr)          # data wrangling (use dplyr:: where ambiguous)
library(stringr)        # filename parsing
library(purrr)          # iteration

# ---- User paths (edit these to your setup) ----
baseline_path <- "raw/clim_data_CZ_reference_period/tas_1991-2020__04-09_mean.tif"
annual_dir    <- "raw/clim_data_CZ_annual/output_annual_means"   # folder with annual rasters, e.g., temp_mean_2000.tif
out_dir       <- "outData/anomalies_2018_2024"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# (Optional) Terra settings for progress/temp files
terraOptions(progress = 1)

# ---- Load baseline ----
baseline <- rast(baseline_path)

if (!file.exists(baseline_path)) {
  stop("Baseline raster not found at: ", baseline_path)
}

# Ensure baseline has a defined CRS
if (is.na(crs(baseline))) {
  stop("Baseline raster has no CRS. Please define its CRS before proceeding.")
}

# ---- Discover annual rasters and extract years ----
annual_files <- list.files(annual_dir, pattern = "^tas.*\\.tif(f)?$", full.names = TRUE)  # only tas rasters, 

if (length(annual_files) == 0) {
  stop("No annual GeoTIFFs found in: ", annual_dir)
}

# Extract a 4-digit year from filenames
annual_tbl <- tibble::tibble(
  file = annual_files,
  year = str_extract(basename(annual_files), "(?<!\\d)(19|20)\\d{2}(?!\\d)") |> as.integer()
) |>
  dplyr::filter(!is.na(year)) |>
  dplyr::arrange(year)

# Keep only 2018–2024
target_years <- 2018:2024
annual_tbl <- dplyr::filter(annual_tbl, year %in% target_years)

if (nrow(annual_tbl) == 0) {
  stop("No annual rasters for years 2018–2024 were found in: ", annual_dir)
}

# Warn about any missing years in the range
missing_years <- setdiff(target_years, annual_tbl$year)
if (length(missing_years) > 0) {
  warning("Missing rasters for years: ", paste(missing_years, collapse = ", "))
}

# ---- Helper: align a raster to the baseline grid ----
# This handles CRS differences and grid alignment.
align_to_baseline <- function(x, template) {
  # If CRS differs, project first (bilinear for continuous temperature)
  if (!terra::same.crs(x, template)) {
    x <- terra::project(x, template, method = "bilinear")
  }
  # If geometry differs (extent/resolution/grid), resample to template
  if (!terra::compareGeom(x, template, stopOnError = FALSE, rowcol = TRUE, crs = TRUE, ext = TRUE, res = TRUE, orig = TRUE, rotation = TRUE)) {
    x <- terra::resample(x, template, method = "bilinear")
  }
  # If extents don't perfectly overlap, crop to common area to avoid NA expansion
  # (Optional) You can mask to the template’s footprint if desired:
  # x <- terra::mask(x, template)
  x
}

# ---- Helper: align a raster to the baseline grid (terra-safe) ----
align_to_baseline <- function(x, template) {
  # 1) CRS
  if (!terra::same.crs(x, template)) {
    x <- terra::project(x, template, method = "bilinear")
  }
  # 2) Grid (nrow/ncol/res/origin/rotation). compareGeom with stopOnError only.
  #    If FALSE, resample to the template grid.
  if (!isTRUE(terra::compareGeom(x, template, stopOnError = FALSE))) {
    x <- terra::resample(x, template, method = "bilinear")
  }
  # 3) Optional: clip/mask to template footprint to avoid edge NAs growing
  # x <- terra::mask(x, template)
  x
}

# ---- Compute anomalies and write outputs ----
# ---- Compute anomalies, write outputs, and collect min/max stats ----

# helper to compute global min/max (NA-safe) and return a tibble row
# safer helper: compute min / max / median via separate global() calls
# Robust helper: min, max, median, mean
.compute_stats <- function(r, year, what) {
  get1 <- function(fn, ...) {
    out <- try(terra::global(r, fun = fn, na.rm = TRUE, ...), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out) || nrow(out) == 0) return(NA_real_)
    as.numeric(out[[1]][1])
  }
  
  # basic stats
  v_min  <- get1(min)
  v_max  <- get1(max)
  v_mean <- get1(mean)
  
  # median with fallback to quantile(0.5)
  v_med <- get1(median)
  if (is.na(v_med)) {
    v_med <- get1(quantile, probs = 0.5)
  }
  
  tibble::tibble(
    year   = year,
    what   = what,      # "baseline", "annual", "anomaly"
    min    = v_min,
    max    = v_max,
    median = v_med,
    mean   = v_mean
  )
}



# stats container
stats_list <- list()

# baseline stats (once)
stats_list[[length(stats_list)+1]] <- .compute_stats(baseline, year = NA_integer_, what = "baseline")

out_files <- vector("character", nrow(annual_tbl))

for (i in seq_len(nrow(annual_tbl))) {
  yr   <- annual_tbl$year[i]
  f_in <- annual_tbl$file[i]
  
  message("Processing year: ", yr, " | file: ", basename(f_in))
  ann <- rast(f_in)
  
  # Align annual raster to baseline grid/CRS
  ann_aligned <- align_to_baseline(ann, baseline)
  
  # Annual stats (aligned)
  stats_list[[length(stats_list)+1]] <- .compute_minmax(ann_aligned, year = yr, what = "annual")
  
  # Simple anomaly (annual - 1991–2010 mean); NA propagates automatically
  anom <- ann_aligned - baseline
  names(anom) <- paste0("anom_", yr)
  
  # Anomaly stats:
  stats_list[[length(stats_list)+1]] <- .compute_stats(ann_aligned, year = yr, what = "annual")
  stats_list[[length(stats_list)+1]] <- .compute_stats(anom,        year = yr, what = "anomaly")
  
  # Write anomaly GeoTIFF
  f_out <- file.path(out_dir, paste0("temp_anomaly_", yr, ".tif"))
  writeRaster(anom, filename = f_out, overwrite = TRUE)
  out_files[i] <- f_out
}

# bind all stats rows
stats_tbl <- dplyr::bind_rows(stats_list) |>
  dplyr::arrange(what, year) %>% 
  na.omit()

# write the table next to rasters
stats_csv <- file.path(out_dir, "temp_anomaly_minmax_2018_2024.csv")
readr::write_csv(stats_tbl, stats_csv)

# quick peek
print(stats_tbl)
message("Saved stats table: ", stats_csv)

# ---- (Optional) Build a multiband stack of all anomalies ----
# Only for the files that were successfully created
anoms_exist <- out_files[file.exists(out_files)]
if (length(anoms_exist) > 0) {
  anom_stack <- rast(anoms_exist)
  # Save as a single multi-layer file
  stack_out <- file.path(out_dir, "temp_anomalies_2018_2024_stack.tif")
  writeRaster(anom_stack, filename = stack_out, overwrite = TRUE)
}

# ---- (Optional) Quick check plot ----
# Plot a single year to visually confirm (change year as needed)
check_year <- 2019
check_file <- file.path(out_dir, paste0("temp_anomaly_", check_year, ".tif"))
if (file.exists(check_file)) {
  plot(rast(check_file), main = paste0("Temperature anomaly (", check_year, ") | annual - 1991–2010 mean"))
}

