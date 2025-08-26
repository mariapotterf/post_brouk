# ===============================
# Multi-baseline anomalies for tas or pr (terra)
# - Absolute anomaly: year - baseline  (°C or mm)
# - Relative anomaly (optional): (year - baseline) / baseline * 100  (%)
# ===============================

# ---- Packages ----
library(terra)
library(dplyr)     # use dplyr::<fn> if ambiguity
library(stringr)
library(purrr)
library(readr)
library(tibble)

terraOptions(progress = 1)

# ---- User paths / params ----
var_prefix <- "tas"   # set to "tas" or "pr"
annual_dir <- "raw/clim_data_CZ_annual"
baseline_dir <- "raw/clim_data_CZ_reference_period"
out_dir <- "outData/anomalies_2018_2024"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Baseline files (pattern like "tas_1981-2010__04-09_mean.tif" / "pr_1981-2010__04-09_mean.tif")
baseline_periods <- c("1961-1990","1981-2010","1991-2020")
baseline_files <- setNames(
  file.path(baseline_dir, paste0(var_prefix, "_", baseline_periods, "__04-09_mean.tif")),
  baseline_periods
)

# Years to process
target_years <- 2018:2024
# Days in Apr–Sep (04–09). Needed because pr baseline is a *daily mean*.
season_days <- 183L

# Relative anomalies:
# - Typically **only meaningful for precipitation**. Default: on for pr, off for tas.
compute_relative <- (var_prefix == "pr")
# To avoid insane % where the baseline is ~0, mask pixels with very small baseline.
# Choose a sensible threshold for your units (e.g. 1 mm). Not used for tas.
rel_min_baseline <- if (var_prefix == "pr") 1 else 0

# ---- Discover annual rasters for the selected variable ----
annual_files <- list.files(
  annual_dir,
  pattern = paste0("^", var_prefix, ".*\\.tif(f)?$"),
  full.names = TRUE
)
if (length(annual_files) == 0) stop("No annual GeoTIFFs found for '", var_prefix, "' in: ", annual_dir)

# Extract the first 4-digit year found in the filename
annual_tbl <- tibble(
  file = annual_files,
  year = str_extract(basename(annual_files), "(?<!\\d)(19|20)\\d{2}(?!\\d)") |> as.integer()
) |>
  dplyr::filter(!is.na(year), year %in% target_years) |>
  dplyr::arrange(year)

if (nrow(annual_tbl) == 0) {
  stop("No ", var_prefix, " annual rasters for years ", paste(range(target_years), collapse = "-"),
       " found in: ", annual_dir)
}
missing_years <- setdiff(target_years, annual_tbl$year)
if (length(missing_years) > 0) warning("Missing rasters for years: ", paste(missing_years, collapse = ", "))

# ---- Helpers ----
align_to_baseline <- function(x, template) {
  if (!terra::same.crs(x, template)) {
    x <- terra::project(x, template, method = "bilinear")
  }
  if (!isTRUE(terra::compareGeom(x, template, stopOnError = FALSE))) {
    x <- terra::resample(x, template, method = "bilinear")
  }
  x
}

# robust global stats (min, max, median, mean)
.compute_stats <- function(r, year, what, baseline_period) {
  get1 <- function(fn, ...) {
    out <- try(terra::global(r, fun = fn, na.rm = TRUE, ...), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out) || nrow(out) == 0) return(NA_real_)
    as.numeric(out[[1]][1])
  }
  v_min  <- get1(min)
  v_max  <- get1(max)
  v_mean <- get1(mean)
  v_med  <- get1(median)
  if (is.na(v_med)) v_med <- get1(quantile, probs = 0.5)
  
  tibble(
    baseline_period = baseline_period,
    year            = year,
    what            = what,      # "baseline","annual","anomaly_abs","anomaly_rel_pct"
    min             = v_min,
    max             = v_max,
    median          = v_med,
    mean            = v_mean
  )
}

# ===============================
# Main: loop over baselines
# ===============================
stats_list  <- list()

for (bp in names(baseline_files)) {
  bpath <- baseline_files[[bp]]
  if (!file.exists(bpath)) stop("Baseline raster not found: ", bpath)
  
  message("=== Baseline: ", bp, " ===")
  baseline <- rast(bpath)
  if (is.na(crs(baseline))) stop("Baseline '", bp, "' has no CRS. Please define it.")
  
  # for precipitation, baseline is per-day mean; scale to Apr–Sep totals (183 days)
  if (var_prefix == "pr") {
    message("Scaling precipitation baseline from daily mean to Apr–Sep seasonal totals (x ", season_days, ").")
    baseline <- baseline * season_days
  }
  
  # Baseline stats (once per baseline)
  stats_list[[length(stats_list)+1]] <- .compute_stats(baseline, year = NA_integer_,
                                                       what = "baseline", baseline_period = bp)
 
  abs_out_files <- vector("character", nrow(annual_tbl))
  rel_out_files <- vector("character", nrow(annual_tbl))
  
  for (i in seq_len(nrow(annual_tbl))) {
    yr   <- annual_tbl$year[i]
    f_in <- annual_tbl$file[i]
    message("Processing year ", yr, " for baseline ", bp, " | ", basename(f_in))
    
    ann <- rast(f_in)
    ann_aligned <- align_to_baseline(ann, baseline)
    
    # Annual stats (aligned)
    stats_list[[length(stats_list)+1]] <- .compute_stats(ann_aligned, year = yr,
                                                         what = "annual", baseline_period = bp)
    
    # ---- Absolute anomaly (same units) ----
    anom_abs <- ann_aligned - baseline
    names(anom_abs) <- paste0(var_prefix, "_anom_abs_", yr, "_ref", bp)
    f_abs <- file.path(out_dir, paste0(var_prefix, "_anom_abs__ref", bp, "__", yr, ".tif"))
    writeRaster(anom_abs, filename = f_abs, overwrite = TRUE)
    abs_out_files[i] <- f_abs
    
    stats_list[[length(stats_list)+1]] <- .compute_stats(anom_abs, year = yr,
                                                         what = "anomaly_abs", baseline_period = bp)
    
    # ---- Relative anomaly (%), optional (recommended for precipitation) ----
    if (isTRUE(compute_relative)) {
      # mask where baseline <= threshold to avoid division by (near) zero
      # NOTE: units of rel_min_baseline should match the baseline units (e.g., mm)
      valid <- baseline > rel_min_baseline
      anom_rel <- ifel(valid, 100 * (ann_aligned - baseline) / baseline, NA)
      names(anom_rel) <- paste0(var_prefix, "_anom_relpct_", yr, "_ref", bp)
      
      f_rel <- file.path(out_dir, paste0(var_prefix, "_anom_relpct__ref", bp, "__", yr, ".tif"))
      writeRaster(anom_rel, filename = f_rel, overwrite = TRUE)
      rel_out_files[i] <- f_rel
      
      stats_list[[length(stats_list)+1]] <- .compute_stats(anom_rel, year = yr,
                                                           what = "anomaly_rel_pct", baseline_period = bp)
    }
  }
  
  # Optional: multiband stacks (per baseline)
  abs_exist <- abs_out_files[file.exists(abs_out_files)]
  if (length(abs_exist) > 0) {
    writeRaster(rast(abs_exist),
                filename = file.path(out_dir, paste0(var_prefix, "_anomalies_abs_",
                                                     min(target_years), "_", max(target_years),
                                                     "__ref", bp, "_stack.tif")),
                overwrite = TRUE)
  }
  
  if (isTRUE(compute_relative)) {
    rel_exist <- rel_out_files[file.exists(rel_out_files)]
    if (length(rel_exist) > 0) {
      writeRaster(rast(rel_exist),
                  filename = file.path(out_dir, paste0(var_prefix, "_anomalies_relpct_",
                                                       min(target_years), "_", max(target_years),
                                                       "__ref", bp, "_stack.tif")),
                  overwrite = TRUE)
    }
  }
}

# ---- One combined CSV with all baselines & all stats ----
stats_tbl <- dplyr::bind_rows(stats_list) |>
  dplyr::mutate(variable = var_prefix) |>
  dplyr::relocate(variable, .before = baseline_period) |>
  dplyr::mutate(dplyr::across(c(min, max, median, mean), ~ round(.x, 3))) |>
  dplyr::arrange(baseline_period, dplyr::across(c(what, year)))

csv_out <- file.path(out_dir, paste0(var_prefix, "_anomaly_stats_", min(target_years), "_", max(target_years), ".csv"))
readr::write_csv(stats_tbl, csv_out)

print(stats_tbl)
message("Saved combined stats: ", csv_out)
message("Anomaly rasters in: ", out_dir)
