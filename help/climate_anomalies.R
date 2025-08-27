# ===============================
# Multi-variable anomalies for tas, pr, vpd (terra)
# ===============================

# ---- Packages ----
library(terra)
library(dplyr)     # use dplyr::<fn> if ambiguity
library(stringr)
library(purrr)
library(readr)
library(tibble)

terraOptions(progress = 1)

<<<<<<< HEAD
# ---- Constants ----
=======
# ---- User paths / params ----
var_prefix <- "tas"   # set to "tas" or "pr"
>>>>>>> 553bb8a8dcc00ed758727caca4adeab52383f1dd
annual_dir <- "raw/clim_data_CZ_annual"
baseline_dir <- "raw/clim_data_CZ_reference_period"
out_dir <- "outData/anomalies_2000_2024"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_years <- 2000:2024
season_days <- 183L

# ---- Functions ----
align_to_baseline <- function(x, template) {
  if (!terra::same.crs(x, template)) {
    x <- terra::project(x, template, method = "bilinear")
  }
  if (!isTRUE(terra::compareGeom(x, template, stopOnError = FALSE))) {
    x <- terra::resample(x, template, method = "bilinear")
  }
  x
}

.compute_stats <- function(r, year, what, baseline_period) {
  get1 <- function(fn, ...) {
    out <- try(terra::global(r, fun = fn, na.rm = TRUE, ...), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out) || nrow(out) == 0) return(NA_real_)
    as.numeric(out[[1]][1])
  }
  tibble(
    baseline_period = baseline_period,
    year            = year,
    what            = what,
    min             = get1(min),
    max             = get1(max),
    median          = get1(median),
    mean            = get1(mean),
    sd              = get1(sd)
  )
}

process_variable <- function(var_prefix) {
  message("\n>>> Processing variable: ", var_prefix)
  
  baseline_periods <- c("1961-1990","1981-2010","1991-2020")
  baseline_files <- setNames(
    file.path(baseline_dir, paste0(var_prefix, "_", baseline_periods, "__04-09_mean.tif")),
    baseline_periods
  )
  
  compute_relative <- (var_prefix == "pr")
  rel_min_baseline <- if (var_prefix == "pr") 1 else 0
  
  annual_files <- list.files(
    annual_dir,
    pattern = paste0("^", var_prefix, ".*\\.tif(f)?$"),
    full.names = TRUE
  )
  if (length(annual_files) == 0) stop("No annual GeoTIFFs found for '", var_prefix, "'")
  
  annual_tbl <- tibble(
    file = annual_files,
    year = str_extract(basename(annual_files), "(?<!\\d)(19|20)\\d{2}(?!\\d)") |> as.integer()
  ) |>
    dplyr::filter(!is.na(year), year %in% target_years) |>
    dplyr::arrange(year)
  
  if (nrow(annual_tbl) == 0) stop("No ", var_prefix, " rasters for ", paste(target_years, collapse = ", "))
  
  stats_list <- list()
  
  for (bp in names(baseline_files)) {
    bpath <- baseline_files[[bp]]
    if (!file.exists(bpath)) stop("Baseline missing: ", bpath)
    
    baseline <- rast(bpath)
    if (is.na(crs(baseline))) stop("Baseline has no CRS: ", bpath)
    
    if (var_prefix == "pr") baseline <- baseline * season_days
    
    stats_list[[length(stats_list)+1]] <- .compute_stats(baseline, NA, "baseline", bp)
    
    abs_out_files <- rel_out_files <- character(nrow(annual_tbl))
    
    for (i in seq_len(nrow(annual_tbl))) {
      yr <- annual_tbl$year[i]
      ann <- rast(annual_tbl$file[i])
      ann_aligned <- align_to_baseline(ann, baseline)
      
      stats_list[[length(stats_list)+1]] <- .compute_stats(ann_aligned, yr, "annual", bp)
      
      anom_abs <- ann_aligned - baseline
      names(anom_abs) <- paste0(var_prefix, "_anom_abs_", yr, "_ref", bp)
      f_abs <- file.path(out_dir, paste0(var_prefix, "_anom_abs__ref", bp, "__", yr, ".tif"))
      writeRaster(anom_abs, filename = f_abs, overwrite = TRUE)
      abs_out_files[i] <- f_abs
      
      stats_list[[length(stats_list)+1]] <- .compute_stats(anom_abs, yr, "anomaly_abs", bp)
      
      if (compute_relative) {
        valid <- baseline > rel_min_baseline
        anom_rel <- ifel(valid, 100 * (ann_aligned - baseline) / baseline, NA)
        names(anom_rel) <- paste0(var_prefix, "_anom_relpct_", yr, "_ref", bp)
        f_rel <- file.path(out_dir, paste0(var_prefix, "_anom_relpct__ref", bp, "__", yr, ".tif"))
        writeRaster(anom_rel, filename = f_rel, overwrite = TRUE)
        rel_out_files[i] <- f_rel
        
        stats_list[[length(stats_list)+1]] <- .compute_stats(anom_rel, yr, "anomaly_rel_pct", bp)
      }
    }
    
    if (any(file.exists(abs_out_files))) {
      writeRaster(rast(abs_out_files[file.exists(abs_out_files)]),
                  filename = file.path(out_dir, paste0(var_prefix, "_anomalies_abs_",
                                                       min(target_years), "_", max(target_years),
                                                       "__ref", bp, "_stack.tif")),
                  overwrite = TRUE)
    }
    
    if (compute_relative && any(file.exists(rel_out_files))) {
      writeRaster(rast(rel_out_files[file.exists(rel_out_files)]),
                  filename = file.path(out_dir, paste0(var_prefix, "_anomalies_relpct_",
                                                       min(target_years), "_", max(target_years),
                                                       "__ref", bp, "_stack.tif")),
                  overwrite = TRUE)
    }
  }
  
  stats_tbl <- dplyr::bind_rows(stats_list) |>
    dplyr::mutate(variable = var_prefix) |>
    dplyr::relocate(variable, .before = baseline_period) |>
    dplyr::mutate(dplyr::across(c(min, max, median, mean), ~ round(.x, 3))) |>
    dplyr::arrange(baseline_period, dplyr::across(c(what, year)))
  
  csv_out <- file.path(out_dir, paste0(var_prefix, "_anomaly_stats_", min(target_years), "_", max(target_years), ".csv"))
  readr::write_csv(stats_tbl, csv_out)
  message("Saved stats: ", csv_out)
  
  return(stats_tbl)
  message("Stats in memory: ", stats_tbl)
}

# ===============================
# Run for each variable
# ===============================
all_stats <- list()
purrr::walk(c("tas", "pr", "vpd"), function(var) {
  tbl <- process_variable(var)
  all_stats[[var]] <<- tbl  # store it in the list
})
# ---- One combined CSV with all baselines & all stats ----
stats_tbl <- dplyr::bind_rows(stats_list) |>
  dplyr::mutate(variable = var_prefix) |>
  dplyr::relocate(variable, .before = baseline_period) |>
  dplyr::mutate(dplyr::across(c(min, max, median, mean), ~ round(.x, 3))) |>
  dplyr::arrange(baseline_period, dplyr::across(c(what, year)))



# ================================
# Make plots 
# ================================
library(ggplot2)

combined_stats <- dplyr::bind_rows(all_stats) |>
  dplyr::filter(what == "anomaly_abs") |>
  dplyr::select(variable, baseline_period, year, mean, sd)

combined_stats <- combined_stats |>
  dplyr::mutate(
    highlight = year %in% 2018:2020
  )


# Plot with highlighted years
ggplot(combined_stats, aes(x = year, y = mean)) +
  geom_point(aes(color = highlight), size = 2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = highlight), width = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "grey", linewidth = 0.7) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black"), guide = "none") +
  facet_wrap(variable ~ baseline_period, scales = "free_y") +
  labs(
    title = "Absolute Anomalies with Standard Deviation (2000–2024)",
    x = "Year",
    y = "Absolute Anomaly"
  ) +
  theme_bw()
