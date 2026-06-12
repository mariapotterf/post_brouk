
# ==============================================================================
#  Distance to nearest UNDISTURBED forest — CZECHIA ONLY
# ==============================================================================

# get data from Prospers - all types of forests (outside of clear-cuts)
# get disturbance patches
# remove disturbned patches form teh forest
# measure distance from teh subplot to nearest undisturbed forests

#disturbance_year_combined_compressed 


library(terra)
library(sf)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
proj_root   <- getwd() #"C:/Users/potterf/OneDrive - CZU v Praze/Dokumenty/2025_CZU_postbrouk/r_post_brouk/"
dist_path   <- "raw/from_prosper2026/disturbance_year_AEF_prosper"

forest_path <- "raw/from_prosper2026/Species_classification_2017_prosper"
gpkg_sub    <- "outDataShare/Karim_AEF/cleaned/env_chars_subplots.gpkg"

#dist_key <- "czechia"

# ── Load rasters ───────────────────────────────────────────────────────────────
disturb_r <- rast(file.path(dist_path, "disturbance_32633.tif"))
forest_r  <- rast(file.path(forest_path, "species_classification2017_test_prediction_v13_SVM.tif"))

# reproject raster
#disturb_r_proj <- terra::project(disturb_r, crs(forest_r), method = "near")
#writeRaster(disturb_r_proj, "raw/from_prosper2026/disturbance_year_32633_R.tif",
         #   overwrite = TRUE)

# next time just load the already-reprojected version
#disturb_r_proj <- rast("raw/from_prosper2026/disturbance_year_32633_R.tif")




cat("Disturbance CRS:", crs(disturb_r, describe = TRUE)$code, "\n")
cat("Forest CRS:     ", crs(forest_r,  describe = TRUE)$code, "\n")
cat("Disturbance res:", res(disturb_r), "\n")
cat("Forest res:     ", res(forest_r),  "\n")

# ── Load subplots — Czechia only ───────────────────────────────────────────────
all_subplots <- vect(file.path(gpkg_sub))
pts_cz       <- all_subplots[all_subplots$country_name == "Czech Republic",
                             c("subplot", "plot_id", "year", "disturbance_year")]
pts_cz       <- terra::project(pts_cz, disturb_r)

cat("Czechia subplots:", nrow(pts_cz), "\n")

# ==============================================================================
#  TEST — one point visually before full run
# ==============================================================================

test_point   <- pts_cz[10, ]
pt_dist_year <- as.numeric(test_point$disturbance_year) - 2000
buffer_test  <- 150

buff              <- buffer(test_point, buffer_test)
disturb_crop      <- crop(disturb_r, buff)
forest_crop       <- crop(forest_r,  buff)

plot(disturb_crop)
plot(forest_crop)

forest_crop_proj  <- if (!same.crs(forest_crop, disturb_crop)) {
  terra::project(forest_crop, crs(disturb_crop), method = "near")
} else forest_crop

forest_aligned <- terra::resample(forest_crop_proj, disturb_crop, method = "near")

# disturbance mask — adjust years if your CZ plots span wider window
dist_year_min <- 18
dist_year_max <- 25

disturb_mask      <- ifel(disturb_crop >= dist_year_min & disturb_crop <= dist_year_max, 1, NA)
# undisturbed forest = any forest class (0-6) AND not in disturbed window
undist_forest <- ifel(
  forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask),
  1, NA
)


dist_to_forest <- distance(undist_forest)
d <- terra::extract(dist_to_forest, test_point)[1, 2]  # extract distance 
cat("Test point distance to nearest undisturbed forest:", round(d, 0), "m\n")

# visual check
par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))
plot(disturb_crop,          main = paste0("1. Dist years", pt_dist_year),        legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
plot(forest_aligned,   main = "2. Forest (aligned)",         legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
plot(undist_forest,    main = "3. Undisturbed forest only",  legend = TRUE, col = "darkgreen")
plot(test_point, add = TRUE, col = "red",  pch = 19, cex = 1.5)
plot(dist_to_forest,        main = paste0("4. Distance to forest: ", round(d), " m"), legend = TRUE)
plot(test_point, add = TRUE, col = "red",  pch = 19, cex = 1.5)
par(mfrow = c(1, 1))



# ==============================================================================
#  FULL RUN — all Czechia subplots
# ==============================================================================
buffer_dist <- 1500

results_list <- lapply(seq_len(nrow(pts_cz)), function(i) {
  
  if (i %% 100 == 0) cat("  Point", i, "/", nrow(pts_cz), "\n")
  
  tryCatch({
    pt           <- pts_cz[i, ]
    pt_dist_year <- as.numeric(pt$disturbance_year) - 2000
    
    buff      <- buffer(pt, buffer_dist)
    d_crop    <- crop(disturb_r, buff)
    f_crop    <- crop(forest_r,  buff)
    
    f_crop_proj <- if (!same.crs(f_crop, d_crop)) {
      terra::project(f_crop, crs(d_crop), method = "near")
    } else f_crop
    
    forest_aligned <- terra::resample(f_crop_proj, d_crop, method = "near")
    
    # disturbance mask up to this point's disturbance year
    disturb_mask <- ifel(d_crop >= 18 & d_crop <= pt_dist_year, 1, NA)
    
    # undisturbed forest = classes 0-6 (not 9) AND not in disturbed window
    undist_forest <- ifel(
      forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask),
      1, NA
    )
    
    # disturbed forest = was forest (not 9) AND IS in disturbed window
    dist_forest <- ifel(
      forest_aligned != 9 & !is.na(forest_aligned) & !is.na(disturb_mask),
      1, NA
    )
    
    # ── pixel counts for proportions ──────────────────────────────
    total_pixels   <- ncell(forest_aligned)
    n_undist       <- global(undist_forest, "sum", na.rm = TRUE)[[1]]
    n_dist         <- global(dist_forest,   "sum", na.rm = TRUE)[[1]]
    pct_undist     <- round(n_undist / total_pixels * 100, 2)
    pct_dist       <- round(n_dist   / total_pixels * 100, 2)
    
    # ── expand buffer if no undisturbed forest found ───────────────
    if (n_undist == 0) {
      cat("    Point", i, "— expanding to 4000m\n")
      buff2          <- buffer(pt, buffer_dist * 2)
      d_crop         <- crop(disturb_r, buff2)
      f_crop         <- crop(forest_r,  buff2)
      f_crop_proj    <- if (!same.crs(f_crop, d_crop)) {
        terra::project(f_crop, crs(d_crop), method = "near")
      } else f_crop
      forest_aligned <- terra::resample(f_crop_proj, d_crop, method = "near")
      disturb_mask   <- ifel(d_crop >= 18 & d_crop <= pt_dist_year, 1, NA)
      undist_forest  <- ifel(
        forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask),
        1, NA
      )
      dist_forest    <- ifel(
        forest_aligned != 9 & !is.na(forest_aligned) & !is.na(disturb_mask),
        1, NA
      )
      total_pixels   <- ncell(forest_aligned)
      n_undist       <- global(undist_forest, "sum", na.rm = TRUE)[[1]]
      n_dist         <- global(dist_forest,   "sum", na.rm = TRUE)[[1]]
      pct_undist     <- round(n_undist / total_pixels * 100, 2)
      pct_dist       <- round(n_dist   / total_pixels * 100, 2)
    }
    
    # ── distances ──────────────────────────────────────────────────
    d_any <- terra::extract(distance(undist_forest), pt)[1, 2]
    
    data.frame(
      subplot          = as.character(pt$subplot),
      plot_id          = as.character(pt$plot_id),
      year             = as.numeric(pt$year),
      disturbance_year = as.numeric(pt$disturbance_year),
      dist_to_forest_m = round(d_any, 0),
      pct_undist_forest = pct_undist,
      pct_dist_forest   = pct_dist
    )
    
  }, error = function(e) {
    warning("Point ", i, " failed: ", e$message)
    data.frame(
      subplot           = as.character(pts_cz[i, ]$subplot),
      plot_id           = as.character(pts_cz[i, ]$plot_id),
      year              = as.numeric(pts_cz[i, ]$year),
      disturbance_year  = as.numeric(pts_cz[i, ]$disturbance_year),
      dist_to_forest_m  = NA_real_,
      pct_undist_forest = NA_real_,
      pct_dist_forest   = NA_real_
    )
  })
})

# ── combine and export ─────────────────────────────────────────
result_cz <- do.call(rbind, results_list)

cat("\nDone. Points:", nrow(result_cz), "\n")
cat("NAs dist_to_forest:", sum(is.na(result_cz$dist_to_forest_m)), "\n")
summary(result_cz)

data.table::fwrite(result_cz,
                   "outDataShare/Karim_AEF/cleaned/dist_to_forest_czechia.csv")

cat("Exported to outDataShare/Karim_AEF/cleaned/dist_to_forest_czechia.csv\n")