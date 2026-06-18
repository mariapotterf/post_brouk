# ==============================================================================
#  Distance to nearest UNDISTURBED forest — CZECHIA ONLY
# ==============================================================================

library(terra)
library(sf)
library(dplyr)
library(landscapemetrics)
library(stars)        # ── NEW: needed for landscapemetrics

# ── Paths ──────────────────────────────────────────────────────────────────────
proj_root   <- getwd()
dist_path   <- "raw/from_prosper2026/disturbance_year_AEF_prosper"
forest_path <- "raw/from_prosper2026/Species_classification_2017_prosper"
gpkg_sub    <- "outDataShare/Karim_AEF/cleaned/env_chars_subplots.gpkg"

# ── Load rasters ───────────────────────────────────────────────────────────────
disturb_r <- rast(file.path(dist_path, "disturbance_32633.tif"))
forest_r  <- rast(file.path(forest_path, "species_classification2017_test_prediction_v13_SVM.tif"))

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

test_point   <- pts_cz[160, ]

test_point <- pts_cz[pts_cz$subplot == "13_15_126_4", ]
buffer_test  <- 500  # ── NEW: changed to 500m (your actual analysis buffer)

# ── NEW: extract disturbance year from raster at point location
dist_year_18_25 <- terra::extract(disturb_r, test_point)[1, 2]
if (is.na(dist_year_18_25)) {
  dist_year_18_25 <- as.numeric(test_point$disturbance_year) - 2000
  cat("Fallback to table disturbance year:", dist_year_18_25, "\n")
}
cat("Disturbance year from raster:", dist_year_18_25, "\n")
cat("Disturbance year from table: ", as.numeric(test_point$disturbance_year) - 2000, "\n")

buff         <- buffer(test_point, buffer_test)
disturb_crop <- crop(disturb_r, buff)
forest_crop  <- crop(forest_r,  buff)

forest_crop_proj <- if (!same.crs(forest_crop, disturb_crop)) {
  terra::project(forest_crop, crs(disturb_crop), method = "near")
} else forest_crop

forest_aligned <- terra::resample(forest_crop_proj, disturb_crop, method = "near")

# ── NEW: use raster-extracted year as cutoff (not fixed 18-25)
disturb_mask <- ifel(disturb_crop >= 18 & disturb_crop <= dist_year_18_25, 1, NA)

undist_forest <- ifel(
  forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask),
  1, NA
)

# ── NEW: disturbed forest layer
dist_forest <- ifel(
  forest_aligned != 9 & !is.na(forest_aligned) & !is.na(disturb_mask),
  1, NA
)

# ── NEW: pixel counts
total_pixels <- ncell(forest_aligned)
n_undist     <- global(undist_forest, "sum", na.rm = TRUE)[[1]]
n_dist       <- global(dist_forest,   "sum", na.rm = TRUE)[[1]]
cat("Undisturbed forest pixels:", n_undist, "(", round(n_undist/total_pixels*100,1), "%)\n")
cat("Disturbed forest pixels:  ", n_dist,   "(", round(n_dist  /total_pixels*100,1), "%)\n")

# ── NEW: spruce composition of surrounding undisturbed forest
n_undist_spruce <- global(
  ifel(forest_aligned == 3 & is.na(disturb_mask), 1, NA),
  "sum", na.rm = TRUE
)[[1]]
pct_spruce_in_undist <- ifelse(n_undist > 0, round(n_undist_spruce / n_undist * 100, 2), NA)
cat("Spruce in undisturbed forest:", pct_spruce_in_undist, "%\n")

dist_to_forest <- distance(undist_forest)
d <- terra::extract(dist_to_forest, test_point)[1, 2]
cat("Distance to nearest undisturbed forest:", round(d, 0), "m\n")

# ── NEW: patch metrics
patch_r       <- ifel(!is.na(dist_forest),  1L,
                      ifel(!is.na(undist_forest), 2L, NA))
patch_r_stars <- stars::st_as_stars(patch_r)
check_landscape(patch_r_stars)  # check it works before full run

# disturbed area
n_cells_class1 <- global(ifel(patch_r == 1, 1, NA), "sum", na.rm = TRUE)[[1]]
cat("Disturbed forest area in buffer (ha):", round(n_cells_class1 * 100 / 10000, 2), "\n")

# undisturbed area for comparison
n_cells_class2 <- global(ifel(patch_r == 2, 1, NA), "sum", na.rm = TRUE)[[1]]
cat("Undisturbed forest area in buffer (ha):", round(n_cells_class2 * 100 / 10000, 2), "\n")

# total buffer area
total_ha <- ncell(patch_r) * 100 / 10000
cat("Total buffer area (ha):", round(total_ha, 2), "\n")

# ratio
cat("Disturbed / total forest (%):", 
    round(n_cells_class1 / (n_cells_class1 + n_cells_class2) * 100, 1), "\n")


lm_np    <- lsm_c_np(patch_r_stars)
lm_area  <- lsm_c_area_mn(patch_r_stars)
lm_shape <- lsm_c_shape_mn(patch_r_stars)
lm_para <- lsm_p_para(patch_r_stars)
para_disturbed <- mean(lm_para$value[lm_para$class == 1], na.rm = TRUE)
cat("Mean perimeter-area ratio (disturbed):", para_disturbed, "\n")

# or total disturbed area in buffer — simpler and ecologically direct
cat("Disturbed area in buffer (ha):", n_cells_class1 * 100 / 10000, "\n")

cat("N disturbed patches:",     lm_np$value[lm_np$class == 1], "\n")
cat("Mean patch area (ha):",    lm_area$value[lm_area$class == 1], "\n")
cat("Edge density:",            lm_ed$value[lm_ed$class == 1], "\n")
cat("Mean shape index:",        lm_shape$value[lm_shape$class == 1], "\n")

# visual check
par(mfrow = c(2, 3), mar = c(2, 2, 3, 1))
plot(disturb_crop,     main = paste0("1. Disturbance (cutoff: ", dist_year_18_25, ")"), legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
plot(forest_aligned,   main = "2. Forest species classes", legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
plot(disturb_mask,     main = "3. Masked disturbances", col = "orange", legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
plot(undist_forest,    main = "4. Undisturbed forest", col = "darkgreen", legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
plot(dist_forest,      main = "5. Disturbed forest", col = "firebrick", legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
plot(dist_to_forest,   main = paste0("6. Distance to forest: ", round(d), "m"), legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
par(mfrow = c(1, 1))


# check how patches are actually calculated
# verify patch_r classes are correct
cat("Cells with class 1 (disturbed forest):", 
    global(ifel(patch_r == 1, 1, NA), "sum", na.rm = TRUE)[[1]], "\n")
cat("Cells with class 2 (undisturbed forest):", 
    global(ifel(patch_r == 2, 1, NA), "sum", na.rm = TRUE)[[1]], "\n")

# cross-check: does class 1 in patch_r match dist_forest?
cat("Cells in dist_forest:   ", global(dist_forest,   "sum", na.rm = TRUE)[[1]], "\n")
cat("Cells in undist_forest: ", global(undist_forest, "sum", na.rm = TRUE)[[1]], "\n")

# plot separately to confirm
par(mfrow = c(1, 3))
plot(dist_forest,   main = "dist_forest (source)",  col = "firebrick")
plot(undist_forest, main = "undist_forest (source)", col = "darkgreen")
plot(patch_r,       main = "patch_r (combined)",     col = c("firebrick", "darkgreen"))
par(mfrow = c(1, 1))


# count separate patches
# count separate patches per class
patches <- patches(patch_r, directions = 8)  # 8-connectivity (includes diagonals)
plot(patches, main = "Individual patches (each colour = one patch)")

# how many patches per class?
patch_vals <- as.data.frame(c(patches, patch_r), na.rm = TRUE)
names(patch_vals) <- c("patch_id", "class")

patch_vals %>%
  group_by(class) %>%
  summarise(n_patches = n_distinct(patch_id))



# ==============================================================================
#  FULL RUN — all Czechia subplots
# ==============================================================================
buffer_dist <- 500

results_list <- lapply(seq_len(nrow(pts_cz)), function(i) {
  
  if (i %% 100 == 0) cat("  Point", i, "/", nrow(pts_cz), "\n")
  
  tryCatch({
    pt <- pts_cz[i, ]
    
    buff      <- buffer(pt, buffer_dist)
    d_crop    <- crop(disturb_r, buff)
    f_crop    <- crop(forest_r,  buff)
    
    f_crop_proj    <- if (!same.crs(f_crop, d_crop)) {
      terra::project(f_crop, crs(d_crop), method = "near")
    } else f_crop
    forest_aligned <- terra::resample(f_crop_proj, d_crop, method = "near")
    
    # extract disturbance year from raster at point, fallback to table
    dist_year_18_25 <- terra::extract(d_crop, pt)[1, 2]
    if (is.na(dist_year_18_25)) {
      dist_year_18_25 <- as.numeric(pt$disturbance_year) - 2000
    }
    
    # masks using point-specific disturbance year cutoff
    disturb_mask  <- ifel(d_crop >= 18 & d_crop <= dist_year_18_25, 1, NA)
    undist_forest <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask),  1, NA)
    dist_forest   <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & !is.na(disturb_mask), 1, NA)
    
    # pixel counts
    total_pixels <- ncell(forest_aligned)
    n_undist     <- global(undist_forest, "sum", na.rm = TRUE)[[1]]
    n_dist       <- global(dist_forest,   "sum", na.rm = TRUE)[[1]]
    pct_undist   <- round(n_undist / total_pixels * 100, 2)
    pct_dist     <- round(n_dist   / total_pixels * 100, 2)
    
    # expand buffer if no undisturbed forest found
    if (n_undist == 0) {
      cat("    Point", i, "— expanding to 1000m\n")
      buff2          <- buffer(pt, buffer_dist * 2)
      d_crop         <- crop(disturb_r, buff2)
      f_crop         <- crop(forest_r,  buff2)
      f_crop_proj    <- if (!same.crs(f_crop, d_crop)) {
        terra::project(f_crop, crs(d_crop), method = "near")
      } else f_crop
      forest_aligned <- terra::resample(f_crop_proj, d_crop, method = "near")
      disturb_mask   <- ifel(d_crop >= 18 & d_crop <= dist_year_18_25, 1, NA)
      undist_forest  <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask),  1, NA)
      dist_forest    <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & !is.na(disturb_mask), 1, NA)
      total_pixels   <- ncell(forest_aligned)
      n_undist       <- global(undist_forest, "sum", na.rm = TRUE)[[1]]
      n_dist         <- global(dist_forest,   "sum", na.rm = TRUE)[[1]]
      pct_undist     <- round(n_undist / total_pixels * 100, 2)
      pct_dist       <- round(n_dist   / total_pixels * 100, 2)
    }
    
    # distance to nearest undisturbed forest
    d_any <- terra::extract(distance(undist_forest), pt)[1, 2]
    
    # spruce composition of surrounding undisturbed forest
    n_undist_spruce      <- global(ifel(forest_aligned == 3 & is.na(disturb_mask), 1, NA), 
                                   "sum", na.rm = TRUE)[[1]]
    pct_spruce_in_undist <- ifelse(n_undist > 0, 
                                   round(n_undist_spruce / n_undist * 100, 2), 
                                   NA_real_)
    
    # patch metrics
    get_val <- function(lm_df, cls = 1) {
      v <- lm_df$value[lm_df$class == cls]
      if (length(v) == 0) NA_real_ else v
    }
    
    patch_r       <- ifel(!is.na(dist_forest),  1L,
                          ifel(!is.na(undist_forest), 2L, NA))
    patch_r_stars <- stars::st_as_stars(patch_r)
    
    patch_n       <- get_val(lsm_c_np(patch_r_stars),      1)
    patch_area_mn <- get_val(lsm_c_area_mn(patch_r_stars), 1)
    patch_para    <- mean(lsm_p_para(patch_r_stars)$value[
      lsm_p_para(patch_r_stars)$class == 1], na.rm = TRUE)
    patch_area_ha <- round(n_dist * 100 / 10000, 2)  # ── use n_dist directly
    
    # ── use already-computed n_dist and n_undist — no redundant global() calls
    pct_dist_of_forest <- ifelse(
      (n_dist + n_undist) > 0,
      round(n_dist / (n_dist + n_undist) * 100, 1),
      NA_real_
    )
    
    data.frame(
      subplot              = as.character(pt$subplot),
      plot_id              = as.character(pt$plot_id),
      year                 = as.numeric(pt$year),
      disturbance_year     = as.numeric(pt$disturbance_year),
      dist_year_18_25      = dist_year_18_25,
      dist_to_forest_m     = round(d_any, 0),
      pct_undist_forest    = pct_undist,
      pct_dist_forest      = pct_dist,
      pct_spruce_in_undist = pct_spruce_in_undist,
      patch_n_disturbed    = patch_n,
      patch_area_mn_ha     = patch_area_mn,
      patch_area_ha        = patch_area_ha,
      patch_para           = patch_para,
      pct_dist_of_forest   = pct_dist_of_forest
    )
    
  }, error = function(e) {
    warning("Point ", i, " failed: ", e$message)
    data.frame(
      subplot              = as.character(pts_cz[i, ]$subplot),
      plot_id              = as.character(pts_cz[i, ]$plot_id),
      year                 = as.numeric(pts_cz[i, ]$year),  # year of field works
      disturbance_year     = as.numeric(pts_cz[i, ]$disturbance_year), # from Cornelius
      dist_year_18_25      = NA_real_,  # disturbance year from Prosper raster (2-digit, 18-25)
      dist_to_forest_m     = NA_real_,  # distance to nearest undisturbed forest [m], 0-197m
      pct_undist_forest    = NA_real_,  # % of 100ha buffer (1 km2) = undisturbed forest, 6.8-96.5%
      pct_dist_forest      = NA_real_,  # % of 100ha buffer = disturbed forest, 0.1-78.7%
      pct_spruce_in_undist = NA_real_,  # % spruce within undisturbed forest (seed source), 22-100%
      patch_n_disturbed    = NA_real_,  # n separate disturbed patches in buffer, 1-41
      patch_area_mn_ha     = NA_real_,  # mean disturbed patch size [ha], 0.03-76.6ha
      patch_area_ha        = NA_real_,  # total disturbed forest area in buffer [ha], 0.1-78.7ha
      patch_para           = NA_real_,  # mean patch shape complexity (PARA), 0.02-0.31 - 0 - circular, 1 = compex shape/elongated, ..
      pct_dist_of_forest   = NA_real_   # disturbed % of all forest in buffer, 0.6-86.3%
    )
  })
})

# ── combine and export ─────────────────────────────────────────────────────────
result_cz <- do.call(rbind, results_list)

cat("\nDone. Points:", nrow(result_cz), "\n")
cat("NAs dist_to_forest:", sum(is.na(result_cz$dist_to_forest_m)), "\n")
summary(result_cz)
View(result_cz)

result_cz %>% 
  arrange(desc(patch_n_disturbed)) %>%
  select(subplot, plot_id, disturbance_year, patch_n_disturbed, 
         patch_area_ha, patch_area_mn_ha, dist_to_forest_m) %>%
  head(20)



data.table::fwrite(result_cz,
                   "outData/dist_and_patch_metrics_czechia.csv")
#cat("Exported to outDataShare/Karim_AEF/cleaned/dist_and_patch_metrics_czechia.csv\n")