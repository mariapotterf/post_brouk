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

test_point <- pts_cz[pts_cz$subplot == "13_15_106_4", ]
test_point <- pts_cz[5, ]
buffer_test <- 500

cat("Plot:", test_point$plot_id, "\n")
cat("Disturbance year (table):", test_point$disturbance_year, "\n")
cat("Sampling year:", test_point$year, "\n")

buff         <- buffer(test_point, buffer_test)
disturb_crop <- crop(disturb_r, buff)
forest_crop  <- crop(forest_r,  buff)

forest_crop_proj <- if (!same.crs(forest_crop, disturb_crop)) {
  terra::project(forest_crop, crs(disturb_crop), method = "near")
} else forest_crop

forest_aligned <- terra::resample(forest_crop_proj, disturb_crop, method = "near")

# ── T1: disturbance year cutoff ────────────────────────────────────────────
dist_year_18_25 <- terra::extract(disturb_crop, test_point)[1, 2]
if (is.na(dist_year_18_25)) {
  dist_year_18_25 <- as.numeric(test_point$disturbance_year) - 2000
}
cat("\nT1 cutoff (disturbance year, 2-digit):", dist_year_18_25, "\n")

disturb_mask_t1  <- ifel(disturb_crop >= 18 & disturb_crop <= dist_year_18_25, 1, NA)
undist_forest_t1 <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask_t1),  1, NA)
dist_forest_t1   <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & !is.na(disturb_mask_t1), 1, NA)

n_undist_t1 <- global(undist_forest_t1, "sum", na.rm = TRUE)[[1]]
n_dist_t1   <- global(dist_forest_t1,   "sum", na.rm = TRUE)[[1]]
total_pixels <- ncell(forest_aligned)

dist_to_forest_t1 <- terra::extract(distance(undist_forest_t1), test_point)[1, 2]

cat("T1 — undisturbed forest:", round(n_undist_t1/total_pixels*100, 1), "%\n")
cat("T1 — disturbed forest:  ", round(n_dist_t1/total_pixels*100, 1), "%\n")
cat("T1 — distance to forest:", round(dist_to_forest_t1, 0), "m\n")

# ── T2: sampling year cutoff ───────────────────────────────────────────────
sample_year_18_25 <- as.numeric(test_point$year) - 2000
cat("\nT2 cutoff (sampling year, 2-digit):", sample_year_18_25, "\n")

disturb_mask_t2  <- ifel(disturb_crop >= 18 & disturb_crop <= sample_year_18_25, 1, NA)
undist_forest_t2 <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask_t2),  1, NA)
dist_forest_t2   <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & !is.na(disturb_mask_t2), 1, NA)

n_undist_t2 <- global(undist_forest_t2, "sum", na.rm = TRUE)[[1]]
n_dist_t2   <- global(dist_forest_t2,   "sum", na.rm = TRUE)[[1]]

dist_to_forest_t2 <- terra::extract(distance(undist_forest_t2), test_point)[1, 2]

cat("T2 — undisturbed forest:", round(n_undist_t2/total_pixels*100, 1), "%\n")
cat("T2 — disturbed forest:  ", round(n_dist_t2/total_pixels*100, 1), "%\n")
cat("T2 — distance to forest:", round(dist_to_forest_t2, 0), "m\n")

# ── change summary ─────────────────────────────────────────────────────────
cat("\n=== CHANGE T1 → T2 ===\n")
cat("Years elapsed:", sample_year_18_25 - dist_year_18_25, "\n")
cat("Δ undisturbed forest (%):", round(n_undist_t2/total_pixels*100 - n_undist_t1/total_pixels*100, 1), "\n")
cat("Δ disturbed forest (%):  ", round(n_dist_t2/total_pixels*100 - n_dist_t1/total_pixels*100, 1), "\n")
cat("Δ distance to forest (m):", round(dist_to_forest_t2 - dist_to_forest_t1, 0), "\n")

# ── visual comparison ──────────────────────────────────────────────────────
par(mfrow = c(2, 5), mar = c(2, 2, 3, 1))


plot(disturb_crop, main = paste0("all disturbnaces"),
      legend = TRUE)


plot(dist_forest_t1, main = paste0("T1 distance to forest (", dist_year_18_25, ")"),
     legend = TRUE)
plot(distance(undist_forest_t1), main = paste0("T1 distance to forest (", dist_year_18_25, ")"),
     legend = TRUE)
plot(distance(undist_forest_t2), main = paste0("T1 distance to forest (", sample_year_18_25, ")"),
     legend = TRUE)



plot(dist_forest_t2, main = paste0("T2 distance to forest (", sample_year_18_25, ")"),
     legend = TRUE)

plot(disturb_mask_t1, main = paste0("T1 mask (cutoff ", dist_year_18_25, ")"),
     col = "orange", legend = FALSE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)

plot(undist_forest_t1, main = paste0("T1 undisturbed (", round(n_undist_t1/total_pixels*100,1), "%)"),
     col = "darkgreen", legend = FALSE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)

plot(dist_forest_t1, main = paste0("T1 disturbed (", round(n_dist_t1/total_pixels*100,1), "%)"),
     col = "firebrick", legend = FALSE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)

plot(disturb_mask_t2, main = paste0("T2 mask (cutoff ", sample_year_18_25, ")"),
     col = "orange", legend = FALSE)
plot(test_point, add = TRUE, col = "blue", pch = 19, cex = 1.5)

plot(undist_forest_t2, main = paste0("T2 undisturbed (", round(n_undist_t2/total_pixels*100,1), "%)"),
     col = "darkgreen", legend = FALSE)
plot(test_point, add = TRUE, col = "blue", pch = 19, cex = 1.5)

plot(dist_forest_t2, main = paste0("T2 disturbed (", round(n_dist_t2/total_pixels*100,1), "%)"),
     col = "firebrick", legend = FALSE)
plot(test_point, add = TRUE, col = "blue", pch = 19, cex = 1.5)

par(mfrow = c(1, 1))


# ==============================================================================
#  FULL RUN — all Czechia subplots
#  Computes landscape metrics at TWO time points:
#    T1 = at disturbance (dist_year_18_25)
#    T2 = at field sampling (year)
#  Plus delta (T2 - T1) to quantify how much the patch grew / forest shrank
# ==============================================================================
buffer_dist <- 500

# ── helper: compute one full snapshot given a cutoff year ─────────────────────
process_snapshot <- function(pt, d_crop, f_crop, cutoff_year, buffer_dist) {
  
  f_crop_proj    <- if (!same.crs(f_crop, d_crop)) {
    terra::project(f_crop, crs(d_crop), method = "near")
  } else f_crop
  forest_aligned <- terra::resample(f_crop_proj, d_crop, method = "near")
  
  disturb_mask  <- ifel(d_crop >= 18 & d_crop <= cutoff_year, 1, NA)
  undist_forest <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & is.na(disturb_mask),  1, NA)
  dist_forest   <- ifel(forest_aligned != 9 & !is.na(forest_aligned) & !is.na(disturb_mask), 1, NA)
  
  total_pixels <- ncell(forest_aligned)
  n_undist     <- global(undist_forest, "sum", na.rm = TRUE)[[1]]
  n_dist       <- global(dist_forest,   "sum", na.rm = TRUE)[[1]]
  pct_undist   <- round(n_undist / total_pixels * 100, 2)
  pct_dist     <- round(n_dist   / total_pixels * 100, 2)
  
  # expand buffer if no undisturbed forest found
  if (n_undist == 0) {
    buff2          <- buffer(pt, buffer_dist * 2)
    d_crop         <- crop(disturb_r, buff2)
    f_crop         <- crop(forest_r,  buff2)
    f_crop_proj    <- if (!same.crs(f_crop, d_crop)) {
      terra::project(f_crop, crs(d_crop), method = "near")
    } else f_crop
    forest_aligned <- terra::resample(f_crop_proj, d_crop, method = "near")
    disturb_mask   <- ifel(d_crop >= 18 & d_crop <= cutoff_year, 1, NA)
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
  patch_area_ha <- round(n_dist * 100 / 10000, 2)
  
  pct_dist_of_forest <- ifelse(
    (n_dist + n_undist) > 0,
    round(n_dist / (n_dist + n_undist) * 100, 1),
    NA_real_
  )
  
  list(
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
}

# ── NA template for error fallback (matches column structure below) ──────────
na_result <- function(pt) {
  data.frame(
    subplot                 = as.character(pt$subplot),
    plot_id                 = as.character(pt$plot_id),
    year                    = as.numeric(pt$year),
    disturbance_year        = as.numeric(pt$disturbance_year),
    dist_year_18_25         = NA_real_,
    sample_year_18_25       = NA_real_,
    dist_to_forest_m_t1     = NA_real_,
    pct_undist_forest_t1    = NA_real_,
    pct_dist_forest_t1      = NA_real_,
    pct_spruce_in_undist_t1 = NA_real_,
    patch_n_disturbed_t1    = NA_real_,
    patch_area_mn_ha_t1     = NA_real_,
    patch_area_ha_t1        = NA_real_,
    patch_para_t1           = NA_real_,
    pct_dist_of_forest_t1   = NA_real_,
    dist_to_forest_m_t2     = NA_real_,
    pct_undist_forest_t2    = NA_real_,
    pct_dist_forest_t2      = NA_real_,
    pct_spruce_in_undist_t2 = NA_real_,
    patch_n_disturbed_t2    = NA_real_,
    patch_area_mn_ha_t2     = NA_real_,
    patch_area_ha_t2        = NA_real_,
    patch_para_t2           = NA_real_,
    pct_dist_of_forest_t2   = NA_real_,
    delta_dist_to_forest    = NA_real_,
    delta_pct_dist_forest   = NA_real_,
    delta_patch_area_ha     = NA_real_,
    delta_patch_n           = NA_real_
  )
}

# ── main loop ───────────────────────────────────────────────────────────────
results_list <- lapply(seq_len(nrow(pts_cz)), function(i) {
  
  if (i %% 100 == 0) cat("  Point", i, "/", nrow(pts_cz), "\n")
  
  tryCatch({
    pt <- pts_cz[i, ]
    
    buff   <- buffer(pt, buffer_dist)
    d_crop <- crop(disturb_r, buff)
    f_crop <- crop(forest_r,  buff)
    
    # T1 cutoff — disturbance year (raster, fallback to table)
    dist_year_18_25 <- terra::extract(d_crop, pt)[1, 2]
    if (is.na(dist_year_18_25)) {
      dist_year_18_25 <- as.numeric(pt$disturbance_year) - 2000
    }
    
    # T2 cutoff — field sampling year
    sample_year_18_25 <- as.numeric(pt$year) - 2000
    
    # compute both snapshots
    t1 <- process_snapshot(pt, d_crop, f_crop, dist_year_18_25,   buffer_dist)
    t2 <- process_snapshot(pt, d_crop, f_crop, sample_year_18_25, buffer_dist)
    
    data.frame(
      subplot                 = as.character(pt$subplot),
      plot_id                 = as.character(pt$plot_id),
      year                    = as.numeric(pt$year),
      disturbance_year        = as.numeric(pt$disturbance_year),
      dist_year_18_25         = dist_year_18_25,
      sample_year_18_25       = sample_year_18_25,
      
      # T1 — at disturbance
      dist_to_forest_m_t1     = t1$dist_to_forest_m,
      pct_undist_forest_t1    = t1$pct_undist_forest,
      pct_dist_forest_t1      = t1$pct_dist_forest,
      pct_spruce_in_undist_t1 = t1$pct_spruce_in_undist,
      patch_n_disturbed_t1    = t1$patch_n_disturbed,
      patch_area_mn_ha_t1     = t1$patch_area_mn_ha,
      patch_area_ha_t1        = t1$patch_area_ha,
      patch_para_t1           = t1$patch_para,
      pct_dist_of_forest_t1   = t1$pct_dist_of_forest,
      
      # T2 — at field sampling
      dist_to_forest_m_t2     = t2$dist_to_forest_m,
      pct_undist_forest_t2    = t2$pct_undist_forest,
      pct_dist_forest_t2      = t2$pct_dist_forest,
      pct_spruce_in_undist_t2 = t2$pct_spruce_in_undist,
      patch_n_disturbed_t2    = t2$patch_n_disturbed,
      patch_area_mn_ha_t2     = t2$patch_area_mn_ha,
      patch_area_ha_t2        = t2$patch_area_ha,
      patch_para_t2           = t2$patch_para,
      pct_dist_of_forest_t2   = t2$pct_dist_of_forest,
      
      # change T1 → T2
      delta_dist_to_forest    = t2$dist_to_forest_m   - t1$dist_to_forest_m,
      delta_pct_dist_forest   = t2$pct_dist_forest     - t1$pct_dist_forest,
      delta_patch_area_ha     = t2$patch_area_ha       - t1$patch_area_ha,
      delta_patch_n           = t2$patch_n_disturbed   - t1$patch_n_disturbed
    )
    
  }, error = function(e) {
    warning("Point ", i, " failed: ", e$message)
    na_result(pts_cz[i, ])
  })
})

# ── combine and export ─────────────────────────────────────────────────────────
result_cz <- do.call(rbind, results_list)

cat("\nDone. Points:", nrow(result_cz), "\n")
cat("NAs dist_to_forest_t1:", sum(is.na(result_cz$dist_to_forest_m_t1)), "\n")
cat("NAs dist_to_forest_t2:", sum(is.na(result_cz$dist_to_forest_m_t2)), "\n")
summary(result_cz)

# ── quick check: did things actually change between T1 and T2? ───────────────
cat("\n=== Change summary (T1 → T2) ===\n")
summary(result_cz[, c("delta_dist_to_forest", "delta_pct_dist_forest",
                      "delta_patch_area_ha", "delta_patch_n")])

result_cz %>% 
  arrange(desc(delta_pct_dist_forest)) %>%
  select(subplot, plot_id, dist_year_18_25, sample_year_18_25,
         pct_dist_forest_t1, pct_dist_forest_t2, delta_pct_dist_forest) %>%
  head(20)

data.table::fwrite(result_cz,
                   "outData/dist_and_patch_metrics_czechia_temporal.csv")
cat("Exported to outData/dist_and_patch_metrics_czechia_temporal.csv\n")
