# ------------------------------------------------------------------------------
#   Get terrain characteristics
# ------------------------------------------------------------------------------

# run on subplot data, then calculate averages fro plot data
# merge back to vegetation chars and create a table having both plot and subplot data

# ── Paths ──────────────────────────────────────────────────────
proj1_root <- "C:/Users/potterf/OneDrive - CZU v Praze/Dokumenty/2023_PanEuropean/r_paneurop/"
elev_path  <- file.path(proj1_root, "rawData/dem")
dist_path  <- file.path(proj1_root, "rawData/disturb_data")

gpkg_sub  <- "C:/Users/potterf/OneDrive - CZU v Praze/Dokumenty/2025_CZU_postbrouk/r_post_brouk/outDataShare/Karim_AEF/regeneration_chars_subplot_3035.gpkg"
gpkg_plot  <- "C:/Users/potterf/OneDrive - CZU v Praze/Dokumenty/2025_CZU_postbrouk/r_post_brouk/outDataShare/Karim_AEF/regeneration_chars_plot_3035.gpkg"

# ── Libraries ──────────────────────────────────────────────────
library(terra)
library(dplyr)
library(tidyr)
library(data.table)

# ── Country map: gpkg value → file naming convention ───────────
country_dem_map <- c(
  "Czech Republic" = "Czechia",
  "Germany"        = "Germany",
  "Austria"        = "Austria",
  "Slovakia"       = "Slovakia",
  "Poland"         = "Poland",
  "Italy"          = "Italy",
  "Switzerland"    = "Switzerland",
  "Slovenia"       = "Slovenia",
  "France"         = "France"
)

# disturbance files use lowercase country names
country_dist_map <- c(
  "Czech Republic" = "czechia",
  "Germany"        = "germany",
  "Austria"        = "austria",
  "Slovakia"       = "slovakia",
  "Poland"         = "poland",
  "Italy"          = "italy",
  "Switzerland"    = "switzerland",
  "Slovenia"       = "slovenia",
  "France"         = "france"
)

# ── Load all subplots once ─────────────────────────────────────
all_subplots <- vect(gpkg_sub)  # already EPSG:3035
all_plots   <- vect(gpkg_plot)     # already EPSG:3035


cols_remove <- c("pre_dist_trees_n", "area_m2", "pre_dist_dens_ha",
                 "time_snc_full_disturbance", "disturbance_year",
                 "forest_year", "disturbance_length")

all_subplots <- all_subplots[, !names(all_subplots) %in% cols_remove]

# verify
names(all_subplots)

# ── Helper: load agent raster or return NA dummy ───────────────
read_or_dummy <- function(path, reference) {
  if (file.exists(path)) {
    terra::rast(path)
  } else {
    warning(paste("File not found, using NA dummy:", path))
    dummy <- terra::rast(reference, nlyr = 1)
    values(dummy) <- NA
    dummy
  }
}

# ── Helper: distance to patch edge for one point ───────────────
process_point_edge <- function(point, disturbance, buffer_dist = 1500) {
  
  crs(point) <- crs(disturbance)
  
  # buffer, crop, mask
  buff             <- buffer(point, buffer_dist)
  disturbance_crop <- crop(disturbance, buff)
  disturbance_mask <- mask(disturbance_crop, buff)
  
  # outside buffer → 0 placeholder
  disturbance_mask[is.na(disturbance_mask)] <- 0
  
  # reclassify:
  # recent patch (2018-2020) → NA  (target: distance measured TO these cells)
  # old disturbance + outside → 1  (origin: non-NA cells get distance values)
  reclass_matrix <- matrix(c(0,      0,      1,     # outside buffer → non-patch
                             1986,   2017.5, 1,     # old disturbance → non-patch
                             2017.6, 2020,   NA),   # recent patch → NA (target)
                           ncol = 3, byrow = TRUE)
  
  reclassified <- classify(disturbance_mask, reclass_matrix)
  
  # distance() measures from each non-NA cell to nearest NA (= patch edge)
  dist_rast <- distance(reclassified)
  
  # point sits inside patch (NA cell) → extract may return NA
  point_dist <- terra::extract(dist_rast, point, method = "simple")[1, 2]
  
  # if NA (point inside patch), get value from nearest non-NA cell
  if (is.na(point_dist)) {
    nearest    <- terra::nearby(point, dist_rast, n = 1)
    point_dist <- terra::extract(dist_rast, nearest)[1, 2]
  }
  
  return(round(point_dist, 0))
}

# ── Main function: terrain + disturbance + edge distance ───────
extract_env_info <- function(country_key, points_all, buffer_dist = 500) {
  print(paste("Processing", country_key))
  
  # 1. subset points
  pts <- points_all[points_all$country_name == country_key, ]
  if (nrow(pts) == 0) {
    warning(paste("No points for", country_key, "— skipping"))
    return(NULL)
  }
  
  dist_key <- country_dist_map[country_key]
  dem_key  <- country_dem_map[country_key]
  
  # 2. disturbance year (reference raster)
  disturb_name <- paste0("disturbance_year_1986-2020_", dist_key, ".tif")
  disturbance  <- rast(file.path(dist_path, dist_key, disturb_name))
  
  # 3. severity
  severity_name <- paste0("disturbance_severity_1986-2020_", dist_key, ".tif")
  severity      <- rast(file.path(dist_path, dist_key, severity_name))
  severity_proj <- terra::project(severity, disturbance, method = "near")
  
  # 4. elevation + terrain
  elevation      <- rast(file.path(elev_path, paste0("dem_", dem_key, ".tif")))
  elevation_proj <- terra::project(elevation,       disturbance, method = "near")
  elevation_proj <- terra::resample(elevation_proj, disturbance, method = "near")
  slope          <- terra::terrain(elevation_proj, "slope",  neighbors = 8)
  aspect         <- terra::terrain(elevation_proj, "aspect", neighbors = 8)
  
  # 5. reproject points to disturbance CRS
  pts_proj <- terra::project(pts, disturbance)
  
  # 6. stack and extract
  env_stack <- c(disturbance, severity_proj, elevation_proj, slope, aspect)
  names(env_stack) <- c("disturbance_year", "disturbance_severity",
                        "elevation", "slope", "aspect")
  
  result <- as.data.frame(
    terra::extract(env_stack, pts_proj, method = "simple", bind = TRUE)
  )
  
  # 7. distance to patch edge — point by point
  cat("  Computing distance to edge for", nrow(pts_proj), "points...\n")
  dist_to_edge <- sapply(1:nrow(pts_proj), function(i) {
    tryCatch(
      process_point_edge(pts_proj[i, ], disturbance, buffer_dist),
      error = function(e) {
        warning(paste("Edge distance failed for point", i, ":", e$message))
        NA_real_
      }
    )
  })
  
  result$dist_to_edge <- dist_to_edge
  
  return(result)
}

# ── Test ───────────────────────────────────────────────────────
test <- extract_env_info("Slovenia", all_subplots)
head(test[, c("subplot", "disturbance_year", "disturbance_severity",
              "elevation", "slope", "aspect", "dist_to_edge")])
summary(test$dist_to_edge)

# ── Run all countries ──────────────────────────────────────────
out_ls     <- lapply(names(country_dem_map), extract_env_info,
                     points_all = all_subplots)
final_sub_env  <- do.call(rbind, out_ls)

# ── Fill NAs ───────────────────────────────────────────────────
all_sub_env <- final_sub_env %>%
  mutate(ID_short = gsub('.{2}$', '', subplot)) %>%
  group_by(ID_short) %>%
  arrange(ID_short, disturbance_year) %>%
  fill(disturbance_year,     .direction = "down") %>%
  fill(disturbance_severity, .direction = "down") %>%
  ungroup() %>%
  select(-ID_short) %>%
  mutate(
    slope        = ifelse(is.na(slope),        median(slope,        na.rm = TRUE), slope),
    aspect       = ifelse(is.na(aspect),       median(aspect,       na.rm = TRUE), aspect),
    elevation    = ifelse(is.na(elevation),    median(elevation,    na.rm = TRUE), elevation),
    dist_to_edge = ifelse(is.na(dist_to_edge), median(dist_to_edge, na.rm = TRUE), dist_to_edge),
    disturbance_year = ifelse(is.na(disturbance_year), median(disturbance_year, na.rm = TRUE), disturbance_year),
    disturbance_severity = ifelse(is.na(disturbance_severity), median(disturbance_severity, na.rm = TRUE), disturbance_severity)
  )

# Get value as average per plot for disturbnace chars
plot_means <- all_sub_env %>%
  group_by(plot_id, year) %>%
  summarise(
    x                    = mean(x,                    na.rm = TRUE),
    y                    = mean(y,                    na.rm = TRUE),
    disturbance_year     = mean(disturbance_year,     na.rm = TRUE),
    disturbance_severity = mean(disturbance_severity, na.rm = TRUE),
    elevation            = mean(elevation,            na.rm = TRUE),
    slope                = mean(slope,                na.rm = TRUE),
    aspect               = mean(aspect,               na.rm = TRUE),
    dist_to_edge         = mean(dist_to_edge,         na.rm = TRUE),
    .groups = "drop"
  )

# add averaged distrubance data back to vegetation chars 
all_plots_env <- all_plots %>%
  as.data.frame() %>% 
  select(-pre_dist_trees_n, -area_m2, -pre_dist_dens_ha, 
         -time_snc_full_disturbance, -disturbance_year, 
         -forest_year, -disturbance_length,
         -x, -y) %>%
  left_join(plot_means, by = c("plot_id", "year"))

names(all_plots_env)
names(all_sub_env)

str(all_plots_env)
str(all_sub_env)


# merge into the same table
# ── 1. add level indicator & subplot ID to plot-level ──────────
all_sub_env_out <- all_sub_env %>%
  mutate(level = "subplot") 

all_plots_env_out <- all_plots_env %>%
  mutate(level = "plot",
         subplot = NA_character_) %>% # subplots don't exist at plot level
  select(all_of(names(all_sub_env_out))) # match column order from subplot table


# get both levels
df_both <- bind_rows(all_sub_env_out, all_plots_env_out) %>%
  arrange(plot_id, year, level, subplot) %>% 
  mutate(time_since = year - disturbance_year)


# ── 3. export ──────────────────────────────────────────────────
fwrite(df_both, "EU_scripts_help/both_levels_EU_comb.csv")


# ── Export -------------------------------------------
fwrite(final_sub_env, "EU_scripts_help/env_chars_subplots.csv")
fwrite(all_plots_env,  "EU_scripts_help/env_chars_plots.csv")


# check distance calculation -------------------------------------------
# ── Visual check for one point before running all ──────────────
# ── Visual check ───────────────────────────────────────────────
# load Slovenia disturbance as test case
disturb_slov  <- rast(file.path(dist_path, "slovenia",
                                "disturbance_year_1986-2020_slovenia.tif"))

# pick first Slovenia point
pts_slov      <- all_subplots[all_subplots$country_name == "Slovenia", ]
pts_slov_proj <- terra::project(pts_slov, disturb_slov)
test_point    <- pts_slov_proj[50, ]
crs(test_point) <- crs(disturb_slov)

# buffer, crop, mask
buff             <- buffer(test_point, 90)
disturbance_crop <- crop(disturb_slov, buff)
disturbance_mask <- mask(disturbance_crop, buff)

cat("Unique values in buffer:\n")
print(sort(unique(values(disturbance_mask), na.rm = TRUE)))

# reclassify: patch (2018-2020) = NA, everything else = 1
disturbance_mask[is.na(disturbance_mask)] <- 0

reclass_matrix <- matrix(c(0,    0,      1,     # outside buffer → non-patch
                           1986, 2017.5, 1,     # old disturbance → non-patch  
                           2017.6, 2020,   NA),   # recent patch → NA (target)
                         ncol = 3, byrow = TRUE)

reclassified <- classify(disturbance_mask, reclass_matrix)

# distance() measures from each non-NA cell to nearest NA (= patch edge)
dist_rast  <- distance(reclassified)

# extract at point — point sits in the patch (NA), so we need
# to extract from the surrounding distance surface
# → use the nearest non-NA cell value instead
point_dist <- terra::extract(dist_rast, test_point, method = "simple")

cat("Distance to patch edge:", round(point_dist[1, 2], 0), "m\n")

# visual
par(mfrow = c(1, 3))

plot(disturbance_mask, main = "1. Raw disturbance years", legend = TRUE)
plot(test_point, add = TRUE, col = "red", pch = 19, cex = 1.5)
legend("bottomleft", legend = "subplot", pch = 19, col = "red", bty = "n", cex = 0.8)

plot(reclassified, main = "2. Reclassified\n(NA = patch, 1 = non-patch)", legend = TRUE)
plot(test_point, add = TRUE, col = "blue", pch = 19, cex = 1.5)
legend("bottomleft", legend = "subplot", pch = 19, col = "blue", bty = "n", cex = 0.8)

plot(dist_rast, main = paste0("3. Distance to patch edge: ",
                              round(point_dist[1, 2], 0), " m"), legend = TRUE)
plot(test_point, add = TRUE, col = "blue", pch = 19, cex = 1.5)
legend("bottomleft", legend = "subplot", pch = 19, col = "blue", bty = "n", cex = 0.8)

par(mfrow = c(1, 1))




# Get climate data --------------------------------------------------------------

dat <- fread("EU_scripts_help/both_levels_EU_comb.csv")

# extract unique plot coordinates
plots <- dat %>%
  group_by(plot_id) %>%
  summarise(x = mean(x),
            y = mean(y),
            disturbance_year = first(disturbance_year))

plots_sf <- st_as_sf(
  plots,
  coords = c("x","y"),
  crs = 3035
)


plots_wgs <- st_transform(plots_sf, 4326)

library(terra)
library(geodata)

tavg <- geodata::worldclim_global(
  var = "tavg",
  res = 2.5,        # 2.5 arcmin ≈ 5 km
  path = "climate"
)

prec <- geodata::worldclim_global(
  var = "prec",
  res = 2.5,
  path = "climate"
)

tavg_vals <- terra::extract(tavg, vect(plots_wgs))
prec_vals <- terra::extract(prec, vect(plots_wgs))
