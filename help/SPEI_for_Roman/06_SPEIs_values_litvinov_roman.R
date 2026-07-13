

# Calculate SPEI characteristics for individual trees in 
# Litvinov (Roman, 07/13/2026)
# from SPEI12_harg1_split, spei12 from Ledvinka, generated in 2025


## ============================================================
## Extract yearly SPEI12 (mean + SD) per tree from monthly rasters
## Strategy: polygons -> centroids -> point extraction
## Rationale: raster res = 500 m, trees = 5-10 m -> polygon zonal
##            stats and point extraction are numerically equivalent
##            here, so centroids avoid unnecessary overhead.
## ============================================================

library(terra)
library(dplyr)   # only used for the final summarise; drop if you'd rather stay base-R

## ---- 1. paths ------------------------------------------------
## relative to getwd() == ".../Dokumenty/2025_CZU_postbrouk/r_post_brouk"
poly_path <- "../../2026_roman/2026_litvinov_stromy/raw/LitvinovSelectedTrees/SelectedTreesLitvPoly.shp"

# raw_SPEI12 folder - just this one
raster_dir <- "raw/SPEI12_new/SPEI12_harg1_split"

raster_patt <- "^spei12_[0-9]{4}_[0-9]{2}_01\\.tif$"   # adjust to your actual naming

id_field <- "Nazev"   # <- tree ID column in the attribute table

## ---- 2. read trees, make centroids ----------------------------
trees_poly <- vect(poly_path)
trees_pts  <- centroids(trees_poly, inside = TRUE)  # inside=TRUE keeps centroid within odd-shaped polygons

## ---- 3. build the raster stack --------------------------------
tif_files <- list.files(raster_dir, pattern = raster_patt, full.names = TRUE)
stopifnot(length(tif_files) > 0)

spei_stack <- rast(tif_files)

## Parse date from filename -> layer names, e.g. spei12_2007_09_01.tif -> 2007-09
dates <- sub("^spei12_([0-9]{4})_([0-9]{2})_01\\.tif$", "\\1-\\2", basename(tif_files))
names(spei_stack) <- dates

## sanity check: one layer per file, chronological
spei_stack <- spei_stack[[order(names(spei_stack))]]

## ---- 4. CRS alignment ------------------------------------------
# IMPORTANT: these SPEI12 rasters come from CHMI (Czech Hydrometeorological
# Institute), not a standard downloaded product, and their embedded CRS is
# broken: datum is tagged "unknown" and there's no EPSG authority code, so
# GDAL reports it as an "unnamed" CRS (sometimes shown as IGNF:ETRS89LAEA).
# The actual projection parameters, however, are unambiguously EPSG:3035
# (ETRS89-LAEA Europe: lat_0=52, lon_0=10, x_0=4321000, y_0=3210000, in
# metres). Relying on terra to auto-match an "unknown" datum CRS via
# same.crs()/project() is risky -- it can warn, silently skip a required
# shift, or mismatch. Safer to force the known-correct CRS explicitly:
# name the projection
crs(spei_stack) <- "EPSG:3035"

# Reproject points  to match the raster CRS 
if (!same.crs(trees_pts, spei_stack)) {
  trees_pts <- project(trees_pts, crs(spei_stack))
}

# sanity check: points should now fall within the raster extent
stopifnot(all(relate(trees_pts, as.polygons(ext(spei_stack), crs = crs(spei_stack)), "intersects")))

## ---- 5. extract --------------------------------------------------
# method = "bilinear" gives smoother, tree-differentiated values;
# use "simple" (nearest cell) if you want the literal cell value instead.
# bilinear represent the average values of teh four neighbouring cells - hence, can make some variation if trees are loated within the opoosite parts of teh pixel
vals <- extract(spei_stack, trees_pts, method = "bilinear", ID = FALSE)

# attach tree ID
vals[[id_field]] <- values(trees_pts)[[id_field]]

## ---- 6. reshape to long format -----------------------------------
long_df <- vals %>%
  tidyr::pivot_longer(
    cols = -all_of(id_field),
    names_to = "yearmonth",
    values_to = "spei12"
  ) %>%
  mutate(
    year  = as.integer(substr(yearmonth, 1, 4)),
    month = as.integer(substr(yearmonth, 6, 7))
  )

## ---- 7. yearly mean + sd per tree ---------------------------------
yearly_stats <- long_df %>%
  group_by(.data[[id_field]], year) %>%
  summarise(
    spei12_mean = mean(spei12, na.rm = TRUE),
    spei12_sd   = sd(spei12,   na.rm = TRUE),
    n_months    = sum(!is.na(spei12)),
    .groups = "drop"
  ) %>%
  arrange(.data[[id_field]], year) %>% 
  select(-n_months)

#View(yearly_stats)

yearly_stats_90 <- yearly_stats %>% 
  filter(year >= 1990)

## ---- 8. save -------------------------------------------------------
#write.csv(yearly_stats, "tree_spei12_yearly_stats.csv", row.names = FALSE)
out_path <- "../../2026_roman/2026_litvinov_stromy/outTable/litvinov_tree_spei12_1990_2024.csv"

write.csv(yearly_stats_90, out_path, row.names = FALSE)

## quick look
head(yearly_stats)


library(ggplot2)

ggplot(yearly_stats_90, aes(x = year, 
                            y = spei12_mean, 
                            color = factor(Nazev))) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  labs(
    x = "Year",
    y = "Mean SPEI12",
    color = "Tree (Nazev)",
    title = "Yearly mean SPEI12 per tree"
  ) +
  theme_minimal()
