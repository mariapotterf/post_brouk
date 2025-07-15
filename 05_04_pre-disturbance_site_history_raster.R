
# ----------------------------------------------------------------------------
# pre-disturbnace site history: from raster data
# -----------------------------------------------------------------------------

# - get sites history - from tree density cover, dominant leaf type
# compare stem density, dominant tree species, species diversity
# vertical and horizontal structure across sscales: fiedl vs drone


gc()

library(terra)
library(sf)
#library(DBI)
library(ggplot2)
#library(dbscan)
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)


# Read files --------------------------------------

# 2023 
dat23_subplot <- fread('outData/subplot_full_2023.csv')
dat23_cluster <- fread('outData/df_cluster_2023.csv')

# Save spatial subplot data with cluster IDs
dat23_sf <- st_read("outData/sf_context_2023.gpkg")


# 2025
dat25_subplot <- fread('outData/subplot_full_2025.csv')
dat25_cluster <- fread('outData/df_cluster_2025.csv')

# Save spatial subplot data with cluster IDs
dat25_sf <- st_read("outData/subplot_with_clusters_2025.gpkg")


# Pre-disturbance history: automatic
#  - get forest pre-disturbance characteristics: tree density, dominant leaf type
tree_cover_dens15    <- rast("raw/forest_heights_composition/TCD_2015_020m_CR.tif") # tree cover density
dominant_leaf_type15 <- rast("raw/forest_heights_composition/DLT_2015_020m_CR.tif") # dominant leaf type


# recode 23 values if neede:

str(dat23_subplot)


# tree cover density: 0-100%
# leaf type: 
# - 1 - deciduous, 
# - 2 - coniferous 


# Get pre-disturbance forest characteristics in point sdurrounding ----------------

# Get target CRS from one of the rasters
target_crs <- crs(tree_cover_dens15)

# Reproject both sf objects to match raster CRS
dat23_sf <- st_transform(dat23_sf, target_crs)
dat25_sf <- st_transform(dat25_sf, target_crs)

# Step 4: Add year and combine
dat23_sf$year <- 2023
dat25_sf$year <- 2025

# Calculate mean coords per cluster and keep one point per cluster
# Add x/y coordinates as new columns
dat23_sf_mean <- dat23_sf %>%
  mutate(x = st_coordinates(.)[, 1],
         y = st_coordinates(.)[, 2])


dat25_sf_mean <- dat25_sf %>%
  mutate(x = st_coordinates(.)[, 1],
         y = st_coordinates(.)[, 2])



# Unify number of columns between two years
dat23_clusters <- dat23_sf_mean %>%
  group_by(cluster) %>%
  summarise(
    x = mean(x),
    y = mean(y),
    year = first(year)
  ) %>%
  st_as_sf(coords = c("x", "y"), crs = target_crs)

dat25_clusters <- dat25_sf_mean %>%
  group_by(cluster) %>%
  summarise(
    x = mean(x),
    y = mean(y),
    year = first(year)
  ) %>%
  st_as_sf(coords = c("x", "y"), crs = target_crs) %>% 
  mutate(cluster = as.character(cluster))



# Combine both years
combined_clusters <- bind_rows(dat23_clusters, dat25_clusters)

# Convert to terra vect
cluster_vect <- vect(combined_clusters)

# Step 3: Buffer XXm around points
my_buffer <- buffer(cluster_vect, width = 60)

# Step 1: Make lookup table
buffer_lookup <- as.data.frame(my_buffer) %>%
  mutate(ID = row_number()) %>%
  dplyr::select(ID, cluster, year)


# Step 4: Extract raster values for both rasters
vals_cover <- terra::extract(tree_cover_dens15, my_buffer, ID = TRUE)
vals_leaf  <- terra::extract(dominant_leaf_type15, my_buffer, ID = TRUE)

# Step 2: Join to extracted values
vals_cover <- left_join(vals_cover, buffer_lookup, by = "ID")
vals_leaf  <- left_join(vals_leaf,  buffer_lookup, by = "ID")


# Recode leaf type to factor for clarity
vals_leaf_recode <- vals_leaf %>%
  mutate(leaf_type = factor(DLT_2015_020m_CR,
                            levels = c(1, 2),
                            labels = c("deciduous", "coniferous")))

# Compute Shannon diversity index per buffer
leaf_diversity <- vals_leaf_recode %>%
  dplyr::filter(!is.na(DLT_2015_020m_CR)) %>%
  group_by(cluster, leaf_type) %>%
  tally() %>%
  group_by(cluster) %>%
  mutate(p = n / sum(n)) %>%
  summarise(shannon = -sum(p * log(p)), .groups = "drop")

hist(leaf_diversity$shannon)


# Calculate coniferous share per cluster buffer !!! need to finalize! I think i have only one value per 
# cluster! aslo, get some sort of aggregation/dispersion of coniferous forests
leaf_stats <- vals_leaf %>%
  group_by(cluster) %>%
  summarise(
    total_cells = n(),  # all cells, including NA
    n_coniferous = sum(DLT_2015_020m_CR == 2, na.rm = TRUE),
    n_deciduous  = sum(DLT_2015_020m_CR == 1, na.rm = TRUE),
    forest_cells = n_coniferous + n_deciduous,
    share_coniferous = ifelse(forest_cells > 0, n_coniferous / forest_cells, NA_real_),
    .groups = "drop"
  ) %>% 
  left_join(leaf_diversity)



# Get summarize statistics for continuous distribution: density for each point ID
summary_stats <- function(df, value_col) {
  df %>%
    group_by(cluster) %>%
    summarise(
      mean   = mean(.data[[value_col]], na.rm = TRUE),
      median = median(.data[[value_col]], na.rm = TRUE),
      cv     = sd(.data[[value_col]], na.rm = TRUE) / mean(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

cover_stats <- summary_stats(vals_cover, "TCD_2015_020m_CR")

# Rename columns before joining
cover_stats_renamed <- cover_stats %>%
  rename_with(~ paste0(.x, "_cover_dens"), .cols = c("mean", "median", "cv"))


# Full join on cluster and year
combined_stats <- full_join(cover_stats_renamed, leaf_stats, by = c("cluster"))

combined_stats %>% 
  ggplot(aes(x = shannon,
             y = share_coniferous)) + 
  geom_point() +
  geom_smooth()


# Reshape to long format
cv_long <- combined_stats %>%
  select(cluster, year, cv_cover, cv_leaf) %>%
  pivot_longer(cols = c(cv_cover, cv_leaf),
               names_to = "variable", values_to = "cv") %>%
  mutate(variable = recode(variable,
                           cv_cover = "Tree Cover Density",
                           cv_leaf  = "Dominant Leaf Type"))

# Plot
ggplot(cv_long, aes(x = variable, y = cv)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.shape = 21) +
  labs(x = "", y = "Coefficient of Variation (CV)",
       title = "Variation in 2015 Pre-Disturbance Forest Structure") +
  theme_minimal(base_size = 14)

hist(combined_stats$median_cover)
hist(combined_stats$median_leaf)


# save data ---------------------------------------------------------------------

fwrite(combined_stats, "outTable/pre_disturb_history_raster.csv")

