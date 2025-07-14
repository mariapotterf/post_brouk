
# ------------------------------------------
# Compare 2025 with 2023 datasets
# -----------------------------------------

# read data: spatial and vegetation from 2025
# read data: spatial and vegetation from 2023

# identify the overlapping clusters -> make common_clust_ID
# - changes in stem density: oevrall, between species
# - get sites history - from tree density cover, dominant leaf type
# compare stem density, dominant tree species, species diversity
# vertical and horizontal structure among 2023-2025

# identify damaged terminal: from raw data czechia 20203 - export gpkg
# identify types of damages: 2025 - export gpkg
# Drones: 
# - get vertical strucure

gc()

library(terra)
library(sf)
library(DBI)
library(ggplot2)
library(dbscan)
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
dat23_sf <- st_read("raw/collected_2023/dat_czechia_2023.gpkg")


# 2025
dat25_subplot <- fread('outData/subplot_full_2025.csv')
dat25_cluster <- fread('outData/df_cluster_2025.csv')

# Save spatial subplot data with cluster IDs
dat25_sf <- st_read("outData/subplot_with_clusters_2025.gpkg")


# get forest pre-disturbance characteristics: tree density, dominant leaft type
tree_cover_dens15    <- rast("raw/forest_heights_composition/TCD_2015_020m_CR.tif") # tree cover density
dominant_leaf_type15 <- rast("raw/forest_heights_composition/DLT_2015_020m_CR.tif") # dominant leaf type

# tree cover density: 0-100%
# leaf type: 1 - deciduous, 2 - coniferous 

# Find overlapping clusters -----------------------------------
st_geometry(dat23_sf) <- "geom"
st_geometry(dat25_sf) <- "geom"


# Summarize to cluster centroids -. get proximity only for centroids
centroids_2023 <- dat23_sf %>%
  group_by(cluster) %>%
  summarise(do_union = TRUE) %>%
  st_centroid() %>%
  rename(cluster_2023 = cluster)

centroids_2025 <- dat25_sf %>%
  group_by(cluster) %>%
  summarise(do_union = TRUE) %>%
  st_centroid() %>%
  rename(cluster_2025 = cluster)


# Match CRS
centroids_2023 <- st_transform(centroids_2023, st_crs(centroids_2025))

# Find clusters within 50 m
overlaps <- st_join(centroids_2025, centroids_2023, join = st_is_within_distance, dist = 25) %>%
  dplyr::filter(!is.na(cluster_2023)) %>%
  mutate(common_cluster_ID = paste0("c_", cluster_2023, "_", cluster_2025))

overlaps

# export overlaps -----------------------------------

# Save spatial subplot data with cluster IDs
st_write(overlaps, "outData/overlaps.gpkg", delete_layer = TRUE)


# Get pre-disturbance forest characteristics ----------------
# Step 2: Get target CRS from one of the rasters
target_crs <- crs(tree_cover_dens15)

# Step 3: Reproject both sf objects to match raster CRS
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
buffer_lookup <- as.data.frame(buffer_60m) %>%
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
  rename_with(~ paste0(.x, "_cover"), .cols = c("mean", "median", "cv"))


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


# get stem density and shannon for field data -----------------------------------


# make working example: 

dd <- data.frame(
  cluster = c(1, 1,   
              2, 2, 2,  
              3,   
              4, 
              5),
  species = c("a", "b",
              "a", "b", "d",
              NA, 
              "d",
              "d"),
  basal_area = c(5, 5,  
                 1, 1, 8,  
                 0,  
                 10, 
                 0.5)
)

# Define Shannon + Effective Species Number
shannon_stats <- function(ba) {
  total <- sum(ba)
  if (total == 0) return(c(H = 0, ESN = 0))
  p <- ba / total
  H <- -sum(p * log(p), na.rm = TRUE) # shannon index
  ESN <- exp(H)  # effective number of species
  c(H = H, ESN = ESN)
}

# All clusters
all_clusters <- dd %>% distinct(cluster)

# Compute diversity only for rows with species info
shannon_summary <- dd %>%
  filter(!is.na(species)) %>%
  group_by(cluster) %>%
  reframe(
    shannon = shannon_stats(basal_area)[["H"]],
    effective_species = shannon_stats(basal_area)[["ESN"]]
  )

# Merge back to include missing clusters
shannon_per_cluster <- all_clusters %>%
  left_join(shannon_summary, by = "cluster") %>%
  mutate(
    shannon = replace_na(shannon, 0),
    effective_species = replace_na(effective_species, 0)
  )

print(shannon_per_cluster)

# summarize across species: get shannon per subplot, per plot
dat23_subplot_sum_cl <- dat23_subplot %>% 
  group_by( cluster, species ) %>% 
  summarise(stem_density = sum(stem_density, na.rm = T))
  

# calculate shannon index: --------------------------
# 1. Total stem density per cluster
total_stems <- dat23_subplot_sum_cl %>%
  group_by(cluster) %>%
  summarise(total_density = sum(stem_density), .groups = "drop")

# 2. Shannon index only for clusters with regeneration
shannon_nonzero <- dat23_subplot_sum_cl %>%
  left_join(total_stems, by = "cluster") %>%
  dplyr::filter(total_density > 0) %>%
  mutate(p = stem_density / total_density) %>%
  group_by(cluster) %>%
  summarise(shannon_index = -sum(p * log(p)), .groups = "drop")

# 3. Join with all clusters, assigning 0 to those with no regeneration
shannon_all <- total_stems %>%
  left_join(shannon_nonzero, by = "cluster") %>%
  mutate(shannon_index = if_else(is.na(shannon_index), 0, shannon_index))

# View result
head(shannon_all)
hist(shannon_all$shannon_index)
