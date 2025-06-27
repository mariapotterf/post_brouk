
# ------------------------------------------
# Compare 2025 with 2023 datasets
# -----------------------------------------

# read data: spatial and vegetation from 2025
# read data: spatial and vegetation from 2023

# identify the overlapping clusters -> make common_clust_ID
# changes in stem density: oevrall, between species

# identify demaged terminal: from raw data czechia 20203

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



# Find overlapping clusters -----------------------------------
st_geometry(dat23_sf) <- "geom"
st_geometry(dat25_sf) <- "geom"


# Summarize to cluster centroids
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