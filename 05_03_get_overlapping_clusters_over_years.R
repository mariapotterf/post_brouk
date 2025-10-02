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
dat23_sf <- st_read("outData/sf_context_2023.gpkg")


# 2025
dat25_subplot <- fread('outData/subplot_full_2025.csv')
dat25_cluster <- fread('outData/df_cluster_2025.csv')

# Save spatial subplot data with cluster IDs
dat25_sf <- st_read("outData/subplot_with_clusters_2025.gpkg")


# recode 23 values if neede:

str(dat23_subplot)


# tree cover density: 0-100%
# leaf type: 
# - 1 - deciduous, 
# - 2 - coniferous 

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
st_write(overlaps, "outData/google_my_map/overlaps.kml", driver = "KML", delete_dsn = TRUE)

# Save spatial subplot data with cluster IDs
st_write(centroids_2025, "outData/centroids_2025.gpkg", delete_layer = TRUE)
st_write(centroids_2023, "outData/centroids_2023.gpkg", delete_layer = TRUE)
