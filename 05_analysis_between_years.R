
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


# save files: -------------------
dat25_subplot <- fread('outData/subplot_full_2025.csv')
dat25_cluster <- fread('outData/df_cluster_2025.csv')

# Save spatial subplot data with cluster IDs
dat25_sf <- st_read("outData/subplot_with_clusters_2025.gpkg")
