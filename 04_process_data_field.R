

# decode field data
# what species?
# what stem density?
# heights?
# damages?
# make a cluster based on 

# read gpkg - has location as point, and additional vegetation classes,
# damages, context
# identify unique name by cluster
# make clusters automatically


library(terra)
library(DBI)
library(ggplot2)
library(dbscan)


# paths -----------------
raw_path   <- "raw/outField"
tablet     <- "T2_AH_20250528"
file_name  <- paste0("forest_structure_survey_v2[1]_", tablet, ".gpkg")

# Load data: spatial geometry
subplot <- st_read(gpkg_path, layer = "subplot")

# Read non spatial data
con <- dbConnect(RSQLite::SQLite(), gpkg_path)

# liat all tables
dbListTables(con)

non_spatial_tables <- c("context", "damage_list", "mature_test", 
                        "regeneration_adv2", "regeneration_small", 
                        "species_list2", "species_list")

# Read into named list
tables <- lapply(non_spatial_tables, function(tbl) dbReadTable(con, tbl))
names(tables) <- non_spatial_tables

dbDisconnect(con)


# test the plotting
ggplot(subplot) +
  geom_sf() +
  ggtitle("Spatial distribution of subplots") +
  theme_minimal()


# create unique cluster ID based on distance -----------------

# Ensure geometry is in meters (projected CRS like EPSG:3035)
subplot_proj <- st_transform(subplot, 3035)

# Extract coordinates
coords <- st_coordinates(subplot_proj)

# Run DBSCAN: eps = max intra-cluster distance (e.g., 20 m)
db <- dbscan::dbscan(coords, eps = 30, minPts = 2)

# Add cluster ID to dataframe
subplot_proj$cluster_id <- db$cluster  # 0 means noise

# Convert cluster 0 (noise) to NA or a separate group if needed
#do it using subplot_proj$cluster_id <- ifelse(subplot_proj$cluster_id == 0, NA, subplot_proj$cluster_id)
# Export to GeoPackage
st_write(subplot_proj, "outData/subplot_with_clusters.gpkg", delete_layer= TRUE)
