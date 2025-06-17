

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

# look up tables
species <- fread(paste0(raw_path, "/look_up_tables/full_sp_list.csv"))

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


# get look up tables ------------------------------



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


# Interpret table:  ---------------------------------------------
# get species, counts, vertical classes per plot and cluster
# investigate individual outputs
# overall, get the list of species and aboundaces per plot
tables$damage_list
tables$mature_test
tables$regeneration_adv2$height

# investigate if my species list are the same? 
identical(tables$species_list2$species, tables$species_list$species)
# YES!!!

# correctly interpret: 
#  - species type
#  - counts: properly listed: 1=1, 2=2, ..17 is >17
# 



# Create whole database: -------------------------------------------------------------

# Intrepret values based on look up tables: species and accronymes



### merge to species table to get the Value -----------------
tables$species_list <- merge(
  tables$species_list,
  species_codes,
  by = "species",
  all.x = TRUE,
  sort = FALSE
)

# ideal strcuture

cluster_id  subplot VegTYpe (mature, advanced, small) acc (species_acronym) 