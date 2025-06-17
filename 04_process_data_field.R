

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

# Load data field data: spatial geometry
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



# Create whole database: -------------------------------------------------------------

# correctly interpret: 
#  - species type
#  - counts: properly listed: 1=1, 2=2, ..17 is >17
# 


# Intrepret values based on look up tables: species and accronymes



### merge to species table to get the Value -----------------
tables$species_list <- merge(
  tables$species_list,
  species,
  by = c("Value", "species"),
  all.x = TRUE,
  sort = FALSE
) 


# Helper: get species acronym
get_species_code <- function(df, species_list) {
  df <- merge(df, species_list[, c("Value", "acc")], by.x = "species_id", by.y = "Value", all.x = TRUE)
  return(df)
}
# ideal strcuture

# cluster_id  subplot VegTYpe (mature, advanced, small) acc (species_acronym) height dbh dmg_terminal dmg_foliage ..


# Add type column
tables$mature_test$VegType <- "mature"
tables$regeneration_adv2$VegType <- "advanced"
tables$regeneration_small$VegType <- "small"


# Combine and annotate species
mature <- get_species_code(tables$mature_test, tables$species_list)
adv    <- get_species_code(tables$regeneration_adv2, tables$species_list)
small  <- get_species_code(tables$regeneration_small, tables$species_list)

# Step 1: Join cluster_id
subplot_cols <- subplot_proj[, c("plot_id", "cluster_id")]

cols_needed <- c(
  "height", "count", "dbh",
  "dmg_terminal", "dmg_terminal_photo",      
  "dmg_foliage", "dmg_stem_horizont", "dmg_stem_horizont_photo",  
  "dmg_root_stem", "dmg_root_stem_horiz_photo", "dmg_stem",                
  "dmg_stem_cause", "dmg_root_stem_horizon", "dmg_root_stem_vert",       
  "dmg_root_stem_cause", "dmg_foliage_leaf_int", "dmg_stem_vert_photo",   
  "plot_id", "dmg_stem_vertical", "dmg_root_stem_vert_photo",
  "dmg_stem_grass", "dmg_term_sample", "dmg_stem_sample",       
  "dmg_foliage_sample", "dmg_term_similar", "dmg_foliage_similar",    
  "dmg_foliage_detail_photo", "dmg_foliage_overall_photo", "dmg_root_stem_sample",    
  "dmg_type", "dmg_bool", "dmg_jedinec"
)


# Step 2: Harmonize each table
# Harmonization function
standardize_columns <- function(df, vegtype) {
  df$VegType <- vegtype
  
  # Ensure all needed columns exist
  for (col in cols_needed) {
    if (!col %in% colnames(df)) df[[col]] <- NA
  }
  
  # Add species and VegType
  df <- df[, c("plot_id", "species_id", "acc", "VegType", cols_needed)]
  
  # Join with subplot info to get cluster_id
  df <- merge(df, subplot_cols, by = "plot_id", all.x = TRUE)
  
  # Track damage occurrence class
  df$dmg_type_terminal <- ifelse(!is.na(df$dmg_terminal), vegtype, NA)
  df$dmg_type_foliage  <- ifelse(!is.na(df$dmg_foliage),  vegtype, NA)
  
  return(df)
}

# Step 3: Process each vegetation layer
mature_h <- standardize_columns(mature, "mature")
adv_h    <- standardize_columns(adv, "advanced")
small_h  <- standardize_columns(small, "small")

# Step 4: Combine all
combined <- dplyr::bind_rows(mature_h, adv_h, small_h)
