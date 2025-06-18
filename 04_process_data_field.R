

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

gc()

library(terra)
library(sf)
library(DBI)
library(ggplot2)
library(dbscan)
library(data.table)


# paths -----------------
raw_path   <- "raw/outField"
tablet     <- "T2_AH_20250528"
file_name  <- paste0("forest_structure_survey_v2[1]_", tablet, ".gpkg")
gpkg_path  <- paste(raw_path, tablet, file_name, sep = "/") 

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




# LOOP ----TESt START ----------------------------------------------------

library(sf)
library(DBI)
library(data.table)
library(dplyr)
library(dbscan)

# ---- Setup ----
species <- fread("raw/outField/look_up_tables/full_sp_list.csv")

# List relevant .gpkg files
gpkg_files <- list.files(
  "raw/outField",
  pattern = "forest_structure_survey_v2.*\\.gpkg$",
  full.names = TRUE,
  recursive = TRUE
)

gpkg_files

# ---- PHASE 1: Collect all subplot geometries ----
subplot_all <- lapply(gpkg_files, function(path) {
  print(path)
  print(basename(dirname(path)))
  tryCatch({
    df <- st_read(path, layer = "subplot", quiet = TRUE)
    df$source_file <- basename(path)
    df$source_folder <- basename(dirname(path))
    df
  }, error = function(e) NULL)
})

subplot_all <- do.call(rbind, subplot_all)
subplot_all <- st_transform(subplot_all, 3035)
subplot_all$plot_key <- paste(subplot_all$plot_id, subplot_all$source_folder, sep = "_")

# ---- Cluster globally ----
coords <- st_coordinates(subplot_all)
db <- dbscan::dbscan(coords, eps = 30, minPts = 2)
subplot_all$cluster_id <- ifelse(db$cluster == 0, NA, db$cluster)

# For merging later
cluster_lookup <- subplot_all |> 
  st_drop_geometry() |> 
  dplyr::select(plot_id, source_folder, cluster_id) |>
  mutate(plot_key = paste(plot_id, source_folder, sep = "__"))

# Save spatial subplot data with cluster IDs
st_write(subplot_all, "outData/subplot_with_clusters_global.gpkg", delete_layer = TRUE)

# ---- PHASE 2: Process vegetation data ----

# Harmonization function
standardize_columns <- function(df, vegtype, source_folder) {
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
  
  df$VegType <- vegtype
  for (col in cols_needed) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df <- df[, c("plot_id", "species_id", "acc", "VegType", cols_needed)]
  
  # Add plot_key using known folder
  df$source_folder <- source_folder
  df$plot_key <- paste(df$plot_id, df$source_folder, sep = "_")
  
  # DEBUG: Check if plot_key exists in both sides
  if (!"plot_key" %in% names(cluster_lookup)) stop("cluster_lookup missing plot_key column")
  if (!"plot_key" %in% names(df)) stop("df missing plot_key column")
  
  
  # Merge cluster ID
  df <- dplyr::left_join(df, cluster_lookup, by = "plot_key")
  
  # Clean up — keep key columns only
  df <- df %>%
    dplyr::select(
      all_of(c("plot_key", "source_folder", "species_id", "acc", "VegType", "cluster_id")),
      everything(), 
      -plot_id
    )
  
  
  df$dmg_type_terminal <- ifelse(!is.na(df$dmg_terminal), vegtype, NA)
  df$dmg_type_foliage  <- ifelse(!is.na(df$dmg_foliage),  vegtype, NA)
  return(df)
}


# ---- Process all gpkg vegetation data ----
all_combined <- list()

for (gpkg_path in gpkg_files) {
  message("Processing vegetation: ", gpkg_path)
  
  try({
    con <- dbConnect(RSQLite::SQLite(), gpkg_path)
    tabs <- dbListTables(con)
    
    non_spatial_tables <- c("mature_test", "regeneration_adv2", "regeneration_small")
    tables <- lapply(non_spatial_tables, function(tbl) {
      if (tbl %in% tabs) dbReadTable(con, tbl) else NULL
    })
    names(tables) <- non_spatial_tables
  
    dbDisconnect(con)
    
    # Skip file if no tables present
    if (all(sapply(tables, is.null))) {
      message("No vegetation tables found in: ", gpkg_path)
      next
    }
    # Extract source_folder
    source_folder <- basename(dirname(gpkg_path))
    
    # replace numberic species_id by acronyms
    get_species_code <- function(df) {
      merge(df, species[, c("Value", "acc")], by.x = "species_id", by.y = "Value", all.x = TRUE)
    }
    
    tables$mature_test$VegType         <- "mature"
    tables$regeneration_adv2$VegType   <- "advanced"
    tables$regeneration_small$VegType  <- "small"
    
    mature <- get_species_code(tables$mature_test)
    adv    <- get_species_code(tables$regeneration_adv2)
    small  <- get_species_code(tables$regeneration_small)
    
    mature_h <- standardize_columns(mature, "mature", source_folder)
    adv_h    <- standardize_columns(adv, "advanced", source_folder)
    small_h  <- standardize_columns(small, "small", source_folder)
    
    combined <- bind_rows(mature_h, adv_h, small_h)
    combined$source_file <- basename(gpkg_path)
    all_combined[[length(all_combined) + 1]] <- combined
  }, silent = TRUE)
}

combined_vegetation_data <- bind_rows(all_combined)

head(combined_vegetation_data)
fwrite(combined_vegetation_data, "outData/combined_vegetation_data.csv")


      