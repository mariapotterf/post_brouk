

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

# make summary stat5istics: 
# how many cluysters
# species composition
# stem density
# damage characteristics
# context information

gc()

library(terra)
library(sf)
library(DBI)
library(ggplot2)
library(dbscan)
library(data.table)
library(dplyr)
library(stringr)

# look up tables
species <- fread(paste0(raw_path, "/look_up_tables/full_sp_list.csv"))

# replace numberic species_id by acronyms
get_species_code <- function(df) {
  df <- dplyr::mutate(df, species_id = suppressWarnings(as.integer(species_id)))
  dplyr::left_join(df, species[, c("Value", "acc")], by = c("species_id" = "Value"))
}

# Test on single file -----------------------------------------------------
# paths
raw_path   <- "raw/outField"
tablet     <- "T2_AH_20250528"
file_name  <- paste0("forest_structure_survey_v2[1]_", tablet, ".gpkg")
gpkg_path  <- paste(raw_path, tablet, file_name, sep = "/") 


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


## create unique cluster ID based on distance -----------------

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


## Interpret table:  ---------------------------------------------
# get species, counts, vertical classes per plot and cluster
# investigate individual outputs
# overall, get the list of species and aboundaces per plot
tables$damage_list
tables$mature_test
tables$regeneration_adv2$height

# investigate if my species list are the same? 
identical(tables$species_list2$species, tables$species_list$species)
# YES!!!



## Create whole database: -------------------------------------------------------------

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


head(combined)

combined$cluster_id




# LOOP ----TESt START ----------------------------------------------------
# some plot_id can be the same across tables and recordiong dates:
# instead, prepare unique plot_key - a combination betwen plot_id and a recording date
# to link between spatial data and no geometry tables



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
 # print(path)
#  print(basename(dirname(path)))
  tryCatch({
    df <- st_read(path, layer = "subplot", quiet = TRUE)
    df$source_file <- basename(path)
    df$source_folder <- basename(dirname(path))
    df
  }, error = function(e) NULL)
})


# make a unique plot key: combination of plot it and redording date (folder name)
subplot_all <- do.call(rbind, subplot_all)
subplot_all <- st_transform(subplot_all, 3035)
subplot_all$plot_key <- paste(subplot_all$plot_id, subplot_all$source_folder, sep = "_")


# test the plotting
ggplot(subplot_all) +
  geom_sf() +
  ggtitle("Spatial distribution of subplots") +
  theme_minimal()


# ---- Cluster globally ----
coords <- st_coordinates(subplot_all)
db <- dbscan::dbscan(coords, eps = 30, minPts = 2)
subplot_all$cluster_id <- ifelse(db$cluster == 0, NA, db$cluster)

# For merging later
cluster_lookup <- subplot_all |> 
  st_drop_geometry() |> 
  dplyr::select(plot_id, source_folder, cluster_id, plot_key)# |>
#  mutate(plot_key = paste(plot_id, source_folder, sep = "__"))

# Save spatial subplot data with cluster IDs
st_write(subplot_all, "outData/subplot_with_clusters_global.gpkg", delete_layer = TRUE)

# ---- PHASE 2: Process vegetation data ----

# defien columsn to keep 
cols_needed <- c(# "plot_id", 
  "height", "count", "dbh",
  "dmg_terminal", "dmg_terminal_photo",      
  "dmg_foliage", "dmg_stem_horizont", "dmg_stem_horizont_photo",  
  "dmg_root_stem", "dmg_root_stem_horiz_photo", "dmg_stem",                
  "dmg_stem_cause", "dmg_root_stem_horizon", "dmg_root_stem_vert",       
  "dmg_root_stem_cause", "dmg_foliage_leaf_int", "dmg_stem_vert_photo",   
  "dmg_stem_vertical", "dmg_root_stem_vert_photo",
  "dmg_stem_grass", "dmg_term_sample", "dmg_stem_sample",       
  "dmg_foliage_sample", "dmg_term_similar", "dmg_foliage_similar",    
  "dmg_foliage_detail_photo", "dmg_foliage_overall_photo", "dmg_root_stem_sample",    
  "dmg_type", "dmg_bool", "dmg_jedinec"
)



# Harmonization function:
# read data in, fill in vegetattion type: mature, advanced, small
# make all columsn to allow merge the data
# create unique plot_key: combines plot_id and collection date
standardize_columns <- function(df, vegtype, source_folder) {
  
  df$VegType <- vegtype
  
  for (col in cols_needed) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  
  df <- df[, c("plot_id", "species_id", "acc", "VegType", cols_needed)]
  
  # Add plot_key using known folder
  df$source_folder <- source_folder
  df$plot_key <- paste(df$plot_id, df$source_folder, sep = "_")
  
  # DEBUG: Check if plot_key exists in both sides
  #if (!"plot_key" %in% names(cluster_lookup)) stop("cluster_lookup missing plot_key column")
  #if (!"plot_key" %in% names(df)) stop("df missing plot_key column")
  
  #df$dmg_type_terminal <- ifelse(!is.na(df$dmg_terminal), vegtype, NA)
 # df$dmg_type_foliage  <- ifelse(!is.na(df$dmg_foliage),  vegtype, NA)
  return(df)
}


# ---- Process all gpkg vegetation data ----
all_combined <- list()

safe_get_species_table <- function(df, vegtype, source_folder) {
  if (!is.null(df) && nrow(df) > 0) {
    df$VegType <- vegtype
    df <- get_species_code(df)
  } else {
    df <- data.frame(
      species_id = NA_integer_,
      acc = NA_character_,
      VegType = vegtype,
      height = NA, count = NA, dbh = NA,
      dmg_terminal = NA, dmg_foliage = NA,
      plot_id = NA_integer_,
      stringsAsFactors = FALSE
    )
  }
  df$source_folder <- source_folder
  return(df)
}

for (gpkg_path in gpkg_files) {
  message("Processing vegetation: ", gpkg_path)
  #gpkg_path = "raw/outField/T1_Jitka_20250514/T1_forest_structure_survey_v2.gpkg"
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
    print(source_folder)
    
    mature <- safe_get_species_table(tables$mature_test, "mature", source_folder)
    adv    <- safe_get_species_table(tables$regeneration_adv2, "advanced", source_folder)
    small  <- safe_get_species_table(tables$regeneration_small, "small", source_folder)
    
    mature_h <- standardize_columns(mature, "mature", source_folder)
    adv_h    <- standardize_columns(adv, "advanced", source_folder)
    small_h  <- standardize_columns(small, "small", source_folder)
    
    combined <- bind_rows(mature_h, adv_h, small_h)
   # combined$source_file <- basename(gpkg_path)
    #combined$plot_key <- paste(combined$plot_id, combined$source_folder, sep = "_")
    
    combined2<- combined %>% 
      dplyr::select(-plot_id, -source_folder )
    
    all_combined[[length(all_combined) + 1]] <- combined2
  }, silent = FALSE)
}

combined_vegetation_data <- bind_rows(all_combined)
str(combined_vegetation_data)

# add cluster indication from spatial data
subplot_all_df <- subplot_all %>% 
  st_drop_geometry() %>%
  # your dplyr operations continue here
  dplyr::select(-plot_id, -source_folder, -source_file ,
                -time, -plot_uuid                )

# remove folder and plot_id
combined_vegetation_data2 <- combined_vegetation_data %>% 
  full_join(subplot_all_df)

head(combined_vegetation_data2)
fwrite(combined_vegetation_data2, "outData/combined_vegetation_data.csv")


# Summary: -----------------------------------------------------------------------
# file contains aslo empty and erroneous plots
# calculate how many plot_key I have per cluster_id? - keep only oones with proper counts
# Calculate: stem density 
# dominant species
# vertical categories
# keep only correct number of plots
# analyze on level of subplot, not yet on level of clusters
dat_subplot <- combined_vegetation_data2 %>% 
  group_by(cluster_id) %>% 
  mutate(n_plots = dplyr::n_distinct(plot_key)) %>% 
  dplyr::filter(n_plots %in% c( 4:6))


# how many clusters?
n_clusters <- length(unique(dat_subplot$cluster_id))  # 35-38
n_samples_terminal <- dat_subplot %>%
  dplyr::filter(!is.na(dmg_term_sample)) %>%
  dplyr::filter(str_starts(dmg_term_sample, "T")) %>%
  dplyr::pull(dmg_term_sample) %>%
  unique() 


n_samples_foliage <- dat_subplot %>%
  dplyr::filter(!is.na(dmg_foliage_sample )) %>%
  dplyr::filter(str_starts(dmg_foliage_sample, "T")) %>%
  dplyr::pull(dmg_foliage_sample)%>%
  unique() 


n_samples_root_stem <- dat_subplot %>%
  dplyr::filter(!is.na(dmg_root_stem_sample  )) %>%
  dplyr::filter(str_starts(dmg_root_stem_sample, "T")) %>%
  dplyr::pull(dmg_root_stem_sample ) %>%
  unique() 

n_samples_stem <- dat_subplot %>%
  dplyr::filter(!is.na(dmg_stem_sample   )) %>%
  dplyr::filter(str_starts(dmg_stem_sample, "T")) %>%
  dplyr::pull(dmg_stem_sample  )%>%
  unique() 


n_samples_terminal
n_samples_foliage
n_samples_root_stem
n_samples_stem


# -----------------------------------------------------------
head(dat_subplot)

# check vzroky kmen? 
dat_subplot %>% 
  dplyr::filter(!is.na(dmg_stem_sample   )) %>%
  View()



# how many stems/plot/species? 
dat_subplot %>% 
  ggplot(aes())
      