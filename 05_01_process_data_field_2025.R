

# decode field data
# what species?
# what stem density?
# heights?
# damages?
# make a cluster based on  proximity

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

# plot_key = surrender for plot_id - note, that they can be duplicated!!! need to figure this out

# read throught photos - copy them all in a single 'damage_photo' folder
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


# Test on single file -----------------------------------------------------
# paths
raw_path   <- "raw/collected_2025"
#gpkg_path  <- paste(raw_path, tablet, file_name, sep = "/") 

# look up tables
species <- fread(paste0(raw_path, "/look_up_tables/full_sp_list.csv"))


# replace numberic species_id by acronyms
get_species_code <- function(df) {
  df <- dplyr::mutate(df, species_id = suppressWarnings(as.integer(species_id)))
  dplyr::left_join(df, species[, c("Value", "acc")], by = c("species_id" = "Value"))
}


# define columns to keep 
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


# get car parking -------------------------------

# List all car_parking.gpkg files
car_parking_files <- list.files(
  "raw/collected_2025",
  pattern = "^car_parking\\.gpkg$",
  full.names = TRUE,
  recursive = TRUE
)

print(car_parking_files)

# Read all car_parking layers into a list
car_parking_list <- lapply(car_parking_files, function(path) {
  tryCatch({
    sf_obj <- st_read(path, quiet = TRUE)
    sf_obj$source_file <- basename(path)
    sf_obj$source_folder <- basename(dirname(path))
    sf_obj
  }, error = function(e) NULL)
})

# Merge into a single sf object
car_parking_all <- do.call(rbind, car_parking_list)

# Export merged result as KML
#st_write(car_parking_all, "outData/google_my_map/merged_car_parking.kml", driver = "KML", delete_dsn = TRUE)



# Process in field data ----------------------------------------------------
# some plot_id can be the same across tables and recordiong dates:
# instead, prepare unique plot_key - a combination betwen plot_id and a recording date
# to link between spatial data and no geometry tables



# List relevant .gpkg files
gpkg_files <- list.files(
  "raw/collected_2025",
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

length(unique(subplot_all$plot_key))
# test the plotting
ggplot(subplot_all) +
  geom_sf() +
  ggtitle("Spatial distribution of subplots") +
  theme_minimal()

table(subplot_all$cluster)


# ---- Cluster globally ----
coords <- st_coordinates(subplot_all)
db <- dbscan::dbscan(coords, eps = 30, minPts = 2)
subplot_all$cluster <- ifelse(db$cluster == 0, NA, db$cluster)

# For merging later
cluster_lookup <- subplot_all |> 
  st_drop_geometry() |> 
  dplyr::select(plot_id, source_folder, cluster, plot_key) |>
  dplyr::distinct()  # keep the latest version to avoid duplicated records
#  mutate(plot_key = paste(plot_id, source_folder, sep = "__"))

#  subplots 518_T4, 519_T4 had cluster NA -> changed manually to cluster 143


# ---- PHASE 2: Process vegetation data ----


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

# keep the table even is no species is present - fill in with NAs
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

# run for loop to prcess the all gpkgs
for (gpkg_path in gpkg_files) {
  message("Processing vegetation: ", gpkg_path)
  #gpkg_path = "raw/collected_2025/T1_Jitka_20250514/T1_forest_structure_survey_v2.gpkg"
  try({
    con <- dbConnect(RSQLite::SQLite(), gpkg_path)
    tabs <- dbListTables(con)
    
    non_spatial_tables <- c("mature_test", "regeneration_adv2", "regeneration_small", 
                            "context")
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
    
    # combine tree layeras: mature, adv, small
    combined_veg <- bind_rows(mature_h, adv_h, small_h)
  
      
    combined_veg<- combined_veg %>% 
      dplyr::select(-plot_id, -source_folder )
    
    
    # add context information per subplot (4m2)
    context <- tables$context %>%
      dplyr::select(-context_uuid, -plot_uuid) %>%
      dplyr::mutate(source_folder = source_folder,
                    plot_key = paste(plot_id, source_folder, sep = "_"))
    
    # combine if teh information is split between two columns (eg is NA in one row but filled in another)
    # keeps simply teh first value that is not NA - eg, if I have true/false, it keeeps TRUE (fiurst value)
    context_dedup <- context %>%
      group_by(plot_key) %>%
      summarise(across(
        .cols = -c(plot_id, source_folder),
        .fns = ~ reduce(.x, coalesce),
        .names = "{.col}"
      )) %>%
      ungroup()
    
    # Left join context info to vegetation
    combined2 <- left_join(combined_veg, context_dedup, by = c("plot_key"))
    
    all_combined[[length(all_combined) + 1]] <- combined2
  }, silent = FALSE)
}

combined_vegetation_data <- bind_rows(all_combined)
str(combined_vegetation_data)
#View(combined_vegetation_data)

# add cluster indication from spatial data
subplot_all_df <- subplot_all %>% 
  st_drop_geometry() %>%
  # your dplyr operations continue here
  dplyr::select(-plot_id, -source_folder, -source_file ,
                -time, -plot_uuid                )

# remove folder and plot_id
combined_vegetation_data2 <- combined_vegetation_data %>% 
  full_join(subplot_all_df)

# add damage cause infomration
combined_vegetation_data_recode <- 
  combined_vegetation_data2 %>%
  mutate(
    # Height labels (used only when VegType is small/advanced)
    height = case_when(
      VegType == "small" ~ recode(as.character(height),
                                  "1" = "0.2–0.4",
                                  "2" = "0.4–0.6",
                                  "3" = "0.6–0.8",
                                  "4" = "0.8–1.0",
                                  "5" = "1.0–1.3",
                                  "6" = "1.3–2.0",
                                  .default = NA_character_
      ),
      VegType == "advanced" ~ recode(as.character(height),
                                     "1" = "2–4",
                                     "2" = ">4",
                                     .default = NA_character_
      ),
      TRUE ~ NA_character_
    ),
    
    # DBH labels (used only when VegType is mature)
    dbh = if_else(VegType == "mature",
                        recode(as.character(dbh),
                               "1" = "10–20cm",
                               "2" = "20–40cm",
                               "3" = "40–60cm",
                               "4" = ">60cm",
                               .default = NA_character_
                        ),
                        NA_character_
    )
  ) %>% 
  mutate(
    dmg_terminal = recode(
      dmg_terminal,
      "1" = "zivy",
      "2" = "chybajuci",
      "3" = "odumrety",
    )
  ) %>% 
  mutate(
    dmg_type = recode(
      dmg_type,
      "1" = "terminal",
      "2" = "kmen",
      "3" = "baze_kmene",
      "4" = "olistení"
    )
  ) %>% 
  mutate(across(
    .cols = contains("_cause"),  # for cause of damage
    .fns = ~ dplyr::recode(
      .,
      "1" = "zver",
      "2" = "mechanizace",
      "3" = "mysovite",
      "4" = "ine bioticke",
      "5" = "nejasna"
    )
  )) %>% 
    mutate(
      dmg_stem = recode(
        dmg_stem,
        "1" = "cerstve",
        "2" = "stare"
      )
    ) %>% 
  mutate(
    dmg_stem_horizont = recode(
      dmg_stem_horizont,
      "1" = "0-20",
      "2" = "20-40",
      "3" = "40-60",
      "4" = "60-80",
      "5" = "80-100"
    )
  ) %>% 
  mutate(
    dmg_stem_vertical = recode(
      dmg_stem_vertical,
      "1" = "0-10",
      "2" = "10-20",
      "3" = "20-30",
      "4" = ">30"
    )
  ) %>% 
  mutate(
    dmg_root_stem_horizon = recode(
      dmg_root_stem_horizon,
      "1" = "0-20",
      "2" = "20-40",
      "3" = "40-60",
      "4" = "60-80",
      "5" = "80-100"
    )
  ) %>% 
  mutate(
    dmg_root_stem_vert = recode(
      dmg_root_stem_vert,
      "1" = "0-5",
      "2" = "5-10",
      "3" = ">10"
    )
  ) %>%
    mutate(
      dmg_root_stem = recode(
        dmg_root_stem,
        "1" = "cerstve",
        "2" = "stare"
      )
    ) %>% 
    mutate(
      dmg_foliage = recode(
        dmg_foliage,
        "1" = "0-20",
        "2" = "20-40",
        "3" = "40-60",
        "4" = "60-80",
        "5" = "80-100"
      )
    ) %>% 
    mutate(
      dmg_foliage_leaf_int = recode(
        dmg_foliage_leaf_int,
        "1" = "0-20",
        "2" = "20-40",
        "3" = "40-60",
        "4" = "60-80",
        "5" = "80-100"
      )
    ) %>% 
  mutate(count = as.integer(count))

dat_subplot <- combined_vegetation_data_recode %>% 
  group_by(cluster) %>% 
  mutate(n_plots = dplyr::n_distinct(plot_key)) #%>% 
  #dplyr::filter(n_plots >4 ) #%>% 
  
length(unique(dat_subplot$cluster))

table(dat_subplot$cluster)

# create only a database to chcek for poltID, damage type, sample and photo: damage indication -----------------
dat_dmg <-  dat_subplot %>%
  dplyr::select(
    species_id,
    acc,
    VegType,
    starts_with("dmg_"),
    plot_key,
    comments,
    cluster,
    ends_with("_sample"),
  )

# convert from wide format to long format: keep record of all samples and photos 
# filter only record that have sample, and also a photo (many have instead of sample name 'mraz')
dat_dmg_filtered <- dat_dmg %>%
  dplyr::filter(
    if_any(ends_with("_sample"), ~ !is.na(.)) |
      if_any(ends_with("_photo"), ~ !is.na(.) & . != "")
  )

dat_dmg_filtered_min <- dat_dmg_filtered %>%
  dplyr::select(
    acc, plot_key, cluster,VegType,
    matches("(_sample|_photo)$")
  )


# Define a mapping of sample-photo pairs: which sample belong to wchih photo (as I have several categorises for photos)
sample_photo_map <- tibble::tribble(
  ~part,         ~sample_col,              ~photo_col,                  ~photo_type,
  "terminal",    "dmg_term_sample",        "dmg_terminal_photo",        NA,
  "stem",        "dmg_stem_sample",        "dmg_stem_horizont_photo",   "horizontal",
  "stem",        "dmg_stem_sample",        "dmg_stem_vert_photo",       "vertical",
  "foliage",     "dmg_foliage_sample",     "dmg_foliage_detail_photo",  "detail",
  "foliage",     "dmg_foliage_sample",     "dmg_foliage_overall_photo", "overall",
  "root_stem",   "dmg_root_stem_sample",   "dmg_root_stem_horiz_photo", "horizontal",
  "root_stem",   "dmg_root_stem_sample",   "dmg_root_stem_vert_photo",  "vertical"
)

# Pivot to long based on the mapping
dat_long <- sample_photo_map %>%
  pmap_dfr(function(part, sample_col, photo_col, photo_type) {
    dat_dmg_filtered_min %>%
      select(plot_key, acc, cluster, 
             sample = all_of(sample_col),
             photo = all_of(photo_col)) %>%
      mutate(part = part,
             photo_type = photo_type)
  }) %>%
  dplyr::filter(!is.na(sample) & sample != "")

# filter only proper samples, not descripitons (Okus, mraz)
dat_long_T <- dat_long %>%
  dplyr::filter(str_starts(sample, "T")) %>% 
  mutate(photo = str_remove(photo, "^DCIM/")) %>% 
  arrange(sample) %>% # remove DCIM/ from photo name
  group_by(photo, sample) %>%
  slice_max(order_by = plot_key, n = 1, with_ties = FALSE) %>% # keep only later data if recorded deveral times
  ungroup()
#df_unique <- df %>%
length(unique(dat_long_T$sample))

#View(dat_long_T)
# export final table as csv
#fwrite(dat_long_T, 'outTable/samples_list.csv')
# export final table as csv
fwrite(dat_long_T, 'outShare/samples_list.csv')

# Summary: -----------------------------------------------------------------------
# file contains aslo empty and erroneous plots
# calculate how many plot_key I have per cluster? - keep only oones with proper counts
# Calculate: stem density 
# dominant species
# vertical categories
# keep only correct number of plots
# analyze on level of subplot
head(dat_subplot)
#View(dat_subplot)





# save files: -------------------
fwrite(dat_subplot, 'outData/subplot_full_2025.csv')

# Save spatial subplot data with cluster IDs
st_write(subplot_all, "outData/subplot_with_clusters_2025.gpkg", delete_layer = TRUE)
st_write(subplot_all, "outData/google_my_map/subplot_with_clusters_2025.kml", driver = "KML", delete_dsn = TRUE)



# Clean upda data for Michal: -------------------------------------------------------
# read all 2025 data, sf
# clean up and recode
# keep empty subplots & plots
# export final veg data and gpkg
# keep all records (even with missing information) as the damage can still be recorded!

gc()

library(terra)
library(sf)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggpubr)

## read data from 2025 ----------------
dat25_subplot    <- data.table::fread("outData/subplot_full_2025.csv")   # subplot-level table
dat25_sf         <- sf::st_read("outData/subplot_with_clusters_2025.gpkg")          # subplot spatial data

# Select and rename
dat25_sf_min <- dat25_sf %>%
  dplyr::select(plot_key, cluster,plot_id)

length(unique(dat25_subplot$plot_key))  # 1009


dat25_subplot <- dat25_subplot %>% 
  rename(
    subplot = plot_key,
    plot = cluster,
    vegtype = VegType,
    hgt = height,
    n = count,
    species = acc,
    clear = clearing,
    grndwrk = site_prep
  ) 


# filter bad plots
dat25_subplot <- dat25_subplot[!(subplot %in% dat25_subplot)]



n25_subplots <- dat25_subplot %>%
  #filter(!is.na(plot), !is.na(subplot)) %>%
  distinct(plot, subplot) %>%       # drop species/vegtype duplicates
  count(plot, name = "n_subplots") #%>%
#arrange(plot)


# recode teh data
## Recode field data -----------------------------------
# guestimate dbh and ba per individual based on height distribution 
dat25_subplot_recode <- dat25_subplot %>% 
  # remove if no species is defined
  #filter(!is.na(species)) %>% # yes, I can safely exclude those, they are not representing 'empty plots' 
  # adjust tree species naming
  # dplyr::mutate(
  #   species = str_trim(species),
  #   species = dplyr::case_when(
  #     species == "abal" ~ "absp",
  #     species %in% c("quro", "quru") ~ "qusp",
  #     TRUE ~ species
  #   )
  # ) %>% 
  mutate(
    # Create a numeric height estimate (keep original height class string)
    hgt_est = case_when(
      vegtype == "mature" & dbh == "10–20cm" ~ 10.0,
      vegtype == "mature" & dbh == "20–40cm" ~ 20.0,
      vegtype == "mature" & dbh == "40–60cm" ~ 30.0,
      hgt == "0.2–0.4"                       ~ 0.3,
      hgt == "0.4–0.6"                       ~ 0.5,
      hgt == "0.6–0.8"                       ~ 0.7,
      hgt == "0.8–1.0"                       ~ 0.9,
      hgt == "1.0–1.3"                       ~ 1.2,
      hgt == "1.3–2.0"                       ~ 1.7,
      hgt == "2–4"                           ~ 3.0,
      hgt == ">4"                            ~ 5.0,
      TRUE                                   ~ NA_real_
    ),
    
    # Estimate DBH (already numeric)
    dbh_est = case_when(
      vegtype == "mature" & dbh == "10–20cm" ~ 15.0,
      vegtype == "mature" & dbh == "20–40cm" ~ 30.0,
      vegtype == "mature" & dbh == "40–60cm" ~ 50.0,
      vegtype == "mature" & dbh == ">60cm"   ~ 70.0,
      hgt == "0.2–0.4"                       ~ 0.3,
      hgt == "0.4–0.6"                       ~ 0.5,
      hgt == "0.6–0.8"                       ~ 0.7,
      hgt == "0.8–1.0"                       ~ 0.9,
      hgt == "1.0–1.3"                       ~ 1.2,
      hgt == "1.3–2.0"                       ~ 1.7,
      hgt == "2–4"                           ~ 3.0,
      hgt == ">4"                            ~ 5.0,
      TRUE                                   ~ NA_real_
    )
  )# %>% 
# # Calculate basal area
# mutate(
#   basal_area_cm2 = pi * (dbh_est / 2)^2,
#   ba_total_cm2   = basal_area_cm2 * n,
#   ba_total_m2    = ba_total_cm2 / 10000,
#   #ba_ha_m2       = ba_total_m2 * scaling_factor
# ) #%>% 
#left_join(traits_full)  # add trait table to identify early vs late seral

length(unique(dat25_subplot_recode$subplot))

# 

fwrite(dat25_subplot_recode, "outDataShare/dat25_subplot_recode.csv")
sf::st_write(dat25_sf_min, "outDataShare/dat25_sf_min.gpkg", delete_dsn = TRUE)



# Summarize infor for overall presentation ------------------------------------
# Make sure the data is a data.table (should already be based on your structure)
dat <- dat25_subplot_recode

# Filter only valid tree records (non-NA and n > 0)
dat_trees <- dat %>% 
  filter(!is.na(n) & n > 0)

# 1. How many trees do I have?
total_trees <- sum(dat_trees$n, na.rm = TRUE)

# 2. How many trees per species?
trees_per_species <- dat_trees %>%
  group_by(species) %>%
  summarise(total_n = sum(n, na.rm = TRUE)) %>%
  arrange(desc(total_n)) %>% 
  mutate(share = round(total_n/total_trees*100,2))

# 3. How many plots do not have any n > 0 (i.e., no trees present)?
plots_no_trees <- dat %>%
  group_by(plot) %>%
  summarise(total_n = sum(n, na.rm = TRUE)) %>%
  filter(total_n == 0) %>%
  nrow()

# 4. How many unique subplots?
n_subplots <- dat %>%
  distinct(subplot) %>%
  nrow()

# 5. How many unique plots?
n_plots <- dat %>%
  distinct(plot) %>%
  nrow()

# 6. Species richness per subplot (number of unique species per subplot)
richness_per_subplot <- dat_trees %>%
  group_by(subplot) %>%
  summarise(species_richness = n_distinct(species))

summary(richness_per_subplot$species_richness)


# 7. Species richness per plot
richness_per_plot <- dat_trees %>%
  group_by(plot) %>%
  summarise(species_richness = n_distinct(species))

summary(richness_per_plot$species_richness)

# Count number of subplots by richness category
subplot_richness_freq <- richness_per_subplot %>%
  count(species_richness) %>% 
  mutate(share = round(n/n_subplots*100,1))

# Count number of plots by richness category
plot_richness_freq <- richness_per_plot %>%
  count(species_richness) %>% 
  mutate(share = round(n/n_plots*100,1))


# 8. management summary

# Select management-related columns
mgmt_cols <- c("clear", "grndwrk", "logging_trail", "windthrow", "standing_deadwood", "planting", "anti_browsing")

# Summarize presence/absence/NA of management types per subplot
management_summary <- dat25_subplot_recode %>%
  distinct(subplot, .keep_all = TRUE) %>%  # One row per subplot
  select(subplot, all_of(mgmt_cols)) %>%
  pivot_longer(cols = all_of(mgmt_cols), names_to = "management_type", values_to = "value") %>%
  mutate(value = replace_na(value, 0),  # Treat NA as 0
         status = ifelse(value == 1, "present", "absent")) %>%
  group_by(management_type, status) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(management_type, status) %>% 
  mutate(share = round(n/n_subplots*100, 1))


management_summary


# 9. species per plot 

species_plot_counts <- dat25_subplot_recode %>%
  filter(!is.na(n) & n > 0) %>%
  distinct(plot, species) %>%
  group_by(species) %>%
  summarise(n = n(), .groups = "drop") %>% 
  mutate(share = round(n/n_plots*100)) %>% 
  arrange(desc(share)) %>% 
  slice_max(share, n = 10) %>% 
  mutate(species = factor(species, levels = rev(species))) 








## Make plots --------------------------------------------

### Tree species composition ---------------------
# Czech labels for species
species_labels_cz <- c(
  "piab" = "Smrk ztepilý",
  "besp" = "Bříza",
  "pisy" = "Borovice lesní",
  "quru" = "Dub",
  "fasy" = "Buk lesní",
  "lade" = "Modřín opadavý",
  "saca" = "Vrba jíva", # mléč
  "acps" = "Javor klen",
  "soau" = "Jeřáb ptačí",
  "potr" = "Topol osika"
)

# Take top 10 species
top_species <- trees_per_species %>%
  slice_max(share, n = 10) %>%
  mutate(species = factor(species, levels = rev(species)))# %>%   # Reverse for top-down order
  
# Plot

my_greens <- colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(nrow(top_species))

p_species_share_stems <- ggplot(top_species, aes(x = share, y = species, fill = species)) +
  geom_col(color = "black", width = 0.6) +
  scale_fill_manual(values = my_greens) +
  scale_y_discrete(labels = species_labels_cz) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(
    x = "Podíl jedinců [%]",
    y = "",
    title = "Dřeviny podle zastoupení"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10, face = "italic"),
    text = element_text(size = 11)
  )


p_species_share_plots <- species_plot_counts %>% 
  ggplot(aes(x = share, y = species, fill = species)) +
  geom_col(color = "black", width = 0.6) +
  scale_fill_manual(values = my_greens) +
  scale_y_discrete(labels = species_labels_cz) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(
    x = "Podíl ploch [%]",
    y = "",
    title = "Výskyt na plochách"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10, face = "italic"),
    text = element_text(size = 11)
  )

ggarrange(p_species_share_stems, p_species_share_plots,
          labels = "auto")



### Species richness ---------------------------------------

p_rich_plot <- ggplot(subplot_richness_freq, aes(x = factor(species_richness), y = share)) +
  geom_col(fill = "darkolivegreen3", color = "black") +
  labs(
    x = "Počet druhů \n(ploška)",
    y = "Podíl plošek [%]",
    title = ""
  ) +
  theme_classic()


p_rich_sub <- ggplot(plot_richness_freq, aes(x = factor(species_richness), y = share)) +
  geom_col(fill = "darkseagreen4", color = "black") +
  labs(
    x = "Počet druhů \n(plocha)",
    y = "Podíl ploch [%]",
    title = ""
  ) +
  theme_classic()
ggarrange(p_rich_plot,p_rich_sub)


### Management ---------------------------------------------


# Czech labels for activities
activity_labels <- c(
  "clear" = "Sanitárni těžba",
  "grndwrk" = "Příprava půdy",
  "planting" = "Výsadba",
  "anti_browsing" = "Ochrana proti zvěři",
  "logging_trail" = "Vyklizovací linka",
  "windthrow" = "Vývrat",
  "standing_deadwood" = "Stojící mrtvé dřevo"
)

# Convert management_summary to plotting format
mng_sub_conv <- management_summary %>%
  rename(activity = management_type) %>%
  mutate(applied = ifelse(status == "present", 1, 0),
         proportion = ifelse(applied == 1, share, -share),
         applied = ifelse(applied == 1, "Presence", "Absence")) %>%
  mutate(activity = factor(activity, levels = names(activity_labels))) #%>%
  #mutate(activity = fct_relabel(activity, ~ activity_labels[.]))  # relabel to Czech

# Order activities by proportion of presence
activity_order <- mng_sub_conv %>%
  filter(applied == "Presence") %>%
  arrange(desc(proportion)) %>%
  pull(activity)

mng_sub_conv <- mng_sub_conv %>%
  mutate(activity = factor(activity, levels = rev(activity_order)))

# Plot
ggplot(mng_sub_conv, aes(x = proportion, y = activity, fill = applied)) +
  geom_col(width = 0.2, color = 'black') +
  scale_x_continuous(
    breaks = seq(-100, 100, 50),
    limits = c(-100, 100),
    labels = function(x) paste0(abs(x), "%")
  ) +
  scale_fill_manual(values = c("Presence" = "red", "Absence" = "darkgreen")) +
  geom_vline(xintercept = 0, color = 'darkgrey', lty = 'dashed') +
  labs(
    x = "Podíl plošek [%]",
    y = "",
    title = "",
    fill = ""
  ) +
  scale_y_discrete(labels = activity_labels) +
  annotate("text", x = -80, y = length(activity_order) + 0.5, label = "Nevyskytuje se", hjust = 0, size = 3, fontface = "bold") +
  annotate("text", x =  80, y = length(activity_order) + 0.5, label = "Vyskytuje se", hjust = 1, size = 3, fontface = "bold") +
  theme_classic() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10, face = "italic"),
    legend.position = "none",
    text = element_text(size = 10)
  )




## Heigh distribution -------------------------
hist(dat25_subplot_recode$hgt_est, breaks = 60 )

# Expand dataset by n (each row = 1 individual)
hgt_expanded <- dat25_subplot_recode %>%
  filter(!is.na(n) & n > 0, !is.na(hgt_est)) %>%
  select(hgt_est, n) %>%
  tidyr::uncount(n)

summary(hgt_expanded)
# Calculate 90% quantile
q90 <- quantile(hgt_expanded$hgt_est, 0.9)

# Plot CDF with 90% quantile line
ggplot(hgt_expanded, aes(x = hgt_est)) +
  stat_ecdf(geom = "step", color = "steelblue", size = 1) +
  geom_vline(xintercept = q90, linetype = "dashed", color = "red") +
  annotate("text", x = q90, y = 0.92, label = paste0("90 % < ", round(q90, 2), " m"),
           hjust = 0, vjust = 0, color = "red", size = 4) +
  labs(
    title = "Kumulativní distribuce výšek stromků",
    x = "Výška (m)",
    y = "Kumulativní podíl jedinců"
  ) +
  theme_classic()

dat25_subplot_recode %>%
  filter(!is.na(n) & n > 0) %>%
  ggplot(aes(x = hgt_est)) +
  geom_histogram(binwidth = 0.3, fill = "grey", color = "black") +
  labs(
    title = "Distribuce výšek",
    x = "Odhadovaná výška (m)",
    y = "Počet jedinců"
  ) +
  geom_vline(xintercept = q90+0.2, linetype = "dashed", 
             color = "red") +
  annotate("text", x = q90+2, y = 400, label = paste0("90 % < ", round(q90, 2), " m"),
           hjust = 0, vjust = 0, color = "red", size = 4) +
  theme_classic2()


dat25_subplot_recode %>%
  filter(!is.na(n) & n > 0, !is.na(hgt_est), hgt_est > 5) %>%
  summarise(total_individuals = sum(n))


# 2. Boxplot of heights per species
library(forcats)
dat25_subplot_recode %>%
  filter(!is.na(n) & n > 0) %>%
  filter(species %in% top_species$species) %>% 
  ggplot(aes(x = fct_reorder(species, -hgt_est, .fun = median, na.rm = TRUE), y = hgt_est,
             fill = species)) +
  scale_x_discrete(labels = species_labels_cz) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_hline(yintercept = 1, lty = 'dashed', col = 'grey') +
  #geom_violin() +
  labs(
    title = "Boxplot výšek stromků podle druhu",
    x = "Druh (seřazeno podle výšky)",
    y = "Odhadovaná výška (m)"
  ) +
  scale_fill_manual(values = my_greens) +
  coord_cartesian(y = c(0,10))+
  theme_classic2() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))
  

# Summary table of median heights by species (only top species)
summary_medians <- dat25_subplot_recode %>%
  filter(!is.na(n) & n > 0) %>%
  filter(species %in% top_species$species) %>%
  group_by(species) %>%
  summarise(
    median_hgt = median(hgt_est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_hgt)) %>%
  mutate(species_cz = species_labels_cz[species]) %>%
  select(species, species_cz, median_hgt)

summary_medians
  