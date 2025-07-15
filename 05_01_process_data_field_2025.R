

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

# read throught photos - copy them all iin a single 'damage_photo' folder
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

# # test counts: compare counts for czechia: if my new data and old data show the same numbers
# # new one
# new_counts  <- fread("raw/collected_2023/cleaned_data_frame_2023.csv")
# old_counts  <- fread("raw/collected_2023/sub_plots_counts_DBH_2023.csv")
# 
# # filter only Czechia:
# #new_counts_cz <- 
#   new_counts %>% 
#   dplyr::filter(country == 13) %>% 
#   mutate(group = group + 100) %>% 
#   mutate(ID = paste(country, region, group, point, sep = "_"),
#          cluster = paste(region, group, sep = "_")) %>% 
#    # dplyr::filter(r_piab_n == 17)#
#     dplyr::filter(cluster == "15_102") %>%
#     View()
# 
# #old_counts_cz <- 
#   old_counts %>% 
#   dplyr::filter(country == "CZ") %>% 
#   dplyr::filter(count > 0) #%>% 
#     dplyr::filter(cluster == "15_102") %>% 
#     arrange(-count)
# #  filter(ID == "13_15_102_4")
# 


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




# LOOP ----TESt START ----------------------------------------------------
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
  dplyr::distinct()  # how to filter dusplicated values????
# need to figure it out!!
#  mutate(plot_key = paste(plot_id, source_folder, sep = "__"))




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
      `1` = "zver",
      `2` = "mechanizace",
      `3` = "mysovite",
      `4` = "ine bioticke",
      `5` = "nejasna"
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
 

View(dat_long_T)
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
# analyze on level of subplot, not yet on level of clusters

## on cluster level ----------------------------------------------------

df_cluster <- dat_subplot %>% 
  group_by(cluster, acc) %>% 
  mutate(area = n_plots *4,
         scaling_factor = 10000/area,
         stem_density = count*scaling_factor) %>%   # study site area
  summarize(stem_density = sum(stem_density, na.rm = T))


df_cluster %>% 
  ggplot(aes(x = acc, y = stem_density, fill = acc)) + 
  geom_boxplot() + geom_jitter()
  

# save files: -------------------
fwrite(dat_subplot, 'outData/subplot_full_2025.csv')
fwrite(df_cluster, 'outData/df_cluster_2025.csv')

# Save spatial subplot data with cluster IDs
st_write(subplot_all, "outData/subplot_with_clusters_2025.gpkg", delete_layer = TRUE)
st_write(subplot_all, "outData/google_my_map/subplot_with_clusters_2025.kml", driver = "KML", delete_dsn = TRUE)
