

#             Temporal change in species composition 
#                          and structure


# read data: 
# - field data from 2023 & 2025
# - get sites history - from tree density (manually mapped)

# get summary statistics
# merge data from 2023 & 2025 - compare mean heights, CVs, richness, shannon between years
# cross scale interaction: 

#  - vertical structure
# effect of pre-disturbance condistions


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
library(corrplot)
library(GGally)
library(vegan) # for diversity indices

library(mgcv)
library(ggeffects)

library(RColorBrewer)

library(ggridges)
library(scales)


theme_set(theme_classic2(base_size = 10) +
            theme(axis.title = element_text(size = 10),
                  axis.text  = element_text(size = 10)))


# Read files: --------------
# get ecological traits
traits_full <- fread('outData/my_species_traits.csv')

traits_full<- traits_full %>% 
  mutate(
  seral_stage = if_else(Shade_tolerance < 2.5, "early", "late") # low value are less totalera = sun loving, high values are 
)

## --- Field data: 2023 
dat23_subplot    <- data.table::fread("outData/subplot_full_2023.csv")   # subplot-level table
dat23_sf         <- sf::st_read("outData/sf_context_2023.gpkg")          # subplot spatial data

# Select and rename
dat23_sf_min <- dat23_sf %>%
  dplyr::select(subplot = ID, plot = cluster)

## read data from 2025 ----------------
dat25_subplot    <- data.table::fread("outData/subplot_full_2025.csv")   # subplot-level table
dat25_sf         <- sf::st_read("outData/subplot_with_clusters_2025.gpkg")          # subplot spatial data

# Select and rename
dat25_sf_min <- dat25_sf %>%
  dplyr::select(plot_key, cluster,plot_id)

length(unique(dat25_subplot$plot_key))  # 1009

# needed columns to merge field data from 2023&2025:
management_types_v <- c("clear", "grndwrk", "logging_trail", "planting", "anti_browsing")

target_cols <- c("plot",
                 "subplot",
                 "species",
                 "vegtype",
                 "hgt",
                 "n",
                 "dbh",
                 "year",
                 management_types_v)


# I have few other species - keep simply as 'other' to not misrepresent their occerence elsewhere
# just remove


#ots1
#138_T1_JC_20250717 - not identified
#358_T2_AH_20250827 - Sambucus nigra
#396_T2_AH_20250827 - Frangula alnus
#495_T4_TP_20250827 - jablon malus
#619_T2_AH_20250827 - Pinus strobus]



# filter only relevant columns: 
dat25_subplot_sub <- dat25_subplot %>% 
  filter(!is.na(cluster)) %>% # filter empty clusters - happend in 3 subplots that they do not exist in points
  dplyr::filter(naruseni == TRUE & nzchar(trimws(photo_e))) %>% # filter out the errorroneous records
  dplyr::select(acc, VegType, height, count, dbh, plot_key, cluster, 
                clearing, site_prep, logging_trail, windthrow, standing_deadwood, planting, anti_browsing) %>% 
  rename(
    subplot = plot_key,
    plot = cluster,
    vegtype = VegType,
    hgt = height,
    n = count,
    species = acc,
    clear = clearing,
    grndwrk = site_prep
         ) %>% 
  mutate(year = "2025") %>% 
  # target column order for the bind
  dplyr::select(all_of(target_cols))

nrow(dat25_subplot_sub)
head(dat25_subplot_sub)
unique(dat25_subplot_sub$species)

subplot_to_plot25 <- dat25_subplot_sub %>%
  distinct(subplot, plot)

# Expand dataset to fit 2023 dataset - have all species, but filled with NA instead of species
# having markjed as NA (in 2025)
species_vec25 <- unique(dat25_subplot_sub$species)

# remove otsp and empty species - now emply subplots will be still recorded
species_vec25 <- species_vec25[species_vec25 != ""]
species_vec25 <- species_vec25[species_vec25 != "otsp1"]


dat25_expanded <-
  dat25_subplot_sub %>%
  as_tibble() %>%
    dplyr::select(-plot) %>% 
  tidyr::complete(
    subplot,
    vegtype = c("small", "advanced", "mature"),
    species = species_vec25,
    fill = list(
      n = NA_integer_,
      dbh = NA_character_,
      hgt = NA_character_
    )
  ) %>%
  # Reattach plot info
  left_join(subplot_to_plot25, by = "subplot") %>%
  mutate(year = "2025") %>%
  # Reorder columns if needed
  dplyr::select(all_of(target_cols)) %>% 
  dplyr::filter(!species %in% c("ots1", "")) 

# How many plots/subplots do not have any stems present? 
dat25_expanded %>%
  group_by(plot) %>%
  summarise(all_stems_missing = all(is.na(n) | n == 0)) %>%
  filter(all_stems_missing) #%>%  # just 1???
  #nrow()

dat25_expanded %>%
  group_by(subplot) %>%
  summarise(all_stems_missing = all(is.na(n) | n == 0)) %>%
  filter(all_stems_missing) #%>%
# 196 ~ 20%

length(unique(dat25_expanded$species))
length(unique(dat25_expanded$subplot))
# 375_T2_AH_20250827


## Field data 2023: subplot level --------------------------------------
dat23_subplot_sub <- 
  dat23_subplot %>% 
  dplyr::select(vegtype, species, ID, cluster,n, hgt, dbh,
                clear, grndwrk, logging_trail, planting, anti_browsing) %>% 
  rename(subplot = ID,
         plot = cluster) %>% 
  mutate(year = "2023") %>% 
  mutate(
    species = as.character(species),
    vegtype = as.character(vegtype),
    hgt     = as.character(hgt),
    n       = as.integer(n),
    dbh     = as.character(dbh),
    subplot = as.character(subplot),
    plot    = as.character(plot),
    year    = as.factor(year)
  ) %>%
  # assure teh same order
  dplyr::select(all_of(target_cols))

dat25_subplot_sub <- dat25_expanded %>%
  mutate(
    species = as.character(species),
    vegtype = as.character(vegtype),
    hgt     = as.character(hgt),
    n       = as.integer(n),
    dbh     = as.character(dbh),
    subplot = as.character(subplot),
    plot    = as.character(plot),
    year    = as.factor(year)
  ) %>%
  select(all_of(target_cols))

#  remove the erroneous subplots, if i have 6 subplots in cluster: 

# plots to check :
# 15_102 ???  - kept 6
# 15_124
# 15_145
# 26_101
# 26_142
# --- 2025: 
# 48 - 6 -> remove 375_T2_AH_20250827
# 117 - 6 -> 306_T3_JL_20250717
# 125 - 6 -> 424_T4_TP_20250827            
# 147 - 6 -> 539_T4_TP_20250827            
# 167 - 6 -> 644_T4_TP_20250827
# 184 - 6 -> 741_T4_TP_20250827
# NA - has only 3

# subplots: 
# 13_15_102_5
# 13_15_124_2
# 13_15_145_4
# 13_26_101_2
# 13_26_142_1

# bind field data from both years:
dat_subplots <- bind_rows(dat23_subplot_sub, dat25_subplot_sub)


bad_subplots <- c("306_T3_JL_20250717",
                  "424_T4_TP_20250827", 
                  "539_T4_TP_20250827",
                  "644_T4_TP_20250827",
                  "741_T4_TP_20250827",
                  "13_15_102_5","13_15_124_2",
                  "13_15_145_4","13_26_101_2","13_26_142_1") 

dat_subplots <- dat_subplots[!(subplot %in% bad_subplots)]

# remove whole plots (robust to numeric/character mix)
bad_plots <- c("15_104", "26_134")

dat_subplots <- dat_subplots %>%
  dplyr::filter(!as.character(plot) %in% bad_plots)

n_subplots <- dat_subplots %>%
  #filter(!is.na(plot), !is.na(subplot)) %>%
  distinct(plot, subplot) %>%       # drop species/vegtype duplicates
  count(plot, name = "n_subplots") #%>%
#arrange(plot)

# combine years 2023 and 2025
dat_subplots <- dat_subplots %>% 
  left_join(n_subplots) %>% 
  dplyr::filter(n_subplots == 5)


## Tree-based pre-disturbance history (vector layers) ------------------------
convex_hull      <- terra::vect("raw/pre-disturb_history_trees/cvx_hull_completed_3035.gpkg")
pre_trees        <- terra::vect("raw/pre-disturb_history_trees/pre-disturbance_trees_3035.gpkg")

## Data clean up 

# Get polygon area (in m²) & perimeter - this is convex hull + buffer around it to avoid edge effects
convex_hull$area_m2     <- expanse(convex_hull, unit = "m")
convex_hull$perimeter_m <- perim(convex_hull)

# clean up tree characteristics: get species, ..
pre_trees_df <- as.data.frame(pre_trees) %>%
  mutate(
    original_species = species,
    species = case_when(
      is.na(species) ~ "piab",
      species == "l" ~ "deciduous",
      species == "d" ~ "deciduous",
      species == "s" ~ "piab",
      TRUE ~ species
    ),
    state = case_when(
      original_species == "s" ~ "dry",
      TRUE ~ "living"
    )) %>%
  dplyr::select(-original_species) 

# Reattach cleaned attributes back to geometry
pre_trees_3035_clean          <- pre_trees
values(pre_trees_3035_clean)  <- pre_trees_df

## Clean up convex hull characteristics 
cvx_clean <- convex_hull %>% 
  as.data.frame() %>%   # attributes only (no geometry)
  mutate(
    # rok_disturbancia like "2019-…", or "nie"
    disturbance_year = case_when(
      rok_disturbancia == "nie" ~ 2023L,
      TRUE ~ suppressWarnings(as.integer(str_extract(rok_disturbancia, "^\\d{4}")))
    ),
    disturbance_note = case_when(
      rok_disturbancia == "nie" ~ "nie",
      TRUE ~ str_trim(str_remove(rok_disturbancia, "^\\d{4}[- ]*"))
    ),
    forest_year = suppressWarnings(as.integer(rok_les)),
    disturbance_length = disturbance_year - forest_year
  ) %>% 
  #rename(plot = cluster) %>% 
  dplyr::select(-rok_les,-rok_disturbancia, 
                - disturbance_note
                ) %>% 
  mutate(plot_comb = dplyr::coalesce(as.character(cluster_2023),
                                     as.character(cluster_2025)))


# Write cleaned attributes back to the terra object
convex_hull_3035_clean  <-convex_hull[, c("cluster_2023","cluster_2025")]
values(convex_hull_3035_clean)  <- cvx_clean


# split the dataset to isolate only shared plots between 2 years
# make oone table for year 2023
cvx_both = cvx_clean %>% 
  dplyr::filter(status == "both")


# kolko ploch mame? 
table(convex_hull_3035_clean$status)

# both only_2023 only_2025 
# 130         9        74 

# Add plot (convex hull) info to each tree  (spatial join)
pre_trees_cvx_joined <- terra::intersect(pre_trees_3035_clean, convex_hull_3035_clean)

### Get pre-disturbnace stem density per plot 
cvx_df <- as.data.frame(pre_trees_cvx_joined)

####  Plot level = CVX: stem density 
cvx_stem_density <- cvx_df %>%
#  ungroup(.) %>% 
  group_by(common_cluster_ID, plot_comb, status, cluster_2023, cluster_2025,
           disturbance_year, forest_year, disturbance_length) %>%  #, area_m2
  dplyr::reframe(
    pre_dist_trees_n = n(),
    area_m2 = mean(area_m2, na.rm  = T),  # keep are here instead of grouping
    pre_dist_dens_ha = pre_dist_trees_n / area_m2 * 10000
  )

# filter pre-disturbnace trees per year
cvx_stem_density23 <- cvx_stem_density %>% 
  filter(!is.na(cluster_2023) )

cvx_stem_density25 <- cvx_stem_density %>% 
  filter(!is.na(cluster_2025) ) %>% 
  mutate(cluster_2025 = as.factor(cluster_2025))

# add to field data year by year
dat_subplots23 <- dat_subplots %>% 
  filter(year == "2023") %>% 
  left_join(cvx_stem_density23, by = c("plot" = "cluster_2023" )) %>% 
  dplyr::select(-cluster_2025, -common_cluster_ID, -plot_comb ) %>% 
  mutate(time_snc_full_disturbance = 2023 - disturbance_year,
         time_snc_part_disturbance = 2023 - forest_year + 1)

dat_subplots25 <- 
  dat_subplots %>% 
  filter(year == "2025") %>% 
  left_join(cvx_stem_density25, by = c("plot" = "cluster_2025" )) %>% 
  select(-plot, -cluster_2023, -common_cluster_ID ) %>% 
  rename(plot = plot_comb) %>% 
  dplyr::select(plot, subplot, species,  vegtype,    hgt,     n,   
         dbh,   year, 
         clear,     grndwrk,  logging_trail,  planting,
         anti_browsing,  
         n_subplots,status, disturbance_year, 
         forest_year, disturbance_length,
         pre_dist_trees_n,  
         area_m2, pre_dist_dens_ha
         ) %>% 
  mutate(time_snc_full_disturbance = 2025 - disturbance_year,
         time_snc_part_disturbance = 2025 - forest_year + 1)

names(dat_subplots23)
names(dat_subplots25)


# merge data with cleaned up naming
dat_subplots_merged <- rbind(dat_subplots23, dat_subplots25) %>% 
  mutate(plot    = as.factor(plot),
         subplot = as.factor(subplot)) 

## Recode field data -----------------------------------
# guestimate dbh and ba per individual based on height distribution 
dat_subplot_recode <- dat_subplots_merged %>% 
  # remove if no species is defined
  #filter(!is.na(species)) %>% # yes, I can safely exclude those, they are not representing 'empty plots' 
  # adjust tree species naming
  dplyr::mutate(
    species = str_trim(species),
    species = dplyr::case_when(
      species == "abal" ~ "absp",
      species %in% c("quro", "quru") ~ "qusp",
      TRUE ~ species
    )
  ) %>% 
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
  ) %>% 
  # Calculate basal area
  mutate(
    basal_area_cm2 = pi * (dbh_est / 2)^2,
    ba_total_cm2   = basal_area_cm2 * n,
    ba_total_m2    = ba_total_cm2 / 10000,
    #ba_ha_m2       = ba_total_m2 * scaling_factor
  ) %>% 
  left_join(traits_full)  # add trait table to identify early vs late seral



## Clean up management type ------------------
dat_subplot_mng <- dat_subplot_recode %>%
  mutate(across(all_of(management_types_v), ~ ifelse(is.na(.), 0, .)))

# --- 2) Subplot-level scores (use max within subplot to avoid duplicates)
mng_subplot_scores <- dat_subplot_mng %>%
  group_by(year, plot, subplot) %>%
  summarise(
    clear          = max(clear),
    grndwrk        = max(grndwrk),
    logging_trail  = max(logging_trail),
    planting       = max(planting),
    anti_browsing  = max(anti_browsing),
    .groups = "drop"
  ) %>%
  mutate(
    salvage_sub    = clear + grndwrk + logging_trail,          # 0..3 - 
    protection_sub = planting + anti_browsing,                  # 0..2
    management_sub = salvage_sub + protection_sub               # 0..5
  )

# --- 3) Plot-level intensities (scaled 0–1)
# continue here to get proper intensities!!! note that i have duplicated values in dat-overlap
mng_plot_intensity <- mng_subplot_scores %>%
  group_by(plot) %>%
  summarise(
    n_subplots           = n_distinct(subplot),
    salvage_sum          = sum(salvage_sub),
    protection_sum       = sum(protection_sub),
    
    # get sums
    clear_sum            = sum(clear),
    grndwrk_sum          = sum(grndwrk),
    logging_trail_sum    = sum(logging_trail),
    planting_sum         = sum(planting),
    anti_browsing_sum    = sum(anti_browsing),
    
    management_sum       = sum(management_sub),
    
    clear_intensity           = clear_sum / n_subplots,
    grndwrk_intensity         = grndwrk_sum / n_subplots,
    logging_trail_intensity   = logging_trail_sum / n_subplots,
    planting_intensity        = planting_sum / n_subplots,
    anti_browsing_intensity   = anti_browsing_sum / n_subplots,
    
    salvage_intensity    = salvage_sum    / (3 * n_subplots), # 3 types
    protection_intensity = protection_sum / (2 * n_subplots), # 2 types
    management_intensity = management_sum / (5 * n_subplots), # 5 types
    .groups = "drop"
  ) %>%
  # make sure the ranges are teh same
  mutate(
    clear_intensity           = pmin(pmax(clear_intensity, 0), 1),
    grndwrk_intensity         = pmin(pmax(grndwrk_intensity, 0), 1),
    logging_trail_intensity   = pmin(pmax(logging_trail_intensity, 0), 1),
    planting_intensity        = pmin(pmax(planting_intensity, 0), 1),
    anti_browsing_intensity   = pmin(pmax(anti_browsing_intensity, 0), 1),
   
    salvage_intensity     = pmin(pmax(salvage_intensity, 0), 1),
    protection_intensity  = pmin(pmax(protection_intensity, 0), 1),
    management_intensity  = pmin(pmax(management_intensity, 0), 1)
  ) %>% 
  select(-n_subplots)


dat_subplot_mng2 <- dat_subplot_mng %>% 
  left_join(mng_plot_intensity, by = c('plot'))

# keep only overlapping plots
dat_overlap <- dat_subplot_mng2 %>% 
  filter(status == 'both')
# filter(status == "both") %>% # keep only overlapping sites
length(unique(dat_overlap$plot))

hist(traits_full$Shade_tolerance)


### export important tables ---------------

fwrite(dat_subplot_mng2, 'outData/full_table_23_25.csv')
fwrite(dat_overlap, 'outData/full_table_overlap_23_25.csv')





