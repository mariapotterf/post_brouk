

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



# Read files: --------------
# --- Field data: 2023 
dat23_subplot    <- data.table::fread("outData/subplot_full_2023.csv")   # subplot-level table
dat23_sf         <- sf::st_read("outData/sf_context_2023.gpkg")          # subplot spatial data

# Select and rename
dat23_sf_min <- dat23_sf %>%
  dplyr::select(subplot = ID, plot = cluster)

# read data from 2025
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

subplot_to_plot <- dat25_subplot_sub %>%
  distinct(subplot, plot)

# Expand dataset
dat25_expanded <-
  dat25_subplot_sub %>%
  as_tibble() %>%
    dplyr::select(-plot) %>% 
  tidyr::complete(
    subplot,
    vegtype = c("small", "advanced", "mature"),
    species = unique(dat25_subplot_sub$species),
    fill = list(
      n = NA_integer_,
      dbh = NA_character_,
      hgt = NA_character_
    )
  ) %>%
  # Reattach plot info
  left_join(subplot_to_plot, by = "subplot") %>%
  mutate(year = "2025") %>%
  # Reorder columns if needed
  dplyr::select(all_of(target_cols))


#length(unique(dat25_subplot_sub$species))
# 375_T2_AH_20250827


### Field data 2023: subplot level --------------------------------------
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

# bind field data from both years:
dat_subplots <- bind_rows(dat23_subplot_sub, dat25_subplot_sub)

#  remove the erroneous subplots, if i have 6 subplots in cluster: 

# plots to check :
  # 15_102 ???  - kept 6
  # 15_124
  # 15_145
  # 184
  # 26_101
  # 26_142
  # 48
  
# subplots: 
# 424_T4_TP_20250827 - to remove, empty
# 539_T4_TP_20250827 - removed, empty
# 375_T2_AH_20250827 - to remove, empty
# 744_T4_TP_20250827
# 13_15_102_5
# 13_15_124_2
# 13_15_145_4
# 13_26_101_2
# 13_26_142_1


bad_subplots <- c("424_T4_TP_20250827","539_T4_TP_20250827","375_T2_AH_20250827",
                  "744_T4_TP_20250827","13_15_102_5","13_15_124_2",
                  "13_15_145_4","13_26_101_2","13_26_142_1",
                  "506_T4_TP_20250827") # has missing tree species

dat_subplots <- dat_subplots[!(subplot %in% bad_subplots)]

# remove whole plots (robust to numeric/character mix)
bad_plots <- c("15_104", "143", "26_134")

dat_subplots <- dat_subplots %>%
  dplyr::filter(!as.character(plot) %in% bad_plots)

n_subplots <- dat_subplots %>%
  #filter(!is.na(plot), !is.na(subplot)) %>%
  distinct(plot, subplot) %>%       # drop species/vegtype duplicates
  count(plot, name = "n_subplots") #%>%
#arrange(plot)

dat_subplots <- dat_subplots %>% 
  left_join(n_subplots) %>% 
  dplyr::filter(n_subplots == 5)


# --- Tree-based history (vector layers)
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

## Clean up convex hull characteristics -------------------------
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

# Add subplot (square)/plot (convex hull) info to each tree  (spatial join)
pre_trees_cvx_joined <- terra::intersect(pre_trees_3035_clean, convex_hull_3035_clean)

### Get pre-disturbnace stem density per plot 
cvx_df <- as.data.frame(pre_trees_cvx_joined)

####  Plot = CVX: stem density ------------------------------
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

### Clean up field data -----------------------------------
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
  ) 


# define recovery type: based on species that are likely planted/pioneer
dat_subplot_recode <- dat_subplot_recode %>%
  mutate(recovery_type = case_when(
    # --- Planted / late-successional / non-native ---
    species %in% c("piab","pisy","absp","lade","psme","taba",
                   "fasy","qusp","acca","acpl","acps","frex",
                   "casa","aehi","saca","rops") ~ "late",
    
    # --- Pioneer / early successional ---
    species %in% c("besp","alin","algl","alvi","potr","posp","prav",
                   "tisp","soau","soto","soar","cabe","ulsp","aial",
                   "fror","juni","jure","qusp","sasp","osca") ~ "early",
    
    # --- Everything else / not clearly one of the two ---
    TRUE ~ "other"
  ))


# clean up management type ------------------
dat_subplot_mng <- dat_subplot_recode %>%
  mutate(across(all_of(management_types_v), ~ ifelse(is.na(.), 0, .)))

# --- 2) Subplot-level scores (use max within subplot to avoid duplicates)
mng_subplot_scores <- dat_subplot_mng %>%
  group_by(plot, subplot) %>%
  summarise(
    clear          = max(clear),
    grndwrk        = max(grndwrk),
    logging_trail  = max(logging_trail),
    planting       = max(planting),
    anti_browsing  = max(anti_browsing),
    .groups = "drop"
  ) %>%
  mutate(
    salvage_sub    = clear + grndwrk + logging_trail,          # 0..3
    protection_sub = planting + anti_browsing,                  # 0..2
    management_sub = salvage_sub + protection_sub               # 0..5
  )

# --- 3) Plot-level intensities (scaled 0–1)
mng_plot_intensity <- mng_subplot_scores %>%
  group_by(plot) %>%
  summarise(
    n_subplots          = n_distinct(subplot),
    salvage_sum         = sum(salvage_sub),
    protection_sum      = sum(protection_sub),
    management_sum      = sum(management_sub),
    salvage_intensity   = salvage_sum    / (3 * n_subplots),
    protection_intensity= protection_sum / (2 * n_subplots),
    management_intensity= management_sum / (5 * n_subplots),
    .groups = "drop"
  ) %>%
  mutate(
    salvage_intensity    = pmin(pmax(salvage_intensity, 0), 1),
    protection_intensity = pmin(pmax(protection_intensity, 0), 1),
    management_intensity = pmin(pmax(management_intensity, 0), 1)
  )


dat_subplot_mng <- dat_subplot_mng %>% 
  left_join(mng_plot_intensity, by = c('plot', 'n_subplots'))

fwrite(dat_subplot_mng, 'outData/full_table_23_25.csv')

# keep only overlapping plots
dat_overlap <- dat_subplot_mng %>% 
  filter(status == 'both')
# filter(status == "both") %>% # keep only overlapping sites
length(unique(dat_overlap$plot))


# get master table, having all unique plots and subplots
dat_master_subplot <- dat_overlap %>% 
  dplyr::select(plot, subplot, year) %>% 
  distinct()

table(dat_overlap$year)  

n_plots_total <- length(unique(dat_master_subplot$plot))     # 126
n_plots_total

n_subplots_total <-length(unique(dat_master_subplot$subplot))  # 1250
n_subplots_total

# histogram of stem denisty per vertcal class:
dat_overlap %>% 
  dplyr::filter(n > 0) %>% 
  ggplot(aes(n , fill = year)) +
  geom_histogram() + 
  facet_grid(year~vegtype, scales = 'free')


# get disturbance characteristics on plot level
plot_disturb_chars <- dat_overlap %>% 
  dplyr::select(plot, year, disturbance_year, 
                forest_year, disturbance_length, 
                time_snc_full_disturbance, 
                time_snc_part_disturbance #,
                #clear, grndwrk, logging_trail, planting, anti_browsing
                )  %>%  
  distinct()

plot_disturb_chars %>% 
  ggplot(aes(time_snc_full_disturbance)) + 
  geom_histogram() + 
  facet_grid(.~year)


# make a master table having all (even empty subplots and plots)
df_master_overlap <- dat_overlap %>% 
  distinct(plot, subplot)




df_master_mng <- dat_overlap %>% 
  dplyr::filter(year == "2023") %>% # keep management oionly frm 2023 for consistency
  distinct(plot, subplot,year,
           clear,
           grndwrk,
           logging_trail,
           planting,
           anti_browsing)

prop.table(table(df_master_mng$clear ))
prop.table(table(df_master_mng$grndwrk))
prop.table(table(df_master_mng$logging_trail))
prop.table(table(df_master_mng$planting))
prop.table(table(df_master_mng$anti_browsing))

# > table(df_master_mng$clear )
# 0   1 
# 7 618 
# > table(df_master_mng$grndwrk)
# 0   1 
# 79 546 

# > table(df_master_mng$logging_trail)
# 0   1 
# 522 103 

# > table(df_master_mng$planting)
# 0   1 
# 297 328 

# > table(df_master_mng$anti_browsing)
# 0   1 
# 377 248

# get management characteristics only from 2023
#   clear,
# grndwrk,
# logging_trail,
# planting,
# anti_browsing


nrow(mng_context)


# these are doibles because i have records for 2023 and 2025
hist(mng_context$salvage_intensity)
hist(mng_context$protection_intensity)
hist(mng_context$management_intensity)




## get share between pioneers vs planted on landscape level ---------------
# Summarize total number of trees per year and recovery type
total_per_year <- dat_overlap %>%
  group_by(year) %>%
  summarise(
    total_trees = sum(n, na.rm = TRUE),
    .groups = "drop"
  )

tree_summary <- dat_overlap %>%
  group_by(year, recovery_type) %>%
  summarise(
    n_trees_recovery = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(total_per_year) %>% 
  mutate(share = n_trees_recovery /total_trees * 100)

print(tree_summary)

# Stacked barplot
tree_summary %>% 
  dplyr::filter(recovery_type != 'other') %>% 
  ggplot( aes(x = year, y = share, fill = recovery_type)) +
  geom_bar(stat = "identity", color = "black") +
  labs(
    title = "Tree Counts by Recovery Type per Year",
    x = "Year",
    y = "Share [%]",
    fill = "Recovery Type"
  ) +
  scale_fill_manual(values = c("early" = "#66c2a5", "late" = "#fc8d62")) +
  theme_classic2(base_size = 8)




dat_sum_recovery <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>% 
  # mutate(year = as.factor(year)) %>% 
  group_by(recovery_type, year) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE))

dat_sum_recovery


# chnage of shares betwen early vs. late species?

dat_overlap %>% 
  group_by(plot, recovery_type) %>%
  summarise(n_stems = n(), .groups = "drop") %>%
  pivot_wider(names_from = recovery_type,
              values_from = n_stems,
              values_fill = 0)  # fill missing types with 0


# check up development of early vs ;late shares givet time since disturbance

# Calculate stem counts by recovery type at the plot level
share_early_vs_late <- 
  dat_overlap %>%
  filter(!is.na(n)) %>%  # Optional: remove NAs if present
  group_by(plot, recovery_type, time_snc_full_disturbance) %>%
  summarise(n_stems = n(), .groups = "drop") %>%
  pivot_wider(names_from = recovery_type,
              values_from = n_stems,
              values_fill = 0) %>% 
  mutate(total = early + late,
         share_early = early/total*100,
         share_late = late/total*100) %>% 
  select(plot, time_snc_full_disturbance, share_early, share_late) %>%
  pivot_longer(cols = starts_with("share_"),
               names_to = "recovery_type",
               values_to = "share")# %>%

### Early vs late: plot level --------------------------
df_plot_share_early <-  
  dat_overlap %>%
  filter(!is.na(n)) %>%  # Optional: remove NAs if present
  group_by(plot,year, recovery_type) %>%
  summarise(n_stems = n(), .groups = "drop") %>%
  pivot_wider(names_from = recovery_type,
              values_from = n_stems,
              values_fill = 0) %>% 
  mutate(total = early + late,
         share_early = early/total*100,
         share_late = late/total*100) %>% 
  select(plot, year, share_early, share_late) 

### early vs late : subplot level ---------------------
df_sub_share_early <-  
  dat_overlap %>%
  filter(!is.na(n)) %>%  # Optional: remove NAs if present
  group_by(subplot, plot,year, recovery_type) %>%
  summarise(n_stems = n(), .groups = "drop") %>%
  pivot_wider(names_from = recovery_type,
              values_from = n_stems,
              values_fill = 0) %>% 
  mutate(total = early + late,
         share_early = early/total*100,
         share_late = late/total*100) %>% 
  select(subplot, plot, year, share_early, share_late) 

   
# Plot
ggplot(share_early_vs_late, 
       aes(x = factor(time_snc_full_disturbance),
           y = share,
           fill = recovery_type)) +
  geom_boxplot() +
  #geom_jitter() +
  #geom_bar(stat = "identity", position = "stack") +
  labs(x = "Time since full disturbance (years)",
       y = "Share of stems (%)",
       fill = "Recovery type") 



# --- Field data: Subplot-level metrics 
field_sub_summ <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>% 
 # mutate(year = as.factor(year)) %>% 
  group_by(plot, subplot, year, 
           time_snc_full_disturbance, time_snc_part_disturbance,disturbance_year, 
           forest_year, disturbance_length,
           clear,
           grndwrk,
           logging_trail,
           planting,
           anti_browsing,
           management_intensity) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE),
    sp_richness = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp  = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
    # weighted mean height using only present stems
    mean_hgt    = if (stems_total > 0) weighted.mean(hgt_est[n > 0], n[n > 0]) else NA_real_,
    
    # compute weighted variance whenever we have ≥2 stems and any finite heights
    var_hgt = {
      if (stems_total > 1) {
        sel <- (n > 0) & is.finite(hgt_est)
        h   <- hgt_est[sel]; ww <- n[sel]
        if (length(h) >= 1 && sum(ww) > 1) {
          mu <- weighted.mean(h, ww)
          v  <- sum(ww * (h - mu)^2) / (sum(ww) - 1)  # freq-weighted, Bessel corrected
          if (is.nan(v)) NA_real_ else v              # will be 0 if all h are equal
        } else NA_real_
      } else NA_real_
    },
    
    cv_hgt = if (is.finite(mean_hgt) && mean_hgt > 0 && !is.na(var_hgt))
      sqrt(var_hgt) / mean_hgt else
        if (stems_total > 1 && is.finite(mean_hgt) && mean_hgt > 0) 0 else NA_real_,
    .groups = "drop"  #,
    #range_hgt = if (sum(n > 0) >= 2) diff(range(hgt_est[n > 0], na.rm = TRUE)) else NA_real_,
   
  )  %>% 
  left_join(cwm_subplot)
  
 # mutate(cv_hgt = ifelse(is.na(cv_hgt), 0L, cv_hgt),
  #       mean_hgt = ifelse(is.na(mean_hgt), 0L, mean_hgt)) # replace NA by 0 if stems are missing


# get stem densiity by management ---------------------
df_long <- field_sub_summ %>%
  dplyr::filter(stems_total > 0) %>% 
  filter(cv_hgt >0) %>% 
  #ungroup() %>%
  select(year,
         clear,
         grndwrk,
         logging_trail,
         planting,
         anti_browsing,
         time_snc_full_disturbance, time_snc_part_disturbance, 
         management_intensity,
         mean_hgt, cv_hgt, shannon_sp, sp_richness,
         CWM_shade ,
         CWM_drought ) %>%
  pivot_longer(-c(year,
                    time_snc_full_disturbance,
                    time_snc_part_disturbance,
                  CWM_shade ,
                  CWM_drought,
                    clear,
                    grndwrk,
                    logging_trail,
                    planting,
                    anti_browsing,
                  management_intensity),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value))# %>%
  #filter(cv_hgt >0)
  # mutate(metric = recode(metric,
  #                        mean_hgt    = "Mean height (m)",
  #                        cv_hgt      = "CV of height",
  #                        shannon_sp  = "Shannon (H')",
  #                        sp_richness = "Species richness"
  # )) #%>%



## Time since disturbnace ----------------------------------------
df_long %>% 
  ggplot(aes(x = clear, y = value)) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, , title = "Clearing, subplot") +
  theme_classic2()

df_long %>% 
  ggplot(aes(x = planting, y = value)) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = "Planting, subplot") +
  theme_classic2()





# see CV with time since disturbnace : poartial disturbance
p_partial_disturbance <- df_long %>% 
  ggplot(aes(x = time_snc_part_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_part_disturbance ), outlier.shape = NA) +
  
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Partial disturbance') +
  theme_classic2()


# see CV with time since disturbnace : full disturbance
p_full_disturbance <- field_sub_summ %>%
  #ungroup() %>%
  select(year, time_snc_full_disturbance, time_snc_part_disturbance, mean_hgt, cv_hgt, shannon_sp, sp_richness) %>%
  pivot_longer(c(-year, 
                 - time_snc_full_disturbance, 
                 - time_snc_part_disturbance),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Full disturbance') +
  theme_classic2()


# traits analysis: subplot  -------------------
field_sub_summ %>%
  #ungroup() %>%
  select(year, time_snc_full_disturbance, time_snc_part_disturbance, 
         mean_hgt, cv_hgt, shannon_sp, sp_richness,  CWM_shade ,
         CWM_drought) %>%
  pivot_longer(c(-year, 
                 - time_snc_full_disturbance, 
                 - time_snc_part_disturbance),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Full disturbance') +
  theme_classic2()



df_long %>% 
  ggplot(aes(x = time_snc_part_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_part_disturbance ), outlier.shape = NA) +
  
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Partial disturbance') +
  theme_classic2()






subplot_summary_tbl <- field_sub_summ %>%
  ungroup() %>%
  select(year, mean_hgt, cv_hgt, shannon_sp, sp_richness) %>%
  pivot_longer(-year, names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  group_by(metric, year) %>%
  summarise(
    n_total  = n(),
    n_non_na = sum(!is.na(value)),
    n_na     = sum(is.na(value)),
    mean     = mean(value, na.rm = TRUE),
    sd       = sd(value, na.rm = TRUE),
    se       = sd/sqrt(n_non_na),
    median   = median(value, na.rm = TRUE),
    p25      = quantile(value, 0.25, na.rm = TRUE),
    p75      = quantile(value, 0.75, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(metric, year)

subplot_summary_tbl


# summarize subplot information 
# how many trees?
sum(field_sub_summ$stems_total)            # 3658 stems
length(unique(field_sub_summ$subplot))     # 1250
sum(field_sub_summ$stems_total == 0)       # 271 
sum(field_sub_summ$cv_hgt > 1, na.rm = T)  # 20


field_sub_summ_filt <- field_sub_summ %>%
  filter(stems_total > 0 & cv_hgt > 0)

p1 <- ggplot(field_sub_summ_filt, aes(x = stems_total, y = cv_hgt, color = year,
                                      fill = year)) +
  geom_point(alpha = 0.4, size = 2 ) +
 # geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3),
              se = TRUE,  linetype = "dashed") +
  labs(
    x = "Stem count per subplot",
    y = "CV of tree height",
    title = "Stem density vs. Vertical Structural Variation",
    subtitle = "Empty subplots excluded"
  ) +
  theme_classic2(base_size = 8)


p2 <- ggplot(field_sub_summ_filt, aes(x = shannon_sp, y = cv_hgt, color = year,
                                      fill = year)) +
  geom_point(alpha = 0.4, size = 2) +
 # geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3),
              se = TRUE, linetype = "dashed") +
  labs(
    x = "Shannon per subplot",
    y = "CV of tree height",
    title = "Shannon vs. Vertical Structural Variation",
    subtitle = "Empty subplots excluded"
  ) +
  theme_classic2(base_size = 8)

ggarrange(p1, p2)
# ssame analysis on both scales? 

## get summary acroos all trees and study sites --------------------------
# find species with the highest share of stems overall
# 0) Safe counts (treat NA counts as 0)
df <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>%
  filter(!is.na(species) & species != "")

# 1) Totals
overall_totals <- df %>%
  summarise(total_trees = sum(n),
            n_plots = n_distinct(plot),
            n_subplots = n_distinct(subplot))

by_year_totals <- df %>%
  group_by(year) %>%
  summarise(total_trees = sum(n),
            #n_plots = n_distinct(plot),
           # n_subplots = n_distinct(subplot),
            .groups = "drop")

n_trees23 <- by_year_totals %>% 
  filter(year == 2023) %>% 
  pull()

n_trees25 <- by_year_totals %>% 
  filter(year == 2025) %>% 
  pull()

# 2) Species × year counts and presence
species_stem_share_year <- 
  df %>%
  group_by(year, species) %>%
  summarise(
    stems = sum(n),                                   # number of trees
    #plots_present = n_distinct(plot[n > 0]),          # plots where species occurs
    .groups = "drop"
  ) %>% 
  tidyr::pivot_wider(
    names_from  = year,
    values_from = c(stems),
    names_glue  = "{.value}_{year}",
    values_fill = list(stems = 0L)
  ) %>%
  dplyr::arrange(species) %>% 
  mutate(trees23 = n_trees23,
         trees25 = n_trees25,
         share23 = round(stems_2023/trees23*100,2),
         share25 = round(stems_2025/trees25*100,2),
         total_stems = stems_2023 + stems_2025,
         total_trees = n_trees23 + n_trees25,
         total_share = round(total_stems/total_trees * 100,2)) 

# merge counst across 2 years:     
top_overall_stem_share <- 
  species_stem_share_year %>%
  dplyr::select(species, total_stems, total_share) %>%
  dplyr::slice_max(order_by = total_share, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup() 

# get top 10 per each year (the order changes a bit)
top_stems_by_year <- species_stem_share_year %>%
  dplyr::select(species, share23, share25) %>%
  tidyr::pivot_longer(
    dplyr::starts_with("share"),
    names_to = "year",
    names_prefix = "share",
    values_to = "share"
  ) %>%
  dplyr::mutate(year = factor(paste0("20", year), levels = c("2023", "2025"))) %>%
  dplyr::group_by(year) %>%
  dplyr::slice_max(order_by = share, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

top_2023 <- top_stems_by_year %>% dplyr::filter(year == "2023") %>% dplyr::pull(species)
top_2025 <- top_stems_by_year %>% dplyr::filter(year == "2025") %>% dplyr::pull(species)

length(union(top_2023, top_2025))      # 11 (your result)
length(intersect(top_2023, top_2025))  # 9  (overlap size)

setdiff(top_2025, top_2023)  # species only in 2025's top10
setdiff(top_2023, top_2025)  # species only in 2023's top10

unique(top_stems_by_year$species)

v_top_species_overall <- top_overall_stem_share %>% pull(species)


# Reverse the color palette and map to the species in the desired order
n_colors  <- length(v_top_species_overall)  # Number of species
my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n_colors)  # Generate colors
# make / use a palette of exactly the needed length
pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(length(v_top_species_overall))
pal <- rev(pal)  # start with dark green

# map colors to species automatically
species_colors <- setNames(pal, v_top_species_overall)
species_colors
#

# species_colors
# piab      besp      pisy      qusp      fasy      lade      saca      soau      acps      potr      absp      sasp 
# "#006837" "#17934D" "#58B65F" "#94D168" "#C6E77F" "#EDF7A7" "#FEF0A7" "#FDCD7B" "#FA9C58" "#EE613D" "#D22B26" "#A50026" 


# Print the color assignments for confirmation
print(species_colors)


# update species labels
species_labels <- c(
  piab = "Picea abies",
  besp = "Betula sp.",
  pisy = "Pinus sylvestris",
  qusp = "Quercus sp.",
  fasy = "Fagus sylvatica",
  lade = "Larix decidua",
  saca = "Salix caprea",
  soau = "Sorbus aucuparia",
  acps = "Acer pseudoplatanus",
  potr = "Populus tremula",
  absp = "Abies sp.",
  sasp = "Salix sp."
)


# Identify plots without any stems present ---------------------------------------------------------------------------------

# Step 1: Summarise total stems per plot and year
plot_year_summ <- df %>%
  group_by(plot, year) %>%
  summarise(total_stems = sum(n, na.rm = TRUE), .groups = "drop")

# Step 2: Reshape to wide format
plot_year_wide <- plot_year_summ %>%
  tidyr::pivot_wider(names_from = year, values_from = total_stems, values_fill = 0)

# Step 3: Filter based on presence in each year
plots_empty23    <- plot_year_wide %>% filter(`2023` == 0) %>% pull(plot) # 6 ~ 4.7%
plots_empty25    <- plot_year_wide %>% filter(`2025` == 0) %>% pull(plot) # 2 ~ 1.5%
plots_empty_both <- plot_year_wide %>% filter(`2023` == 0 & `2025` == 0)  # zero

plots_empty23
plots_empty25



## get average stem density per species per top 10 species --------------------------------
df_stem_dens_species <- df %>% 
  group_by(plot, species, year, n_subplots ) %>%
  summarize(sum_n = sum(n, na.rm =T)) %>% 
  mutate(scaling_factor = 10000/(n_subplots * 4),
         stem_dens = sum_n*scaling_factor) %>% 
  mutate(log_sum_stem_density = log10(stem_dens + 1)) #%>%  # Adding 1 to avoid log(0)
  #ungroup()

# get total sum and calculate as average value over all sites 
df_stem_dens_species_sum <- 
  df_stem_dens_species %>% 
  group_by(species, year) %>% 
  summarise(stem_dens = sum(stem_dens, na.rm = T),
            log_sum_stem_density = sum(log_sum_stem_density, na.rm = T)) %>%
    mutate(stem_dens_avg = stem_dens/n_plots_total,
           log_sum_stem_density_avg = log_sum_stem_density/n_plots_total)
  

df_stem_dens_species_year <- df_stem_dens_species %>% 
  ungroup(.) %>% 
  filter(sum_n >0) %>% 
  filter(species %in% v_top_species_overall) %>% 
  dplyr::group_by(species, year) %>%
  dplyr::mutate(median_stem_density = median(stem_dens, na.rm = TRUE)) %>% 
  dplyr::ungroup(.) %>%
  mutate(species = factor(species, levels = rev(v_top_species_overall))) # Set custom order



# Boxplot for stem density -------------
df_stem_dens_species_year2 <- df_stem_dens_species_year %>%
  dplyr::filter(!is.na(log_sum_stem_density) & sum_n > 0) %>%
  dplyr::mutate(year = factor(year, levels = c("2023","2025")),
                # order by mean log density (ascending → highest ends up at the TOP after coord_flip)
                species = forcats::fct_reorder(species, log_sum_stem_density, .fun = mean, na.rm = TRUE))

# 
p_density<-df_stem_dens_species_year2 %>% 
  filter(!is.na(species)) %>% 
  ggplot(aes(x = log_sum_stem_density, y = species,
             fill = year)) +
  geom_boxplot(
    #aes(group = interaction(species, year), 
     #   alpha = factor(year)),
    position = position_dodge(width = 0.6),
    outlier.shape = NA,
    width = 0.45#,
   # color = "black"
  ) +
  
 # coord_flip() +
  labs(
    x = "log(sum stem density)",
    y = "",  
    fill = "Year"
  ) +
  # scale_fill_manual(values= species_colors) +
  # scale_alpha_manual(
  #   values = c("2023" = 0.5, "2025" = 1)#,  # 2023 = lighter
  #   #guide = "none"  # hides alpha legend
  # ) +
  theme_classic(base_size = 10) +
  scale_y_discrete(labels = species_labels) +
  theme(axis.text.y = element_text(face = "italic", size = 8))


p_density


# barplot of change over years

# Order species by their maximum share across years (nice stable ordering)
order_levels <- top_stems_by_year %>%
  group_by(species) %>%
  summarise(max_share = max(share, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_share)) %>%
  pull(species)

top_stems_by_year <- top_stems_by_year %>%
  mutate(
    species = factor(species, levels = order_levels),
    year    = factor(year, levels = c("2023","2025"))
  )

# master mapping: code -> Latin name
species_labels_all <- c(
  piab = "Picea abies",
  besp = "Betula sp.",
  pisy = "Pinus sylvestris",
  qusp = "Quercus sp.",
  fasy = "Fagus sylvatica",
  saca = "Salix caprea",
  lade = "Larix decidua",
  soau = "Sorbus aucuparia",
  acps = "Acer pseudoplatanus",
  potr = "Populus tremula",
  absp = "Abies sp.",
  sasp = "Salix sp.",
  cabe = "Carpinus betulus"
)

p_bar <- ggplot(top_stems_by_year, aes(x = share, y = species, fill = year)) +
  geom_col(aes(#group = interaction(species, year)#,
               #alpha = factor(year)
               ), 
           position = position_dodge(width = 0.7), width = 0.6) +
  # if 'share' is 0–100 already, just add a % suffix:
  scale_x_continuous(labels = label_number(accuracy = 0.1, suffix = "")) +
  scale_y_discrete(
    limits = rev(names(species_labels)),
    labels = species_labels,
    drop = FALSE
  ) +
  # scale_alpha_manual(
  #   values = c("2023" = 0.5, "2025" = 1.0)#,  # 2023 = lighter
  #   #guide = "none"  # hides alpha legend
  # ) +
  # if you prefer proportions (0–1), use:
  # scale_x_continuous(labels = label_percent()) 
  labs(
    x = "Share of stems",
    y = "Species",
    fill = "Year"#,
    #title = "Top 12 species by share, by year"
  ) +
 # scale_fill_manual(values = species_colors) +
  theme_classic2(base_size = 10) +
  theme(axis.text.y = element_text(size = 8, face = "italic"))

p_bar

# Get species occurence from total number of plots 
# Total number of unique plots
total_plots <- df_stem_dens_species %>%
  pull(plot) %>%
  n_distinct()

# Share of plots per species (where species has non-zero stems)
species_occurence <- 
  df_stem_dens_species %>%
  ungroup(.) %>% 
  dplyr::filter(sum_n > 0) %>%                 # Only where species occurred
  distinct(year, species, plot) %>%           # Unique species × plot combos
  count(year, species, name = "n_plots") %>%  # Count number of plots per species
  mutate(share_of_plots = n_plots / total_plots*100) %>% 
  arrange()

species_occurence

# Optional: order species by max share across years
species_order <- species_occurence %>%
  group_by(species) %>%
  summarise(max_share = max(share_of_plots)) %>%
  arrange(desc(max_share)) %>%
  pull(species)

species_plot_share <- species_occurence %>%
  mutate(species = factor(species, levels = rev(species_order)))

# Plot
p_occurence <- species_plot_share %>% 
  filter(species %in% v_top_species_overall ) %>% 
  ggplot(aes(x = share_of_plots, y = species, fill = year)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Species occurence over plots (%)",
    y = "Species",
    fill = "Year"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "right"
  )

p_occurence



ggarrange(p_bar, p_occurence, p_density,  
          ncol = 3, common.legend = T)


# Get functional traits database ---------------------------------------

# tolerance scales range from 0 (no tolerance) to 5 (maximal tolerance)

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

#  Input
niin_path <- "raw/traits_database/Niinemets_2006.xls"  # adjust if needed

# Your species codes from dat_subplot_mng$species
sp_codes <- unique(dat_overlap$species)
sp_codes <- sp_codes[!sp_codes %in% c("", "ots1")]

# Read Niinemets + keep essential columns 
eco_traits_raw <- read_excel(
  path  = niin_path,
  sheet = "Niinemets_2006_appendix",
  skip  = 3,
  .name_repair = function(x) gsub("\\s+", "_", x)
)

# filter out only teh important columns
eco_traits <- eco_traits_raw %>%
  select(Species, Shade_tolerance, Drought_tolerance) %>%
  mutate(Species = str_squish(Species))


# Map my acronyms to latin ones, indicate reason fror merge
# 'species_latin' must either be an exact species in Niinemets OR a genus label
# for which we will compute a genus mean (e.g., "Quercus", "Salix", "Betula").
sp_map <- tribble(
  ~species_code, ~species_latin,            ~notes,
  "piab",        "Picea abies",             "",
  "pisy",        "Pinus sylvestris",        "",
  "lade",        "Larix decidua",           "",
  "absp",        "Abies alba",              "Abies sp. -> proxy A. alba",
  "psme",        "Pseudotsuga menziesii",   "Fix spelling if needed (not 'mensiesii')",
  "taba",        "Taxus baccata",           "",
  "fasy",        "Fagus sylvatica",         "",
  "qusp",        "Quercus",                 "Quercus spp. -> genus mean (Q. robur & Q. petraea)",
  "acca",        "Acer campestre",          "",
  "acpl",        "Acer platanoides",        "",
  "acps",        "Acer pseudoplatanus",     "",
  "algl",        "Alnus glutinosa",         "",
  "alin",        "Alnus incana",            "",
  "alvi",        "Alnus viridis",           "Close to A. alnobetula complex",
  "potr",        "Populus tremula",         "",
  "posp",        "Populus tremula",         "Populus spp. -> use P. tremula proxy (keep consistent)",
  "besp",        "Betula",                  "Betula spp. -> genus mean (B. pendula & B. pubescens)",
  "frex",        "Fraxinus excelsior",      "",
  "tisp",        "Tilia cordata",           "Tilia spp. -> proxy T. cordata",
  "prav",        "Prunus avium",            "",
  "soau",        "Sorbus aucuparia",        "",
  "soto",        "Sorbus torminalis",       "",
  "soar",        "Sorbus aria",             "",
  "casa",        "Castanea sativa",         "",
  "aehi",        "Aesculus hippocastanum",  "",
  "cabe",        "Carpinus betulus",        "",
  "ulsp",        "Ulmus glabra",            "Ulmus spp. -> proxy U. glabra (switch to U. minor if that’s your system)",
  "rops",        "Robinia pseudoacacia",    "",
  "saca",        "Salix caprea",            "",
  "juni",        "Juniperus communis",      "",
  "jure",        "Juglans regia",           "",
  "sasp",        "Salix",                   "Salix spp. -> genus mean (S. caprea & S. alba)",
  "aial",        "Alnus alnobetula",        "Map to A. viridis in Niinemets if needed",
  "osca",        "Ostrya carpinifolia",     "If missing in Niinemets, proxy with Carpinus betulus",
  "fror",        "Fraxinus ornus",          ""
 
)

# get compleyte binomial latin names for my acronyms
df_filter_binomial <- sp_map %>% 
  filter(!species_latin %in% c("Quercus", "Salix", "Betula")) %>%  # filter out only Genus names
  select(-notes)

# filter complete records 
eco_traits_binomial <- eco_traits %>% 
  filter(Species %in% df_filter_binomial$species_latin) %>% 
  rename(species = Species)



#  3) Create genus means for Quercus / Salix / Betula 
# Helper function: mean traits for a set of species -> single row with genus name
genus_mean <- function(species_vec, genus_label) {
  eco_traits %>%
    filter(Species %in% species_vec) %>%
    summarise(
      species = genus_label,
      Shade_tolerance   = mean(Shade_tolerance, na.rm = TRUE),
      Drought_tolerance = mean(Drought_tolerance, na.rm = TRUE)#,
      #n_species_used    = dplyr::n()
    )
}

# Define which species to use per genus (adjust if your region differs)
q_species <- c("Quercus robur", "Quercus petraea")
s_species <- c("Salix caprea", "Salix alba")
b_species <- c("Betula pendula", "Betula pubescens") # if only one is present, mean() still works

traits_quercus <- genus_mean(q_species, "Quercus")
traits_salix   <- genus_mean(s_species, "Salix")
traits_betula  <- genus_mean(b_species, "Betula")

# Bind genus means
traits_full <- bind_rows(
  eco_traits_binomial,
  traits_quercus,
  traits_salix,
  traits_betula
) 


# add back acronyms
traits_full <- traits_full %>% 
  left_join(sp_map, by = c("species" = "species_latin")) %>% 
  select(-species, -notes) %>% 
  rename(species = species_code)


# add traits into full data
dat_overlap <-  dat_overlap %>% 
  left_join(traits_full)


dat_subplot_mng <-  dat_subplot_mng %>% 
  left_join(traits_full)


## get Weighted community mean per subplot and Plot ------------------

# choose  weighting : by stems (for early communities), by structure?
# Option A (default): weight by stem counts
wvar <- "n"

# Option B: weight by structure (uncomment ONE)
# wvar <- "basal_area_cm2"
# wvar <- "hgt_est"

# helper to pull a numeric weight safely
wfun <- function(x) ifelse(is.na(x) | x < 0, 0, x)

# Subplot × year CWMs 
cwm_subplot <- dat_overlap %>%
  filter(species != 'ots1') %>% 
  mutate(w = wfun(.data[[wvar]])) %>%
  # keep only rows that contribute weight and have trait scores
  filter(w > 0) %>%
  group_by(plot, subplot, year) %>%
  summarise(
    stems_with_traits = sum(w[!is.na(Shade_tolerance) & !is.na(Drought_tolerance)], na.rm = TRUE),
    stems_total       = sum(w, na.rm = TRUE),
    CWM_shade   = ifelse(stems_with_traits > 0,
                         sum(w * Shade_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    CWM_drought = ifelse(stems_with_traits > 0,
                         sum(w * Drought_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    trait_coverage = stems_with_traits / pmax(stems_total, 1e-9),
    .groups = "drop"
  )

# 2) Plot × year CWMs
cwm_plot <- dat_overlap %>%
  filter(species != 'ots1') %>% 
  mutate(w = wfun(.data[[wvar]])) %>%
  filter(w > 0) %>%
  group_by(plot, year) %>%
  summarise(
    stems_with_traits = sum(w[!is.na(Shade_tolerance) & !is.na(Drought_tolerance)], na.rm = TRUE),
    stems_total       = sum(w, na.rm = TRUE),
    CWM_shade   = ifelse(stems_with_traits > 0,
                         sum(w * Shade_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    CWM_drought = ifelse(stems_with_traits > 0,
                         sum(w * Drought_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    trait_coverage = stems_with_traits / pmax(stems_total, 1e-9),
    .groups = "drop"
  )

### quick plots on subplot and plot level 
ord_year <- function(x) factor(as.character(x), levels = c("2023","2025"))


# Long format for subplot
subplot_long <- cwm_subplot %>%
  transmute(level = "subplot",
            plot, subplot,
            year = ord_year(year),
            CWM_shade, CWM_drought) %>%
  pivot_longer(c(CWM_shade, CWM_drought),
               names_to = "trait", values_to = "CWM")

# Long format for plot
plot_long <- cwm_plot %>%
  transmute(level = "plot",
            plot,
            year = ord_year(year),
            CWM_shade, CWM_drought) %>%
  pivot_longer(c(CWM_shade, CWM_drought),
               names_to = "trait", values_to = "CWM")

# Combine
all_long <- bind_rows(subplot_long, plot_long) %>%
  filter(!is.na(CWM))

# Nice facet labels
trait_labs <- c(CWM_shade = "Shade tolerance", CWM_drought = "Drought tolerance")

ggplot(all_long, aes(x = year, y = CWM, fill = year)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.15) +
  geom_jitter(width = 0.1, alpha = 0.25, size = 0.7) +
  facet_grid(level ~ trait, labeller = labeller(trait = trait_labs)) +
  labs(x = "Year", y = "Community-weighted mean (CWM)",
       title = "Trait CWMs by level and year") +
  theme_classic2(base_size = 10)



# Summary stats on plot level stem density per m² --------------------
area_subplot_m2 <- 4      # 4 m²
area_plot_m2    <- 5*4    # 20 m²



# --- Plot-level metrics (aggregate over subplots) ---
plot_metrics_mean <- field_sub_summ %>%
  group_by(plot, year) %>%
  summarise(
    mean_sp_richness = mean(sp_richness, na.rm = TRUE),
    #var_sp_richness  = var(sp_richness,  na.rm = TRUE),
    mean_shannon_sp  = mean(shannon_sp,  na.rm = TRUE),
    mean_evenness_sp = mean(evenness_sp, na.rm = TRUE),
    mean_mean_hgt    = mean(mean_hgt,    na.rm = TRUE),
    mean_cv_hgt      = mean(cv_hgt,      na.rm = TRUE),
    mean_eff_numb    = mean(effective_numbers, na.rm = TRUE),    
    #var_cv_hgt       = var(cv_hgt,       na.rm = TRUE),
    #mean_range_hgt   = mean(range_hgt,   na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(cwm_plot)


# --- pooled CV directly from dat23_subplot_recode ---
plot_metrics_pooled  <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>%                # no NA counts
  group_by(plot, year) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE),
    sp_richness = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp  = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    mean_hgt    = if (stems_total > 0) weighted.mean(hgt_est, n, na.rm = TRUE) else NA_real_,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
    # compute weighted variance whenever we have ≥2 stems and any finite heights
    var_hgt = {
      if (stems_total > 1) {
        sel <- (n > 0) & is.finite(hgt_est)
        h   <- hgt_est[sel]; ww <- n[sel]
        if (length(h) >= 1 && sum(ww) > 1) {
          mu <- weighted.mean(h, ww)
          v  <- sum(ww * (h - mu)^2) / (sum(ww) - 1)  # freq-weighted, Bessel corrected
          if (is.nan(v)) NA_real_ else v              # will be 0 if all h are equal
        } else NA_real_
      } else NA_real_
    },
    
    cv_hgt = if (is.finite(mean_hgt) && mean_hgt > 0 && !is.na(var_hgt))
      sqrt(var_hgt) / mean_hgt else
        if (stems_total > 1 && is.finite(mean_hgt) && mean_hgt > 0) 0 else NA_real_,
    .groups = "drop" 
  ) %>%
  left_join(cwm_plot)
#mutate(cv_hgt = ifelse(is.na(cv_hgt), 0L, cv_hgt)) # replace NA by 0 if stems are missing

# get only context and disturbance information on plot level
df_plot_context <- dat_overlap %>%
  dplyr::select(plot,year,  
                pre_dist_trees_n,  area_m2, 
                pre_dist_dens_ha, 
                time_snc_full_disturbance,
                time_snc_part_disturbance,
                disturbance_year, forest_year, disturbance_length,
                protection_intensity, management_intensity) %>% 
  distinct() 
  



# create final table for both levels -----------------------------------------------------
# 1) Subplot table (has subplot mean_hgt and stems_total as weights)
sub_df <- field_sub_summ %>%
  filter(stems_total > 0) %>% #, cv_hgt > 0
  transmute(
    ID       = subplot,
    plot_id  = plot, #str_replace(subplot, "^[^_]+_([^_]+_[^_]+)_.*$", "\\1"),
    year     = year,
    level    = "subplot",
    dens_m2  = stems_total / area_subplot_m2,     # 4 m² subplot
    cv_hgt   = cv_hgt,
    mean_hgt = mean_hgt,
    sp_richness = sp_richness,
    shannon_sp  = shannon_sp ,
    evenness_sp = evenness_sp ,
    effective_numbers = effective_numbers,
    w        = stems_total
  ) %>% 
  mutate(dens_ha = dens_m2*10000) %>% 
  left_join(df_plot_context, by = c("plot_id" = "plot",
                                    "year" = "year")) %>% 
  left_join(cwm_subplot, by = c("plot_id" = "plot",
                                "year" = "year",
                                "ID" = 'subplot') )

# 2) Plot table (pooled metrics already computed)
plot_df <- plot_metrics_pooled %>%
  transmute(
    ID       = plot,
    plot_id  = plot,
    year     = year,
    level    = "plot",
    dens_m2  = stems_total / area_plot_m2,    # 5×4 m² = 20 m²
    cv_hgt   = cv_hgt,
    mean_hgt = mean_hgt,
    sp_richness = sp_richness,
    shannon_sp  = shannon_sp ,
    evenness_sp = evenness_sp ,
    effective_numbers = effective_numbers,
    w        = stems_total
  ) %>%
  #filter(!is.na(cv_hgt), cv_hgt > 0) %>% 
  mutate(dens_ha = dens_m2*10000,
         mean_hgt = replace_na(mean_hgt, 0)) %>% # Replace NA with 0)
  left_join(df_plot_context, by = c("plot_id" = "plot",
                                  "year" = "year")) %>% 
  left_join(cwm_plot, by = c("plot_id" = "plot",
                                "year" = "year") )
  
#  3) Bind & clean final table with both levels
both_levels_re2 <- bind_rows(sub_df, plot_df) %>%
  mutate(
    level   = factor(level, levels = c("subplot","plot")),
    plot_id = factor(plot_id),
    w       = pmin(pmax(w, 1), 50)   # cap weights so a few dense plots don't dominate
  ) %>% 
  left_join(df_plot_share_early, by = c('plot_id' = 'plot',
                                        "year" = "year"))

p_plot <- plot_df %>%
  #ungroup() %>%
  select(year, time_snc_full_disturbance, time_snc_part_disturbance, 
         mean_hgt, cv_hgt, shannon_sp, sp_richness,  CWM_shade ,
         CWM_drought) %>%
  pivot_longer(c(-year, 
                 - time_snc_full_disturbance, 
                 - time_snc_part_disturbance),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Plot level') +
  theme_classic2()


p_subplot <- sub_df %>%
  #ungroup() %>%
  select(year, time_snc_full_disturbance, time_snc_part_disturbance, 
         mean_hgt, cv_hgt, shannon_sp, sp_richness,  CWM_shade ,
         CWM_drought) %>%
  pivot_longer(c(-year, 
                 - time_snc_full_disturbance, 
                 - time_snc_part_disturbance),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Subplot level') +
  theme_classic2()
  
ggarrange(p_plot, p_subplot)
# add comparison between years at plot and subplot levels

make_violin_per_year <- function(df, y,
                                 ylim = NULL,
                                 drop_zeros = FALSE,
                                # p_y = NULL,           # e.g., p_y = 5  (data units)
                                 p_y_npc = 0.9,       # or p_y_npc = 0.9 (90% up)
                                 p_size = 3,
                                 p_method = "wilcox.test") {
  
  pd <- position_dodge(0.9)
  
  d <- df %>% filter(!is.na({{y}}))
  if (drop_zeros) d <- d %>% filter({{y}} > 0)
  
  p <- ggplot(d, aes(x = year, y = {{y}}, fill = year, color = year)) +
    geom_violin(alpha = 0.5, trim = TRUE, width = 0.8, position = pd) +
    geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.5,
                 width = 0.2, position = pd) +
    theme_grey(base_size = 8)
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  
  p + ggpubr::stat_compare_means(method = p_method, label = "p.format",
                                 size = p_size,
                                 #label.y = p_y,       # use either…
                                 label.y.npc = p_y_npc)  # …or this (0–1)
}


# porovnanie cez subplot
p.dens     <- make_violin_per_year(plot_df, dens_ha, ylim = c(0, 45000), p_y_npc = 0.80) #,   
p.dens
p.height   <- make_violin_per_year(plot_df, mean_hgt, drop_zeros = TRUE, ylim = c(0, 6), 
                                   p_y_npc = 0.15) #,  
p.cv       <- make_violin_per_year(plot_df, cv_hgt)
p.shannon  <- make_violin_per_year(plot_df, shannon_sp)
p.richness <- make_violin_per_year(plot_df, sp_richness)
p.eveness  <- make_violin_per_year(plot_df, evenness_sp)
p.eff      <- make_violin_per_year(plot_df, effective_numbers) #, ylim = c(0, 10)

# Arrange with a shared legend
out_plot <- ggarrange(p.dens, p.height, p.cv, p.shannon, p.richness,p.eveness,p.eff,
          common.legend = TRUE, legend = "bottom")

annotate_figure(out_plot, top = text_grob("Plot level", 
                                      color = "black", face = "bold", size = 14))


# get summary statistics
out_summary_full <- plot_df %>%
  ungroup() %>%
  select(year, level, mean_hgt, cv_hgt, shannon_sp, sp_richness, evenness_sp,effective_numbers) %>%
  pivot_longer(-c(year, level), names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  group_by(metric, year, level) %>%
  summarise(
    n_total  = n(),
    n_non_na = sum(!is.na(value)),
    n_na     = sum(is.na(value)),
    mean     = mean(value, na.rm = TRUE),
    sd       = sd(value, na.rm = TRUE),
    se       = sd/sqrt(n_non_na),
    median   = median(value, na.rm = TRUE),
    p25      = quantile(value, 0.25, na.rm = TRUE),
    p75      = quantile(value, 0.75, na.rm = TRUE),
    .groups  = "drop"
  ) 

out_summary_full

out_summary_years <- both_levels_re2 %>%
  ungroup() %>%
  select(year,  mean_hgt, cv_hgt, shannon_sp, sp_richness, evenness_sp,effective_numbers) %>%
  pivot_longer(-c(year), names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  group_by(metric, year) %>%
  summarise(
    n_total  = n(),
    n_non_na = sum(!is.na(value)),
    n_na     = sum(is.na(value)),
    mean     = mean(value, na.rm = TRUE),
    sd       = sd(value, na.rm = TRUE),
    se       = sd/sqrt(n_non_na),
    median   = median(value, na.rm = TRUE),
    p25      = quantile(value, 0.25, na.rm = TRUE),
    p75      = quantile(value, 0.75, na.rm = TRUE),
    .groups  = "drop"
  ) 

out_summary_years


# how does the diversity and composition changes over years?

ggplot(plot_df, aes(x = mean_hgt, y = shannon_sp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~year) 

ggplot(plot_df, aes(x = cv_hgt, y = shannon_sp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~year) 


ggplot(plot_df, aes(x = time_snc_full_disturbance, y =cv_hgt )) +
  geom_point() +
  geom_smooth(method = "lm") 

ggplot(plot_df, aes(x = time_snc_full_disturbance, y =mean_hgt )) +
  geom_point() +
  geom_smooth(method = "lm") 

summary(lm(sp_richness ~ time_snc_full_disturbance, data = plot_df))


# change over time




# Step 1: Reduce to one row per plot × year
plot_wide <- plot_df %>%
  select(plot_id, year, mean_hgt, shannon_sp, cv_hgt) %>%
  distinct() %>%  # make sure it's one row per plot-year
  pivot_wider(
    names_from = year,
    values_from = c(mean_hgt, shannon_sp,cv_hgt),
    names_sep = "_"
  ) %>%
  mutate(
    delta_hgt = mean_hgt_2025 - mean_hgt_2023,
    delta_div = shannon_sp_2025 - shannon_sp_2023,
    delta_cv = cv_hgt_2025 - cv_hgt_2023
  )

# View result
head(plot_wide)

ggplot(plot_wide, aes(x = delta_hgt, y = delta_div)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = "Δ Mean Height", y = "Δ Shannon Diversity") +
  theme_classic()

ggplot(plot_wide, aes(x = delta_hgt, y = delta_cv)) +
  geom_point() +
  geom_smooth(method = "gam", se = TRUE) +
  labs(x = "Δ Mean CV", y = "Δ Shannon Diversity") +
  theme_classic()


# Make a threshold for legacy effects: > 4 m

# ---- 1) Tall-stem legacy (change threshold to 4 for sensitivity) ----
tall_thresh <- 3.5   # meters

tall_sub <- dat_overlap %>%
  filter(!is.na(hgt_est), n > 0) %>%
  mutate(is_tall = hgt_est >= tall_thresh) %>%
  group_by(plot, subplot, year) %>%
  summarise(
    tall_n      = sum(n[is_tall], na.rm = TRUE),
    stems_total = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(tall_density_m2 = tall_n / 4)  # 4 m² subplot

legacy_plot <- tall_sub %>%
  group_by(plot, year) %>%
  summarise(
    tall_presence        = as.integer(any(tall_n > 0)),
    tall_share_subplots  = mean(tall_n > 0),
    tall_density_m2      = sum(tall_n) / 20,     # 20 m² per plot (5×4 m²)
    .groups = "drop"
  )

# graded class; collapse to binary later if sample sizes are tiny
nz <- legacy_plot %>% filter(tall_density_m2 > 0)
cuts <- if (nrow(nz) >= 3) quantile(nz$tall_density_m2, probs = c(1/3, 2/3), na.rm = TRUE) else c(0, 0)

legacy_plot <- legacy_plot %>%
  mutate(legacy_class = case_when(
    tall_density_m2 == 0 ~ "none",
    tall_density_m2 <= cuts[1] ~ "low",
    tall_density_m2 <= cuts[2] ~ "mid",
    TRUE ~ "high"
  )) %>%
  mutate(legacy_class = factor(legacy_class, 
                               levels = c("none","low","mid","high")))

# ---- 2) Join onto your modelling frame (both_levels_re2 already has plot_id) ----
both_levels_re4 <- both_levels_re2 %>%
  left_join(legacy_plot, by = c("plot_id" = "plot")) %>%
  left_join(cvx_stem_density, by = c("plot_id" = "plot")) %>% # year 2025 is automatically filtered
  mutate(
    legacy_class = as.character(legacy_class),
    legacy_class = ifelse(is.na(legacy_class), "none", legacy_class),
    legacy_class = factor(legacy_class, levels = c("none","low","mid","high")),
    tall_presence = ifelse(is.na(tall_presence), 0L, tall_presence)
  )


# If any legacy_class × level cell has < 10 rows, collapse to binary
tbl <- table(both_levels_re4$legacy_class, both_levels_re4$level)
if (length(tbl) > 0 && min(tbl) < 10) {
  both_levels_re4 <- both_levels_re4 %>%
    mutate(legacy_class = factor(ifelse(tall_presence == 1, "present", "absent"),
                                 levels = c("absent","present")))
}

# make sure types are right and there are no NAs for covariates/weights
# build the interaction factor explicitly
both_levels_re4 <- both_levels_re4 %>%
  mutate(
    plot_id      = factor(plot_id),
    level        = factor(level, levels = c("subplot","plot")),
    legacy_class = factor(legacy_class, levels = c("absent","present")),
    lev_legacy   = interaction(level, legacy_class, drop = TRUE),
    w            = pmin(pmax(w, 1), 50)
  ) %>% 
  filter(!is.na(mean_hgt), !is.na(w), !is.na(cv_hgt), !is.na(dens_m2))




# ---- 3) Fit GAM with level × legacy-specific smooths (controls: mean_hgt; RE: plot) ----
m_cv_legacy <- gam(
  cv_hgt ~ level + legacy_class +
    s(dens_m2, by = interaction(level, legacy_class), k = 5) +
    s(mean_hgt, k = 5) +
    s(plot_id, bs = "re"),
  data    = both_levels_re4,
  weights = w,
  family  = Gamma(link = "log"),
  method  = "REML"
)
summary(m_cv_legacy)
library(gratia)
appraise(m_cv_legacy)
plot.gam(m_cv_legacy, page = 1)

# test pre-disturbance trees:
both_levels_re4 %>% 
  ggplot(aes(x = pre_dist_dens_ha   ,
             y = dens_m2            )) +
  geom_point() + 
  geom_smooth()


both_levels_re4 %>% 
  ggplot(aes(x = pre_dist_dens_ha   ,
             y = cv_hgt            )) +
  geom_point() + 
  geom_smooth()


both_levels_re4 %>% 
  ggplot(aes(x = pre_dist_dens_ha   ,
             y = shannon_sp            )) +
  geom_point() + 
  geom_smooth()



m_tw_ML <- mgcv::gam(
  cv_hgt ~ level + legacy_class +
    s(dens_m2, by = interaction(level, legacy_class), k = 5) +
    s(mean_hgt, k = 5) + 
    s(plot_id, bs = "re"),
  data = both_levels_re4, 
  weights = w,
  family = mgcv::tw(link = "log"), method = "ML", select = TRUE
)

fin.m <- m_tw_ML
plot.gam(fin.m, page = 1, shade = T)


# addd pre disturbance conditison:
both_levels_re4 <- both_levels_re4 %>%
  mutate(
    pre_dens_m2 = pre_dist_dens_ha / 1e4,                   # trees / m^2
    pre_sc      = as.numeric(scale(log1p(pre_dens_m2)))     # stable, centered
  )


m_tw_pre1 <- gam(
  cv_hgt ~ level + legacy_class +
    s(dens_m2, by = interaction(level, legacy_class), k = 5) + #, bs = "cs"
    s(mean_hgt, k = 5, bs = "cs") +
    s(pre_sc,  by = level, k = 4, bs = "cs") +       # <-- history term
    s(plot_id, bs = "re"),
  data = both_levels_re4,
  weights = w,
  family = mgcv::tw(link = "log"),
  method = "ML", select = TRUE
)

AIC(m_tw_pre1, m_tw_ML)

plot.gam(m_tw_pre1, page = 1)

mh <- median(both_levels_re4$mean_hgt, na.rm = TRUE)
pd <- median(both_levels_re4$pre_sc,   na.rm = TRUE)



p <- predict_response(m_tw_pre1, terms = c("dens_m2 [all]"))

plot(p, one_plot = TRUE)

# Continue from HERE!!!!!! 10022025




# anchors for conditioning
mh <- median(both_levels_re4$mean_hgt, na.rm = TRUE)
md <- median(both_levels_re4$dens_m2,  na.rm = TRUE)



## ------------ B) CV ~ dens_m2 (two lines), legacy fixed ----------------------------
g_dens <- ggpredict(
  m_tw2,
  terms          = c("dens_m2 [all]", paste0("lev_legacy [", paste(levels_abs, collapse=","), "]")),
  condition      = list(mean_hgt = mh),
  type           = "fixed",
  back.transform = TRUE
) |> as.data.frame() |>
  mutate(level = ifelse(grepl("^subplot", group), "subplot", "plot"))

p.density <- ggplot(g_dens, aes(x, predicted, color = level, fill = level)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, color = NA) +
  geom_line(size = 1.1) +
  labs(x = "Stem density (per m²)", y = "CV of tree height",
       title = "CV ~ stem density",
       subtitle = sprintf("mean_hgt fixed at %.2f m", mh)) +
  theme_classic2(base_size = 8)


# test for Shannon - does shannon predicts vertical diversity? -------------

# Tweedie (log link)
m_tw_shan <- gam(
  cv_hgt ~ level + legacy_class +
    s(shannon_sp, by = lev_legacy, k = 3) +
    s(mean_hgt, k = 3) +
    #s(dens_m2, k = 3 ) +
    s(plot_id, bs = "re"),
  data = both_levels_re4, weights = w,
  family = mgcv::tw(link = "log"), method = "ML", select = TRUE
)

m_tw_shan2 <- gam(
  cv_hgt ~ level + legacy_class +
    s(shannon_sp, by = lev_legacy, k = 3) +
    s(mean_hgt, k = 3) +
    s(dens_m2, k = 3 ) +
    s(plot_id, bs = "re"),
  data = both_levels_re4, weights = w,
  family = mgcv::tw(link = "log"), method = "ML", select = TRUE
)


m_cs_shan <- gam(
  cv_hgt ~ level + legacy_class +
    s(shannon_sp, by = lev_legacy, k = 3, bs = "cs") +  # ‘cs’ lets the wiggle shrink away
    s(dens_m2,  k = 3, bs = "cs") +
    s(mean_hgt, k = 3, bs = "cs") +
    s(plot_id, bs = "re"),
  data = both_levels_re4, weights = w,
  family = mgcv::tw(link = "log"), method = "ML", select = TRUE
)


m_lin_shan <- gam(
  cv_hgt ~ legacy_class + level * shannon_sp +                 # <- linear slopes by level
    s(dens_m2,  k = 3, bs = "cs") +
    s(mean_hgt, k = 3, bs = "cs") +
    s(plot_id,  bs = "re"),
  data = both_levels_re4, weights = w,
  family = tw(link = "log"), method = "ML", select = TRUE
)

# lev_legacy has levels like "subplot.absent", "plot.absent", ...
m_lin_4 <- gam(
  cv_hgt ~ 0 + lev_legacy + shannon_sp:lev_legacy +         # separate intercepts & slopes per group
    s(dens_m2,  k = 3, bs = "cs") +
    s(mean_hgt, k = 3, bs = "cs") +
    s(plot_id,  bs = "re"),
  data = both_levels_re4, weights = w,
  family = tw(link = "log"), method = "ML", select = TRUE
)

p4 <- predict_response(m_lin_4, terms = c("shannon_sp [all]", "lev_legacy"))
plot(p4, one_plot = TRUE)

AIC(m_lin_4, m_lin_shan)

both_levels_re4 %>% 
  ggplot(aes(x = shannon_sp,
             y = cv_hgt,
             color = level)) + 
  geom_point() + 
  geom_smooth()


# !!!!
AIC(m_tw_shan, m_tw_shan2,m_cs_shan,m_lin_shan)
summary(m_tw_shan)
summary(m_tw_shan2)
summary(m_cs_shan)

plot.gam(m_lin_shan, page = 1)

# Gamma (log link) with shrinkage
m_gam_shan <- gam(
  cv_hgt ~ level + legacy_class +
    s(shannon_sp, by = lev_legacy, k = 3, bs = "cs") +
    s(mean_hgt, k = 3, bs = "cs") +
    s(plot_id, bs = "re"),
  data = both_levels_re4, weights = w,
  family = Gamma(link = "log"), method = "ML", select = TRUE
)

AIC(m_tw_shan, m_gam_shan)
gratia::appraise(m_tw_shan); 
gratia::appraise(m_gam_shan)

p <- ggpredict(m_cs_shan)

plot(p, one_plot = TRUE )

p <- predict_response(m_lin_shan, terms = c("shannon_sp [all]" , "level [all]"
                                           ))
fin.m.shannon <- m_lin_shan
p
plot(p, one_plot = TRUE)

# anchors for conditioning
mh <- median(both_levels_re4$mean_hgt, na.rm = TRUE)
md <- median(both_levels_re4$dens_m2,  na.rm = TRUE)

# 1) Observed weights of legacy within each scale (subplot/plot)
w_legacy <- both_levels_re4 %>%
  count(level, legacy_class, name = "n") %>%
  group_by(level) %>%
  mutate(w = n / sum(n)) %>%
  ungroup()
# w_legacy has columns: level, legacy_class, w

# 2) Get predictions for all level×legacy combos across Shannon
p_all <- ggpredict(
  fin.m.shannon,
  terms = c(
    "shannon_sp [all]",
    "level [all]"
  ),
  condition      = list(mean_hgt = mh, dens_m2 = md),
  type           = "fixed",
  back.transform = TRUE
) %>% 
  as.data.frame() %>%
  mutate(
    level  = ifelse(grepl("^subplot", group), "subplot", "plot"),
    legacy = ifelse(grepl("present$", group), "present", "absent"),
    se     = (conf.high - conf.low) / (2 * 1.96)
  ) %>%
  left_join(w_legacy, by = c("level" = "level", "legacy" = "legacy_class"))

# 3) Legacy-averaged (marginal) curves per scale
p_marg <- p_all %>%
  group_by(level, x) %>%
  summarise(
    predicted = sum(w * predicted, na.rm = TRUE),
    # simple pooled SE under independence
    se        = sqrt(sum((w * se)^2, na.rm = TRUE)),
    .groups   = "drop"
  ) %>%
  mutate(
    conf.low  = predicted - 1.96 * se,
    conf.high = predicted + 1.96 * se
  )

# 4) Plot: two lines (subplot vs plot), legacy averaged
p.shannon <- ggplot(p_marg, aes(x, predicted, color = level, fill = level)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, color = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "Species Shannon (subplot/plot scale-consistent)",
    y = "Predicted CV of tree height",
    title = "CV ~ Shannon, legacy-marginalized",
    subtitle = sprintf("Density fixed at %.2f per m²; mean height fixed at %.2f m", md, mh)
  ) +
  theme_classic2(base_size = 8)


ggarrange(p.density, p.shannon, common.legend = TRUE)





# how many tree species? 
length(unique(dat23_subplot_recode$species[dat23_subplot_recode$n > 0]))

sort(unique(dat23_subplot_recode$species[dat23_subplot_recode$n > 0]))

# which species has teh most stems?
dat23_subplot_recode %>%
  group_by(species) %>%
  summarise(total_stems = sum(n, na.rm = TRUE)) %>%
  mutate(share = total_stems / sum(total_stems)) %>%
  arrange(desc(total_stems))# %>%
  #slice(1)

# which species occured in most plots?
# Total number of unique subplots in the whole dataset (including empty ones)
total_subplots <- length(unique(dat23_subplot_recode$subplot))

# Species occurrence with subplot share
dat23_subplot_recode %>%
  filter(n > 0) %>%
  distinct(species, subplot) %>%
  count(species, name = "n_subplots") %>%
  mutate(share = n_subplots / total_subplots) %>%
  arrange(desc(n_subplots))


# get height structure
summary(field_sub_summ$median_hgt)
hist(dat23_subplot_recode$hgt_est)
summary(dat23_subplot_recode$hgt_est)

mean(dat23_subplot_recode$hgt_est < 1, na.rm = TRUE) * 100 # % of trees less than 1 m tall?

# merge subplot levels summaries with sf data to see it on the map
# 

field_sub_summ_sf <- left_join(dat23_sf_min,
                               field_sub_summ)
# Select variables to map
field_sub_summ_sf %>%
  dplyr::select(geom, 
               # stems_total, 
                #mean_hgt,
               cv_hgt 
               # sp_richness 
               ) %>%
  pivot_longer(cols = -geom, names_to = "variable", values_to = "value") %>%
  ggplot() +
  geom_sf(aes(color = value), size = 2, alpha = 0.5) +
  facet_wrap(~variable) +
  scale_color_viridis_c(option = "C", na.value = "grey80", guide = "colorbar") +
  theme_minimal() +
  theme(legend.position = "right")

summary(field_sub_summ$cv_hgt)




### Merge data: Correlation of height structure between field& drone data ----------
# rename cols for correlation

field_sub_summ_sub <- 
  field_sub_summ %>% 
  dplyr::select(subplot, median_hgt, mean_hgt, cv_hgt) %>% 
  rename(
    field_median_hgt = median_hgt,
    field_mean_hgt   = mean_hgt,
    field_cv_hgt     = cv_hgt)

drone_sub_summ_sub <- drone_sub_summ %>% 
  dplyr::select(subplot, median_hgt, mean_hgt, cv_hgt) %>% 
  rename(
    drone_median_hgt = median_hgt,
    drone_mean_hgt = mean_hgt,
    drone_cv_hgt = cv_hgt)
  

# join both tables to run corelations
df_cor_comb <- field_sub_summ_sub %>% 
  inner_join(drone_sub_summ_sub) %>% 
  na.omit() %>% 
  # remove extreme values
  dplyr::filter(field_median_hgt<5 & field_median_hgt>-0.5) %>%
  dplyr::filter(drone_median_hgt<5 & drone_median_hgt>-0.5) #%>%
  

plot(df_cor_comb$field_median_hgt, df_cor_comb$drone_median_hgt)
plot(df_cor_comb$field_cv_hgt, df_cor_comb$drone_cv_hgt)

cor(df_cor_comb$field_cv_hgt, df_cor_comb$drone_cv_hgt)

plot(df_cor_comb$field_cv_hgt, df_cor_comb$drone_cv_hgt)

# very low correlation between height structure between drones and field data on subplot level!!! 

# maybe try again on mean plot/pooled plot level?

# test for Spatial autocorrelation ----------------------------------------------
# summarize first on plot level (as subplots arre not independent): stem density, species composition (Hill number)
str(dat23_sf)

# Step 1: rename columns
dat23_sf <- dat23_sf %>%
  rename(
    subplot = ID,
    plot = cluster
  )

# Step 2: calculate mean coordinates per plot
# First extract coordinates
coords <- st_coordinates(dat23_sf)

dat23_plot_sf <- dat23_sf %>%
  mutate(
    X = coords[,1],
    Y = coords[,2]
  ) %>%
  st_drop_geometry() %>%       # drop geometry so we can summarise easily
  group_by(plot) %>%
  summarise(
    mean_X = mean(X, na.rm=TRUE),
    mean_Y = mean(Y, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("mean_X","mean_Y"), crs = st_crs(dat23_sf))

# Check result
plot(dat23_plot_sf)

# join attributes to spatial object
dat23_plot_merged <- dat23_plot_sf %>%
  left_join(plot_density, by = "plot")

# check result
dat23_plot_merged

library(spdep)

# extract coordinates
coords <- st_coordinates(dat23_plot_merged)

# neighbor list: here I set max distance = 50 km (adjust to your context)
nb <- dnearneigh(coords, d1 = 0, d2 = 60000)   # neighbors within 0–50 km
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)                 # row-standardized weights

# run Moran's I test
moran_res <- moran.test(dat23_plot_merged$sum_density, lw, zero.policy=TRUE)

moran_res


# use k-nearest neightbour - make suer each plot has a neighbours

nb <- knn2nb(knearneigh(coords, k = 5))  # 5 nearest neighbors for each plot
lw <- nb2listw(nb, style="W")
moran.test(dat23_plot_merged$sum_density, lw)

# strong spatial correlation up to 60 km!!!!

# see correlogram

library(ncf)

coords <- st_coordinates(dat23_plot_merged)

# correlogram with 10 km bins
correlog_res <- correlog(
  x = coords[,1],
  y = coords[,2],
  z = dat23_plot_merged$sum_density,
  increment = 10000,  # 10 km bin width
  resamp = 1000       # permutation test
)

# inspect results
head(correlog_res)
# keep only bins up to 60 km (60000 m)
correlog_trim <- correlog_res
keep <- correlog_res$mean.of.class <= 60000
correlog_trim$mean.of.class <- correlog_res$mean.of.class[keep]
correlog_trim$cor           <- correlog_res$cor[keep]
correlog_trim$p             <- correlog_res$p[keep]

# plot
plot(correlog_trim$mean.of.class, correlog_trim$cor,
     type="b", pch=16,
     xlab="Distance class (m)", ylab="Moran's I",
     main="Spatial correlogram up to 60 km")
abline(h=0, lty=2, col="red")




# Decompose the heterogeneity across scales -----------------------

## TEST START -----------------------------------------

trees <- dat23_subplot_recode

# total species per plot
species_totals_plot <- trees %>%
  group_by(plot, species) %>%
  summarise(
    stems_n      = sum(n, na.rm = TRUE),
    stems_density= sum(stem_density, na.rm = TRUE),
    .groups = "drop"
  )


# per subplot height stats 
sub_height <- trees %>%
  group_by(plot, subplot) %>%
  summarise(
    N_sub   = sum(n, na.rm = TRUE),                              # total trees in this subplot
    mean_h  = ifelse(N_sub > 0, stats::weighted.mean(hgt_est, w = n, na.rm = TRUE), NA_real_),
    # weighted "within-subplot" SS and variance (unbiased if N_sub>1)
    SS_within_sub = ifelse(N_sub > 1, sum(n * (hgt_est - mean_h)^2, na.rm = TRUE), 0),
    var_h   = ifelse(N_sub > 1, SS_within_sub / (N_sub - 1), NA_real_),
    .groups = "drop"
  )

# decompose height variance within a plot: split variability within subplots and 
# between subplots using sum of squares

comp_plot <- sub_height %>%
  group_by(plot) %>%
  summarise(
    N_tot     = sum(N_sub, na.rm = TRUE),
    # Weighted mean height across *trees* in the plot
    mu_plot   = ifelse(N_tot > 0, stats::weighted.mean(mean_h, w = N_sub, na.rm = TRUE), NA_real_),
    
    # Sum of within-subplot SS (already computed per subplot)
    SS_within = sum(SS_within_sub, na.rm = TRUE),
    
    # Between-subplots SS: variation of subplot means around the plot mean, weighted by N_sub
    SS_between_sub = sum(N_sub * (mean_h - mu_plot)^2, na.rm = TRUE),
    
    # Total SS and convert to variances with common denominator (N_tot - 1)
    var_total_plot   = ifelse(N_tot > 1, (SS_within + SS_between_sub) / (N_tot - 1), NA_real_),
    var_within_plot  = ifelse(N_tot > 1,  SS_within                 / (N_tot - 1), NA_real_),
    var_bsub_plot    = ifelse(N_tot > 1,  SS_between_sub            / (N_tot - 1), NA_real_),
    
    # Convert components to CVs (unitless). Use the *plot* mean as the scale.
    cv_total_plot   = ifelse(!is.na(mu_plot) && mu_plot > 0, sqrt(var_total_plot)  / mu_plot, NA_real_),
    cv_within_plot  = ifelse(!is.na(mu_plot) && mu_plot > 0, sqrt(var_within_plot) / mu_plot, NA_real_),
    cv_bsub_plot    = ifelse(!is.na(mu_plot) && mu_plot > 0, sqrt(var_bsub_plot)   / mu_plot, NA_real_),
    .groups = "drop"
  )


# aggregate on landscape and get between plots

# Keep SS to aggregate
plot_ss <- sub_height %>%
  group_by(plot) %>%
  summarise(
    N_tot   = sum(N_sub, na.rm = TRUE),
    mu_plot = ifelse(N_tot > 0, stats::weighted.mean(mean_h, w = N_sub, na.rm = TRUE), NA_real_),
    SS_within = sum(SS_within_sub, na.rm = TRUE),
    .groups = "drop"
  )

# SS between subplots summed across plots
SS_between_sub_all <- sub_height %>%
  group_by(plot) %>%
  summarise(mu_plot = stats::weighted.mean(mean_h, w = N_sub, na.rm = TRUE),
            SS_between_sub = sum(N_sub * (mean_h - mu_plot)^2, na.rm = TRUE),
            N_tot = sum(N_sub, na.rm = TRUE), .groups="drop") %>%
  summarise(SS_between_sub_all = sum(SS_between_sub, na.rm=TRUE),
            N_all = sum(N_tot, na.rm=TRUE), .groups="drop")

# Grand mean across all trees
mu_grand <- with(plot_ss, ifelse(sum(N_tot, na.rm=TRUE) > 0,
                                 stats::weighted.mean(mu_plot, w = N_tot, na.rm = TRUE), NA_real_))

# SS between plots
SS_between_plots <- plot_ss %>%
  summarise(SS = sum(N_tot * (mu_plot - mu_grand)^2, na.rm = TRUE),
            N_all = sum(N_tot, na.rm=TRUE), .groups="drop")

# Components on a common denominator (N_all - 1)
SS_within_all  <- sum(plot_ss$SS_within, na.rm=TRUE)
SS_bsub_all    <- SS_between_sub_all$SS_between_sub_all
SS_bplots_all  <- SS_between_plots$SS
N_all          <- SS_between_plots$N_all

var_within_all <- SS_within_all / (N_all - 1)
var_bsub_all   <- SS_bsub_all   / (N_all - 1)
var_bplots_all <- SS_bplots_all / (N_all - 1)

cv_within_all  <- sqrt(var_within_all) / mu_grand
cv_bsub_all    <- sqrt(var_bsub_all)   / mu_grand
cv_bplots_all  <- sqrt(var_bplots_all) / mu_grand

# species composition ---------------------

# Shannon entropy and Hill number q = 1
H <- function(counts) {
  p <- counts / sum(counts)
  p <- p[p > 0]   # drop zeros
  -sum(p * log(p))
}
D1 <- function(counts) exp(H(counts))   # effective number of species




species_sub <- trees %>%
  group_by(plot, subplot, species) %>%
  summarise(stems = sum(n, na.rm=TRUE), .groups="drop")


species_plot <- trees %>%
  group_by(plot, species) %>%
  summarise(stems = sum(n, na.rm=TRUE), .groups="drop")

# Recompute LANDSCAPE species vector correctly: pool across PLOTS by SPECIES
spp_land <- species_plot %>%
  group_by(species) %>%
  summarise(stems = sum(stems, na.rm = TRUE), .groups = "drop")

# subplot diversity 

sub_D1 <- species_sub %>%
  group_by(plot, subplot) %>%
  summarise(
    stems_total = sum(stems, na.rm = TRUE),
    D1_sub = ifelse(stems_total > 0, D1(stems), NA_real_),
    .groups = "drop"
  )

head(sub_D1, 20)

# get alpha, and gamma diversity and betta as its difference 
# γ diversity per plot (all stems pooled)
plot_gamma <- species_plot %>%
  group_by(plot) %>%
  summarise(
    stems_plot = sum(stems, na.rm = TRUE),
    D1_gamma_plot = ifelse(stems_plot > 0, D1(stems), NA_real_),
    .groups = "drop"
  )

# α diversity per plot (mean of subplot diversities, weighted by stems)
plot_alpha <- sub_D1 %>%
  group_by(plot) %>%
  summarise(
    D1_alpha_plot = weighted.mean(D1_sub, w = stems_total, na.rm = TRUE),
    stems_plot    = sum(stems_total, na.rm = TRUE),
    .groups = "drop"
  )

# β diversity between subplots → plot
plot_beta <- left_join(plot_gamma, plot_alpha, by = c("plot","stems_plot")) %>%
  mutate(D1_beta_subplots = D1_gamma_plot / D1_alpha_plot)

head(plot_beta, 20)  

# alpha = 1 - monoculture, NA - empty plots, 


# Landscape γ (all stems pooled across plots)
# Landscape γ (pooled across all plots, by species)
land_gamma <- spp_land %>%
  summarise(D1_gamma_land = ifelse(sum(stems, na.rm=TRUE) > 0,
                                   D1(stems),
                                   NA_real_)) %>%
  pull(D1_gamma_land)

# Mean γ across plots (weighted by stems)
mean_gamma_plot <- weighted.mean(plot_beta$D1_gamma_plot,
                                 w = plot_beta$stems_plot,
                                 na.rm = TRUE)

# β diversity among plots (plot → landscape)
D1_beta_plots <- land_gamma / mean_gamma_plot

# Summarise
landscape_div <- tibble::tibble(
  D1_alpha_sub_mean  = weighted.mean(sub_D1$D1_sub, 
                                     w = sub_D1$stems_total, 
                                     na.rm = TRUE),
  D1_beta_subplots_mean = weighted.mean(plot_beta$D1_beta_subplots, w = plot_beta$stems_plot, na.rm = TRUE),
  D1_gamma_land = land_gamma,
  D1_beta_plots = D1_beta_plots
)

landscape_div

# make a tidy table with 3 rows (scales)
cross_scale <- tibble::tibble(
  scale = c("within_subplot", "between_subplots", "between_plots"),
  x_height = c(cv_within_all, cv_bsub_all, cv_bplots_all),
  y_composition = c(landscape_div$D1_alpha_sub_mean,
                    landscape_div$D1_beta_subplots_mean,
                    landscape_div$D1_beta_plots)
)

cross_scale

cross_scale$scale <- factor(
  cross_scale$scale,
  levels = c("within_subplot", "between_subplots", "between_plots")
)

# make a plot cross-scale post-disturbance
ggplot(cross_scale, aes(x = x_height, 
                        y = y_composition,
                        color = scale)) +
  geom_point(aes(shape = scale,
                 color = scale), size = 5, fill = "white") +
  #geom_text(aes(label = scale), vjust = -1, size = 4.5) +
  # scale_shape_manual(values = c(
  #   within_subplot = 21,   # circle with fill
  #   between_subplots = 21,
  #   between_plots = 21
  # )) +
  labs(
    x = "Height heterogeneity (CV)",
    y = "Compositional heterogeneity (Hill number, q=1)",
    title = "Cross-scale heterogeneity (post-disturbance)"
  ) +
  theme_classic(base_size = 8) + 
  lims(x = c(0,1),
       y = c(0,4))




# Cross scale variation ----------------------------------------------
##  variation of height classes across scales: subplot vs plot -----------------------

# vertical structure per subplot& plot:
# plot values are generated as 
# 1. average across plots 
# (compare how different the subplot is compared to plit mean )
# 2. as pooled all values (no subplots boundaries)


# merge plot versions: mean vs pool -----------------

plot_compare <- plot_metrics_mean %>%
  left_join(plot_metrics_pooled %>% 
              dplyr::select(
                     plot,
                     stems_total,
                     pooled_sp_richness = sp_richness,
                     pooled_shannon_sp  = shannon_sp,
                     pooled_evenness_sp = evenness_sp,
                     pooled_mean_hgt    = mean_hgt,
                     pooled_cv_hgt      = cv_hgt#,
                     #pooled_range_hgt   = range_hgt
                     ),
            by = c("plot"))


# calculate differene between two metrics
plot_compare <- plot_compare %>%
  mutate(
    delta_cv_hgt   = pooled_cv_hgt - mean_cv_hgt,
    delta_richness = pooled_sp_richness - mean_sp_richness
  ) %>% 
  mutate(
    height_type = case_when(
      delta_cv_hgt < -0.1 ~ "Cancellation",
      delta_cv_hgt >  0.1 ~ "Reinforcement",
      TRUE                ~ "Persistence"
    ),
    richness_type = case_when(
      delta_richness < 0  ~ "Cancellation",
      delta_richness == 0 ~ "Persistence",
      delta_richness > 0  ~ "Reinforcement"
    )
  )

# SKIP ----------------
hist(plot_compare$delta_richness)
hist(plot_compare$delta_cv_hgt)

# op <- par(mfrow = c(1, 2))  # two plots on one page
# #windows()
# hist(na.omit(plot_compare$delta_richness),
#      main = "Δ Richness (pooled − mean)", xlab = "Difference",
#      col = "grey85", border = "grey40")
# abline(v = 0, col = "red", lty = 2, lwd = 2)
# 
# hist(na.omit(plot_compare$delta_cv_hgt),
#      main = "Δ CV Height (pooled − mean)", xlab = "Difference",
#      col = "grey85", border = "grey40")
# abline(v = 0, col = "red", lty = 2, lwd = 2)
# 
# par(op)

# -----------------------------------------
# Linear model
# -----------------------------------------7
# x - mean of subplot values 
# y - pooled plot value
# what the subplot tells me on average? e
# what the whole plot looks like if i pool everything?














# prepare table:
df_fin <- plot_compare %>% 
 # left_join(pre_dist_wide) %>%  # add rasters for rought estimation: density_cover and % coniferous (leaf type)
 # left_join(pre_disturb_stem_density_filt) %>%  # stem density based on tree calculation
  left_join(plot_density) %>%   # get stem density by vertical classes
  left_join(mng_plot_intensity ) %>% # add context information
  #left_join(subplot_group_density_wide) %>%  # get stem density by trees likely planted/pioneer - need to group on plot level!!
  mutate(
    log_mean = log1p(mean_cv_hgt),
    log_resp = log1p(pooled_cv_hgt)
  ) 


# merge data based on subplot and plot level to make a scatter plot ----------------
# select from teh subplot data
df_subplot_select <- field_sub_summ %>% 
  dplyr::select(subplot, 
                shannon_sp ,
                cv_hgt) %>% 
  mutate(type = 'subplot') %>% 
  rename(ID = subplot,
         shannon = shannon_sp,
         cv = cv_hgt
         )


# select from teh plot data
df_plot_select <- df_fin %>% 
  dplyr::select(plot, 
                mean_shannon_sp, 
               # mean_mean_hgt, 
                mean_cv_hgt, 
                pooled_shannon_sp,
                #pooled_mean_hgt, 
                pooled_cv_hgt )  %>% 
  # Step 2: Pivot longer to get variable and value columns
  pivot_longer(cols = -plot, names_to = "variable", values_to = "value") %>%
  
  # Step 3: Separate the variable into type and metric
  separate(variable, into = c("scale", "metric"), sep = "_", extra = "merge") %>%
  
  # Step 4: Clean up metric names
  mutate(metric = recode(metric,
                         "shannon_sp" = "shannon",
                         "cv_hgt" = "cv")) %>%
  
  # Step 5: Pivot wider to get shannon and cv as separate columns
  pivot_wider(names_from = metric, values_from = value) %>% 
  rename(ID = plot)


# merge data from subplot and plot level 
df_merge <- rbind(df_subplot_select, df_plot_select)

df_merge <- df_merge %>%
  mutate(type = case_when(
    type == "mean"   ~ "mean_plot",
    type == "pooled" ~ "pooled_plot",
    TRUE             ~ type  # keep "subplot" as is
  ))

df_merge %>% 
  ggplot(aes(x = cv,
             y = shannon,
             color = type,
             fill = type)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.68, alpha = 0.5) +  # 68% confidence ellipse
   geom_point(alpha = 0.3) +
 theme_bw()



# subset only important variables
vars <- c(#"mean_sp_richness" ,     
          # "mean_shannon_sp",         
          #"mean_mean_hgt" ,  
          "delta_cv_hgt",                 # pooled_cv_hgt - mean_cv_hgt: if + => larger variation at large scale
          "mean_cv_hgt" ,          
        #"mean_range_hgt",
       # "pooled_sp_richness" ,   
         # "pooled_shannon_sp" ,
        #  "pooled_mean_hgt"  ,
       "stems_total",
       "stem_regeneration",
       "pooled_cv_hgt" ,        
       #  "pooled_range_hgt",
       # "pre_dst_median_cover_dens_b25", 
        "pre_dst_median_cover_dens_b100",
       "pre_dst_cv_cover_dens_b100",
       # "pre_dst_share_coniferous_b25",
        "pre_dst_share_coniferous_b100",
       # "disturbance_year" ,     
      #  "forest_year" ,           
       # "disturbance_length" ,   
        # "n_trees_total",          
        "pre_dst_density_ha" ,           
         "pre_dst_share_spruce",
      "management_intensity",
      "salvage_intensity",
      "protection_intensity")

# select jus selected variables
df_fin_sub <- df_fin %>% 
  dplyr::select(all_of(vars))

pairs(df_fin_sub)

#library(corrplot)

# compute correlation matrix (excluding NAs)
cor_mat <- cor(df_fin_sub, use = "pairwise.complete.obs")

# plot correlogram
corrplot(cor_mat, method = "color", type = "upper",
         addCoef.col = "black",#, # add correlation values
         tl.col = "black", 
        tl.srt = 45,
         #diag = FALSE
         )

# link current structure with pre-disturbance history
plot(df_fin$delta_cv_hgt      ~ df_fin$salvage_intensity)
plot(df_fin$mean_cv_hgt      ~ df_fin$stems_total)
plot(df_fin$stems_total ~ df_fin$share_spruce)
plot(df_fin$stems_total ~ df_fin$cv_cover_dens_b100)

plot(df_fin$pooled_cv_hgt      ~ df_fin$share_spruce)

plot(df_fin$disturbance_length ~ df_fin$share_spruce)

df_fin %>% 
  ggplot(aes(x = salvage_intensity, y = delta_cv_hgt)) +
  geom_smooth(method = 'lm') +
  geom_point()

df_fin %>% 
  ggplot(aes(x = mean_cv_hgt, y = share_spruce)) +
  geom_smooth(method = 'lm') +
  geom_point()



# test: delta vs managemnet ------------------------------------------------------------------------
library(mgcv)
library(gratia)

m1 <- gam(delta_cv_hgt ~ s(protection_intensity, k = 5), data = df_fin)

summary(m1)
plot.gam(m1)
# # plots --------------------
# 
# 
# op <- par(mfrow = c(1, 2))  # two plots on one page
# 
# hist(na.omit(df_fin$share_spruce),
#      main = "Spruce share", xlab = "Share [%]",
#      col = "grey85", border = "grey40")
# 
# hist(na.omit(df_fin$density_ha),
#      main = "Mature trees stem density", xlab = "#/ha",
#      col = "grey85", border = "grey40")
# 
# par(op)
# 

# drivers ---------------------------------------------------------
# model testing -> best model is with logged values
model <- lm(pooled_cv_hgt ~ mean_cv_hgt, data = df_fin)
summary(model)

model_poly <- lm(pooled_cv_hgt ~ poly(mean_cv_hgt, 2), data = df_fin)
summary(model_poly)

model_log <- lm(log1p(pooled_cv_hgt) ~ log1p(mean_cv_hgt), data = df_fin)
summary(model_log)


library(MASS)
model_robust <- rlm(pooled_cv_hgt ~ mean_cv_hgt, data = df_fin)
summary(model_robust)

model_gam <- gam(pooled_cv_hgt ~ s(mean_cv_hgt, k = 4), data = df_fin)
summary(model_gam)
gam.check(model_gam)
appraise(model_gam)
plot(model_gam)
draw(model_gam, residuals = TRUE, rug = TRUE)


eps <- max(1e-6, min(plot_compare$pooled_cv_hgt[plot_compare$pooled_cv_hgt > 0], na.rm=TRUE) / 10)
m_gam_gamma <- mgcv::gam(I(pooled_cv_hgt + eps) ~ s(mean_cv_hgt),
                         family = Gamma(link="log"),
                         data = df_fin, method="REML")

m_gam_tw <- mgcv::gam(pooled_cv_hgt ~ s(mean_cv_hgt),
                      family = mgcv::tw(link="log"),  # p will be estimated
                      data = df_fin, method="REML")

m_gam_t <- mgcv::gam(pooled_cv_hgt ~ s(mean_cv_hgt),
                     family = mgcv::scat(), data = df_fin, method="REML")

# log transformed is teh best!!
m_gam_log <- gam(log_resp ~ s(log_mean), data = df_fin, method = "REML")

AIC(m_gam_gamma,model_gam, m_gam_tw , m_gam_t, m_gam_log)

fin.m.hgt_cv <- m_gam_log


# add pre-disturbance condistion
m_gam_log1 <- gam(log_resp ~ s(log_mean) + s(share_coniferous),
                  data = df_fin, method = "REML")
m_gam_log2 <- gam(log_resp ~ s(log_mean) + s(shannon),
                  data = df_fin, method = "REML")
m_gam_log3 <- gam(log_resp ~ s(log_mean) + s(cv_cover_dens),
                  data = df_fin, method = "REML")
# indluclude tree based metrics
m_gam_log4 <- gam(log_resp ~ s(log_mean) + s(cv_cover_dens) +
                    s(share_spruce) +s(density_ha) + s(disturbance_length, k = 3 ) ,
                  data = df_fin, method = "REML")

m_gam_log4 <- gam(log_resp ~ s(log_mean) + s(cv_cover_dens) +
                    s(share_spruce) +s(density_ha) + s(disturbance_length, k = 3 ) ,
                  data = df_fin, method = "REML")

# CV of heighst ~ stem density 
df_fin %>% 
  ggplot(aes(y = pooled_mean_hgt, x = mean_mean_hgt )) +
  geom_point() +
  geom_smooth()

cor(df_fin$pooled_mean_hgt, df_fin$mean_mean_hgt, use = "complete.obs")

hist(df_fin$mean_cv_hgt)

AIC(m_gam_log2, m_gam_log1, m_gam_log3, m_gam_log, m_gam_log4)

m<-m_gam_log4

summary(m)
appraise(m)
plot.gam(m, page = 1, shade = T)


# backtransform teh values and plot
library(ggeffects)

p_log_mean  <- ggpredict(m, terms = "log_mean [all]")
p_spruce    <- ggpredict(m, terms = "share_spruce [all]")

# --- Back-transform to original scale (pooled_cv_hgt) ---
bt <- function(df) {
  df %>%
    mutate(predicted = expm1(predicted),
           conf.low  = expm1(conf.low),
           conf.high = expm1(conf.high))
}

p_log_mean_bt <- bt(p_log_mean)
p_spruce_bt   <- bt(p_spruce)

# --- 2) Simple plots on the log1p scale (quick + minimal) ---
pl1 <- plot(p_log_mean_bt) +
  labs(x = "log(mean CV of height)", y = "log1p(pooled CV of height)") +
  theme_bw()

pl2 <- plot(p_spruce_bt) +
  labs(x = "Spruce share (proportion)", y = "log1p(pooled CV of height)") +
  theme_bw()

pl1 + pl2




gp <- ggpredict(m_gam_log, terms = "log_mean [all]")  # on log-scale
df <- as.data.frame(gp) %>%
  transmute(
    mean_cv_hgt = expm1(x),                 # back-transform x
    fit        = expm1(predicted),          # back-transform y
    lo         = expm1(conf.low),
    hi         = expm1(conf.high)
  )


ggplot() +
  geom_point(aes(mean_cv_hgt, pooled_cv_hgt),
             data = df_fin, color = "grey60", size = 1.8) +
  geom_ribbon(aes(mean_cv_hgt, ymin = lo, ymax = hi), data = df, alpha = 0.2) +
  geom_line(aes(mean_cv_hgt, fit), data = df, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +  # 1:1 reference
  labs(x = "Mean subplot CV (heights)",
       y = "Pooled plot CV (heights)",
       title = "GAM (log–log) predictions back-transformed") +
  theme_minimal()



AIC(model, model_poly, model_log, model_robust, model_gam)


# run derivatives - gradia: to understand at what 
# numbers the where the curve is increasing (reinforcement)
# flast (persistence) or decreasing (cancellation)

# Derivatives of smooth (with 95% CI)
d1 <- derivatives(model_gam, select = "s(mean_cv_hgt)")

# Quick look at where slope is significantly >0 or <0
d1 <- d1 %>%
  mutate(cs_type = case_when(
    .derivative > 0 ~ "Reinforcement",   # strictly positive slope
    .derivative < 0 ~ "Cancellation",    # strictly negative slope
    TRUE          ~ "Persistence"      # includes 0
  ))


# Plot derivative curve
ggplot(d1, aes(x = mean_cv_hgt, y = .derivative)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  labs(x = "Mean subplot CV (heights)",
       y = "d(Pooled CV)/d(Mean subplot CV)",
       title = "GAM first derivative with 95% CI") +
  theme_minimal()






# Basic scatterplot
plot(plot_compare$mean_cv_hgt, plot_compare$pooled_cv_hgt,
     xlab = "Mean subplot CV (heights)",
     ylab = "Pooled plot CV (heights)",
     main = "Cross-scale regression: subplot → plot CV",
     pch = 19, col = "grey50")

# Add fitted regression line
fit <- lm(pooled_cv_hgt ~ mean_cv_hgt, data = plot_compare)
abline(fit, col = "red", lwd = 2)

# Add 1:1 persistence line
abline(0, 1, lty = 2, col = "blue")

# Add legend
legend("topleft", legend = c("Plots", "Fitted regression", "1:1 line"),
       col = c("grey50", "red", "blue"), pch = c(19, NA, NA),
       lty = c(NA, 1, 2), lwd = c(NA, 2, 1))


### get stem density and shannon for field data --------------
df_cluster <- dat23_subplot_recode %>% 
  group_by(cluster, manag_intensity, salvage_intensity, protection_intensity) %>% 
  summarise(stem_density = sum(stem_density, na.rm = T),
            basal_area_ha_m2 = sum(ba_ha_m2, na.rm = TRUE))#,
           # mean_hgt  = mean)



# 

# Ensure height is treated as a factor, ignore missing/empty strings
dat_clean <- dat23_subplot_recode %>%
  filter(hgt_est != "" & !is.na(hgt_est))

# 1. Number of unique height classes per plot
height_div_plot <- dat_clean %>%
  group_by(ID, cluster) %>%
  summarise(n_height_classes = n_distinct(hgt_est)) %>%
  ungroup()

# 2. Number of unique height classes per cluster
height_div_cluster <- dat_clean %>%
  group_by(cluster) %>%
  summarise(n_height_classes = n_distinct(hgt_est)) %>%
  ungroup()

# 3. Number of unique height classes across landscape
n_height_classes_landscape <- dat_clean %>%
  summarise(n_height_classes = n_distinct(hgt_est)) %>%
  pull(n_height_classes)



# Compute medians and IQRs
summary_stats <- tibble(
  spat_scale = factor(c("Plot", "Cluster", "Landscape"),
                    level = c("Plot", "Cluster", "Landscape")),
  median = c(
    median(height_div_plot$n_height_classes),
    median(height_div_cluster$n_height_classes),
    n_height_classes_landscape
  ),
  IQR_low = c(
    quantile(height_div_plot$n_height_classes, 0.25),
    quantile(height_div_cluster$n_height_classes, 0.25),
    n_height_classes_landscape
  ),
  IQR_high = c(
    quantile(height_div_plot$n_height_classes, 0.75),
    quantile(height_div_cluster$n_height_classes, 0.75),
    n_height_classes_landscape
  )
)


ggplot(summary_stats, aes(y = spat_scale, x = median)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = IQR_low, xmax = IQR_high), 
                 height = 0) +
  labs(
    x = "Number of Height Classes",
    y = "Spatial Scale",
    title = "Height Class Diversity Across Spatial Scales"
  ) +
  theme_classic2()

## Plots vs cluster --------------------------------------------------------------
## HOw value in one plot affect clusre level mean? sd?
# split table in towo group: 
# select one random subplot
# subset remaining subplots - calculate means and sd

# One random row per cluster
height_div_subset <- height_div_plot %>%
  group_by(cluster) %>%
  slice_sample(n = 1) %>%
  ungroup() %>% 
  rename(plt_n_hgt = n_height_classes)

# Remaining rows
height_div_remaining <- anti_join(height_div_plot, height_div_subset, by = c("ID", "cluster"))

height_div_remaining_summary <- height_div_remaining %>% 
  group_by(cluster) %>% 
  summarize(cl_mean_hgt_n = mean(n_height_classes, na.rm = T),
            cl_median_hgt_n = median(n_height_classes, na.rm = T),
            cl_sd_hgt_n = sd(n_height_classes, na.rm = T),
            cl_cv_hgt_n = cl_sd_hgt_n/cl_mean_hgt_n) %>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  left_join(height_div_subset) #%>% # replace NA by 0s
  
# how does the plot value predict the cluster value? 
# predcit average value per plot
height_div_remaining_summary %>% 
  ggplot(aes(x = plt_n_hgt,
             y = cl_mean_hgt_n )) + 
  geom_point() +
  geom_smooth(method = 'loess') +
  theme_classic()

### predict coeffcient of variation -
height_div_remaining_summary %>% 
  ggplot(aes(x = plt_n_hgt,
             y = cl_cv_hgt_n )) + 
  geom_point() +
  geom_smooth(method = 'loess') +
  theme_classic()



### bootstrap analysis of how the random selection of one point represent the values per ----
# cluster

# Set number of bootstrap iterations
n_boot <- 1000

# Store results
bootstrap_results <- map_dfr(1:n_boot, function(i) {
  
  # 1. Random plot per cluster
  height_div_subset <- height_div_plot %>%
    group_by(cluster) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    rename(plot_ID = ID,
           plt_n_hgt = n_height_classes)
  
  # 2. Remaining plots
  height_div_remaining <- anti_join(
    height_div_plot,
    height_div_subset,
    by = c("ID" = "plot_ID", "cluster", "n_height_classes" = "plt_n_hgt")
  )
  
  # 3. Mean height richness of remaining plots
  cluster_remaining_summary <- height_div_remaining %>%
    group_by(cluster) %>%
    summarize(cl_mean_hgt_n = mean(n_height_classes, na.rm = TRUE)) %>%
    ungroup()
  
  # 4. Merge back to random plots
  paired <- height_div_subset %>%
    left_join(cluster_remaining_summary, by = "cluster") %>%
    mutate(iteration = i)
  
  return(paired)
})


### calculate average LOESS fit -----------------
# Define prediction grid (plot-level richness: 1, 2, 3)
pred_grid <- tibble(plt_n_hgt = seq(1, 3, by = 0.1))

# Apply LOESS fit to each iteration and predict over the grid
loess_fits <- bootstrap_results %>%
  filter(!is.na(cl_mean_hgt_n)) %>%
  group_by(iteration) %>%
  group_split() %>%
  map_dfr(function(df) {
    loess_model <- tryCatch(loess(cl_mean_hgt_n ~ plt_n_hgt, data = df, span = 0.75), error = function(e) NULL)
    if (!is.null(loess_model)) {
      pred <- predict(loess_model, newdata = pred_grid)
      tibble(plt_n_hgt = pred_grid$plt_n_hgt, cl_mean_pred = pred, iteration = unique(df$iteration))
    } else {
      NULL
    }
  })

# Aggregate: mean and 95% CI of LOESS fits
loess_summary <- loess_fits %>%
  group_by(plt_n_hgt) %>%
  summarize(
    mean_fit = mean(cl_mean_pred, na.rm = TRUE),
    lower_CI = quantile(cl_mean_pred, 0.025, na.rm = TRUE),
    upper_CI = quantile(cl_mean_pred, 0.975, na.rm = TRUE)
  )


ggplot(loess_summary, aes(x = plt_n_hgt, y = mean_fit)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray") +
  labs(
    x = "Height Class Richness in Random Plot",
    y = "Mean Height Class Richness in Remaining Plots",
    title = "Average LOESS Fit Across Bootstraps"
  ) +
  theme_classic2()



### calculate pseudo R2 ----------------------- 
# Function to compute pseudo-R² for each LOESS fit
compute_loess_r2 <- function(df) {
  loess_model <- tryCatch(loess(cl_mean_hgt_n ~ plt_n_hgt, data = df, span = 0.75), error = function(e) NULL)
  if (!is.null(loess_model)) {
    preds <- predict(loess_model, newdata = df)
    y <- df$cl_mean_hgt_n
    rss <- sum((y - preds)^2, na.rm = TRUE)
    tss <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- 1 - (rss / tss)
    return(tibble(iteration = unique(df$iteration), r2 = r2))
  } else {
    return(tibble(iteration = unique(df$iteration), r2 = NA_real_))
  }
}

# Apply to each bootstrap iteration
r2_results <- bootstrap_results %>%
  filter(!is.na(cl_mean_hgt_n)) %>%
  group_by(iteration) %>%
  group_split() %>%
  map_dfr(compute_loess_r2)

# Summarize R² distribution
r2_summary <- r2_results %>%
  summarize(
    mean_r2 = mean(r2, na.rm = TRUE),
    median_r2 = median(r2, na.rm = TRUE),
    r2_2.5 = quantile(r2, 0.025, na.rm = TRUE),
    r2_97.5 = quantile(r2, 0.975, na.rm = TRUE)
  )


ggplot(r2_results, aes(x = r2)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  labs(
    x = expression(R^2 ~ " of LOESS fits across bootstraps"),
    y = "Frequency",
    title = "Distribution of Pseudo-" ~ R^2 ~ " Values"
  ) +
  theme_minimal()





# make working example: 

dd <- data.frame(
  cluster = c(1, 1,   
              2, 2, 2,  
              3,   
              4, 
              5),
  species = c("a", "b",
              "a", "b", "d",
              NA, 
              "d",
              "d"),
  basal_area = c(5, 5,  
                 1, 1, 8,  
                 0,  
                 10, 
                 0.5)
)

## merge field data with drone estimation ------------------------------

### Sensitivity analysis ----------

# Function to generate boxplot for a given variable
plot_buffer_variation <- function(data, var_name, y_label) {
  ggplot(data, aes(x = factor(buffer_size), y = .data[[var_name]])) +
    geom_boxplot(fill = "skyblue") +
    labs(
      title = paste("", y_label, ""),
      x = "Buffer size (m)",
      y = y_label
    ) +
    theme_minimal()
}

# Generate plots
p_mean   <- plot_buffer_variation(drone_cv, "mean_height", "Mean Height")
p_median <- plot_buffer_variation(drone_cv, "median_height", "Mean Median")
p_sd     <- plot_buffer_variation(drone_cv, "sd_height", "SD of Height")
p_cv     <- plot_buffer_variation(drone_cv, "cv_height", "CV of Height")
p_max    <- plot_buffer_variation(drone_cv, "max_height", "Max Height")
p_min    <- plot_buffer_variation(drone_cv, "min_height", "Min Height")

# Combine using ggarrange
ggarrange(p_mean, p_median, p_sd, p_cv, p_max, p_min, 
          ncol = 2, nrow = 3, 
          labels = "AUTO")

# convert to wide format as i have 4 different buffer sizes - need to do this for cluster level as well!
drone_cv_wide <- drone_cv %>%
  dplyr::filter(buffer_size == 20) %>%  # filter to single buffer, representing our sample (increase to 22-25 m!!)
  pivot_wider(
    id_cols = c(cluster, drone),
    names_from = buffer_size,
    values_from = c(mean_height, sd_height, cv_height, max_height, min_height),
    names_glue = "{.value}_{buffer_size}"
  )


# merge with structural data ---------------------------------------
# Add prefix to drone data
drone_cv_wide_renamed <- drone_cv_wide %>%
  rename_with(~ paste0("drone_", .), -c(cluster, drone))

# Add prefix to pre-disturbance data
pre_dist_renamed <- pre_dist_history_2023 %>%
  rename_with(~ paste0("pre_", .), -cluster)

# Now join
df_fin <- df_field_diversity %>%
  right_join(pre_dist_renamed, by = "cluster") %>% 
 # right_join(drone_cv_wide_renamed, by = "cluster") %>%
  right_join(df_similarity, by = "cluster") %>% 
  right_join(height_div_remaining_summary)
  

# any correlation between variables?
# split between field based vs drone based. pre-dciturbance are drivers. 
# Step 1: Select numeric columns only
# Step 1: Select numeric variables
df_numeric <- df_fin %>%
  select(where(is.numeric)) %>%
  drop_na()  # drop rows with NAs for clean correlation and VIF

# Step 2: Correlation matrix
cor_mat <- cor(df_numeric)

# Step 3: Visualize with corrplot
corrplot(cor_mat,
         method = "circle",
         type = "lower",
         tl.col = "black",
         tl.cex = 0.8,
         number.cex = 0.7,
         title = "Correlation Matrix (Field & Drone Metrics)",
         mar = c(0,0,1,0),
         addCoef.col = "black")


# is tehre a correlation betwen field and drone derived vertical strcuture?
df_model %>% 
  ggplot(aes(x = shannon_height,
             y = drone_cv_height_20)) + 
  geom_point() + 
  geom_smooth(method = "loess")


#s does vertical diversity between plots predict stem density? -------------

df_fin %>% 
  ggplot(aes(x = shannon_height     ,
             y = stem_density )) + 
  geom_point() + 
  geom_smooth(method = "lm")




# Fit the linear model --------------------------------
df_model <- df_fin %>%
  dplyr::filter(drone_cv_height_20 < 10)  # example threshold


model <- lm(shannon_height ~ drone_cv_height_20, data = df_model)

# Summary output
summary(model)


# compare field vs drone metrics: ------------------------------------------------

# Define columns
field_metrics <- df_fin %>%
  select(stem_density, basal_area_ha_m2, shannon_species, shannon_height)

drone_metrics <- df_fin %>%
  select(drone_mean_height_20, drone_sd_height_20, drone_cv_height_20,
         drone_max_height_20, drone_min_height_20)

# Combine and track origin
combined_df <- bind_cols(field_metrics, drone_metrics)

# Label columns for plot clarity
colnames(combined_df) <- c(
  "StemDensity", "BA_ha", "Shannon_Species", "Shannon_Height",
  "Drone_Mean", "Drone_SD", "Drone_CV", "Drone_Max", "Drone_Min"
)

# Custom ggpairs plot: compare only Field vs Drone (upper part blank)
GGally::ggpairs(
  combined_df,
  columns = 1:9,
  lower = list(
    continuous = wrap("smooth", alpha = 0.6, size = 0.5)
  ),
  upper = list(continuous = "blank"),
  diag = list(continuous = "blankDiag"),
  mapping = aes(color = NULL)
) +
  ggtitle("Field vs Drone Metrics: Pairwise Comparison")



# how does cluster variability affects foresrt resilience? ---------------
head(df_fin)


cor(df_fin$cl_cv_hgt_n , df_fin$stem_density)

hist(df_fin$stem_density)

library(mgcv)
# library(tweedie)
# library(statmod)

# Start with null model
m_null <- gam(stem_density ~ 1, data = df_fin, family = tw(link = "log"))

# Try candidate models with individual smooth terms
m1 <- gam(stem_density ~ s(pre_mean_cover_dens, k = 3), data = df_fin, family = tw(link = "log"))
m2 <- gam(stem_density ~ s(pre_share_coniferous, k = 3), data = df_fin, family = tw(link = "log"))
m3 <- gam(stem_density ~ s(pre_shannon, k = 3), data = df_fin, family = tw(link = "log"))
m4 <- gam(stem_density ~ s(mean_jaccard, k = 3), data = df_fin, family = tw(link = "log"))
m5 <- gam(stem_density ~ s(alpha_mean, k = 3), data = df_fin, family = tw(link = "log"))
m6 <- gam(stem_density ~ s(beta, k = 3), data = df_fin, family = tw(link = "log"))
m7 <- gam(stem_density ~ s(manag_intensity, k = 3), data = df_fin, family = tw(link = "log"))
m8 <- gam(stem_density ~ s(salvage_intensity, k = 3), data = df_fin, family = tw(link = "log"))
m9 <- gam(stem_density ~ s(protection_intensity, k = 3), data = df_fin, family = tw(link = "log"))

# Compare AICs
AIC(m_null, m1, m2, m3, m4, m5, m6, m7, m8, m9)

df_fin %>%
  ungroup(.) %>% 
  dplyr::select(beta, mean_jaccard, salvage_intensity, protection_intensity) %>%
  cor(use = "complete.obs")



best_model <- gam(
  stem_density ~ 
  #  s(alpha_mean, k = 3) +
    s(cl_cv_hgt_n , k = 3) +
    #s(mean_jaccard, k = 3) +
    #s(manag_intensity, k = 3) ,#+
    s(salvage_intensity, k = 3) +
    s(protection_intensity, k = 3) +
    s(pre_share_coniferous, k = 3) +
    s(pre_mean_cover_dens, k = 3),
    
  data = df_fin,
  family = tw(link = "log")
)
summary(best_model)
plot(best_model,page = 1)

library(gratia)
appraise(best_model)
