

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
df_sub_long <- field_sub_summ %>%
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
         stems_total,
         mean_hgt, cv_hgt, shannon_sp, sp_richness,
         CWM_shade ,
         CWM_drought ) %>%
  pivot_longer(-c(year,
                    time_snc_full_disturbance,
                    time_snc_part_disturbance,
                    #CWM_shade ,
                    #CWM_drought,
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
p_partial_disturbance <- df_sub_long %>% 
  ggplot(aes(x = time_snc_part_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_part_disturbance ), outlier.shape = NA) +
  
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Subplot: Partial disturbance') +
  theme_classic2()


# see CV with time since disturbnace : full disturbance
p_full_disturbance <- df_sub_long %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Subplot: Full disturbance') +
  theme_classic2()

ggarrange(p_partial_disturbance, p_full_disturbance, ncol = 2)



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

fwrite(traits_full, 'outData/my_species_traits.csv')


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
trait_labs <- c(CWM_shade = "Shade tolerance", 
                CWM_drought = "Drought tolerance")

ggplot(all_long, aes(x = year, y = CWM, fill = year)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.15) +
  geom_jitter(width = 0.1, alpha = 0.25, size = 0.7) +
  facet_grid(level ~ trait, labeller = labeller(trait = trait_labs)) +
  labs(x = "Year", y = "Community-weighted mean (CWM)",
       title = "Trait CWMs by level and year") +
  theme_classic2(base_size = 10)



# is teh change in community shading/drought tolerance driven by planting????
p_shade_planting <- field_sub_summ %>% 
  ggplot(aes(x = as.factor(time_snc_full_disturbance),
             y = CWM_shade,
             fill = factor(planting))) +
  geom_boxplot()

p_shade_drought <- field_sub_summ %>% 
  ggplot(aes(x = as.factor(time_snc_full_disturbance),
             y = CWM_drought ,
             fill = factor(planting))) +
  geom_boxplot()

ggarrange(p_shade_planting, p_shade_drought, ncol = 1,
          nrow = 2)

# get heights by species
library(forcats)
dat_overlap %>% 
  dplyr::filter(hgt_est>0) %>% 
  mutate(species = fct_reorder(species, hgt_est, .fun = median, .desc = TRUE)) %>%
  ggplot(aes(x = species,
             y = hgt_est,
             fill = year)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,7.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# compare change in feights pioneer vs late
dat_overlap %>% 
  dplyr::filter(recovery_type != "other") %>%
  dplyr::filter(hgt_est>0) %>% 
  #mutate(species = fct_reorder(species, hgt_est, .fun = median, .desc = TRUE)) %>%
  ggplot(aes(x = recovery_type,
             y = hgt_est,
             fill = year)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,6)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


dat_overlap %>% 
  dplyr::filter(recovery_type != "other") %>%
  dplyr::filter(hgt_est>0) %>% 
  #mutate(species = fct_reorder(species, hgt_est, .fun = median, .desc = TRUE)) %>%
  ggplot(aes(x = recovery_type,
             y = n,
             fill = year)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,6)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




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
                                        "year" = "year")) %>% 
  mutate(early_class = case_when(
    share_early >= 60 ~ "early_dom",
    share_early <= 20 ~ "late_dom",
    TRUE              ~ "mixed"
  )) |> 
  mutate(early_class = factor(early_class, 
                              levels = c("late_dom",
                                         "mixed",
                                         "early_dom")))

table(both_levels_re2$early_class)

both_levels_re2 %>% 
  filter(level == 'subplot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance),
             y = cv_hgt,
             fill = early_class)) + 
  geom_boxplot() + 
  facet_grid(.~early_class)

p_early_cls_plot <- both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance),
             y = cv_hgt,
             color = early_class,
             fill = early_class,
             group = early_class)) +
  geom_point(alpha = 0.1) +
  ggtitle('plot') +
  geom_smooth(method = 'lm') +
  coord_cartesian(y = c(0,1))
  #facet_grid(.~early_class)

p_early_cls_subplot <- both_levels_re2 %>% 
  filter(level == 'subplot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance),
             y = cv_hgt,
             color = early_class,
             fill = early_class,
             group = early_class)) +
  geom_point(alpha = 0.1) +
  ggtitle('subplot') +
  geom_smooth(method = 'lm') +
  coord_cartesian(y = c(0,1))
#facet_grid(.~early_class)
ggarrange(p_early_cls_plot, p_early_cls_subplot, common.legend = T)

# test effect of early seral dominance
m_cv_cont <- gam(
  cv_hgt ~ s(share_early, k=5) + s(time_snc_full_disturbance, k=5) +
    te(share_early, time_snc_full_disturbance, k=c(4,4)) +
    level + s(plot_id, bs="re"),
  data = both_levels_re2, 
  family = gaussian()
)

m_cv_cont1 <- gam(
  cv_hgt ~ s(share_early, k=5) + 
    s(time_snc_full_disturbance, k=5) +
    te(share_early, time_snc_full_disturbance, k=c(4,4)) +
    level + s(plot_id, bs="re"),
  data = both_levels_re2, 
  family = gaussian()
)
summary(m_cv_cont1)
appraise(m_cv_cont1)
plot.gam(m_cv_cont1, page = 1)



summary(m_cv_cont)
draw(m_cv_cont)

m_cv <- gam(
  cv_hgt ~ early_class + level +
    s(time_snc_full_disturbance, by = early_class, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  data = both_levels_re2, family = gaussian()
)
summary(m_cv)

m_hgt <- gam(
  mean_hgt ~ early_class + level +
    s(time_snc_full_disturbance, by = early_class, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  data = both_levels_re2, family = tw(link = "log")
)

summary(m_hgt)

boxplot(cv_hgt ~ early_class, 
        data = both_levels_re2, 
        col = "grey90")
boxplot(mean_hgt ~ early_class, 
        data = both_levels_re2, 
        col = "grey90")

# CV smooths by class
draw(m_cv, select = 1:3)

# Mean height smooths by class
draw(m_hgt, select = 1:3)


# test ---------------------
library(mgcv)
library(gratia)
m_hgt_tw <- gam(
  mean_hgt ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2
)
summary(m_hgt_tw)
appraise(m_hgt_tw)
draw(m_hgt_tw, select = c(1,2))
plot(m_hgt_tw, page = 1)

m_stems_tw <- gam(
  stems_total ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2,
)
summary(m_stems_tw)
appraise(m_stems_tw)
draw(m_stems_tw, select = c(1,2))



m_share_early_tw <- gam(
  share_early ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2,
)
summary(m_share_early_tw)
appraise(m_share_early_tw)
draw(m_share_early_tw, select = c(1,2))



m_shade_tw <- gam(
  CWM_shade ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2,
)
summary(m_shade_tw)
appraise(m_shade_tw)
draw(m_shade_tw, select = c(1,2))


hist(both_levels_re2$share_early)

m_cv_tw <- gam(
  cv_hgt ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  #family = tw(link = "log"),  # Tweedie with log-link
  data = filter(both_levels_re2, cv_hgt >0)
)
summary(m_cv_tw)
appraise(m_cv_tw)
draw(m_cv_tw, select = c(1,2))

m_shan_tw <- gam(
  shannon_sp ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  #family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2
)
summary(m_shan_tw)
appraise(m_shan_tw)
draw(m_shan_tw, select = c(1,2))



# get variables over plot level !!!
df_plot_long <- plot_df %>%
  dplyr::filter(stems_total > 0) %>% 
  filter(cv_hgt >0) %>% 
  #ungroup() %>%
  select(year,
         time_snc_full_disturbance, 
         time_snc_part_disturbance, 
         dens_ha, 
         # management_intensity,
         mean_hgt, cv_hgt, shannon_sp, sp_richness,
         CWM_shade ,
         CWM_drought ) %>%
  pivot_longer(-c(year,
                  time_snc_full_disturbance,
                  time_snc_part_disturbance#,
                  #CWM_shade ,
                  #CWM_drought,
                 ),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value))# %>%





df_plot_long %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Plot: Full disturbance') +
  theme_classic2()

#!!!!



# check my pioneer vs early species separation in drought/shade tolerance?
p_traits_development1 <- dat_overlap %>% 
  filter(recovery_type != 'other') %>% 
  ggplot(aes(x = recovery_type, 
             y = Drought_tolerance)) +
  geom_boxplot()

p_traits_development2 <- dat_overlap %>% 
  filter(recovery_type != 'other') %>% 
  ggplot(aes(x = recovery_type, 
             y = Shade_tolerance )) +
  geom_boxplot()


ggarrange(p_traits_development1, p_traits_development2)




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




p <- predict_response(m_tw_pre1, terms = c("dens_m2 [all]"))

plot(p, one_plot = TRUE)

#