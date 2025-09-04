

# --------------------------------------------------------
#             Cross-scale interaction: 
#    field & drone variation in vertical structures 
# --------------------------------------------------------


# read data: 
# - field data from 2023 (2025?)
# - drone 2023&2024

# - get sites history - from tree density cover, dominant leaf type
#                     - from tree density (manually mapped)

# correlate: vertical and horizontal structure across scales: field vs drone
# cross scale interaction: 
#  - vertical structure
# effect of pre-disturbance condistions

# extras: 
# - identify damaged terminal: from raw data czechia 20203 - export gpkg
# - identify types of damages: 2025 - export gpkg

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


# Read files --------------------------------------

## Field data 
# 2023 
dat23_subplot <- fread('outData/subplot_full_2023.csv')
#dat23_cluster <- fread('outData/df_cluster_2023.csv')
# Save spatial subplot data with cluster IDs
dat23_sf <- st_read("outData/sf_context_2023.gpkg")

# 2025
#dat25_subplot <- fread('outData/subplot_full_2025.csv')
#dat25_cluster <- fread('outData/df_cluster_2025.csv')

# Save spatial subplot data with cluster IDs
#dat25_sf <- st_read("outData/subplot_with_clusters_2025.gpkg")

# Read drone CHM per pixel
drone_cv <- fread("outTable/chm_buff_raw.csv")  # get values per pixel

## read pre-disturbance site history: rasters
pre_dist_history <- fread("outTable/pre_disturb_history_raster_ALL_buffers.csv")


# pre-disturbance history: individual trees 
convex_hull <- vect("raw/pre-disturb_history_trees/convex_hull_full_5514_buffer.gpkg")  # replace with correct path
pre_trees  <- vect("raw/pre-disturb_history_trees/pre-disturbance_trees.gpkg")     # replace with correct path


## data clean up -------------------------------------------------


# rename to 'subplots'
drone_cv <- drone_cv %>% 
  rename(subplot = plot_ID) %>% 
  mutate(plot = str_replace(subplot, "^[^_]+_([^_]+_[^_]+)_.*$", "\\1"),
         pixel_value = pmax(pixel_value, 0)   # clamp negatives to 0)
  )


### pre-disturbance cluster: tree density, % share of spruce  
# Reproject to EPSG:3035 (ETRS89-LAEA)
convex_hull_3035 <- project(convex_hull, "EPSG:3035")
pre_trees_3035   <- project(pre_trees, "EPSG:3035")

# Get polygon area (in m²) & perimeter - this is convex hull + buffer around it to avoid edge effects
convex_hull_3035$cvx_area_m2 <- expanse(convex_hull_3035, unit = "m")
convex_hull_3035$cvx_perimeter_m <- perim(convex_hull_3035)

# Add plot info to each tree:  (spatial join)
pre_trees_3035_joined <- terra::intersect(pre_trees_3035, convex_hull_3035)

#  Clean up tree- and convex hull input data
clean_trees <- pre_trees_3035_joined %>%
  as.data.frame() %>%
  mutate(
    original_species = species,
    species = case_when(
      is.na(species) ~ "piab",
      species == "l" ~ "deciduous",
      species == "s" ~ "piab",
      TRUE ~ species
    ),
    status = case_when(
      original_species == "s" ~ "dry",
      TRUE ~ "living"
    ),
    # Split year and note from `rok_disturbancia`
    disturbance_year = case_when(
      rok_disturbancia == "nie" ~ 2024L,
      TRUE ~ as.integer(str_extract(rok_disturbancia, "^\\d{4}"))
    ),
    disturbance_note = case_when(
      rok_disturbancia == "nie" ~ "nie",
      TRUE ~ str_trim(str_remove(rok_disturbancia, "^\\d{4}[- ]*"))
    ),
    # Parse rok_les safely and calculate disturbance length
    disturbance_length = disturbance_year - as.integer(rok_les)
  ) %>%
  dplyr::select(-original_species) %>%
  rename(
    forest_year = rok_les,           
    plot = cluster,
    buff_area = area,
    buff_perimeter = perimeter
  )

# calculate stem density from tree counst (based on mapped area - convex hull)
pre_disturb_stem_density_full <- clean_trees %>% 
  group_by(plot, species, status, 
           cvx_area_m2, 
           disturbance_year,
           forest_year,
           disturbance_length) %>%
  summarise(n_trees_sp_vitality = n(), .groups = "drop") %>%
  group_by(plot) %>%
  mutate(
    n_trees_total = sum(n_trees_sp_vitality),
    n_trees_species = ave(n_trees_sp_vitality, species, FUN = sum)
  ) %>%
  ungroup() %>% 
  # calculate stem density
  mutate(
    cvx_area_ha = cvx_area_m2 / 10000,
    density_ha = round(n_trees_total / cvx_area_ha, 1),
    density_species_ha = round(n_trees_species / cvx_area_ha, 1),
    share_spruce = ifelse(species == "piab", round(n_trees_species / n_trees_total, 2), NA_real_)
  )

# filter only relevant information: per plot level
pre_disturb_stem_density_filt <- pre_disturb_stem_density_full %>% 
  dplyr::filter(species == 'piab') %>% 
  dplyr::select(plot,cvx_area_m2, 
                disturbance_year, 
                forest_year,
                disturbance_length,
                n_trees_total, density_ha,share_spruce) %>% 
  distinct() %>% 
  rename_with(~ paste0("pre_dst_", .x), -c(plot))


### pre-disturbance raster -----------------

pre_dist_history <- pre_dist_history %>%
  mutate(field_year = if_else(str_detect(cluster, "_"), 2023, 2025))

pre_dist_history_2023_rst <- pre_dist_history %>%   # filed only pre dicturbancs history for sites collected in 2023
  dplyr::filter(field_year == "2023",
                buffer_m %in% c(25, 100)) %>% # filter only two buffer widths 
  rename(plot = cluster)

# convert to wide format 
pre_dist_wide <- pre_dist_history_2023_rst %>%
  dplyr::select(-ID, -field_year) %>%                 # drop unneeded cols
  mutate(buff = paste0("b", buffer_m)) %>%     # label buffers b25/b100
  pivot_wider(
    id_cols   = c(plot, year),
    names_from = buff,
    values_from = c(
      mean_cover_dens, median_cover_dens, cv_cover_dens,
      total_cells, n_coniferous, n_deciduous, forest_cells,
      share_coniferous, share_deciduous, shannon
    ),
    names_glue = "{.value}_{buff}",            # e.g., mean_cover_dens_b25
    values_fn = dplyr::first                   # in case of duplicates
  ) %>%
  arrange(plot) %>% 
  rename_with(~ paste0("pre_dst_", .x), -c(plot, year))


### Field data: subplot level --------------------------------------



# guestimate dbh and ba per individual based on height distribution 
dat23_subplot_recode <- dat23_subplot %>% 
  dplyr::select(-country) %>% 
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
    ba_ha_m2       = ba_total_m2 * scaling_factor
  ) %>% 
  rename(plot = cluster, 
         subplot = ID) %>% 
  mutate(n = ifelse(is.na(n), 0L, n),
         stem_density = ifelse(is.na(stem_density), 0L, stem_density))

  

# define recovery type: based on species that are likely planted/pioneer
dat23_subplot_recode <- dat23_subplot_recode %>%
  mutate(recovery_type = case_when(
    # --- Planted / late-successional / non-native ---
    species %in% c("piab","pisy","abal","lade","psme","taba",
                   "fasy","quro","acca","acpl","acps","frex",
                   "casa","aehi","saca","rops") ~ "planted",
    
    # --- Pioneer / early successional ---
    species %in% c("besp","alin","algl","alvi","potr","posp","prav",
                   "tisp","soau","soto","soar","cabe","ulsp","aial",
                   "fror","juni","jure","qusp","sasp","osca") ~ "pioneer",
    
    # --- Everything else / not clearly one of the two ---
    TRUE ~ "other"
  ))

# get stem density by vertical class
plot_density <- dat23_subplot_recode %>%
  group_by(plot, vegtype) %>%
  summarise(stem_density_sum = sum(stem_density, na.rm = TRUE), 
            .groups = "drop") %>% 
  pivot_wider(
    names_from = vegtype,
    values_from = stem_density_sum,
    values_fill = 0
  ) %>% 
  mutate(stem_regeneration = advanced + small)

# get information about context: salvage and protection intensity
str(dat23_subplot_recode)


hist(dat23_subplot_recode$clear)


# Get context information: management intensity at plot level 


# --- 1) Columns to use & NA -> 0 ---
management_types_v <- c("clear", "grndwrk", "logging_trail", "planting", "anti_browsing")

dat23_subplot_recode <- dat23_subplot_recode %>%
  mutate(across(all_of(management_types_v), ~ ifelse(is.na(.), 0, .))) 

# --- 2) Subplot-level scores (use max within subplot to avoid duplicates) ---
mng_subplot_scores <- dat23_subplot_recode %>%
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

# --- 3) Plot-level intensities (scaled 0–1) ---
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






# stem density by species group --------------
# need to get plot level average! 
subplot_group_density_wide <- dat23_subplot_recode %>%
  filter(vegtype %in% c("small", "advanced")) %>%
  group_by(subplot, plot, recovery_type, vegtype) %>%
  summarise(stem_density_sum = sum(stem_density, na.rm = TRUE),
            .groups = "drop") %>%
  unite("group_veg", recovery_type, vegtype) %>%   # combine recovery_type + vegtype
  pivot_wider(
    names_from = group_veg,
    values_from = stem_density_sum,
    values_fill = 0
  ) %>%
  mutate(
    pioneer_total = coalesce(pioneer_small, 0) + coalesce(pioneer_advanced, 0),
    planted_total = coalesce(planted_small, 0) + coalesce(planted_advanced, 0)
  )
head(subplot_group_density_wide)
View(subplot_group_density_wide)


# does the planted vs pioneer trees correlate?
subplot_group_density_wide %>% 
  ggplot(aes(x = pioneer_small , y = planted_small)) +
  geom_point() + 
  geom_smooth(method = "lm")

## summary metrics on subplot level ----------------------------
# drone
# field data

# test
drone_cv_sub <- drone_cv %>% 
  dplyr::filter(subplot == '13_26_105_5')

my_mean = mean(drone_cv_sub$pixel_value)
my_sd = sd(drone_cv_sub$pixel_value)

# Drone: subplot vs plot level 

# Per-subplot
drone_sub_summ <- drone_cv %>%
  mutate(pixel_value = pixel_value +0.001) %>% # add a small values to allow CV calculation
  group_by(plot, subplot, drone) %>%
  summarise(
    n_pix      = sum(!is.na(pixel_value)),
    mean_hgt   = mean(pixel_value, na.rm = TRUE),
    median_hgt = median(pixel_value, na.rm = TRUE),
    sd_hgt     = sd(pixel_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cv_hgt = sd_hgt / mean_hgt)

# Plot-level
drone_plot_summ <- drone_cv %>%
  group_by(plot) %>%
  summarise(
    n_pix      = sum(!is.na(pixel_value)),
    mean_hgt   = mean(pixel_value, na.rm = TRUE),
    median_hgt = median(pixel_value, na.rm = TRUE),
    sd_hgt     = sd(pixel_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cv_hgt =  sd_hgt / mean_hgt)



# --- Field data: Subplot-level metrics ---
field_sub_summ <- dat23_subplot_recode %>%
  group_by(plot, subplot) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE),
    sp_richness = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp  = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    mean_hgt    = if (stems_total > 0) weighted.mean(hgt_est, n, na.rm = TRUE) else NA_real_,
    median_hgt    = if (stems_total > 0) median(hgt_est, na.rm = TRUE) else NA_real_,
    
    # fix CV - if all subplots have not stems, CV is 0
    cv_hgt = if (stems_total > 0 && mean(hgt_est, na.rm = TRUE) > 0) {
      sd(hgt_est, na.rm = TRUE) / mean(hgt_est, na.rm = TRUE)
    } else 0,
    range_hgt   = if (stems_total > 0 && !all(is.na(hgt_est)))
      max(hgt_est, na.rm = TRUE) - min(hgt_est, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>% 
  mutate(cv_hgt = ifelse(is.na(cv_hgt), 0L, cv_hgt),
         mean_hgt = ifelse(is.na(mean_hgt), 0L, mean_hgt)) # replace NA by 0 if stems are missing



## Correlation of height structure between field& drone data ----------
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

# very low correlation between height structure between drones and field data!


# Cross scale variation ----------------------------------------------
##  variation of height classes across scales: subplot vs plot -----------------------

# vertical structure per subplot& plot:
# plot values are generated as 
# 1. average across plots 
# (compare how different the subplot is compared to plit mean )
# 2. as pooled all values (no subplots boundaries)


# --- Plot-level metrics (aggregate over subplots) ---
plot_metrics_mean <- field_sub_summ %>%
  group_by(plot) %>%
  summarise(
    mean_sp_richness = mean(sp_richness, na.rm = TRUE),
    #var_sp_richness  = var(sp_richness,  na.rm = TRUE),
    mean_shannon_sp  = mean(shannon_sp,  na.rm = TRUE),
    mean_evenness_sp = mean(evenness_sp, na.rm = TRUE),
    mean_mean_hgt    = mean(mean_hgt,    na.rm = TRUE),
    mean_cv_hgt      = mean(cv_hgt,      na.rm = TRUE),
    #var_cv_hgt       = var(cv_hgt,       na.rm = TRUE),
    mean_range_hgt   = mean(range_hgt,   na.rm = TRUE),
    .groups = "drop"
  )


# --- pooled CV directly from dat23_subplot_recode ---
plot_metrics_pooled  <- dat23_subplot_recode %>%
  group_by(plot) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE),
    sp_richness = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp  = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    mean_hgt    = if (stems_total > 0) weighted.mean(hgt_est, n, na.rm = TRUE) else NA_real_,
    cv_hgt      = if (stems_total > 0 && mean(hgt_est, na.rm = TRUE) > 0)
      sd(hgt_est, na.rm = TRUE) / mean(hgt_est, na.rm = TRUE) else NA_real_,
    range_hgt   = if (stems_total > 0 && !all(is.na(hgt_est)))
      max(hgt_est, na.rm = TRUE) - min(hgt_est, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>% 
  mutate(cv_hgt = ifelse(is.na(cv_hgt), 0L, cv_hgt)) # replace NA by 0 if stems are missing




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
                     pooled_cv_hgt      = cv_hgt,
                     pooled_range_hgt   = range_hgt),
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
  left_join(pre_dist_wide) %>%  # add rasters for rought estimation: density_cover and % coniferous (leaf type)
  left_join(pre_disturb_stem_density_filt) %>%  # stem density based on tree calculation
  left_join(plot_density) %>%   # get stem density by vertical classes
  left_join(mng_plot_intensity ) %>% # add context information
  #left_join(subplot_group_density_wide) %>%  # get stem density by trees likely planted/pioneer - need to group on plot level!!
  mutate(
    log_mean = log1p(mean_cv_hgt),
    log_resp = log1p(pooled_cv_hgt)
  ) 

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
