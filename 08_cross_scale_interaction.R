

#             Cross-scale interaction: 
#    field & drone variation in vertical structures 


# read data: 
# - field data from 2023 (2025?)
# - drone 2023&2024

# - get sites history - from tree density cover, dominant leaf type
#                     - from tree density (manually mapped)

# correlate: vertical structure and composition across scales: field vs drone
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


# Read files: only field 2023!  --------------
# --- Field data: 2023 
dat23_subplot    <- data.table::fread("outData/subplot_full_2023.csv")   # subplot-level table
dat23_sf         <- sf::st_read("outData/sf_context_2023.gpkg")          # subplot spatial data

# --- Drone CHM data 
drone_cv         <- data.table::fread("outTable/chm_buff_raw.csv")            # pixel-level values

# --- Pre-disturbance history - Raster-based history
pre_dist_history <- data.table::fread("outTable/pre_disturb_history_raster_ALL_buffers.csv")

# --- Tree-based history (vector layers)
convex_hull      <- terra::vect("raw/pre-disturb_history_trees/convex_hull_full_5514_buffer.gpkg")
pre_trees        <- terra::vect("raw/pre-disturb_history_trees/pre-disturbance_trees.gpkg")

# --- Buffer polygons for study sites 
buff_square      <- terra::vect("outData/square_buffers_2m_all.gpkg")

## Data clean up -------------------------------------------------

# rename to 'subplots'
drone_cv <- drone_cv %>% 
  rename(subplot = plot_ID) %>% 
  mutate(plot = str_replace(subplot, "^[^_]+_([^_]+_[^_]+)_.*$", "\\1"),
         pixel_value = pmax(pixel_value, 0)   # clamp negatives to 0)
  )


### pre-disturbance cluster: tree density, % share of spruce, Reproject to EPSG:3035 (ETRS89-LAEA)
convex_hull_3035 <- project(convex_hull, "EPSG:3035")
pre_trees_3035   <- project(pre_trees, "EPSG:3035")
buff_square_3035 <- project(buff_square, "EPSG:3035")

# Get polygon area (in m²) & perimeter - this is convex hull + buffer around it to avoid edge effects
convex_hull_3035$area_m2     <- expanse(convex_hull_3035, unit = "m")
convex_hull_3035$perimeter_m <- perim(convex_hull_3035)

buff_square_3035$area_m2      <- expanse(buff_square_3035, unit = "m")
buff_square_3035$perimeter_m  <- perim(buff_square_3035)


# clean up tree characteristics: get species, ..
pre_trees_df <- as.data.frame(pre_trees_3035) %>%
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
    )) %>%
  dplyr::select(-original_species) 

# Reattach cleaned attributes back to geometry
pre_trees_3035_clean          <- pre_trees_3035
values(pre_trees_3035_clean)  <- pre_trees_df

## Clean up convex hull characteristics -------------------------
cvx_clean <- convex_hull_3035 %>% 
  as.data.frame() %>%   # attributes only (no geometry)
  mutate(
    # rok_disturbancia like "2019-…", or "nie"
    disturbance_year = case_when(
      rok_disturbancia == "nie" ~ 2024L,
      TRUE ~ suppressWarnings(as.integer(str_extract(rok_disturbancia, "^\\d{4}")))
    ),
    disturbance_note = case_when(
      rok_disturbancia == "nie" ~ "nie",
      TRUE ~ str_trim(str_remove(rok_disturbancia, "^\\d{4}[- ]*"))
    ),
    forest_year = suppressWarnings(as.integer(rok_les)),
    disturbance_length = disturbance_year - forest_year
  ) %>% 
  rename(plot = cluster) %>% 
  dplyr::select(-rok_les,-rok_disturbancia, -area, -perimeter,
                - disturbance_note)

# Write cleaned attributes back to the terra object
convex_hull_3035_clean  <-convex_hull_3035[, c("cluster")]
values(convex_hull_3035_clean)  <- cvx_clean


# Add subplot (square)/plot (convex hull) info to each tree  (spatial join)
pre_trees_cvx_joined <- terra::intersect(pre_trees_3035_clean, convex_hull_3035_clean)
# add convex hull information also to subplots
#buff_squared_joined  <- terra::intersect(buff_square_3035, convex_hull_3035_clean)
pre_trees_sq_joined  <- terra::intersect(pre_trees_3035_clean, buff_square_3035)

### Get stem density per Run once for convex hulls, once for squares  -------------------
cvx_df <- as.data.frame(pre_trees_cvx_joined) 
sq_df  <- as.data.frame(pre_trees_sq_joined) 

####  Plot = CVX: stem density ------------------------------
cvx_stem_density <- cvx_df %>%
#  ungroup(.) %>% 
  group_by(plot) %>%  #, area_m2
  dplyr::reframe(
    n_trees = n(),
    area_m2 = mean(area_m2, na.rm  = T),  # keep are here instead of grouping
    density_ha = n_trees / area_m2 * 10000
  )

cxv_stem_density_species <- cvx_df %>%
  group_by(plot, species) %>%
  dplyr::reframe(
    n_trees = n(),
    area_m2 = mean(area_m2, na.rm  = T),  # keep are here instead of grouping
    density_ha = n_trees / area_m2 * 10000
  )

cvx_stem_density_status <- cvx_df %>%
  group_by(plot, species, status) %>%
  summarise(
    n_trees = n(),
    area_m2 = mean(area_m2, na.rm  = T),  # keep are here instead of grouping
    density_sp_status = n_trees / area_m2 * 10000,
    .groups = "drop"
  )



#### Density per square ------------------------------

# Plot = CVX: stem density
sq_stem_density <- sq_df %>%
  #  ungroup(.) %>% 
  group_by(plot_ID) %>%  #, area_m2
  dplyr::reframe(
    n_trees = n(),
    area_m2 = mean(area_m2, na.rm  = T),  # keep are here instead of grouping
    density_ha = n_trees / area_m2 * 10000
  ) %>% 
  rename(subplot = plot_ID)

sq_stem_density_species <- sq_df %>%
  group_by(plot_ID, species) %>%
  dplyr::reframe(
    n_trees = n(),
    area_m2 = mean(area_m2, na.rm  = T),  # keep are here instead of grouping
    density_ha = n_trees / area_m2 * 10000
  ) %>% 
  rename(subplot = plot_ID)

sq_stem_density_status <- sq_df %>%
  group_by(plot_ID, species, status) %>%
  summarise(
    n_trees = n(),
    area_m2 = mean(area_m2, na.rm  = T),  # keep are here instead of grouping
    density_sp_status = n_trees / area_m2 * 10000,
    .groups = "drop"
  ) %>% 
  rename(subplot = plot_ID)


head(sq_stem_density_species)
head(cxv_stem_density_species)

### Pre-disturbance raster -----------------

pre_dist_history <- pre_dist_history %>%
  mutate(field_year = if_else(str_detect(cluster, "_"), 2023, 2025))

pre_dist_history_2023_rst <- pre_dist_history %>%   # filed only pre dicturbancs history for sites collected in 2023
  dplyr::filter(field_year == "2023",
                buffer_m %in% c(25, 100)) %>% # filter only two buffer widths 
  rename(plot = cluster)

# convert to wide format 
pre_dist_2023_wide <- pre_dist_history_2023_rst %>%
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
  mutate(stem_regeneration = advanced + small,
         sum_density = advanced + small + mature)

# get information about context: salvage and protection intensity
str(dat23_subplot_recode)


#hist(dat23_subplot_recode$clear)


### Get context information: management intensity at plot level ----------------------------
# --- 1) Columns to use & NA -> 0 
management_types_v <- c("clear", "grndwrk", "logging_trail", "planting", "anti_browsing")

dat23_subplot_recode <- dat23_subplot_recode %>%
  mutate(across(all_of(management_types_v), ~ ifelse(is.na(.), 0, .))) 

# --- 2) Subplot-level scores (use max within subplot to avoid duplicates)
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

###  Get Regeneration stem density based on height group: small, advanced, mature ------------- 
# need to get plot level average! 
subplot_group_density_wide <- dat23_subplot_recode %>%
  dplyr::filter(vegtype %in% c("small", "advanced")) %>%
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
#head(subplot_group_density_wide)
#View(subplot_group_density_wide)



## Correlation: field vs drone on subplot level ---------------------------
# drone
# field data


drone_cv_sub <- drone_cv %>% 
  dplyr::filter(subplot == '13_26_105_5')

#my_mean = mean(drone_cv_sub$pixel_value)
#my_sd = sd(drone_cv_sub$pixel_value)

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



# --- Field data: Subplot-level metrics 
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



# summarize subplot information 
# how many trees?
sum(field_sub_summ$stems_total)  # 2125 stems/7816 stems/ha
sum(field_sub_summ$stems_total == 0)

# how many tree species? 
length(unique(dat23_subplot_recode$species[dat23_subplot_recode$n > 0]))

sort(unique(dat23_subplot_recode$species[dat23_subplot_recode$n > 0]))




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



# TOY example !!! -------------------
#library(dplyr)
# library(tidyr)

# Variance decomposition: within-subplot, between subplots (within plots), between-plots (landscape)

#library(dplyr)

# --- Toy dataset: two plots -----------------------------------
dd <- tibble::tibble(
  plot     = c("p1","p1","p1","p1","p1",  "p2","p2","p2"),
  subplot  = c("s1","s1","s2","s2","s3",  "s1","s2","s3"),
  species  = c("a","b","c","a", NA,       "a","b","c"),
  counts   = c(5,2,1,1,NA,                 3,2,1),
  height_m = c(5,0.2,3,1,NA,               4,2,6)
)

# --- Step 1: subplot means and within-subplot SS --------------
sub_height <- dd %>%
  group_by(plot, subplot) %>%
  summarise(
    N_sub = sum(counts, na.rm=TRUE),
    mean_h = ifelse(N_sub > 0,
                    stats::weighted.mean(height_m, w=counts, na.rm=TRUE),
                    NA_real_),
    SS_within_sub = ifelse(N_sub > 1,
                           sum(counts * (height_m - mean_h)^2, na.rm=TRUE),
                           0),
    .groups="drop"
  )

# --- Step 2: plot-level decomposition -------------------------
plot_stats <- sub_height %>%
  group_by(plot) %>%
  summarise(
    N_tot = sum(N_sub),
    mu_plot = ifelse(N_tot > 0,
                     stats::weighted.mean(mean_h, w=N_sub, na.rm=TRUE),
                     NA_real_),
    SS_within = sum(SS_within_sub),
    SS_between_sub = sum(N_sub * (mean_h - mu_plot)^2, na.rm=TRUE),
    .groups="drop"
  )

# --- Step 3: landscape-level decomposition -------------------
# grand mean across all trees in all plots
mu_grand <- weighted.mean(plot_stats$mu_plot, 
                          w=plot_stats$N_tot, 
                          na.rm=TRUE)

SS_within_all  <- sum(plot_stats$SS_within)
SS_bsub_all    <- sum(plot_stats$SS_between_sub)
SS_bplots_all  <- sum(plot_stats$N_tot * (plot_stats$mu_plot - mu_grand)^2)

N_all <- sum(plot_stats$N_tot)

var_within_all <- SS_within_all  / (N_all - 1)
var_bsub_all   <- SS_bsub_all    / (N_all - 1)
var_bplots_all <- SS_bplots_all  / (N_all - 1)
var_total_all  <- (SS_within_all + SS_bsub_all + SS_bplots_all) / (N_all - 1)

cv_within_all <- sqrt(var_within_all)/mu_grand
cv_bsub_all   <- sqrt(var_bsub_all)/mu_grand
cv_bplots_all <- sqrt(var_bplots_all)/mu_grand
cv_total_all  <- sqrt(var_total_all)/mu_grand

landscape_stats <- tibble::tibble(
  mu_grand,
  cv_within_all,
  cv_bsub_all,
  cv_bplots_all,
  cv_total_all
)

landscape_stats







# TEST END -----------------------------------











## relationship between Shannon vs CV on subplot ------------------------------------------------------------
# CV of heights vs shannon
field_sub_summ %>% 
  dplyr::filter(shannon_sp >0) %>% 
  ggplot(aes(x = cv_hgt,
             y = shannon_sp)) +
  geom_jitter() + 
  geom_smooth(method = 'gam')





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
