
# ----------------------------------------------------------------------------
# Cross-scale interaction: field vs drone variation in vertical structures
# -----------------------------------------------------------------------------

# read data: 
# - field data from 2023 (2025)
# - drone 2023 (2024?)

# - get sites history - from tree density cover, dominant leaf type
# compare stem density, dominant tree species, species diversity
# vertical and horizontal structure across sscales: fiedl vs drone

# identify damaged terminal: from raw data czechia 20203 - export gpkg
# identify types of damages: 2025 - export gpkg
# Drones: 
# - get vertical strucure

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
library(tidyr)


# Read files --------------------------------------

# 2023 
dat23_subplot <- fread('outData/subplot_full_2023.csv')
dat23_cluster <- fread('outData/df_cluster_2023.csv')

# Save spatial subplot data with cluster IDs
dat23_sf <- st_read("outData/sf_context_2023.gpkg")


# 2025
#dat25_subplot <- fread('outData/subplot_full_2025.csv')
#dat25_cluster <- fread('outData/df_cluster_2025.csv')

# Save spatial subplot data with cluster IDs
#dat25_sf <- st_read("outData/subplot_with_clusters_2025.gpkg")

# Read drone CHM summary
drone_cv <- fread("outTable/chm_buff_summary.csv")

# read pre-disturbance site history
pre_dist_history <- fread("outTable/pre_disturb_history_raster.csv")

pre_dist_history <- pre_dist_history %>%
  mutate(field_year = if_else(str_detect(cluster, "_"), 2023, 2025))

pre_dist_history_2023 <- pre_dist_history %>%   # filed only pre dicturbancs history for sites collected in 2023
  dplyr::filter(field_year == "2023")


# Sumarize field observation per cluster  --------------------------------------

# guestimate dbh and ba per individual based on height distribution -------------------------
dat23_subplot2 <- dat23_subplot %>% 
  mutate(
    dbh_est = case_when(
      vegtype == "mature"           ~ case_when(
        dbh == "10–20cm" ~ 15.0,
        dbh == "20–40cm" ~ 30.0,
        dbh == "40–60cm" ~ 50.0,
        TRUE             ~ NA_real_
      ),
      hgt == "0.2–0.4"              ~ 0.2,
      hgt == "0.4–0.6"              ~ 0.4,
      hgt == "0.6–0.8"              ~ 0.7,
      hgt == "0.8–1.0"              ~ 1.0,
      hgt == "1.0–1.3"              ~ 1.5,
      hgt == "1.3–2.0"              ~ 2.5,
      hgt == "2–4"                  ~ 4.5,
      hgt == ">4"                   ~ 7.0,
      TRUE                          ~ NA_real_
    ),
    basal_area_cm2 = pi * (dbh_est / 2)^2
  )




# get stem density and shannon for field data 
df_stem_density <- dat23_subplot %>% 
  group_by(cluster) %>% 
  summarise(stem_density = sum(stem_density, na.rm = T))

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

# Define Shannon + Effective Species Number
shannon_stats <- function(ba) {
  total <- sum(ba)
  if (total == 0) return(c(H = 0, ESN = 0))
  p <- ba / total
  H <- -sum(p * log(p), na.rm = TRUE) # shannon index
  ESN <- exp(H)  # effective number of species
  c(H = H, ESN = ESN)
}

# All clusters
all_clusters <- dd %>% distinct(cluster)

# Compute diversity only for rows with species info
shannon_summary <- dd %>%
  filter(!is.na(species)) %>%
  group_by(cluster) %>%
  reframe(
    shannon = shannon_stats(basal_area)[["H"]],
    effective_species = shannon_stats(basal_area)[["ESN"]]
  )

# Merge back to include missing clusters
shannon_per_cluster <- all_clusters %>%
  left_join(shannon_summary, by = "cluster") %>%
  mutate(
    shannon = replace_na(shannon, 0),
    effective_species = replace_na(effective_species, 0)
  )

print(shannon_per_cluster)

# summarize across species: get shannon per subplot, per plot
dat23_subplot_sum_cl <- dat23_subplot %>% 
  group_by( cluster, species ) %>% 
  summarise(stem_density = sum(stem_density, na.rm = T))
  

# calculate shannon index: --------------------------
# 1. Total stem density per cluster
total_stems <- dat23_subplot_sum_cl %>%
  group_by(cluster) %>%
  summarise(total_density = sum(stem_density), .groups = "drop")

# 2. Shannon index only for clusters with regeneration
shannon_nonzero <- dat23_subplot_sum_cl %>%
  left_join(total_stems, by = "cluster") %>%
  dplyr::filter(total_density > 0) %>%
  mutate(p = stem_density / total_density) %>%
  group_by(cluster) %>%
  summarise(shannon_index = -sum(p * log(p)), .groups = "drop")

# 3. Join with all clusters, assigning 0 to those with no regeneration
shannon_all <- total_stems %>%
  left_join(shannon_nonzero, by = "cluster") %>%
  mutate(shannon_index = if_else(is.na(shannon_index), 0, shannon_index))

# View result
head(shannon_all)
hist(shannon_all$shannon_index)
