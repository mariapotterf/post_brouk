

# process data: from 2023
# get better inidcation of heights variabuility
# terminal damages
# get spatial attributes

gc()

library(DBI)
library(RSQLite)
library(tidyverse)
library(data.table)


# Helper: reshaping vegetation structure
reshape_layer <- function(df, prefix, traits = c("n", "hgt", "dmg", "dbh")) {
  df %>%
    dplyr::select(fid, group, point, region, country,
                  matches(paste0("^", prefix, "_[a-z0-9]+_(", paste(traits, collapse = "|"), ")$"))) %>%
    mutate(across(matches(paste0("^", prefix, "_[a-z0-9]+_(", paste(traits, collapse = "|"), ")$")), as.character)) %>%
    pivot_longer(
      cols = matches(paste0("^", prefix, "_[a-z0-9]+_(", paste(traits, collapse = "|"), ")$")),
      names_to = "variable",
      values_to = "value"
    ) %>%
    extract(variable, into = c("vegtype", "species", "trait"),
            regex = paste0("^(", prefix, ")_([a-z0-9]+)_(", paste(traits, collapse = "|"), ")$")) %>%
    mutate(value = as.numeric(value))
}

# Helper: full processing pipeline per table
process_table <- function(path, table_name) {
  con <- dbConnect(RSQLite::SQLite(), path)
  df <- dbReadTable(con, table_name)
  dbDisconnect(con)
  
  # Build metadata once
  meta <- df %>%
    dplyr::select(fid, stump, stump_spc, clear, grndwrk, logging_trail,
                  stump_otsp, windthrow, planting, deadwood, anti_browsing)
  
  # Reshape and join
  veg_long <- bind_rows(
    reshape_layer(df, "r", traits = c("n", "hgt", "dmg")) %>% mutate(vegtype = "small"),
    reshape_layer(df, "ar", traits = c("n", "hgt", "dmg")) %>% mutate(vegtype = "advanced"),
    reshape_layer(df, "s", traits = c("n", "dbh")) %>% mutate(vegtype = "mature")
  ) %>%
    left_join(meta, by = "fid") %>%
    mutate(group = group + 100,
           ID = paste(country, region, group, point, sep = "_"),
           cluster = paste(region, group, sep = "_"))
  
  return(veg_long)
}


# Process both GPKG tables
cz1_all <- process_table("raw/collected_2023/cz_final_2023.gpkg", "cz_final")
cz2_all <- process_table("raw/collected_2023/cz4_final_2023.gpkg", "cz4_final")

# Combine into one long dataframe
cz_full <- bind_rows(cz1_all, cz2_all)

# Summary check
length(unique(cz_full$ID))


cz_wide <- cz_full %>%
  pivot_wider(
    names_from = trait,
    values_from = value
  )


head(cz_wide)

# recode the metrix: ---------------------------

# helper to recode 1 = NA, 2 = 1 (TRUE), 3 = 0 (FALSE)
recode_binary <- function(x) {
  x <- as.numeric(x)
  case_when(
    x == 1 ~ NA_real_,
    x == 2 ~ 1,
    x == 3 ~ 0,
    TRUE ~ NA_real_
  )
}


cz_wide2 <- cz_wide %>%
 # dplyr::filter(dist == TRUE) %>%
  mutate(
    planting      = recode_binary(planting),
    anti_browsing = recode_binary(anti_browsing),
    windthrow     = recode_binary(windthrow),
    deadwood      = recode_binary(deadwood),
    logging_trail = ifelse(logging_trail == 1, 1, 0),
    clear         = ifelse(clear == 1, 1, 0),
    grndwrk       = ifelse(grndwrk == 1, 1, 0)
  ) %>%
  mutate(
    manag_intensity      = logging_trail + clear + grndwrk + planting + anti_browsing,
    salvage_intensity    = logging_trail + clear + grndwrk,
    protection_intensity = planting + anti_browsing
  )

# check how counts are looking?
table(cz_wide2$n)

# need to change categorie!::

# because: 1 = >16, 2 = 1, 3 = 2,...17 = 16

cz_wide2 <- cz_wide2 %>%
  mutate(
    n_corr = case_when(
      vegtype %in% c("small", "advanced") & n == 1 ~ 17,
      vegtype %in% c("small", "advanced")         ~ n-1,
     # vegtype == "mature" & n == 5                ~ ">4",
     # vegtype == "mature"                         ~ as.character(n),
      TRUE                                        ~ NA_integer_
    )
  ) %>% 
  group_by(cluster) %>%
  mutate(n_plots = n_distinct(ID)) %>%
  ungroup() %>% 
  mutate(area = n_plots *4,
         scaling_factor = 10000/area,
         stem_density = n_corr*scaling_factor) 


cz_wide2_cluster <- cz_wide2 %>% 
  group_by(species, cluster) %>% 
  summarize(stem_density = sum(stem_density, na.rm = T))

