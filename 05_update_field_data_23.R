

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
