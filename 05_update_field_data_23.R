

# process data: from 2023
# get better inidcation of heights variabuility
# terminal damages
# get spatial attributes


library(DBI)
library(RSQLite)
library(tidyverse)

# Paths to the uploaded GPKG files
gpkg1 <- "raw/collected_2023/cz_final_2023.gpkg"
#gpkg2 <- "/mnt/data/cz4_final_2023.gpkg"

# Connect and list tables
con1 <- dbConnect(RSQLite::SQLite(), gpkg1)
#con2 <- dbConnect(RSQLite::SQLite(), gpkg2)

tables1 <- dbListTables(con1)
#tables2 <- dbListTables(con2)

print(tables1)

#print(tables2)
cz_final_df <- dbReadTable(con1, "cz_final")

# test: convert to a table:

cz_long_r <- cz_final_df %>%
  dplyr::select(fid, group, point, region, country,
                matches("^r_[a-z0-9]+_(n|hgt|dmg)$")) %>%
  mutate(across(matches("^r_[a-z0-9]+_(n|hgt|dmg)$"), as.character)) %>%
  pivot_longer(
    cols = matches("^r_[a-z0-9]+_(n|hgt|dmg)$"),
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  extract(variable, into = c("layer", "species", "trait"),
          regex = "^(r)_([a-z0-9]+)_(n|hgt|dmg)$") %>%
  mutate(value = as.numeric(value))


# make a function:

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
    filter(!is.na(value)) %>%
    extract(variable, into = c("vegtype", "species", "trait"),
            regex = paste0("^(", prefix, ")_([a-z0-9]+)_(", paste(traits, collapse = "|"), ")$")) %>%
    mutate(value = as.numeric(value))
}

# apply for each vegetation layer ----------------------------
cz_long_r <- reshape_layer(cz_final_df, "r", traits = c("n", "hgt", "dmg")) %>%
  mutate(vegtype = "small")

cz_long_ar <- reshape_layer(cz_final_df, "ar", traits = c("n", "hgt","dmg")) %>%
  mutate(vegtype = "advanced")

cz_long_s <- reshape_layer(cz_final_df, "s", traits = c("n", "dbh")) %>%
  mutate(vegtype = "mature")


# merge vegetationb data together

cz_long_all <- bind_rows(cz_long_r, cz_long_ar, cz_long_s)

# get a proper naming:
cz_long_all <- cz_long_all %>% 
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(region, group, point, sep="_")) #%>%  #dplyr::select(6:323) %>% names()
  
