

# process data: from 2023
# get better inidcation of heights variabuility
# terminal damages
# get spatial attributes


library(DBI)
library(RSQLite)
library(tidyverse)
library(data.table)

# # read the cleaned df:
# 
# dat <- fread("raw/collected_2023/cleaned_data_frame.csv")
# 
# unique(dat$country)
# unique(dat$region)
# 
# dat_cz <- dat %>% 
#   dplyr::filter(region %in% c(15, 26)) %>% 
#   mutate(group = group +100) %>%  
#   mutate(cluster_ID_corr =  paste(region, group, '_'))
#  # nrow()
# 
# # 680 rows
# length(unique(dat_cz$region))
# length(unique(dat_cz$ID))
# length(unique(dat_cz$cluster_ID_corr))
# 

# Paths to the uploaded GPKG files---------------------------------------------------
gpkg1 <- "raw/collected_2023/cz_final_2023.gpkg"
gpkg2 <- "raw/collected_2023/cz4_final_2023.gpkg"

# Connect and list tables
con1 <- dbConnect(RSQLite::SQLite(), gpkg1)
con2 <- dbConnect(RSQLite::SQLite(), gpkg2)

tables1 <- dbListTables(con1)
tables2 <- dbListTables(con2)

print(tables1)
print(tables2)

#print(tables2)
cz1 <- dbReadTable(con1, "cz_final")
cz2 <- dbReadTable(con2, "cz4_final")

# test: convert to a table:---------

# cz1_r <- cz1 %>%
#   dplyr::select(fid, group, point, region, country,
#                 matches("^r_[a-z0-9]+_(n|hgt|dmg)$")) %>%
#   mutate(across(matches("^r_[a-z0-9]+_(n|hgt|dmg)$"), as.character)) %>%
#   pivot_longer(
#     cols = matches("^r_[a-z0-9]+_(n|hgt|dmg)$"),
#     names_to = "variable",
#     values_to = "value"
#   ) %>%
#   filter(!is.na(value)) %>%
#   extract(variable, into = c("layer", "species", "trait"),
#           regex = "^(r)_([a-z0-9]+)_(n|hgt|dmg)$") %>%
#   mutate(value = as.numeric(value))
# 

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
  #  filter(!is.na(value)) %>%
    extract(variable, into = c("vegtype", "species", "trait"),
            regex = paste0("^(", prefix, ")_([a-z0-9]+)_(", paste(traits, collapse = "|"), ")$")) %>%
    mutate(value = as.numeric(value))
}

## process atble 1 ----------------------------
cz1_r <- reshape_layer(cz1, "r", traits = c("n", "hgt", "dmg")) %>%
  mutate(vegtype = "small")

cz1_ar <- reshape_layer(cz1, "ar", traits = c("n", "hgt","dmg")) %>%
  mutate(vegtype = "advanced")

cz1_s <- reshape_layer(cz1, "s", traits = c("n", "dbh")) %>%
  mutate(vegtype = "mature")


# merge vegetationb data together
cz1_all <- bind_rows(cz1_r, cz1_ar, cz1_s)

# get a proper naming:
cz1_all <- cz1_all %>% 
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(country, region, group, point, sep = "_"))

## process table 2 ----------------------------  
# Process with the same reshape_layer function
cz2_r <- reshape_layer(cz2, "r", traits = c("n", "hgt", "dmg")) %>%
  mutate(vegtype = "small")
cz2_ar <- reshape_layer(cz2, "ar", traits = c("n", "hgt", "dmg")) %>%
  mutate(vegtype = "advanced")
cz2_s <- reshape_layer(cz2, "s", traits = c("n", "dbh")) %>%
  mutate(vegtype = "mature")

# Merge all layers
cz2_all <- bind_rows(cz2_r, cz2_ar, cz2_s)

# Add unique ID and harmonize group ID
cz2_all <- cz2_all %>% 
  dplyr::mutate(group = group + 100) %>%
  dplyr::mutate(ID = paste(country, region, group, point, sep = "_"))

## Merge both datbles: -----------------------------
head(cz2_all)

cz_full <- bind_rows(cz2_all,
                     cz1_all)

length(unique(cz_full$ID))


# test on different data:  ------------------------------------------

# Apply reshaping for each vegetation layer from dat_cz

cz1_r <- reshape_layer(dat_cz, "r", traits = c("n", "hgt", "dmg")) %>%
  mutate(vegtype = "small")

cz1_ar <- reshape_layer(dat_cz, "ar", traits = c("n", "hgt", "dmg")) %>%
  mutate(vegtype = "advanced")

cz1_s <- reshape_layer(dat_cz, "s", traits = c("n", "dbh")) %>%
  mutate(vegtype = "mature")

# Merge all layers
cz1_all <- bind_rows(cz1_r, cz1_ar, cz1_s)

# Create unique ID
cz1_all <- cz1_all %>%
  dplyr::mutate(group = group + 100) %>%
  dplyr::mutate(ID = paste(country, region, group, point, sep = "_"))
