

# process data: from 2023
# get better inidcation of heights variabuility
# terminal damages
# get spatial attributes
# process raw df tables, and also as sf data

gc()

library(DBI)
library(RSQLite)
library(tidyverse)
library(data.table)
library(sf)

# Read spatial layers
#cz1_sf <- read_sf("raw/collected_2023/cz_final_2023.gpkg", "cz_final")
#cz2_sf <- read_sf("raw/collected_2023/cz4_final_2023.gpkg", "cz4_final")

# ── 0. File inventory ──────────────────────────────────────────────────────────
# All country/region GPKGs (exclude the _2023 suffix ones - those are CZ-specific older files)
gpkg_files <- tibble(
  path = list.files("raw/collected_2023", pattern = "\\.gpkg$", full.names = TRUE)
) %>%
  filter(!str_detect(path, "_2023\\.gpkg$")) %>%          # drop cz_final_2023, cz4_final_2023
  filter(!str_detect(path, "dat_czechia_2023\\.gpkg$")) %>% # drop any other legacy files
  mutate(
    basename   = tools::file_path_sans_ext(basename(path)),
    table_name = basename
  )

print(gpkg_files)  # sanity check — should show 17 files


# sanity check
st_layers("raw/collected_2023/pdb_final2.gpkg")


# TEST START

# ── 1. Helper functions (unchanged from your CZ script) ───────────────────────

reshape_layer <- function(df, prefix, traits = c("n", "hgt", "dmg", "dbh")) {
  df %>%
    dplyr::select(fid, group, point, region, country,
                  matches(paste0("^", prefix, "_[a-z0-9]+_(", paste(traits, collapse = "|"), ")$"))) %>%
    mutate(across(matches(paste0("^", prefix, "_[a-z0-9]+_(", paste(traits, collapse = "|"), ")$")), as.character)) %>%
    pivot_longer(
      cols = matches(paste0("^", prefix, "_[a-z0-9]+_(", paste(traits, collapse = "|"), ")$")),
      names_to  = "variable",
      values_to = "value"
    ) %>%
    extract(variable, into = c("vegtype", "species", "trait"),
            regex = paste0("^(", prefix, ")_([a-z0-9]+)_(", paste(traits, collapse = "|"), ")$")) %>%
    mutate(value = as.numeric(value))
}

process_table <- function(path, table_name) {
  con <- dbConnect(RSQLite::SQLite(), path)
  df  <- dbReadTable(con, table_name)
  dbDisconnect(con)
  
  meta <- df %>%
    dplyr::select(fid, stump, stump_spc, clear, grndwrk, logging_trail,
                  stump_otsp, windthrow, planting, deadwood, anti_browsing)
  
  veg_long <- bind_rows(
    reshape_layer(df, "r",  traits = c("n", "hgt", "dmg")) %>% mutate(vegtype = "small"),
    reshape_layer(df, "ar", traits = c("n", "hgt", "dmg")) %>% mutate(vegtype = "advanced"),
    reshape_layer(df, "s",  traits = c("n", "dbh"))        %>% mutate(vegtype = "mature")
  ) %>%
    left_join(meta, by = "fid") %>%
    mutate(
      group   = group + 100,
      ID      = paste(country, region, group, point, sep = "_"),
      cluster = paste(region, group, sep = "_")
    )
  
  return(veg_long)
}


# ── 2. Recode helpers ─────────────────────────────────────────────────────────

recode_binary <- function(x) {
  x <- as.numeric(x)
  case_when(
    x == 1 ~ NA_real_,
    x == 2 ~ 1,
    x == 3 ~ 0,
    TRUE   ~ NA_real_
  )
}

stump_species_lookup <- c(
  "deciduous", "conifer", "abal", "acca", "acpl", "acps", "aehi", "aial", "algl", "alin",
  "alvi", "besp", "cabe", "casa", "fasy", "frex", "juni", "jure", "lade", "piab",
  "pisy", "potr", "posp", "prav", "psme", "quro", "qusp", "rops", "saca", "sasp",
  "soar", "soau", "soto", "taba", "tisp", "ulsp", "otsp_1"
)

apply_recoding <- function(df) {
  df %>%
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
    ) %>%
    mutate(
      n = case_when(
        vegtype %in% c("small", "advanced") & n == 1 ~ 17,
        vegtype %in% c("small", "advanced")          ~ n - 1,
        TRUE                                         ~ NA_integer_
      ),
      hgt = case_when(
        vegtype == "small"    & hgt == 1 ~ "0.2–0.4",
        vegtype == "small"    & hgt == 2 ~ "0.4–0.6",
        vegtype == "small"    & hgt == 3 ~ "0.6–0.8",
        vegtype == "small"    & hgt == 4 ~ "0.8–1.0",
        vegtype == "small"    & hgt == 5 ~ "1.0–1.3",
        vegtype == "small"    & hgt == 6 ~ "1.3–2.0",
        vegtype == "advanced" & hgt == 1 ~ "2–4",
        vegtype == "advanced" & hgt == 2 ~ ">4",
        TRUE ~ NA_character_
      ),
      dbh = factor(
        recode(dbh,
               "1" = "10–20cm", "2" = "20–40cm",
               "3" = "40–60cm", "4" = ">60cm"),
        levels  = c("10–20cm", "20–40cm", "40–60cm", ">60cm"),
        ordered = TRUE
      ),
      stump_spc = case_when(
        is.na(stump_spc)              ~ NA_character_,
        stump_spc %in% 1:38           ~ stump_species_lookup[stump_spc],
        TRUE                          ~ as.character(stump_spc)
      )
    ) %>%
    group_by(cluster) %>%
    mutate(n_plots = n_distinct(ID)) %>%
    ungroup() %>%
    mutate(
      area           = n_plots * 4,
      scaling_factor = 10000 / area,
      stem_density   = n * scaling_factor
    ) %>%
    dplyr::select(-fid, -group, -point, -region)
}


# ── 3. Main loop ──────────────────────────────────────────────────────────────

all_data   <- vector("list", nrow(gpkg_files))
all_coords <- vector("list", nrow(gpkg_files))

for (i in seq_len(nrow(gpkg_files))) {
  
  path       <- gpkg_files$path[i]
  table_name <- gpkg_files$table_name[i]
  
  message("Processing: ", table_name)
  
  # ── vegetation data
  tryCatch({
    raw <- process_table(path, table_name)
    
    wide <- raw %>%
      pivot_wider(names_from = trait, values_from = value) %>%
      apply_recoding()
    
    all_data[[i]] <- wide
  }, error = function(e) {
    warning("Failed veg data for ", table_name, ": ", e$message)
  })
  
  # ── spatial coords
  tryCatch({
    sf_i <- read_sf(path, table_name) %>%
      mutate(
        ID      = paste(country, region, group + 100, point, sep = "_"),
        cluster = paste(region, group + 100, sep = "_")
      ) %>%
      dplyr::select(ID, cluster)
    
    all_coords[[i]] <- sf_i
  }, error = function(e) {
    warning("Failed spatial data for ", table_name, ": ", e$message)
  })
}


# ── 4. Combine ────────────────────────────────────────────────────────────────

full_data   <- bind_rows(all_data)
full_coords <- bind_rows(all_coords)   # stays as sf

# Quick checks
cat("Total plots (IDs):", length(unique(full_data$ID)), "\n")
cat("Countries:", paste(unique(full_data$country), collapse = ", "), "\n")


# ── 5. Derived outputs ────────────────────────────────────────────────────────

df_context <- full_data %>%
  dplyr::select(ID, cluster, country, clear, grndwrk, logging_trail,
                windthrow, planting, deadwood, anti_browsing) %>%
  distinct()

df_cluster <- full_data %>%
  group_by(species, cluster, country) %>%
  summarise(sum_stem_density = sum(stem_density, na.rm = TRUE), .groups = "drop")

sf_context <- full_coords %>%
  full_join(df_context, by = c("ID", "cluster"))


# ── 6. Save ───────────────────────────────────────────────────────────────────

data.table::fwrite(full_data, "outData/subplot_full_2023_allcountries.csv")
data.table::fwrite(df_cluster,   "outData/df_cluster_2023_allcountries.csv")

sf::st_write(
  sf_context,
  "outData/sf_context_2023_allcountries.gpkg",
  layer = "context",
  delete_layer = TRUE
)


summary(sf_context$x)

message("Done! ", nrow(full_data), " rows across ", 
        length(unique(full_data$country)), " countries.")


# 1. What CRS is actually stored in the sf object?
st_crs(sf_context)

# 2. Check CRS per source file — this is the key question
crs_check <- map_dfr(seq_len(nrow(gpkg_files)), function(i) {
  sf_i <- read_sf(gpkg_files$path[i], gpkg_files$table_name[i])
  tibble(
    file   = gpkg_files$basename[i],
    crs    = st_crs(sf_i)$epsg,
    crs_wkt = st_crs(sf_i)$input
  )
})

print(crs_check)







#### END
# prepare functions to create df -------------------

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


cz_wide_recode <- cz_wide %>%
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
table(cz_wide_recode$n)


# Corrent counts and categories: --------------------------

# because: 1 = >16, 2 = 1, 3 = 2,...17 = 16

# stupm species list: 
# Define lookup vector
stump_species_lookup <- c(
  "deciduous", "conifer", "abal", "acca", "acpl", "acps", "aehi", "aial", "algl", "alin",
  "alvi", "besp", "cabe", "casa", "fasy", "frex", "juni", "jure", "lade", "piab",
  "pisy", "potr", "posp", "prav", "psme", "quro", "qusp", "rops", "saca", "sasp",
  "soar", "soau", "soto", "taba", "tisp", "ulsp", "otsp_1"
)



cz_wide_corrected <- cz_wide_recode %>%
  # rorrect counts
  mutate(
    n = case_when(
      vegtype %in% c("small", "advanced") & n == 1 ~ 17,
      vegtype %in% c("small", "advanced")         ~ n-1,
     # vegtype == "mature" & n == 5                ~ ">4",
     # vegtype == "mature"                         ~ as.character(n),
      TRUE                                        ~ NA_integer_
    )
  ) %>% 
  # height class
  mutate(
    hgt = case_when(
      vegtype == "small" & hgt == 1 ~ "0.2–0.4",
      vegtype == "small" & hgt == 2 ~ "0.4–0.6",
      vegtype == "small" & hgt == 3 ~ "0.6–0.8",
      vegtype == "small" & hgt == 4 ~ "0.8–1.0",
      vegtype == "small" & hgt == 5 ~ "1.0–1.3",
      vegtype == "small" & hgt == 6 ~ "1.3–2.0",
      vegtype == "advanced" & hgt == 1 ~ "2–4",
      vegtype == "advanced" & hgt == 2 ~ ">4",
      TRUE ~ NA_character_  # handles unmatched cases
    )
  )%>%
  # dbh - recorded only for mature
  mutate(
    dbh = factor(
      recode(
        dbh,
        "1" = "10–20cm",
        "2" = "20–40cm",
        "3" = "40–60cm",
        "4" = ">60cm"
      ),
      levels = c("10–20cm", "20–40cm", "40–60cm", ">60cm"),
      ordered = TRUE
    )) %>% 
  # Stump species recoding
  mutate(
    stump_spc = case_when(
      is.na(stump_spc) ~ NA_character_,
      stump_spc %in% 1:38 ~ stump_species_lookup[stump_spc],
      TRUE ~ as.character(stump_spc)  # fallback in case of unexpected values
    )
  ) %>% 
  group_by(cluster) %>%
  mutate(n_plots = n_distinct(ID)) %>%
  ungroup() %>% 
  mutate(area = n_plots *4,
         scaling_factor = 10000/area,
         stem_density = n*scaling_factor) %>% 
  dplyr::select(-fid, -group, -point, -region)


# extract just context situation ----------------------------------------------
df_context <- cz_wide_corrected %>% 
  dplyr::select(ID, cluster, clear, grndwrk, logging_trail, windthrow,planting, deadwood, anti_browsing) %>% 
  distinct()


# summarize on cluster level
cz_wide_corrected_cluster <- cz_wide_corrected %>% 
  group_by(species, cluster) %>% 
  summarize(sum_stem_density = sum(stem_density, na.rm = T))


# clean up spatial data ---------------------------------------------------------
# keep only coordinates
cz1_coords <- cz1_sf %>%
  mutate(ID = paste(country, region, group + 100, point, sep = "_")) %>%
  mutate(cluster = paste(region, group + 100, sep = "_")) %>%
  dplyr::select(ID, cluster)
  
cz2_coords <- cz2_sf %>%
  mutate(ID = paste(country, region, group + 100, point, sep = "_")) %>%
  mutate(cluster = paste(region, group + 100, sep = "_")) %>%
  dplyr::select(ID, cluster)

# Combine and convert to sf
cz_coords_all <- bind_rows(cz1_coords, cz2_coords) 

# add to context information
sf_context <- cz_coords_all %>%
  full_join(df_context)

# add data about missing terminal???
# alsmlno, need to visually correct fencing information! eg if NA, and obvuiously within a fence 


# save files: -------------------
fwrite(cz_wide_corrected, 'outData/subplot_full_2023.csv')
fwrite(cz_wide_corrected_cluster, 'outData/df_cluster_2023.csv')

# Write to new GPKG
st_write(sf_context,
         "outData/sf_context_2023.gpkg",
         layer = "context",
         delete_layer = TRUE)



