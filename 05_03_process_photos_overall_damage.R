

# identify photos of damage - outside of teh samples

# read all folders in
# loop over photos
# filter photos of damage - no patogenes, this is in separate script
# export to final folder to share with Roman and Michal - by July 1st, 2026

# identify photos of samples AND other field-recorded damage
#
# 1) pathogen sample photos  -> matched from outShare/samples_list.csv       -> outShare/sample_photo_renamed/
# 2) other damage photos     -> matched from forest_structure_survey gpkg(s) -> outShare/damage_photo_renamed/
#
# damage is recorded per individual (regeneration_small, regeneration_adv2, mature_test)
# with dmg_bool == 1 flagging that *something* was recorded, and dmg_type coding *what*
# (terminal / stem / root-stem / foliage - see damage_list layer). Each damage "part" can
# have its own photo(s) and its own free-text sample field, which is either a specimen tag
# code (e.g. "T4433KL2") if something was physically collected, or just a Czech cause note
# (e.g. "Okus" = browsing, "Mraz" = frost, "Zver" = wildlife) if not.

library(fs)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(sf)
library(writexl)

base_dir <- "raw/collected_2025"


# ============================================================
# 1) PATHOGEN SAMPLE PHOTOS (unchanged from previous workflow)
# ============================================================

target_dir_sample <- "outShare/sample_photo_renamed"
dir_create(target_dir_sample)

df_sample_photos <- fread('outShare/samples_list.csv') %>%
  dplyr::filter(photo != "")

# all candidate photo files across all collection subfolders
photo_files <- dir_ls(
  path = base_dir,
  recurse = TRUE,
  regexp = "DCIM/(rg|mat).*\\.(jpg|jpeg|png)$",
  type = "file"
)
unique_files <- photo_files[!duplicated(path_file(photo_files))]
photo_lookup <- setNames(photo_files, path_file(photo_files))

matched_photos <- df_sample_photos %>%
  dplyr::filter(photo %in% names(photo_lookup)) %>%
  mutate(original_path = photo_lookup[photo])

matched_photos_renamed <- matched_photos %>%
  mutate(
    sample_clean = word(sample, 1),
    tablet = substr(sample_clean, 1, 2),
    plot_date = str_extract(plot_key, "\\d{8}"),
    photo_date = str_extract(photo, "\\d{8}"),
    photo_description = photo %>%
      str_replace("_\\d{17,}\\.(jpg|jpeg|png)$", "") %>%
      tools::file_path_sans_ext(),
    sort_key = as.integer(photo_date)
  ) %>%
  arrange(sort_key) %>%
  mutate(
    photo_seq = sprintf("%03d", row_number()),
    new_filename = paste0(photo_seq, "_", sample_clean, ".jpg"),
    new_path = file.path(target_dir_sample, new_filename)
  ) %>%
  dplyr::select(-c(sample, photo_seq)) %>%
  rename(sample = sample_clean)

file_copy(
  path = matched_photos_renamed$original_path,
  new_path = matched_photos_renamed$new_path,
  overwrite = TRUE
)

df_sample_photos_updated <- df_sample_photos %>%
  left_join(
    matched_photos_renamed %>% dplyr::select(photo, new_filename),
    by = "photo"
  )

write_xlsx(df_sample_photos_updated, "outShare/samples_list_renamed.xlsx")

cat("Copied", nrow(matched_photos_renamed), "renamed SAMPLE photos to", target_dir_sample, "\n")



# ============================================================
# 2) OTHER DAMAGE PHOTOS (new - from forest structure survey gpkg)
# ============================================================

target_dir_damage <- "outShare/damage_photo_renamed"
dir_create(target_dir_damage)

# there may be one gpkg per transect/date (e.g. .../T1_JC_20250717/forest_structure_survey_v2.gpkg)
# -> pick up all of them automatically so this doesn't need editing as new sessions come in
gpkg_files <- dir_ls(
  path = base_dir,
  recurse = TRUE,
  regexp = "forest_structure_survey.*\\.gpkg$",
  type = "file"
)

if (length(gpkg_files) == 0) {
  stop("No forest_structure_survey*.gpkg files found under ", base_dir)
}

# --- lookup tables (same in every gpkg, just read from the first one) ---
damage_type_lookup <- st_read(gpkg_files[1], layer = "damage_list", quiet = TRUE) %>%
  st_drop_geometry() %>%
  transmute(dmg_type = as.character(Value), damage_type = species)  # 'species' col holds the damage-type label here

species_lookup <- st_read(gpkg_files[1], layer = "species_list", quiet = TRUE) %>%
  st_drop_geometry() %>%
  transmute(species_id = as.character(Value), species_name = species)

# which photo column belongs to which damage "part", and a readable view label
photo_map <- tribble(
  ~photo_col,                   ~damage_part, ~photo_view,
  "dmg_terminal_photo",         "terminal",   "terminal",
  "dmg_stem_horizont_photo",    "stem",       "stem_horizontal",
  "dmg_stem_vert_photo",        "stem",       "stem_vertical",
  "dmg_root_stem_horiz_photo",  "root_stem",  "root_stem_horizontal",
  "dmg_root_stem_vert_photo",   "root_stem",  "root_stem_vertical",
  "dmg_foliage_detail_photo",   "foliage",    "foliage_detail",
  "dmg_foliage_overall_photo",  "foliage",    "foliage_overview"
)

# which free-text "sample/cause" column belongs to which damage part
sample_map <- tribble(
  ~sample_col,             ~damage_part,
  "dmg_term_sample",       "terminal",
  "dmg_stem_sample",       "stem",
  "dmg_root_stem_sample",  "root_stem",
  "dmg_foliage_sample",    "foliage"
)

# pull damage records + their photos (long format) from one survey layer of one gpkg
extract_damage_photos <- function(gpkg_file, layer_name, id_col) {
  
  raw <- st_read(gpkg_file, layer = layer_name, quiet = TRUE) %>%
    st_drop_geometry()
  
  if (!"dmg_bool" %in% names(raw)) return(NULL)
  
  present_photo_cols  <- intersect(photo_map$photo_col, names(raw))
  present_sample_cols <- intersect(sample_map$sample_col, names(raw))
  
  base_tbl <- raw %>%
    dplyr::filter(dmg_bool == 1) %>%
    mutate(
      record_id = .data[[id_col]],
      dmg_type = as.character(dmg_type),
      species_id = as.character(species_id)
    ) %>%
    dplyr::select(record_id, plot_id, species_id, dmg_type,
                  dplyr::all_of(present_photo_cols), dplyr::all_of(present_sample_cols))
  
  if (nrow(base_tbl) == 0) return(NULL)
  
  photos_long <- base_tbl %>%
    dplyr::select(record_id, plot_id, species_id, dmg_type, dplyr::all_of(present_photo_cols)) %>%
    pivot_longer(dplyr::all_of(present_photo_cols), names_to = "photo_col", values_to = "photo_path") %>%
    dplyr::filter(!is.na(photo_path), photo_path != "") %>%
    left_join(photo_map, by = "photo_col")
  
  if (nrow(photos_long) == 0) return(NULL)
  
  samples_long <- base_tbl %>%
    dplyr::select(record_id, dplyr::all_of(present_sample_cols)) %>%
    pivot_longer(dplyr::all_of(present_sample_cols), names_to = "sample_col", values_to = "sample_note") %>%
    dplyr::filter(!is.na(sample_note), str_trim(sample_note) != "") %>%
    left_join(sample_map, by = "sample_col")
  
  photos_long %>%
    left_join(
      samples_long %>% dplyr::select(record_id, damage_part, sample_note),
      by = c("record_id", "damage_part")
    ) %>%
    mutate(
      survey_layer = layer_name,
      gpkg_file = path_file(gpkg_file)
    )
}

layer_id_cols <- list(
  regeneration_small = "reg_uuid",
  regeneration_adv2  = "adv_reg_uuid",
  mature_test        = "mature_uuid"
)

damage_all <- purrr::map_dfr(gpkg_files, function(gf) {
  purrr::map_dfr(names(layer_id_cols), function(ly) {
    extract_damage_photos(gf, ly, layer_id_cols[[ly]])
  })
})

# specimen tag codes look like "T4433KL2" (T + digits); anything else in the sample
# field is a cause note, not a physical sample
damage_all <- damage_all %>%
  mutate(
    sample_note = str_trim(sample_note),
    # NA sample_note (nothing recorded at all) must resolve to FALSE, not NA -
    # otherwise if_else()/paste0() downstream silently prints the string "NA"
    # into the filename instead of treating it as "no sample"
    has_sample  = coalesce(str_detect(sample_note, "^T\\d"), FALSE),
    sample_code = if_else(has_sample, sample_note, NA_character_),
    cause_note  = if_else(!has_sample, sample_note, NA_character_)
  ) %>%
  left_join(damage_type_lookup, by = "dmg_type") %>%
  left_join(species_lookup, by = "species_id")

# match against photo files on disk (reuse the same file index built above)
damage_all <- damage_all %>%
  mutate(photo_basename = path_file(photo_path))

n_missing <- sum(!damage_all$photo_basename %in% names(photo_lookup))
if (n_missing > 0) {
  cat(n_missing, "damage photo(s) referenced in the gpkg were not found under", base_dir,
      "- check they've been copied off the camera/tablet.\n")
}

damage_matched <- damage_all %>%
  dplyr::filter(photo_basename %in% names(photo_lookup)) %>%
  mutate(
    original_path = photo_lookup[photo_basename],
    photo_date = str_extract(photo_basename, "\\d{8}"),
    sort_key = as.integer(photo_date)
  ) %>%
  arrange(sort_key) %>%
  mutate(
    photo_seq = sprintf("%03d", row_number()),
    # kept in the table for filtering, but no longer baked into the filename -
    # since only no-sample photos get copied, every copied file would end in
    # the same "_nosample" suffix, which is redundant
    sample_tag = if_else(has_sample, paste0("sample_", sample_code), "nosample"),
    new_filename = paste0(
      photo_seq, "_", plot_id, "_", damage_type, "_", photo_view,
      ".jpg"
    ),
    new_path = file.path(target_dir_damage, new_filename)
  )

file_copy(
  path = damage_matched %>% dplyr::filter(!has_sample) %>% pull(original_path),
  new_path = damage_matched %>% dplyr::filter(!has_sample) %>% pull(new_path),
  overwrite = TRUE
)

write_xlsx(
  damage_matched %>%
    mutate(copied_to_folder = !has_sample) %>%
    dplyr::select(plot_id, species_id, species_name, dmg_type, damage_type,
                  damage_part, photo_view, has_sample, sample_code, cause_note,
                  copied_to_folder, photo_date, photo_basename, new_filename,
                  survey_layer, gpkg_file),
  "outShare/damage_list_renamed.xlsx"
)

cat("Copied", sum(!damage_matched$has_sample), "renamed DAMAGE photos (no sample) to", target_dir_damage, "\n")
cat("  skipped", sum(damage_matched$has_sample), "photos where a physical sample was collected",
    "(listed in damage_list_renamed.xlsx but not copied).\n")
