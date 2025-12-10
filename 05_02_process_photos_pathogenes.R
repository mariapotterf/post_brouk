
# identify photso of samples

# read all fodlers in
# loop over photos
# filter photos of damage
# export to final folder to share

library(fs)
library(stringr)
library(data.table)
library(dplyr)
library(writexl)  # for Excel writing

# Set base directory and target folder
base_dir <- "raw/collected_2025"
target_dir <- "outShare/sample_photo_renamed"  # You can customize this

# export final table as csv
df_sample_photos <- fread('outShare/samples_list.csv') %>% 
  dplyr::filter(photo != "")  

# writexl::write_xlsx(
#   df_sample_photos,
#   path = "outShare/samples_list.xlsx"
# )

# Create target directory if it doesn't exist
dir_create(target_dir)

# Find all files in subfolders named 'DCIM' that contain '_fol_dtl_'
photo_files <- dir_ls(
  path = base_dir,
  recurse = TRUE,
  #regexp = "DCIM/rg.*\\.(jpg|jpeg|png)$",
  regexp = "DCIM/(rg|mat).*\\.(jpg|jpeg|png)$",
  type = "file"
)

# Deduplicate by file name (assumes duplicates have same filename)
unique_files <- photo_files[!duplicated(path_file(photo_files))]


# filter only photos with sample
# Create a named vector of photo basenames and full paths
photo_lookup <- setNames(photo_files, path_file(photo_files))

# # Match only those in sample_photos$photo
# matched_photos <- df_sample_photos %>%
#   dplyr::filter(photo %in% names(photo_lookup)) %>%
#   mutate(
#     original_path = photo_lookup[photo],
#     new_filename = paste0(sample, "_", photo),
#     new_path = file.path(target_dir, new_filename)
#   )
# 
# # Copy and rename photos
# file_copy(matched_photos$original_path, matched_photos$new_path, overwrite = TRUE)
# 
# # Output result
# cat("Copied", nrow(matched_photos), "renamed sample photos to", target_dir, "\n")


# Join matched photo files - filter only photos that have sample in them
matched_photos <- df_sample_photos %>%
  dplyr::filter(photo %in% names(photo_lookup)) %>%
  mutate(original_path = photo_lookup[photo])

# Update photos names, extract date (YYYYMMDD) from photo name using regex
matched_photos_renamed <- matched_photos %>%
  mutate(
    sample_clean = word(sample, 1),  # get first word before space
    tablet = substr(sample_clean, 1, 2), # get indication of tablet - seems my date on photos is wrong on Tablet1?
    plot_date = str_extract(plot_key, "\\d{8}"),
    photo_date = str_extract(photo, "\\d{8}"), # extract first 8 characters from the time stamp
    photo_description = photo %>%
      str_replace("_\\d{17,}\\.(jpg|jpeg|png)$", "") %>%
      tools::file_path_sans_ext(),
    sort_key = as.integer(photo_date)
  ) %>%
  arrange(sort_key) %>% # arrange by date (as integer)
  mutate(
    photo_seq = sprintf("%03d", row_number()),  # padded sequence
    #new_filename = paste0(plot_date, "_", photo_seq, "_", sample_clean, "_", photo_description, ".jpg"), # , 
    new_filename = paste0(photo_seq, "_", sample_clean, ".jpg"), # , 
    new_path = file.path(target_dir, new_filename)
  ) %>% 
  dplyr::select(-c(sample, photo_seq)) %>% 
  rename(sample = sample_clean)

# Copy and rename photo files
file_copy(
  path = matched_photos_renamed$original_path,  # same as matched_photos$original_path
  new_path = matched_photos_renamed$new_path,
  overwrite = TRUE
)
# Update sample list with new photo names
df_sample_photos_updated <- df_sample_photos %>%
  left_join(
    matched_photos_renamed %>% dplyr::select(photo, new_filename),
    by = "photo"
  )

# Export updated table to Excel
write_xlsx(df_sample_photos_updated, "outShare/samples_list_renamed.xlsx")

# Output summary
cat("Copied", nrow(matched_photos_renamed), "renamed sample photos to", target_dir, "\n")




# read data from Roman - merge orevious and new database -----------------------
# need to marge dataset that was inspected manually with new photos 
# 2025/12/10 
library(readxl)

df_manual <- read_xlsx("outShare/samples_list_zari_maja.xlsx") %>%
  dplyr::filter(photo != "")

df_new <- read_xlsx("outShare/samples_list_renamed.xlsx") %>%
  dplyr::filter(photo != "")

head(df_manual)
head(df_new)

nrow(df_manual)
nrow(df_new)


df_merged <- df_manual %>%
  left_join(
    df_new %>% dplyr::select(plot_key, sample, photo, new_filename),
    by = c("plot_key", "sample", "photo")
  ) %>% 
  select(-c(photo)) %>% 
  rename(photo = new_filename)


# Export updated table to Excel
write_xlsx(df_merged, "outShare/samples_list_merged.xlsx")

