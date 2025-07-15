
# identify photso of samples

# read all fodlers in
# loop over photos
# filter photos of damage
# export to final folder to share

library(fs)
library(stringr)
library(data.table)

# Set base directory and target folder
base_dir <- "raw/collected_2025"
target_dir <- "outShare/sample_photo"  # You can customize this

# export final table as csv
df_sample_photos <- fread('outTable/samples_list.csv') %>% 
  dplyr::filter(photo != "")  

# Create target directory if it doesn't exist
dir_create(target_dir)

# Find all files in subfolders named 'DCIM' that contain '_fol_dtl_'
photo_files <- dir_ls(
  path = base_dir,
  recurse = TRUE,
  regexp = "DCIM/rg.*\\.(jpg|jpeg|png)$",
  type = "file"
)

# Deduplicate by file name (assumes duplicates have same filename)
unique_files <- photo_files[!duplicated(path_file(photo_files))]


# filter only photos with sample
# Create a named vector of photo basenames and full paths
photo_lookup <- setNames(photo_files, path_file(photo_files))

# Match only those in sample_photos$photo
matched_photos <- df_sample_photos %>%
  dplyr::filter(photo %in% names(photo_lookup)) %>%
  mutate(
    original_path = photo_lookup[photo],
    new_filename = paste0(sample, "_", photo),
    new_path = file.path(target_dir, new_filename)
  )

# Copy and rename photos
file_copy(matched_photos$original_path, matched_photos$new_path, overwrite = TRUE)

# Output result
cat("Copied", nrow(matched_photos), "renamed sample photos to", target_dir, "\n")
