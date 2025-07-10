
# identify photso of samples

# read all fodlers in
# loop over photos
# filter photos of damage
# export to final folder to share

library(fs)
library(stringr)

# Set base directory and target folder
base_dir <- "raw/collected_2025"
target_dir <- "processed_photos/foliage_detail"  # You can customize this

# Create target directory if it doesn't exist
dir_create(target_dir)

# Find all files in subfolders named 'DCIM' that contain '_fol_dtl_'
photo_files <- dir_ls(
  path = base_dir,
  recurse = TRUE,
  regexp = "DCIM/.+_fol_dtl_.*\\.(jpg|jpeg|png)$",
  type = "file"
)

# Deduplicate by file name (assumes duplicates have same filename)
unique_files <- photo_files[!duplicated(path_file(photo_files))]

# Copy to target, renaming if necessary (preserve filename)
file_copy(unique_files, file.path(target_dir, path_file(unique_files)), overwrite = FALSE)

# Output
cat("Copied", length(unique_files), "unique _fol_dtl_ photos to", target_dir, "\n")