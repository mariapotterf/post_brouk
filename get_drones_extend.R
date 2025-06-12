
# Get raster extend  from orthophoto

library(terra)

# Set path to the folder with your rasters
raster_folder <- "raw/UAV_Images/2023"

# List all .tif files in the folder
raster_files <- list.files(paste0(getwd(), "/", raster_folder), 
                           pattern = "^CHM_.*\\.tif$", full.names = TRUE)

raster_files

# Function to extract extent polygon and UAV number
get_extent_polygon <- function(file) {
  r <- rast(file)
  e <- ext(r)
  p <- as.polygons(e, crs = crs(r))
  
  # Extract number using regex from filename (e.g., CHM_uav7.tif → 7)
  uav_id <- sub(".*uav(\\d+)\\.tif$", "\\1", basename(file))
  p$uav_id <- uav_id  # Add as attribute
  
  return(p)
}
# Apply to all rasters
extent_list <- lapply(raster_files, get_extent_polygon)

# Combine all polygons into a single SpatVector
extent_polygons <- do.call(rbind, extent_list)

# Preview the result
plot(extent_polygons)

# Export to GeoPackage
output_path <- file.path("outData/drone2023_extents.gpkg")
writeVector(extent_polygons, output_path, filetype = "GPKG", overwrite = TRUE)


# - !!! TEST START -------------------

# Read raster and convert only non-NA areas to polygon
get_extent_polygon <- function(file, year) {
  r <- rast(file)
  
  valid_area <- !is.na(r[[1]])  # safer if multi-band raster
  
  p <- as.polygons(valid_area, dissolve = TRUE)
  crs(p) <- crs(r)
  
  uav_id <- sub(".*uav(\\d+)\\.tif$", "\\1", basename(file))
  p$uav_id <- paste0("uav", uav_id)
  p$year <- year  # this now works
  
  return(p)
}


# Function to extract all drone extents for a given year
extract_drone_extents <- function(year) {
  raster_folder <- file.path("raw/UAV_Images", as.character(year))
  raster_files <- list.files(raster_folder, pattern = "^CHM_.*\\.tif$", full.names = TRUE)
  
  # Apply improved function for each raster file
  extent_list <- lapply(raster_files, get_extent_polygon, year = year)
  
  # Combine into one SpatVector
  extent_polygons <- do.call(rbind, extent_list)
  
  # Export individual year's extents
  output_path <- file.path("outData", paste0("drone", year, "_extents.gpkg"))
  writeVector(extent_polygons, output_path, filetype = "GPKG", overwrite = TRUE)
  
  return(extent_polygons)
}

# Run for 2023 and 2024
ext_2023 <- extract_drone_extents(2023)
ext_2024 <- extract_drone_extents(2024)

# Combine both into one SpatVector
ext_combined <- rbind(ext_2023, ext_2024)

# Export combined dataset
writeVector(ext_combined, "outData/drone_combined_extents.gpkg", filetype = "GPKG", overwrite = TRUE)

# TEST END !!! --------------------
