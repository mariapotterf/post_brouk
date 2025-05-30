
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
