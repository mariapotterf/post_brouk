# process drone imagery

# extract height information
# maybe some additional classification?

library(terra)
library(dplyr)
library(purrr)
library(sf)

# Read subplot data for 2023
dat23_sf <- st_read("raw/collected_2023/dat_czechia_2023.gpkg")

# Read drone CHM rasters from 2023
chm_folder <- "raw/UAV_Images/2023"
chm_files <- list.files(chm_folder, pattern = "^CHM_.*\\.tif$", full.names = TRUE)

chm_rasters <- set_names(chm_files, tools::file_path_sans_ext(basename(chm_files))) %>%
  map(terra::rast)

# Load CHMs as named list
chm_rasters <- map(chm_rasters, function(x) {
  crs(x) <- "EPSG:5514"
  x
})
# ───────────────────────────────────────────────────────────────

# Step 1: Assign proper CRS to rasters (S-JTSK / EPSG:5514)
krovak_crs <- "EPSG:5514"

# Step 2: Reproject vector data to match raster CRS (faster!)
target_crs <- crs(chm_rasters[[1]])
dat23_5514   <- st_transform(dat23_sf, target_crs)

# Step 3: Buffer around each cluster centroid
dat23_5514 <- dat23_5514 %>%
  mutate(year = 2023,
         x = st_coordinates(.)[, 1],
         y = st_coordinates(.)[, 2])

#st_geometry(dat23_5514) <- "geom"

# Start fresh from dat23_5514 — raw points with correct geometry
dat23_clusters <- dat23_5514 %>%
  group_by(cluster) %>%
  summarise(year = first(year), geometry = st_centroid(st_union(st_geometry(.)))) %>%
  ungroup() %>%
  st_as_sf()

# Then convert to terra vector
cluster_vect <- vect(dat23_clusters)


# Create buffer (150m) and add ID
buffer_60 <- buffer(cluster_vect, width = 60)
buffer_60$ID <- seq_len(nrow(buffer_60))  # consistent ID

# ───────────────────────────────────────────────────────────────
# 
r <- chm_rasters[[1]]
plot(r)
points(dat23_clusters, col = "red", pch = 3)

# Convert cluster points to terra and buffer
v <- vect(dat23_clusters)
buf <- terra::buffer(v[1, ], width = 150)

# Check overlap
plot(r)
plot(buf, add = TRUE, border = "red")



# test run: buffer 26_111 overlaps with drone uav6
r <- chm_rasters$CHM_uav6
buf <- buffer_60[buffer_60$cluster == "26_111", ]

# Extract values
vals <- terra::extract(r, dat23_5514, ID = TRUE)


plot(r)
points(buf, col = "red")
points(dat23_5514, col = 'red')
points(dat23_clusters, col = 'green')





plot(chm_rasters$CHM_uav6)
plot(buffer_60, add = T)
plot(dat23_clusters)


# Function to extract CHM data for buffers
extract_chm_for_image <- function(image, image_name, buffer_vect) {
  # Get intersecting buffer indices
  intersecting <- relate(buffer_vect, image, relation = "intersects")[, 1]
  if (length(intersecting) == 0) return(NULL)
  
  buffer_subset <- buffer_vect[intersecting, ]
  
  # Extract values
  vals <- terra::extract(image, buffer_subset, ID = FALSE)
  
  if (nrow(vals) == 0 || all(is.na(vals[[1]]))) {
    message("No valid data in ", image_name, " for overlapping buffers.")
    return(NULL)
  }
  
  vals$ID <- buffer_subset$ID
  
  # Join metadata
  lookup <- as.data.frame(buffer_subset) %>%
    dplyr::select(ID, cluster, year)
  
  vals <- left_join(vals, lookup, by = "ID")
  vals$image <- image_name
  
  return(vals)
}



# Run extraction across all CHM images
chm_extracted <- imap_dfr(chm_rasters, extract_chm_for_image, buffer_vect = buffer_150m)
