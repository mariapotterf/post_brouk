# process drone imagery

# extract height information
# maybe some additional classification?

gc()

library(terra)
library(dplyr)
library(purrr)
library(sf)

# Read subplot data for 2023
dat23_sf   <- st_read("outData/sf_context_2023.gpkg") # data per subplot
drone23_24_sf <- st_read("raw/position_drone.gpkg") # data per subplot

# Read drone CHM rasters from 2023&2024
chm_folder <- "raw/UAV_Images" # /2023
chm_files <- list.files(chm_folder, pattern = "^CHM_.*\\.tif$", full.names = TRUE, recursive = TRUE)


# process data -----------------------
drone23_24_sf <- drone23_24_sf %>% 
  rename(drone_year = year)

# read rasters
chm_rasters <- set_names(chm_files, tools::file_path_sans_ext(basename(chm_files))) %>%
  map(terra::rast)

# Load CHMs as named list
chm_rasters <- map(chm_rasters, function(x) {
  crs(x) <- "EPSG:5514"
  x
})

# Assign proper CRS to rasters (S-JTSK / EPSG:5514)
krovak_crs <- "EPSG:5514"

# Reproject vector data to match raster CRS (faster!)
target_crs <- crs(chm_rasters[[1]])
dat23_5514   <- st_transform(dat23_sf, target_crs)
drone_5514   <- st_transform(drone23_24_sf, target_crs)

# get overlapping clusters with drones
dat23_subset <- st_intersection(dat23_5514, drone_5514)

#Get buffer for cluster (centroid of teh cluster)
dat23_subset <- dat23_subset %>%
  mutate(field_year = 2023,
         x = st_coordinates(.)[, 1],
         y = st_coordinates(.)[, 2])

# Start fresh from dat23_5514 — raw points with correct geometry
dat23_clusters <- dat23_subset %>%
  as.data.frame() %>% 
  group_by(cluster) %>%
  summarise(
    x = mean(x),
    y = mean(y),
    drone_year  = first(drone_year),
    field_year  = first(field_year),
    uav_ID = first(uav_ID),
    .groups = "drop"
  ) %>%
  mutate(uav_ID = paste0("CHM_", uav_ID)) %>% 
  st_as_sf(coords = c("x", "y"), crs = target_crs)  # Use correct CRS

# Convert to terra vector
cluster_vect <- vect(dat23_clusters)

data.frame(cluster_vect)

# Create buffer  and add ID
buff_width = 60
buffs <- buffer(cluster_vect, width = buff_width)
buffs$ID <- seq_len(nrow(buffs))  # consistent ID

#--------------------------------------------------------

# test run: buffer 26_111 overlaps with drone uav6
r <- chm_rasters$CHM_uav6
buf <- buffs[buffs$cluster == "26_111", ]

plot(r)
plot(buf, add = TRUE, border = "red")


# Ensure cluster_vect is a data.frame (not SpatVector or sf)
cluster_lookup <- as.data.frame(cluster_vect)[, c("cluster", "uav_ID")]

# Define buffer widths to test
buffer_sizes <- c(5, 10, 30, 60)

# Prepare output list
all_cluster_heights <- list()

# Loop over buffer sizes
for (buff_width in buffer_sizes) {
  cat("Processing buffer:", buff_width, "\n")
  
  # Create buffer for this size
  buffs <- terra::buffer(cluster_vect, width = buff_width)
  buffs$ID <- seq_len(nrow(buffs))
  
  # Loop through drone rasters
  cluster_heights <- lapply(names(chm_rasters), function(drone_id) {
    r <- chm_rasters[[drone_id]]
    r_name <- names(r)[1]
    
    # Get clusters linked to this drone
    matching_clusters <- cluster_lookup %>%
      filter(uav_ID == drone_id) %>%
      pull(cluster)
    
    buffers_this_drone <- buffs[buffs$cluster %in% matching_clusters, ]
    
    if (nrow(buffers_this_drone) == 0) return(NULL)
    
    # Extract
    vals <- terra::extract(r, buffers_this_drone, ID = TRUE)
    vals$cluster <- buffers_this_drone$cluster[vals$ID]
    
    # Summarize
    vals %>%
      as_tibble() %>%
      group_by(cluster) %>%
      summarize(
        drone       = drone_id,
        buffer_size = buff_width,
        mean_height = mean(.data[[r_name]], na.rm = TRUE),
        sd_height   = sd(.data[[r_name]], na.rm = TRUE),
        cv_height   = sd_height / mean_height,
        max_height  = max(.data[[r_name]], na.rm = TRUE),
        min_height  = min(.data[[r_name]], na.rm = TRUE),
        .groups     = "drop"
      )
  })
  
  all_cluster_heights[[as.character(buff_width)]] <- bind_rows(cluster_heights)
}

# Combine all buffer sizes into one table
chm_summary_multi <- bind_rows(all_cluster_heights)


# save output -----------------------------------------
fwrite(chm_summary_multi, "outTable/chm_buff_summary.csv")
