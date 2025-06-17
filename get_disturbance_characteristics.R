
# 
# Prepare spatial database of the each cluster point
# 

# read field data by level of country
# get raster data:
# elevation
# data from the Cornelius: disturbances
#    - year
#    - agent - NA = no disturbance
# 1 = bark beetle or wind disturbances (both classes had to be grouped due to technical reasons)
# 2 = fire disturbances
# 3 = other disturbances, mostly harvest but might include salvage logging go small-scale natural disturbances and infrequent other natural agents (e.g., defoliation, avalanches, etc.)


#    - distance to nearest edge  (merge 2018-2020 into one!)
# patch zize & patch characteristics


library(terra)
library(raster)
library(dplyr)
library(tidyr)
library(data.table)
library(rasterize)



# read the plots by country
# read respective rasters: DEM, disturbances...
#    = stack them
# extract info by plot level
# move toanother country

# REad input data --------------------------------
# get paths to rasters 
elev_path          <- "raw/dem"
dist_path          <- "raw/disturb_data" # for year and severity

# study area: C4
# read field data
aoi        <- vect("raw/core_4.gpkg")
subplots   <- vect("raw/dat_czechia_2023.gpkg")




# Test distance to edge: prepare for loop over ----------------------------------------

# list all countries
country_names <- list("czechia") 

# distances will contain the distance of each point to the nearest edge of its patch

### TEST example ------------------------------------------------------------------


# Create a small example raster
mat <- matrix(c(3, NA, 2, NA, NA,
                NA, NA, 1, NA, NA,
                NA, NA, 1, 2, NA,
                NA, NA, 2, 1, NA,
                NA, NA, 2, 2, 2), byrow = T, nrow=5, ncol=5)
example_raster <- rast(nrows=5, ncols=5, vals=mat, crs = 'EPSG:3035')
plot(example_raster)

# Step 2: Temporarily replace all NA values with 0 (or another placeholder)
example_raster[is.na(example_raster)] <- 0

plot(example_raster, main = "Raster with NA as 0")

# Reclassify the raster - change values '2' to NA
#reclass_matrix <- matrix(c(0.8,2, NA), ncol=3, byrow=TRUE)  #matrix(c(1, NA), ncol=2, byrow=TRUE)
reclass_matrix <- matrix(c(0.8, 2, NA, 
                           3, 3, 3,     # Keep value 3 as-is
                           0, 0, 0),    # Keep placeholder 0 as-is
                         ncol = 3, byrow = TRUE)
reclassified_raster <- classify(example_raster, reclass_matrix)
plot(reclassified_raster)
# Assuming reclassified_raster is your raster
# Create a data frame with the coordinates of the point
point_df <- data.frame(x = 75, y =0)

# Create a SpatVector from the data frame
point_vector <- vect(point_df, geom = c("x", "y"), crs = crs(reclassified_raster))
plot(point_vector, add = T)

# Calculate the distance to the nearest NA
distance_to_edge <- distance(reclassified_raster) # , point_vector


# Plotting to visualize the results
plot(example_raster, main="example_raster")
plot(point_vector, main="", add = T)
plot(reclassified_raster, main="Reclassified Raster")
plot(point_vector, main="", add = T)
plot(distance_to_edge, main="")
plot(point_vector, main="", add = T)

# Measure the distance from the point to the nearest NA
point_distance <- extract(distance_to_edge, point_vector)
(point_distance)
plot(distance_to_edge, main="Distance to Nearest Edge")
plot(point_vector, main="", add = T)






# distance to edge: test for single country -----------------------------------------------------------
country_name = 'czechia'

# read field data 
#subplots = vect('raw/dat_czechia_2023.gpkg')
subplots_proj <- project(subplots, "EPSG:3035")

# select only one point per plot (as I am aggregation information afterwards anyway)
subplots_proj <- terra::aggregate(subplots_proj, by = "cluster", fun = "first")

# disturbance year
disturb_name = paste0('disturbance_year_1986-2020_', country_name, '.tif')
disturbance  = rast(paste(dist_path, country_name, disturb_name, sep = '/'))

# read raster data
desired_crs <- crs(disturbance)

# Create a SpatVector from the data frame
#point_vector <- subplots_proj[15,]# vect(point_df, geom = c("x", "y"), crs = crs(reclassified_raster))
point_vector <- subplots_proj[subplots_proj$cluster == "14_106", ]
crs(point_vector) <- desired_crs

# get buffer: lower to 1 km diostance
buff <- buffer(point_vector, 100)
crs(buff) <- desired_crs

# crop data and then mask them to have a circle
disturbance_crop <- crop(disturbance, buff)
disturbance_mask <- mask(disturbance_crop, buff)

plot(disturbance_mask)
plot(point_vector, main="", add = T)

# Reclassify the raster - change values '2' to NA

# set NA to values (to handle the empty pixels), my distnace function calculates distance for every 
# cell that is NA
disturbance_mask[is.na(disturbance_mask)] <- 0

# reclassify matrix: disturbances to NAs
reclass_matrix <- matrix(c(#0,0,0,
                           1985, 2017.1, 1,
                           2017.5, 2021, NA), 
                         ncol=3, 
                         byrow=TRUE)
reclassified_raster <- classify(disturbance_mask, reclass_matrix)


plot(disturbance_mask)
plot(reclassified_raster)
plot(point_vector, main="", add = T)

# Calculate the distance to the nearest NA
distance_to_edge <- distance(reclassified_raster)

# Plotting to visualize the results
plot(reclassified_raster, main="Reclassified Raster")
plot(distance_to_edge, main="Distance to Nearest Edge")
plot(point_vector, main="", add = T)

# Measure the distance from the point to the nearest NA
point_distance <- terra::extract(distance_to_edge, point_vector)
(point_distance)
plot(distance_to_edge, main="Distance to Nearest Edge")
plot(disturbance_mask)
plot(point_vector, main="", add = T)


# Make a function distance ti edge -------------------------
process_point <- function(point, disturbance, buffer_dist) {
  # read raster data
  desired_crs <- crs(disturbance)
  crs(point) <- desired_crs
  
  # Create buffer around the point
  buff <- buffer(point, buffer_dist)
  crs(buff) <- desired_crs
  
  # Crop and mask the disturbance raster
  disturbance_crop <- crop(disturbance, buff)
  disturbance_mask <- mask(disturbance_crop, buff)
  
  # set NA to values (to handle the empty pixels)
  disturbance_mask[is.na(disturbance_mask)] <- 0
  
  
  # Reclassify the raster
  reclass_matrix <- matrix(c(0,0,0,
                             1985, 2017.1, 1,
                             2017.5, 2021, NA), 
                           ncol=3, 
                           byrow=TRUE)
  reclassified_raster <- classify(disturbance_mask, reclass_matrix)
  
  # Calculate the distance to the nearest NA
  distance_to_value <- distance(reclassified_raster)
  
  # Extract the distance for the point
  point_distance <- terra::extract(distance_to_value, point)
  
  return(data.frame(ID = point$ID, distance = point_distance))
}


# Loop over all countries
buffer_dist = 1500
extract_distance_to_edge <- function(country_name, buffer_dist=buffer_dist) {
  print(paste("Processing", country_name))
  
  # Read field data
  subplots <- vect(paste0('raw/dat_czechia_2023.gpkg'))
  subplots_proj <- project(country, "EPSG:3035")
  
  # Aggregate points by cluster
  #subplots_proj <- terra::aggregate(subplots_proj, by = "cluster", fun = "first")
  
  
  # Read disturbance raster data
  disturb_name <- paste0('disturbance_year_1986-2020_', country_name, '.tif')
  disturbance <- rast(paste0(dist_path, '/', country_name, '/', disturb_name))
  disturbance <- project(disturbance, "EPSG:3035")
  
  # Process each point and collect results
  results <- lapply(1:nrow(subplots_proj), function(i) {
    point_vector <- subplots_proj[i,]
    process_point(point_vector, disturbance, buffer_dist)
  })
  
  # Combine all results into one dataframe
  final_results <- do.call(rbind, results)
  final_results$country <- country_name
  
  # rename columns to fit
  colnames(final_results) <- c("ID","dist_ID",  "distance", "country")
  return(final_results)
}

# run or all countries
#country_names =
all_results <- lapply(country_names, function(cn) {
  extract_distance_to_edge(cn, buffer_dist)
})

# Combine results from all countries
final_results_distance <- do.call(rbind, all_results)

final_results_distance <- final_results_distance %>% 
  mutate(distance = round(distance,0 ))


# Get landscape patchiness ----------------------------------------
# Ensure AOI and fieold data  are in the same CRS as disturbance raster
aoi_proj <- project(aoi, crs(disturbance))
subplots_proj <- project(subplots, crs(disturbance))


# Crop and mask disturbance raster to AOI extent
disturbance_crop <- crop(disturbance, aoi_proj)
disturbance_mask <- mask(disturbance_crop, aoi_proj)


# reclassify matrix: disturbances to NAs, only 2018:2020 are coded
rcl_disturb_matrix <- matrix(c(
  1986, 2017, NA,   # set all years from 1986 to 2017 to NA
  2018, 2020, 1     # set 2018 to 2020 to 1
), ncol = 3, byrow = TRUE)

disturb_simple <- classify(disturbance_mask, rcl_disturb_matrix)
# convert to binary raster
disturb_bin <- classify(disturb_simple, rcl = matrix(c(0, 0, NA), ncol = 3, byrow = TRUE))

# Detect contiguous patches of 1s
disturb_patches <- patches(disturb_bin, directions = 8)  # 4 = rook, 8 = queen

# Optional: check result
plot(disturb_patches, main = "Disturbance Patches")


# Optional: visualize
plot(disturbance_mask, main = paste("Disturbance in", country_name))
plot(aoi_proj, add = TRUE)

#writeRaster(disturb_patches, "outData/disturb_patches.tif", overwrite = TRUE)

writeRaster(disturb_patches, "outData/patches_compr.tif", overwrite = TRUE,
            wopt = list(gdal = "COMPRESS=LZW", datatype = "INT2U"))

# Get patch sizes
patch_sizes <- freq(disturb_patches)

# Add area if you want size in m² (assuming projected CRS)
patch_sizes_df <- as.data.frame(patch_sizes)
patch_sizes_df$area_ha <- patch_sizes_df$count * res(disturb_simple)[1] * res(disturb_simple)[2] / 10000

# Sort by size if useful
patch_sizes_df <- patch_sizes_df[order(-patch_sizes_df$area_ha), ]

hist(patch_sizes_df$area_ha)

### extract patch information to field observation  --------- : add it as a ratser to disturbances?

# Step 2: Extract patch ID for each point (NA if not in disturbed patch)
subplots_proj$patch_id <- terra::extract(disturb_patches, subplots_proj)[,2]

# add patch informtion: size, area..
subplots_proj_sf <- subplots_proj %>% 
  st_as_sf() %>% 
  dplyr::right_join(patch_sizes_df,by = c("patch_id" = "value"))


# get summary statistics -----------------------------------

# 1. Count how many unique clusters fall into each patch
clusters_per_patch <- subplots_proj_sf %>%
  distinct(patch_id, cluster) %>%
  count(patch_id, name = "n_clusters") #%>% 
  #dplyr::filter(n_clusters > 2)

# 2. Summary statistics for patch area
patch_stats_summary <- subplots_proj_sf %>%
  distinct(patch_id, area_ha) %>%
  summarise(
    mean_area = mean(area_ha, na.rm = TRUE),
    q25 = quantile(area_ha, 0.25, na.rm = TRUE),
    median = median(area_ha, na.rm = TRUE),
    q75 = quantile(area_ha, 0.75, na.rm = TRUE),
    max_area = max(area_ha, na.rm = TRUE),
    min_area = min(area_ha, na.rm = TRUE)
  )
patch_stats_summary


# Extract disturbance and elevation data --------------------------


# Function to read or create a dummy raster
read_or_dummy_raster <- function(path, reference_raster) {
  if (file.exists(path)) {
    rast <- terra::rast(path)
    terra::project(rast, desired_crs)
  } else {
    # Create a dummy raster with NA values
    dummy_rast <- terra::rast(disturbance, nlyr = 1)
    values(dummy_rast) <- NA
    dummy_rast
  }
}


country_name = 'czechia'



# function to extract all data
extract_disturb_info <- function(country_name) {
  # country_name = 'germany' 
  
  print(country_name)
  
  read_or_dummy_raster <- function(path) {
    if (file.exists(path)) {
      rast <- terra::rast(path)
      #terra::project(rast, "EPSG:3035")
    } else {
      # Create a dummy raster with NA values
      dummy_rast <- terra::rast(disturbance, nlyr = 1)
      values(dummy_rast) <- NA
      dummy_rast
    }
  }
  
  # read field data 
  subplots <- vect(paste0('raw/dat_czechia_2023.gpkg'))
  subplots_proj <- project(country, "EPSG:3035")
  
  # disturbance year
  disturb_name = paste0('disturbance_year_1986-2020_', country_name, '.tif')
  disturbance  = rast(paste(dist_path, country_name, disturb_name, sep = '/'))
  
  # read raster data
  desired_crs <- crs(disturbance)
  
  # disturbace severity
  severity_name = paste0('disturbance_severity_1986-2020_', country_name, '.tif')
  severity      = rast(paste(dist_path, country_name, severity_name, sep = '/'))
  severity_proj = terra::project(x = severity, y = disturbance,  method="near")
  
  # disturbance agent
  #read_or_dummy_raster(paste(dist_path, agent_name, sep = '/'), subplots_proj)
  agent_name = paste0('fire_wind_barkbeetle_', country_name, '.tif')
  agent      = read_or_dummy_raster(paste(dist_path, agent_name, sep = '/'))
  crs(agent) <-desired_crs
  agent_proj = terra::project(x = agent, y = disturbance,  method="near")
  
  crs(subplots_proj)
  crs(disturbance)
  crs(severity_proj)
  crs(agent_proj)
  
  
  # elevation
  elev_name      <- paste0(toupper(substr(country_name, 1, 1)), tolower(substr(country_name, 2, nchar(country_name))))
  elevation      <- rast(paste0(elev_path, '/dem_', elev_name, '.tif'))
  elevation_proj <- terra::project(x = elevation, y = disturbance,  method="near")
  elevation_proj <- terra::resample(x = elevation_proj, y = disturbance, method="near")
  
  # create raster stacks
  dist.stack <- c(disturbance, severity_proj, agent_proj, elevation_proj)
  names(dist.stack) <- c("disturbance_year", 
                         "disturbance_severity", 
                         "disturbance_agent", 
                         "elevation")
  
  # extract elevation for every point
  plots.disturbance <- terra::extract(dist.stack, subplots_proj, method = "simple", bind=TRUE)
  
  # export plot
  return(plots.disturbance)
  
}

out <- extract_disturb_info('czechia')
View(as.data.frame(out))
# list all countries
#country_names <- list( "austria", "czechia") #austria" 

out_ls <- lapply(country_names, extract_disturb_info)

# merge them all in one file
disturbance_ls <- do.call("rbind", out_ls)

disturbance_df <- data.frame(disturbance_ls)


# replace NA by the values in the same cluster (~ 20 points in total, lay at the edge of the pixel)
disturbance_df <- disturbance_df %>%
  mutate(ID_short = gsub('.{2}$', '', ID)) %>%
  group_by(ID_short) %>%
  arrange(ID_short, disturbance_year) %>% # Arrange if necessary to ensure the first value is non-NA
  fill(disturbance_year,     .direction = "down") %>% 
  fill(disturbance_severity, .direction = "down") %>% 
  fill(disturbance_agent,    .direction = "down") %>% # Fill NAs downwards 
  ungroup(.) %>% 
  dplyr::select(-ID_short) %>% 
  mutate(disturbance_year = if_else(ID == '11_19_137_4', 2019, disturbance_year))

#11_19_137_4 - fallen outside of teh pixel, checked manually


fwrite(disturbance_df, 'outData/disturbance_chars.csv')
fwrite(final_results_distance, 'outData/distance_to_edge.csv')







