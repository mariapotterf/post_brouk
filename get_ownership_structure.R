

# get ownership class:

#- intersect points from 2023 with ownership structure
#- goal: where to lead the 20205 field work?
  
library(terra)
library(sf)
library(dplyr)
library(ggplot2)


# Input -------------------------------------------------------------------------
# get ownership layers
private    <- vect("raw/Ownership_Czech/Fyzické_osoby.shp")
unknown    <- vect("raw/Ownership_Czech/neznamo.shp")
village    <- vect("raw/Ownership_Czech/Obec.shp")
company    <- vect("raw/Ownership_Czech/Právnické_osoby.shp")
state      <- vect("raw/Ownership_Czech/Stát.shp")

# study area
aoi        <- vect("raw/core_4.gpkg")

subplots <- vect("raw/dat_czechia_2023.gpkg")

# drones extent
drone_ext <- vect("raw/position_drone.gpkg")

# ccheck coordinate system --------
CRS_5514 <- crs(unknown)
crs(private)
crs(unknown)
crs(village)
crs(company)
crs(state)

crs(private) <- CRS_5514

# are the CRS same? 

crs(private) == crs(unknown) &&
  crs(private) == crs(village) &&
  crs(private) == crs(company) &&
  crs(private) == crs(state)

# YES!


# clean up the Ownership data: i need only field: area, (m2), perimeter (m), FID, SpruceArea

private_clean <- private[, c("fid_1", "area", "perimeter")] #, "SpruceArea"
unknown_clean <- unknown[, c("fid_1", "area", "perimeter")]
village_clean <- village[, c("fid_1", "area", "perimeter")]
company_clean <- company[, c("fid_1", "area", "perimeter")] # no SpruceArea - 0
state_clean   <- state[, c("fid_1","area", "perimeter")]

# add types of of ownership
private_clean$owner <- 'private'
unknown_clean$owner <- 'unknown'
village_clean$owner <- 'village'
company_clean$owner <- 'company'
state_clean$owner   <- 'state'

# simplify
my_tolerance =10
private_simplified  <- simplifyGeom(private_clean,  tolerance = my_tolerance)
unknown_simplified  <- simplifyGeom(unknown_clean,  tolerance = my_tolerance)
village_simplified  <- simplifyGeom(village_clean,  tolerance = my_tolerance)
company_simplified  <- simplifyGeom(company_clean,  tolerance = my_tolerance)
state_simplified    <- simplifyGeom(state_clean,    tolerance = my_tolerance)

# merge them together
all_owners <- rbind(private_simplified, 
                    unknown_simplified, 
                    village_simplified, 
                    company_simplified, 
                    state_simplified)

# change projectiuon system 
all_owners_3035 <- project(all_owners, "EPSG:3035")
aoi_3035        <- project(aoi, 'EPSG:3035')  # Transforms coordinates
subplots_3035   <- project(subplots,  'EPSG:3035')
drone_ext_3035   <- project(drone_ext ,  'EPSG:3035')

# check if tehy all have the same projection
# are the CRS same? 

crs(all_owners_3035) == crs(aoi_3035) &&
  crs(subplots_3035) == crs(drone_ext_3035)


crs(aoi_3035)

## PLOT data on map ---------------------------------

# Clean up field data: keep only points within AOI -----------------
subplots_in_aoi     <- crop(subplots_3035, aoi_3035)
# Extract ownership info and preserve point ID
owner_info <- terra::extract(all_owners_3035, subplots_in_aoi)

table(table(owner_info$id.y))
dim(subplots_in_aoi)

# I have one point on twoi wnerships, aslo i have NA values: need to summarize data firt!
owner_summary <- owner_info %>%
  group_by(id.y) %>%
  summarise(property_type = paste(unique(owner), collapse = ";"),
            n_unique_owners = n_distinct(owner)) %>%
  ungroup()

# how many points are problematic?
owner_summary %>% filter(n_unique_owners > 1)

table(owner_summary$id.y)

# Check column names; assume 'owner' is the relevant attribute
# Merge back to point layer using ID
# Initialize and assign ownership info
subplots_in_aoi$property_type <- NA
subplots_in_aoi$property_type[owner_summary$id.y] <- owner_summary$property_type


table(subplots_in_aoi$property_type)/5


# convert to sf for ggplot
owners_sf <- st_as_sf(all_owners_3035)
aoi_sf    <- st_as_sf(aoi_3035)
subplots_sf    <- st_as_sf(subplots_in_aoi)


crs(owners_sf) == crs(subplots_sf)

ggplot() +
  geom_sf(data = owners_sf, aes(fill = owner), color = NA, lwd = 1.2) +
  geom_sf(data = subplots_sf, color = "white", size = 3) +
   geom_sf(data = subplots_sf, aes(fill = property_type), color = 'black', shape = 21,  size = 2, stroke = 0.8) +
  geom_sf(data = aoi_sf, fill =NA, color = "black") +
 # scale_fill_brewer(palette = "Set3") +
  coord_sf(crs = st_crs(3035)) +  # ← forces plot to stay in EPSG:3035
  theme_void()


# export shp with ownership types: 
writeVector(subplots_in_aoi, "outData/subplots_by_owner.gpkg", filetype = "GPKG", overwrite = TRUE)

# some cluysters can have mixed ownership! 
table(subplots_sf$cluster, subplots_sf$property_type )



# intersect drone position with ownership structure ---------

# Create an empty raster with 10 m resolution over the extent of ownership layer
template_raster <- terra::rast(ext(all_owners_3035), resolution = 10, crs = crs(all_owners_3035))

# Rasterize ownership polygons by "owner" field (e.g. 'state', 'private' etc.)
# This creates a categorical raster
ownership_raster <- terra::rasterize(all_owners_3035, 
                                     template_raster, field = "owner")


# loop over drones to clip every drone to ownership structure
# Get unique drone IDs
drone_ids <- unique(drone_ext_3035$uav_ID)

# Initialize results list
ownership_summary_list <- list()

for (id in drone_ids) {
#  id = 'uav5'
  # Subset footprint for current drone
  drone_poly <- drone_ext_3035[drone_ext_3035$uav_ID == id, ]
  
  # Crop and mask ownership raster to drone footprint
  cropped <- terra::crop(ownership_raster, drone_poly)
  masked <- terra::mask(cropped, drone_poly)
  
  # Summarize pixel counts by owner
  freq_table <- terra::freq(masked)
  
  # Add UAV ID and calculate area (10m x 10m = 100 m² = 0.01 ha)
  if (!is.null(freq_table)) {
    freq_table <- as.data.frame(freq_table)
    freq_table$uav_ID <- id
    freq_table$area_ha <- freq_table$count * 0.01
    ownership_summary_list[[id]] <- freq_table
  }
}

ownership_summary <- do.call(rbind, ownership_summary_list)

ownership_summary <- ownership_summary %>% 
  group_by(uav_ID) %>% 
  mutate(extent = sum(area_ha))
