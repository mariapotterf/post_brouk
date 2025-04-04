

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

subplots <- vect("raw/dat_czechia.gpkg")

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


# clean up the data: i need only field: area, (m2), perimeter (m), FID, SpruceArea

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
my_tolerance =5
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
crs(aoi_3035)

## PLOT data on map ---------------------------------

# Clean up field data: keep only points within AOI -----------------


subplots_in_aoi     <- crop(subplots_3035, aoi_3035)
subplots_with_owner <- terra::extract(all_owners_3035, subplots_in_aoi)
# intersect subplots with ownership structure

# 2. Add 'ID' column to subplots (needed for join)
subplots_in_aoi$row <- 1:nrow(subplots_in_aoi)

# 3. Join extracted attributes back to the point layer using ID
subplots_with_owner <- merge(subplots_in_aoi, ownership_info, by.x = "ID", by.y = "ID")

round(table(subplots_with_owner$owner)/5,0)


# convert to sf for ggplot
owners_sf <- st_as_sf(all_owners_3035)
aoi_sf    <- st_as_sf(aoi_3035)
subplots_sf    <- st_as_sf(subplots_in_aoi)

# add ownership data:
subplots_sf <- subplots_sf %>% 
  full_join(subplots_with_owner)

crs(owners_sf) == crs(subplots_sf)

ggplot() +
  geom_sf(data = owners_sf, aes(fill = owner), color = NA, lwd = 1.2) +
  geom_sf(data = subplots_sf, color = "white", size = 3) +
   geom_sf(data = subplots_sf, aes(fill = owner_type), color = 'black', shape = 21,  size = 2, stroke = 0.8) +
  geom_sf(data = aoi_sf, fill =NA, color = "black") +
 # scale_fill_brewer(palette = "Set3") +
  coord_sf(crs = st_crs(3035)) +  # ← forces plot to stay in EPSG:3035
  theme_void()


# export shp with ownership types: 


writeVector(subplots_in_aoi, "outData/subplots_by_owner.gpkg", filetype = "GPKG", overwrite = TRUE)
