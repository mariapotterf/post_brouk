

# get year of disturbance  - from Prosper's maps
# read dat_czechia
library(terra)

dist_2018 <- rast("raw/disturb_recovery_Prosper/2018.tif") 
dist_2019 <- rast("raw/disturb_recovery_Prosper/2019.tif") 
dist_2020 <- rast("raw/disturb_recovery_Prosper/2020.tif") 
dist_2021 <- rast("raw/disturb_recovery_Prosper/2021.tif") 
dist_2022 <- rast("raw/disturb_recovery_Prosper/2022.tif") 
dist_2023 <- rast("raw/disturb_recovery_Prosper/2023.tif") 
dist_2024 <- rast("raw/disturb_recovery_Prosper/2024.tif") 

rasters <- list.files("raw/disturb_recovery_Prosper", pattern = "\\.tif$", full.names = TRUE)

# keep only 2018-2020
rasters <- rasters[1:3]
dat_2023 <- vect('raw/dat_czechia_2023.gpkg')

# Check metadata for all rasters
for (f in rasters) {
  r <- rast(f)
  crs(r) <- "EPSG:5514"
  print(f)
  print(ext(r))
  print(crs(r))
}

# problem is raster 2021

# align raster: 2021 is problematic
aligned_rasters <- list()

for (f in rasters) {
  r <- rast(f)
  
  # If CRS mismatches, reproject (optional, depending on your case)
  if (!compareGeom(r, ref, stopOnError = FALSE)) {
    r <- resample(r, ref, method = "near")  # use nearest neighbor for integer rasters
  }
  crs(r) <- "EPSG:5514"
  aligned_rasters[[f]] <- r
}

#stacked <- rast(aligned_rasters)
stacked <- rast(rasters)

# change vector projection to raster
ref <- rast(aligned_rasters[1])

crs(ref) <- "EPSG:5514"
dat_2023_proj <- project(dat_2023, crs(ref))

vals <- extract(stacked, dat_2023_proj) 
