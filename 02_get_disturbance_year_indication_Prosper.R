

# get year of disturbance  - from Prosper's maps
# read dat_czechia
library(terra)
library(tidyr)
library(dplyr)

dist_2018 <- rast("raw/disturb_recovery_Prosper/2018.tif") 
dist_2019 <- rast("raw/disturb_recovery_Prosper/2019.tif") 
dist_2020 <- rast("raw/disturb_recovery_Prosper/2020.tif") 
# dist_2021 <- rast("raw/disturb_recovery_Prosper/2021.tif") 
# dist_2022 <- rast("raw/disturb_recovery_Prosper/2022.tif") 
# dist_2023 <- rast("raw/disturb_recovery_Prosper/2023.tif") 
# dist_2024 <- rast("raw/disturb_recovery_Prosper/2024.tif") 

crs(dist_2018) <- "EPSG:5514"
crs(dist_2019) <- "EPSG:5514"
crs(dist_2020) <- "EPSG:5514"

str(dist_2018)
str(dist_2019)
str(dist_2020)
#rasters <- list.files("raw/disturb_recovery_Prosper", pattern = "\\.tif$", full.names = TRUE)

# keep only 2018-2020, make raster stack
raster_stack <- terra::rast(list(dist_2018, 
                                 dist_2019,
                                 dist_2020))
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

# reproject vector data into 5514
dat_2023_5514 <- project(dat_2023, crs(ref))

plot(ref)
plot(stacked[1])
plot(rasters[1])
plot(dat_2023_5514, add =T)

vals18 <- terra::extract(dist_2018, dat_2023_5514, bind = TRUE) 
vals19 <- terra::extract(dist_2019, dat_2023_5514, bind = TRUE) 
vals20 <- terra::extract(dist_2020, dat_2023_5514, bind = TRUE) 

# merge data:
vals_temp <- merge(vals18, vals19, by = "ID")
vals <- merge(vals_temp, vals20, by = "ID")

# Assuming your data frame is named `vals`
vals_long <- vals %>%
  data.frame() %>% 
  pivot_longer(cols = c(`2018`, `2019`, `2020`), 
               names_to = "Year", 
               values_to = "Value") %>%
  filter(Value == 1) %>%
  group_by(ID) %>%
  slice_min(order_by = Year, with_ties = FALSE) %>%
  ungroup()

table(vals_long$Year)
