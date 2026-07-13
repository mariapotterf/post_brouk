# --- PACKAGES ---------------------------------------------------------------
library(terra)
library(dplyr)
library(stringr)
library(tibble)
library(readr)
library(fs)
library(tidyr)
library(data.table)
library(writexl)  # export xlsx for Roman
library(rasterVis)
library(RColorBrewer)
library(terra)

#install.packages("RCzechia")
library(RCzechia)    # get spatial data for czechia: Jirka Lacko https://rczechia.jla-data.net/

# Get cities as POINTS
pts <- RCzechia::obce_body()

# Define a vector of major cities to keep
major_cities <- c(
  "Praha",           # Prague
  "Brno",            # Brno
  "Ostrava",         # Ostrava
  "Plzeň",           # Plzen
  # "Liberec",         # Liberec
  # "Olomouc",         # Olomouc
  # "Ústí nad Labem",  # Usti
  # "Hradec Králové",  # Hradec
  "České Budějovice"#,# Ceske Budejovice
  #"Pardubice"        # Pardubice
)

# Filter to just these cities
pts_major <- pts %>%
  filter(NAZ_OBEC %in% major_cities) %>%
  arrange(NAZ_OBEC)


st_write(pts_major, "raw/CR_administrativa/major_czech_cities.gpkg", delete_layer = TRUE)
st_write(pts, "raw/CR_administrativa/cities.gpkg", delete_layer = TRUE)


# Get line geometries of rivers
rivers <- RCzechia::reky()

# Define main rivers of interest
main_rivers <- c("Vltava", "Labe", "Morava")  # Optionally add Morava, Dyje, etc.

# Filter by river name (NAZEV)
rivers_main <- rivers %>%
  filter(NAZEV %in% main_rivers) %>%
  arrange(NAZEV)

# Write to GeoPackage
st_write(rivers_main, "raw/CR_administrativa/main_czech_rivers.gpkg", delete_layer = TRUE)
st_write(rivers, "raw/CR_administrativa/rivers.gpkg", delete_layer = TRUE)

# Get line geometries of roads
roads <- silnice()
st_write(roads,       "raw/CR_administrativa/roads.gpkg",            delete_layer = TRUE)


# 1. Load woodland areas
woodlands <- lesy()
st_write(woodlands,       "raw/CR_administrativa/woodlands.gpkg",            delete_layer = TRUE)


# 2. Load protected natural areas
prot_zones <- chr_uzemi()
st_write(prot_zones,       "raw/CR_administrativa/prot_zones.gpkg",            delete_layer = TRUE)






# --- INPUTS -----------------------------------------------------------------
spei_folder <- "raw/SPEI12_new/spei12_harg1_split"  # 
okres_shp   <- terra::vect("raw/CR_administrativa/OKRESY_P.shp")



# --- READ FILES ---
files <- list.files(spei_folder, pattern = "\\.tif$", full.names = TRUE)
meta <- tibble(file = files) |>
  mutate(
    name = basename(file),
    year = as.integer(str_extract(name, "(?<=_)(\\d{4})(?=_)")),
    month = as.integer(str_extract(name, "(?<=_\\d{4}_)(\\d{2})(?=_)")),
    day = as.integer(str_extract(name, "(?<=_\\d{4}_\\d{2}_)(\\d{2})(?=\\.)"))
  ) |>
  dplyr::filter(!is.na(year), !is.na(month)) |>
  arrange(year, month, day)

# --- PREPARE VECTORS & RASTER TEMPLATE ---
# Convert SpatVector attribute table to tibble
okres_attrib <- as_tibble(okres_shp) |>
  mutate(okres_id = row_number())

# Assign new attribute table back to SpatVector
okres_shp$okres_id <- okres_attrib$okres_id

# Now it's safe to use `okres_shp` with the new field
v_ok <- okres_shp

r_template <- rast(meta$file[1])

if (!terra::same.crs(r_template, v_ok)) {
  v_ok <- terra::project(v_ok, r_template)
}


# --- RASTERIZE OKRES POLYGONS ---
r_okres <- terra::rasterize(v_ok, r_template, field = "okres_id")

# --- STACK RASTERS & EXTRACT VALUES ---
r_stack <- rast(meta$file)
names(r_stack) <- paste0(meta$year, "_", meta$month)

# combine raster stack with rasterized okres
r_all <- c(r_okres, r_stack)
v_all <- values(r_all)

# Convert matrix to data.table directly
dt <- as.data.table(v_all)

# Name columns - columns are properly named
#setnames(dt, 1, "okres_id")
#date_names <- paste0(meta$year, "_", meta$month)
#setnames(dt, 2:ncol(dt), date_names)

# Filter NA and melt to long format
dt_filt <- dt[!is.na(okres_id)]  # remove NA values
dt_long <- melt(dt_filt,
                id.vars = "okres_id",
                variable.name = "date",
                value.name = "value"
)

# Split date into year and month
#dt_long[, c("year", "month") := tstrsplit(date, "_", fixed = TRUE)]
#dt_long[, `:=`(
#  year = as.integer(year),
#  month = as.integer(month)
#)]

# Fast zonal summary
zonal_stats <- dt_long[
  !is.na(value),
  .(
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    min = min(value),
    max = max(value)
  ),
  by = .(okres_id, date)
]

# add new columns
# Split date into components
zonal_stats[, c("year", "month") := tstrsplit(date, "_", fixed = TRUE)]

# Convert to integer first
zonal_stats[, `:=`(
  year = as.integer(year),
  month = as.integer(month)
)]

# Now safely reformat the 'date' column with leading zeros
zonal_stats[, date := sprintf("%04d_%02d", year, month)]

# --- JOIN ATTRIBUTES ---
okres_attrib <- as_tibble(v_ok) |> dplyr::select(okres_id, NAZEV)
df_final <- left_join(zonal_stats, okres_attrib, by = "okres_id") |>
  relocate(NAZEV, .after = okres_id) |>
  arrange(year, month, okres_id)



# get values per year
# Convert to data.table if not already
df_annual <- as.data.table(df_final)[
  , .(
    country_mean = mean(mean, na.rm = TRUE),
    country_sd   = sd(mean, na.rm = TRUE),
    country_median = median(mean, na.rm = TRUE),
    country_min = min(mean, na.rm = TRUE),
    country_max = max(mean, na.rm = TRUE),
    n_okres = .N  # number of contributing okres values
  ),
  by = year
][order(year)]


# Save to CSV and XLSX
write_csv(df_annual, "outTable/country_SPEI12_annual_stats.csv")
write_xlsx(df_annual, "outTable/country_SPEI12_annual_stats.xlsx")


# --- SAVE ---
write_csv(df_final, "outTable/okres_SPEI12_monthly_stats.csv")
write_xlsx(df_final, "outTable/okres_SPEI12_monthly_stats.xlsx")



# --- EXPORT RASTERS OF ZONAL MEANS -------------------------------------------

# average per month ----------------
# Make sure output folder exists
#dir_create("outData/SPEI12_zonal")

# Convert df_final to data.table
dt_zonal <- as.data.table(df_final)

# Loop over each year-month and create a raster
unique_dates <- unique(dt_zonal[, .(year, month)])

unique_dates[, `:=`(
  year = as.integer(year),
  month = as.integer(month)
)]

for (i in seq_len(nrow(unique_dates))) {
  yr <- unique_dates$year[i]
  mo <- unique_dates$month[i]
  date_label <- sprintf("%04d_%02d", as.integer(yr), as.integer(mo))
  #print(date_label)
  # Subset values for this date
  vals <- dt_zonal[year == yr & month == mo, .(okres_id, mean)]
  
  # Create lookup vector: index = okres_id, value = mean
  lookup <- rep(NA_real_, max(r_okres[], na.rm = TRUE))
  lookup[vals$okres_id] <- vals$mean
  
  # Replace cell values: map okres_id to mean value
  r_out <- classify(r_okres, 
                    matrix(c(1:length(lookup), lookup), 
                           ncol = 2))
  
  # Write raster to file
  out_path <- sprintf("outData/SPEI12_zonal/spei12_mean_%s.tif", date_label)
  # print(out_path)  # ccheck if values corresponds with 
  writeRaster(r_out, out_path, overwrite = TRUE)
}

## avg per years -----------------------------------------

# Make sure output folder exists
#dir_create("outData/SPEI12_zonal_year")

# 1. Compute yearly mean SPEI per okres
df_zonal_year <- df_final %>% 
  group_by(okres_id, year) %>% 
  summarise(mean_spei = mean(mean, na.rm = TRUE), .groups = "drop")

# 2. Loop over years to create rasters
unique_years <- sort(unique(df_zonal_year$year))

for (yr in unique_years) {
  message("Processing year: ", yr)
  
  # Subset for the year
  vals <- df_zonal_year %>% 
    filter(year == yr) %>% 
    select(okres_id, mean_spei)
  
  # Create lookup vector for classification
  lookup <- rep(NA_real_, max(r_okres[], na.rm = TRUE))
  lookup[vals$okres_id] <- vals$mean_spei
  
  # Map okres_id values in raster to mean SPEI values
  r_out <- classify(r_okres, 
                    matrix(c(seq_along(lookup), lookup), ncol = 2))
  
  # Output file path
  out_path <- sprintf("outData/SPEI12_zonal/spei12_year_mean_%04d.tif", yr)
  writeRaster(r_out, out_path, overwrite = TRUE)
}







# make a gif from zonal stats ---------------------------


library(terra)
library(magick)

# Diverging color palette (red-white-blue)
spei_palette <- colorRampPalette(c("red", "white", "blue"))

# Normalize to consistent limits
raster_limits <- c(-3, 3)  # typical SPEI range


#  test how to fix limits  on single file ------------------------
r <- rast("outData/SPEI12_zonal/spei12_mean_2006_06.tif")

# using base R: 

#png_file <- tempfile(fileext = ".png")
#png(png_file, width = 600, height = 600)
plot(r, col = spei_palette(100), #zlim = c(-1,1),
     #main = paste("SPEI Mean", date_label),
     axes = FALSE, box = FALSE)

plot(r,
     col = spei_palette(100),
     breaks = seq(-3, 3, length.out = 101),  # exactly 100 intervals
     #main = paste("SPEI Mean", date_label),
     axes = FALSE, box = FALSE)

# using Rastervis - seems working but slow!
levelplot(r,
          margin = FALSE,
          col.regions = spei_palette(100),
          at = seq(-3, 3, length.out = 101),  # fixed color breaks
          scales = list(draw = FALSE),  # <-- Hides x and y axis
          #main = paste("SPEI Mean", date_label),
          colorkey = list(
            space = "right",
            labels = list(at = seq(-3, 3, by = 1), labels = seq(-3, 3, by = 1))
          ))



# Run animation using RasterVis - per months and years, fixed scale ------------

# Define palette
spei_palette <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))  # red = dry, blue = wet

make_plot_image <- function(r_path, date_label) {
  library(terra)
  library(raster)       # required for `raster()` conversion
  library(rasterVis)
  library(magick)
  
  # Read and convert raster
  r_terra <- rast(r_path)
  r <- raster::raster(r_terra)  # proper conversion to RasterLayer
  
  # Temp PNG file
  png_file <- tempfile(fileext = ".png")
  png(png_file, width = 800, height = 800)
  
  # Plot with fixed z-scale and no axes
  plt <- levelplot(
    r,
    margin = FALSE,
    col.regions = spei_palette(100),
    at = seq(-3, 3, length.out = 101),
    scales = list(draw = FALSE),
    main = paste("SPEI Mean", date_label),
    colorkey = list(
      space = "right",
      labels = list(at = seq(-3, 3, by = 1), labels = seq(-3, 3, by = 1))
    )
  )
  
  print(plt)
  dev.off()
  
  magick::image_read(png_file)
}



# List and sort raster files
raster_files <- list.files("outData/SPEI12_zonal", pattern = "^spei12_mean_\\d{4}_\\d{2}\\.tif$", full.names = TRUE)
raster_files <- sort(raster_files)

# Extract date labels
date_labels <- str_extract(basename(raster_files), "\\d{4}_\\d{2}")

# Generate image frames
frames <- mapply(make_plot_image, raster_files, date_labels, SIMPLIFY = FALSE)

# Create animation
animation <- image_animate(image_join(frames), fps = 5)

# Save GIF
image_write(animation, "outData/SPEI12_zonal/spei12_zonal_mean_animation2.gif")



# Export gifs for yearly values ------------------------------------------------


df <- readr::read_csv("outTable/okres_SPEI12_monthly_stats.csv")

df_yearly <- df %>%
  mutate(year = as.integer(year)) %>%
  group_by(okres_id, year) %>%
  summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")

# --- RECREATE YEARLY RASTERS ------------------------------------------------
# Load rasterized okres base
#r_okres <- rast("your_saved_rasterized_okres.tif")  # or use your in-memory `r_okres`

# Get unique year list
years <- sort(unique(df_yearly$year))

# Create raster for each year
r_list <- lapply(years, function(yr) {
  df_yr <- df_yearly %>% filter(year == yr)
  
  # Create full vector of values (fill with NA)
  vals <- rep(NA_real_, ncell(r_okres))
  match_ids <- match(values(r_okres), df_yr$okres_id)
  valid <- !is.na(match_ids)
  vals[valid] <- df_yr$mean[match_ids[valid]]
  
  # Create raster
  r_out <- setValues(r_okres, vals)
  names(r_out) <- as.character(yr)
  return(r_out)
})

r_stack_yearly <- rast(r_list)

make_plot_image <- function(r, year_label) {
  r_raster <- raster::raster(r)
  png_file <- tempfile(fileext = ".png")
  png(png_file, width = 800, height = 800)
  
  plt <- levelplot(
    r_raster,
    margin = FALSE,
    col.regions = spei_palette(100),
    at = seq(-3, 3, length.out = 101),
    scales = list(draw = FALSE),
    main = paste("Yearly SPEI Mean", year_label),
    colorkey = list(
      space = "right",
      labels = list(at = seq(-3, 3, by = 1), labels = seq(-3, 3, by = 1))
    )
  )
  
  print(plt)
  dev.off()
  magick::image_read(png_file)
}

# --- GENERATE PLOTS & ANIMATE -----------------------------------------------
frames <- mapply(make_plot_image, r_list, years, SIMPLIFY = FALSE)
animation <- image_animate(image_join(frames), fps = 2)
image_write(animation, "outData/SPEI12_zonal/spei12_zonal_yearly_mean.gif")

