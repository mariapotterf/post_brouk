
# process drone imagery

# extract height information per square (2x2 m subplot)/circle on teh subplot level
# get CV, height, min-max

gc()

library(terra)
library(dplyr)
library(purrr)
library(sf)
library(data.table)

# Read subplot data for 2023
dat23_sf   <- st_read("outData/sf_context_2023.gpkg") # data per subplot
drone23_24_sf <- st_read("raw/position_drone.gpkg") # data per subplot

# Read drone CHM rasters from 2023&2024
chm_folder <- "raw/UAV_Images" # /2023
chm_files <- list.files(chm_folder, pattern = "^CHM_.*\\.tif$", full.names = TRUE, recursive = TRUE)

# drone heights are in meters!

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
         y = st_coordinates(.)[, 2]) %>% 
  mutate(uav_ID = paste0("CHM_", uav_ID)) %>% 
  rename(plot_ID = ID) %>% 
  st_as_sf(coords = c("x", "y"), crs = target_crs)  # U

# Start fresh from dat23_5514 — raw points with correct geometry
# dat23_clusters <- dat23_subset %>%
#   as.data.frame() %>% 
#   group_by(cluster) %>%
#   summarise(
#     x = mean(x),
#     y = mean(y),
#     drone_year  = first(drone_year),
#     field_year  = first(field_year),
#     uav_ID = first(uav_ID),
#     .groups = "drop"
#   ) %>%
#   mutate(uav_ID = paste0("CHM_", uav_ID)) %>% 
#   st_as_sf(coords = c("x", "y"), crs = target_crs)  # Use correct CRS

# Convert to terra vector
cluster_vect <- vect(ungroup(dat23_subset))

data.frame(cluster_vect)

#  Make a square buffer TEST STARTS ----------------------
#library(terra)

# Define square side length
buffer_width <- 2
half_width <- buffer_width / 2

# Check input
stopifnot(inherits(cluster_vect, "SpatVector"))
stopifnot(geomtype(cluster_vect) == "points")

# Create square buffers as bounding boxes around each point
squares <- vect(lapply(1:nrow(cluster_vect), function(i) {
  pt <- cluster_vect[i, ]
  ext <- ext(
    crds(pt)[1, 1] - half_width,
    crds(pt)[1, 1] + half_width,
    crds(pt)[1, 2] - half_width,
    crds(pt)[1, 2] + half_width
  )
  as.polygons(ext, crs = crs(cluster_vect))
}))

# Add attributes
squares$ID <- seq_len(nrow(squares))
squares$plot_ID <- cluster_vect$plot_ID

# Test run: extract square for cluster "26_111"
r <- chm_rasters$CHM_uav6
buf <- squares[squares$plot_ID == "13_15_104_1", ]

plot(buf)
plot(r, add = T)
plot(buf, add = T, col = "red")
plot(r)


# make square buffer ------------------------------------------------------------

# Lookup from the input points (keeps your original mapping)
cluster_lookup <- as.data.frame(cluster_vect)[, c("plot_ID", "uav_ID")]

# ---- Helper to build square polygons (axis-aligned) around points
make_square_buffers <- function(points, side_m) {
  half <- side_m / 2
  xy <- terra::crds(points)
  # build extents, then polygons
  polys <- vect(lapply(seq_len(nrow(xy)), function(i) {
    e <- terra::ext(xy[i, 1] - half, xy[i, 1] + half,
                    xy[i, 2] - half, xy[i, 2] + half)
    terra::as.polygons(e, crs = terra::crs(points))
  }))
  # carry attributes
  polys$plot_ID <- points$plot_ID
  polys$uav_ID  <- points$uav_ID
  polys$buf_side <- side_m
  polys$poly_row <- seq_len(nrow(polys))  # stable row index for mapping
  polys
}

# if you already have `squares`, skip this:
squares <- make_square_buffers(cluster_vect, side_m = 2)

# convert to sf for robust GPKG writing (multi-layer friendly)
sq_sf <- sf::st_as_sf(squares)

# ------------------------------------------------------------------
# Get all squares together 
# ------------------------------------------------------------------
gpkg_all <- "outData/square_buffers_2m.gpkg"

# overwrite the dataset if it exists
sf::st_write(
  obj   = dplyr::select(sq_sf, plot_ID, uav_ID, buf_side),
  dsn   = gpkg_all,
  layer = "square_buffers_2m",
  delete_dsn = TRUE
)


# --- Build 2 m squares once (you said one size, but vector-ready if you add more)
square_sides <- c(2)  # meters
all_cluster_heights <- list()
all_cluster_raw <- list()

for (side in square_sides) {
  message("Building squares with side = ", side, " m")
  squares <- make_square_buffers(cluster_vect, side_m = side)
  
  drone_stats <- lapply(names(chm_rasters), function(drone_id) {
    r <- chm_rasters[[drone_id]]
    value_col <- names(r)[1]
    
    matching_plots <- cluster_lookup %>%
      dplyr::filter(.data$uav_ID == drone_id) %>%
      dplyr::pull(.data$plot_ID)
    
    squares_this_drone <- squares[squares$plot_ID %in% matching_plots, ]
    if (nrow(squares_this_drone) == 0) return(NULL)
    
    vals <- terra::extract(r, squares_this_drone, ID = TRUE, cells = TRUE)
    if (is.null(vals) || nrow(vals) == 0) return(NULL)
    
    vals$plot_ID     <- squares_this_drone$plot_ID[vals$ID]
    vals$drone       <- drone_id
    vals$buffer_size <- side
    
    vals_tbl <- tibble::as_tibble(vals)
    
    vals_raw <- vals_tbl %>%
      dplyr::mutate(pixel_value = .data[[value_col]]) %>%
      dplyr::select(plot_ID, drone, buffer_size, ID, cell, pixel_value)
    
    vals_stats <- vals_tbl %>%
      dplyr::group_by(drone, buffer_size, plot_ID) %>%
      dplyr::summarize(
        n_cells_total = dplyr::n_distinct(cell[!is.na(cell)]),
        n_cells_nonNA = dplyr::n_distinct(cell[!is.na(.data[[value_col]]) & !is.na(cell)]),
        n_cells_NA    = n_cells_total - n_cells_nonNA,
        mean_height   = mean(.data[[value_col]], na.rm = TRUE),
        median_height = median(.data[[value_col]], na.rm = TRUE),
        sd_height     = sd(.data[[value_col]], na.rm = TRUE),
        max_height    = suppressWarnings(max(.data[[value_col]], na.rm = TRUE)),
        min_height    = suppressWarnings(min(.data[[value_col]], na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_cells_total > 0) %>%
      dplyr::mutate(
        mean_height   = dplyr::if_else(n_cells_nonNA > 0, mean_height,   NA_real_),
        median_height = dplyr::if_else(n_cells_nonNA > 0, median_height, NA_real_),
        sd_height     = dplyr::if_else(n_cells_nonNA > 0, sd_height,     NA_real_),
        max_height    = dplyr::if_else(n_cells_nonNA > 0, max_height,    NA_real_),
        min_height    = dplyr::if_else(n_cells_nonNA > 0, min_height,    NA_real_),
        cv_height     = dplyr::if_else(n_cells_nonNA > 0 & mean_height != 0,
                                       sd_height / mean_height, NA_real_),
        range_height  = dplyr::if_else(n_cells_nonNA > 0 &
                                         is.finite(max_height) & is.finite(min_height),
                                       max_height - min_height, NA_real_)
      )
    
    return(list(stats = vals_stats, raw = vals_raw))
  })
  
  # Bind both raw and summary results
  cluster_stats <- purrr::map(drone_stats, "stats") %>% purrr::compact() %>% dplyr::bind_rows()
  cluster_raw   <- purrr::map(drone_stats, "raw")   %>% purrr::compact() %>% dplyr::bind_rows()
  
  all_cluster_heights[[as.character(side)]] <- cluster_stats
  all_cluster_raw[[as.character(side)]]     <- cluster_raw
}


# --- Final table across (potentially) multiple square sizes
chm_summary_multi <- dplyr::bind_rows(all_cluster_heights)
chm_raw_multi <- dplyr::bind_rows(all_cluster_raw)


# save output -----------------------------------------
fwrite(chm_summary_multi, "outTable/chm_buff_summary.csv")
fwrite(chm_raw_multi, "outTable/chm_buff_raw.csv")


hist(chm_summary_multi$cv_height)
