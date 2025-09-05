# -----------------------------------------
#
#             Process drone imagery 
#
#  -----------------------------------------

rm(list = ls())

# extract height information per square (2x2 m subplot)/circle on teh subplot level
# get CV, height, min-max


library(terra)
library(sf)
library(dplyr)
library(purrr)
library(data.table)

# --- Read ---------------------------------------------------------------------
dat23_sf        <- st_read("outData/sf_context_2023.gpkg")   # subplots (points or polygons)
drone23_24_sf   <- st_read("raw/position_drone.gpkg")        # drone footprints/points
chm_folder      <- "raw/UAV_Images"
chm_files       <- list.files(chm_folder, pattern = "^CHM_.*\\.tif$", full.names = TRUE, recursive = TRUE)
chm_rasters     <- set_names(chm_files, tools::file_path_sans_ext(basename(chm_files))) |> map(terra::rast)
krovak_crs      <- "EPSG:5514"

# Force CRS
chm_rasters <- map(chm_rasters, \(x){ crs(x) <- krovak_crs; x })
target_crs  <- crs(chm_rasters[[1]])
dat23_5514  <- st_transform(dat23_sf, target_crs)
drone_5514  <- st_transform(drone23_24_sf, target_crs) |> rename(drone_year = year)

# If subplots are polygons, use centroids for buffer centers
if (!inherits(dat23_5514$geometry, "sfc_POINT")) {
  dat23_pts <- st_centroid(dat23_5514)
} else {
  dat23_pts <- dat23_5514
}

# Ensure we have an ID column for plotting/grouping
dat23_pts <- dat23_pts |> rename(plot_ID = ID)

# --- Helpers ------------------------------------------------------------------
make_square_buffers <- function(points, side_m) {
  if (inherits(points, "sf")) points <- terra::vect(points)
  
  half <- side_m / 2
  xy <- terra::crds(points)
  polys <- terra::vect(lapply(seq_len(nrow(xy)), function(i) {
    e <- terra::ext(xy[i,1]-half, xy[i,1]+half,
                    xy[i,2]-half, xy[i,2]+half)
    terra::as.polygons(e, crs = terra::crs(points))
  }))
  
  polys$plot_ID  <- points$plot_ID
  if ("uav_ID" %in% names(points)) polys$uav_ID <- points$uav_ID
  polys$buf_side <- side_m
  polys$poly_row <- seq_len(nrow(polys))
  polys
}



# (optional) attach uav_ID to points by intersection with drone footprints/points
attach_uav_id_by_intersection <- function(points_sf, drone_sf, uav_id_col = "uav_ID") {
  inter <- sf::st_intersection(points_sf, drone_sf)
  if (nrow(inter) == 0) return(points_sf |> mutate(!!uav_id_col := NA_character_))
  # keep one (if multiple matches, choose first)
  inter <- inter |> st_drop_geometry() |> group_by(across(everything()), .add = FALSE) |> ungroup()
  # simplest: spatial join & pick first match
  joined <- sf::st_join(points_sf, drone_sf, left = TRUE)
  # assume drone layer has a column 'uav_ID'
  joined <- joined |> mutate(uav_ID = if ("uav_ID" %in% names(joined)) uav_ID else NA_character_)
  joined
}

# --- 1) Buffers for ALL subplots ---------------------------------------------
# (Option A) keep uav_ID as NA for non-overlapping subplots
dat23_pts_all <- dat23_pts

# (Option B, optional) assign uav_ID where available via spatial join
# dat23_pts_all <- attach_uav_id_by_intersection(dat23_pts, drone_5514)

squares_all <- make_square_buffers(terra::vect(dat23_pts_all), side_m = 2)
sq_all_sf   <- sf::st_as_sf(squares_all)

st_write(
  dplyr::select(sq_all_sf, plot_ID, dplyr::any_of("uav_ID"), buf_side),
  "outData/square_buffers_2m_all.gpkg",
  layer = "square_buffers_2m_all",
  delete_dsn = TRUE
)

# --- 2) Buffers for ONLY subplots overlapping drones (if you still want this) -
dat23_overlap <- sf::st_filter(dat23_pts, drone_5514, .predicate = st_intersects)
# If your drone layer has uav_ID, bring it in:
dat23_overlap <- sf::st_join(dat23_overlap, dplyr::select(drone_5514, uav_ID), left = TRUE)

squares_overlap <- make_square_buffers(dat23_overlap, side_m = 2)
sq_ov_sf        <- sf::st_as_sf(squares_overlap)

st_write(
  dplyr::select(sq_ov_sf, plot_ID, dplyr::any_of("uav_ID"), buf_side),
  "outData/square_buffers_2m_drone.gpkg",
  layer = "square_buffers_2m_drone",
  append = TRUE
)





# --- Build squares for multiple sizes -----------------
dat23_subset <- st_intersection(dat23_5514, drone_5514) %>%
  mutate(field_year = 2023,
         x = st_coordinates(.)[, 1],
         y = st_coordinates(.)[, 2]) %>% 
  mutate(uav_ID = paste0("CHM_", uav_ID)) %>% 
  rename(plot_ID = ID) %>% 
  st_as_sf(coords = c("x", "y"), crs = target_crs)

cluster_lookup <- dat23_subset %>%
  st_drop_geometry() %>%
  dplyr::select(plot_ID, uav_ID) %>%
  distinct()




cluster_vect <- terra::vect(dat23_subset)

square_sides <- c(2, 2.5, 5)   # meters; add more if needed

all_cluster_heights <- list()
all_cluster_raw     <- list()

for (side in square_sides) {
  message("Building squares with side = ", side, " m")
  
  squares <- make_square_buffers(cluster_vect, side_m = side)
  
  drone_stats <- lapply(names(chm_rasters), function(drone_id) {
    r <- chm_rasters[[drone_id]]
    value_col <- names(r)[1]
    
    matching_plots <- cluster_lookup %>%
      filter(.data$uav_ID == drone_id) %>%
      pull(.data$plot_ID)
    
    squares_this_drone <- squares[squares$plot_ID %in% matching_plots, ]
    if (nrow(squares_this_drone) == 0) return(NULL)
    
    vals <- terra::extract(r, squares_this_drone, ID = TRUE, cells = TRUE)
    if (is.null(vals) || nrow(vals) == 0) return(NULL)
    
    vals$plot_ID     <- squares_this_drone$plot_ID[vals$ID]
    vals$drone       <- drone_id
    vals$buffer_size <- side
    
    vals_tbl <- tibble::as_tibble(vals)
    
    # --- Raw values
    vals_raw <- vals_tbl %>%
      mutate(pixel_value = .data[[value_col]]) %>%
      select(plot_ID, drone, buffer_size, ID, cell, pixel_value)
    
    # --- Summary stats
    vals_stats <- vals_tbl %>%
      group_by(drone, buffer_size, plot_ID) %>%
      summarize(
        n_cells_total = n_distinct(cell[!is.na(cell)]),
        n_cells_nonNA = n_distinct(cell[!is.na(.data[[value_col]]) & !is.na(cell)]),
        n_cells_NA    = n_cells_total - n_cells_nonNA,
        mean_height   = mean(.data[[value_col]], na.rm = TRUE),
        median_height = median(.data[[value_col]], na.rm = TRUE),
        sd_height     = sd(.data[[value_col]], na.rm = TRUE),
        max_height    = suppressWarnings(max(.data[[value_col]], na.rm = TRUE)),
        min_height    = suppressWarnings(min(.data[[value_col]], na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      filter(n_cells_total > 0) %>%
      mutate(
        mean_height   = if_else(n_cells_nonNA > 0, mean_height,   NA_real_),
        median_height = if_else(n_cells_nonNA > 0, median_height, NA_real_),
        sd_height     = if_else(n_cells_nonNA > 0, sd_height,     NA_real_),
        max_height    = if_else(n_cells_nonNA > 0, max_height,    NA_real_),
        min_height    = if_else(n_cells_nonNA > 0, min_height,    NA_real_),
        cv_height     = if_else(n_cells_nonNA > 0 & mean_height != 0,
                                sd_height / mean_height, NA_real_),
        range_height  = if_else(n_cells_nonNA > 0 &
                                  is.finite(max_height) & is.finite(min_height),
                                max_height - min_height, NA_real_)
      )
    
    return(list(stats = vals_stats, raw = vals_raw))
  })
  
  # Bind results for this buffer size
  cluster_stats <- map(drone_stats, "stats") %>% compact() %>% bind_rows()
  cluster_raw   <- map(drone_stats, "raw")   %>% compact() %>% bind_rows()
  
  all_cluster_heights[[as.character(side)]] <- cluster_stats
  all_cluster_raw[[as.character(side)]]     <- cluster_raw
}

# --- Combine across sizes
chm_summary_multi <- bind_rows(all_cluster_heights)
chm_raw_multi     <- bind_rows(all_cluster_raw)

# save output
fwrite(chm_summary_multi, "outTable/chm_buff_summary.csv")
fwrite(chm_raw_multi,     "outTable/chm_buff_raw.csv")
