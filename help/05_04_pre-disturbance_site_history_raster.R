
# 
#TEST RUN: pre disturbnace history - from rasters estimation
#



# ----------------------------------------------------------------------------
# Pre-disturbance site history from rasters, for multiple buffer widths
# ----------------------------------------------------------------------------
gc()

library(terra)
library(sf)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# -------------------- I/O --------------------
# 2023
dat23_sf <- sf::st_read("outData/sf_context_2023.gpkg")
# 2025
dat25_sf <- sf::st_read("outData/subplot_with_clusters_2025.gpkg")

# Rasters (2015)
tree_cover_dens15    <- rast("raw/forest_composition/TCD_2015_020m_CR.tif") # 0–100 %
dominant_leaf_type15 <- rast("raw/forest_composition/DLT_2015_020m_CR.tif") # 1 decid., 2 conif.

# Output folder
dir.create("outTable", showWarnings = FALSE, recursive = TRUE)

# -------------------- Prep points --------------------
target_crs <- crs(tree_cover_dens15)

dat23_sf <- st_transform(dat23_sf, target_crs) |> mutate(year = 2023)
dat25_sf <- st_transform(dat25_sf, target_crs) |> mutate(year = 2025)

# mean subplot coords per cluster → one point per cluster-year
to_cluster_points <- function(sfobj, target_crs) {
  sfobj <- sf::st_transform(sfobj, target_crs)
  
  # get coords
  coords <- sf::st_coordinates(sfobj)
  df <- sfobj |>
    dplyr::mutate(x = coords[,1], y = coords[,2]) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(cluster, year) |>
    dplyr::summarise(x = mean(x), y = mean(y), .groups = "drop") |>
    dplyr::mutate(cluster = as.character(cluster)) |>
    sf::st_as_sf(coords = c("x","y"), crs = target_crs)  # POINT geometry
  
  df
}

# build both sets with consistent types/geometries
dat23_clusters <- to_cluster_points(dat23_sf, target_crs)
dat25_clusters <- to_cluster_points(dat25_sf, target_crs)

# now this will work
combined_clusters <- dplyr::bind_rows(dat23_clusters, dat25_clusters)
cluster_vect <- vect(combined_clusters)  # terra SpatVector

# -------------------- Helpers --------------------
# Safe CV (avoid division by zero)
cv_safe <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA_real_)
  sd(x, na.rm = TRUE) / m
}

# Summaries for continuous raster (e.g., TCD)
summary_stats <- function(df, value_col) {
  df |>
    group_by(ID, cluster, year) |>
    summarise(
      mean   = mean(.data[[value_col]], na.rm = TRUE),
      median = median(.data[[value_col]], na.rm = TRUE),
      cv     = cv_safe(.data[[value_col]]),
      .groups = "drop"
    )
}

# Leaf-type (categorical) summaries: Shannon + shares
leaf_summaries <- function(vals_leaf) {
  # recode for readability (optional)
  vals_leaf_rec <- vals_leaf |>
    mutate(leaf_type = factor(DLT_2015_020m_CR, levels = c(1, 2),
                              labels = c("deciduous","coniferous")))
  
  # Shannon per buffer ID
  shannon_df <- vals_leaf_rec |>
    filter(!is.na(DLT_2015_020m_CR)) |>
    group_by(ID, cluster, year, leaf_type) |>
    tally(name = "n") |>
    group_by(ID, cluster, year) |>
    mutate(p = n / sum(n)) |>
    summarise(shannon = -sum(p * log(p)), .groups = "drop")
  
  # Shares per buffer ID
  share_df <- vals_leaf |>
    group_by(ID, cluster, year) |>
    summarise(
      total_cells   = n(),
      n_coniferous  = sum(DLT_2015_020m_CR == 2, na.rm = TRUE),
      n_deciduous   = sum(DLT_2015_020m_CR == 1, na.rm = TRUE),
      forest_cells  = n_coniferous + n_deciduous,
      share_coniferous = ifelse(forest_cells > 0, n_coniferous / forest_cells, NA_real_),
      share_deciduous  = ifelse(forest_cells > 0, n_deciduous  / forest_cells, NA_real_),
      .groups = "drop"
    )
  
  left_join(share_df, shannon_df, by = c("ID","cluster","year"))
}

# Core runner for one buffer width
run_for_width <- function(width_m) {
  message("• Processing buffer = ", width_m, " m")
  
  # Build buffers and lookup
  buf <- buffer(cluster_vect, width = width_m)
  lookup <- as.data.frame(buf) |>
    mutate(ID = row_number()) |>
    select(ID, cluster, year)
  
  # Extract rasters
  vals_cover <- terra::extract(tree_cover_dens15, buf, ID = TRUE) |>
    left_join(lookup, by = "ID")
  vals_leaf  <- terra::extract(dominant_leaf_type15, buf, ID = TRUE) |>
    left_join(lookup, by = "ID")
  
  # Summaries
  cover_stats <- summary_stats(vals_cover, "TCD_2015_020m_CR") |>
    rename(
      mean_cover_dens   = mean,
      median_cover_dens = median,
      cv_cover_dens     = cv
    )
  
  leaf_stats <- leaf_summaries(vals_leaf)
  
  combined <- left_join(cover_stats, leaf_stats, by = c("ID","cluster","year")) |>
    mutate(buffer_m = width_m) |>
    relocate(buffer_m, .before = everything())
  
  # Save per-width table
  out_file <- file.path("outTable", sprintf("pre_disturb_history_raster_w%sm.csv", width_m))
  fwrite(combined, out_file)
  message("  └─ saved: ", out_file)
  
  combined
}

# -------------------- Run for multiple radii --------------------
buffer_set <- c(25, 50, 100, 175, 300, 500, 1000)

all_buffers_tbl <- map_dfr(buffer_set, run_for_width)

# Save combined
fwrite(all_buffers_tbl, "outTable/pre_disturb_history_raster_ALL_buffers.csv")



# --------------------- Sensitivity analysis --------------------

# ---- Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(data.table)

all_buffers_tbl <- fread("outTable/pre_disturb_history_raster_ALL_buffers.csv") |> as.data.frame()

# Choose predictors to inspect (names from your combined table)
pred_vars <- c("mean_cover_dens", "cv_cover_dens", "share_coniferous", 
               "shannon") # , "total_cells" 
# cv of tree-cover density, conifer share (~spruce share proxy), leaf-type Shannon

# Long format
sens_long <- all_buffers_tbl %>%
  dplyr::select(cluster, year, buffer_m, all_of(pred_vars)) %>%
  tidyr::pivot_longer(cols = all_of(pred_vars),
                      names_to = "variable", values_to = "value") #%>%
  # mutate(
  #   variable = recode(variable,
  #                     cv_cover_dens   = "Cover density CV",
  #                     share_coniferous = "Coniferous share (%)",
  #                     shannon          = "Leaf-type Shannon"
  #   ),
  #   # Express share as %
  #   value = ifelse(variable == "Coniferous share (%)", value * 100, value)
  # )

# Summary across clusters (median + IQR) by buffer and year
sens_sum <- sens_long %>%
  group_by(variable, buffer_m) %>%
  summarise(
    n = sum(!is.na(value)),
    med = median(value, na.rm = TRUE),
    p25 = quantile(value, 0.25, na.rm = TRUE),
    p75 = quantile(value, 0.75, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

# mean +- SD
sens_sum %>% 
  ggplot( aes(x=buffer_m, y=mean, group=variable)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  geom_line(linewidth = 0.5) + 
  geom_point(size = 1.5) + 
  facet_wrap(.~variable, scale = 'free_y' ) +
  theme_classic()
  


# medians% IQR
sens_sum %>% 
  ggplot( aes(x=buffer_m, y=med, group=variable)) +
  geom_errorbar(aes(ymin=p25, ymax=p75), width=.2,
                position=position_dodge(0.05)) +
  geom_line(linewidth = 0.5, color = 'grey') + 
  geom_point(size = 1.5) + 
  facet_wrap(.~variable, scale = 'free_y' ) +
  theme_classic()

