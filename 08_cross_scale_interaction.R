
# ----------------------------------------------------------------------------
# Cross-scale interaction: field vs drone variation in vertical structures
# -----------------------------------------------------------------------------

# read data: 
# - field data from 2023 (2025)
# - drone 2023 (2024?)

# - get sites history - from tree density cover, dominant leaf type
# compare stem density, dominant tree species, species diversity
# vertical and horizontal structure across sscales: fiedl vs drone

# identify damaged terminal: from raw data czechia 20203 - export gpkg
# identify types of damages: 2025 - export gpkg
# Drones: 
# - get vertical strucure

gc()

library(terra)
library(sf)
#library(DBI)
library(ggplot2)
#library(dbscan)
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggpubr)


library(corrplot)
library(GGally)

# Read files --------------------------------------

# 2023 
dat23_subplot <- fread('outData/subplot_full_2023.csv')
#dat23_cluster <- fread('outData/df_cluster_2023.csv')

# Save spatial subplot data with cluster IDs
dat23_sf <- st_read("outData/sf_context_2023.gpkg")


# 2025
#dat25_subplot <- fread('outData/subplot_full_2025.csv')
#dat25_cluster <- fread('outData/df_cluster_2025.csv')

# Save spatial subplot data with cluster IDs
#dat25_sf <- st_read("outData/subplot_with_clusters_2025.gpkg")

# Read drone CHM summary
drone_cv <- fread("outTable/chm_buff_summary.csv")

# read pre-disturbance site history
pre_dist_history <- fread("outTable/pre_disturb_history_raster.csv")

pre_dist_history <- pre_dist_history %>%
  mutate(field_year = if_else(str_detect(cluster, "_"), 2023, 2025))

pre_dist_history_2023 <- pre_dist_history %>%   # filed only pre dicturbancs history for sites collected in 2023
  dplyr::filter(field_year == "2023")


# Sumarize field observation per cluster  --------------------------------------

# guestimate dbh and ba per individual based on height distribution -------------------------
dat23_subplot_recode <- dat23_subplot %>% 
  mutate(
    dbh_est = case_when(
      vegtype == "mature"           ~ case_when(
        dbh == "10–20cm" ~ 15.0,
        dbh == "20–40cm" ~ 30.0,
        dbh == "40–60cm" ~ 50.0,
        TRUE             ~ NA_real_
      ),
      hgt == "0.2–0.4"              ~ 0.2,
      hgt == "0.4–0.6"              ~ 0.4,
      hgt == "0.6–0.8"              ~ 0.7,
      hgt == "0.8–1.0"              ~ 1.0,
      hgt == "1.0–1.3"              ~ 1.5,
      hgt == "1.3–2.0"              ~ 2.5,
      hgt == "2–4"                  ~ 4.5,
      hgt == ">4"                   ~ 7.0,
      TRUE                          ~ NA_real_
    ),
    basal_area_cm2 = pi * (dbh_est / 2)^2
  ) %>% 
  # calculate basal area from counts per cm2, per m2, and scale up to ha
  mutate(
    ba_total_cm2   = basal_area_cm2 * n,
    ba_total_m2    = ba_total_cm2 / 10000,  # convert to m²
    ba_ha_m2       = ba_total_m2 * scaling_factor  # scale to per hectare
  )




# get stem density and shannon for field data 
df_cluster <- dat23_subplot_recode %>% 
  group_by(cluster) %>% 
  summarise(stem_density = sum(stem_density, na.rm = T),
            basal_area_ha_m2 = sum(ba_ha_m2, na.rm = TRUE))



# shannon estimation:
shannon_species <- dat23_subplot_recode %>%
  dplyr::filter(!is.na(species), ba_ha_m2 > 0) %>%
  group_by(cluster, species) %>%
  summarise(ba = sum(ba_ha_m2), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(p = ba / sum(ba)) %>%
  summarise(shannon_species = -sum(p * log(p)), .groups = "drop")

shannon_height <- dat23_subplot_recode %>%
  dplyr::filter(hgt != "", ba_ha_m2 > 0) %>%
  group_by(cluster, hgt) %>%
  summarise(ba = sum(ba_ha_m2), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(p = ba / sum(ba)) %>%
  summarise(shannon_height = -sum(p * log(p)), .groups = "drop")


df_field_diversity <- df_cluster %>%
  left_join(shannon_species, by = "cluster") %>%
  left_join(shannon_height, by = "cluster") %>% 
  replace_na(list(
    shannon_species = 0,
    shannon_height  = 0
  ))




# make working example: 

dd <- data.frame(
  cluster = c(1, 1,   
              2, 2, 2,  
              3,   
              4, 
              5),
  species = c("a", "b",
              "a", "b", "d",
              NA, 
              "d",
              "d"),
  basal_area = c(5, 5,  
                 1, 1, 8,  
                 0,  
                 10, 
                 0.5)
)

# merge field data with drone estimation ------------------------------

# Sensitivity analysis ----------

# Function to generate boxplot for a given variable
plot_buffer_variation <- function(data, var_name, y_label) {
  ggplot(data, aes(x = factor(buffer_size), y = .data[[var_name]])) +
    geom_boxplot(fill = "skyblue") +
    labs(
      title = paste("", y_label, ""),
      x = "Buffer size (m)",
      y = y_label
    ) +
    theme_minimal()
}

# Generate plots
p_mean   <- plot_buffer_variation(drone_cv, "mean_height", "Mean Height")
p_median <- plot_buffer_variation(drone_cv, "median_height", "Mean Median")
p_sd     <- plot_buffer_variation(drone_cv, "sd_height", "SD of Height")
p_cv     <- plot_buffer_variation(drone_cv, "cv_height", "CV of Height")
p_max    <- plot_buffer_variation(drone_cv, "max_height", "Max Height")
p_min    <- plot_buffer_variation(drone_cv, "min_height", "Min Height")

# Combine using ggarrange
ggarrange(p_mean, p_median, p_sd, p_cv, p_max, p_min, 
          ncol = 2, nrow = 3, 
          labels = "AUTO")

# convert to wide format as i have 4 different buffer sizes
drone_cv_wide <- drone_cv %>%
  dplyr::filter(buffer_size == 20) %>%  # filter to single buffer, representing our sample (increase to 22-25 m!!)
  pivot_wider(
    id_cols = c(cluster, drone),
    names_from = buffer_size,
    values_from = c(mean_height, sd_height, cv_height, max_height, min_height),
    names_glue = "{.value}_{buffer_size}"
  )


# merge with structural data ---------------------------------------
# Add prefix to drone data
drone_cv_wide_renamed <- drone_cv_wide %>%
  rename_with(~ paste0("drone_", .), -c(cluster, drone))

# Add prefix to pre-disturbance data
pre_dist_renamed <- pre_dist_history_2023 %>%
  rename_with(~ paste0("pre_", .), -cluster)

# Now join
df_fin <- df_field_diversity %>%
  right_join(pre_dist_renamed, by = "cluster") %>% 
  right_join(drone_cv_wide_renamed, by = "cluster")# %>%
  

# any correlation between variables?
# split between field based vs drone based. pre-dciturbance are drivers. 
# Step 1: Select numeric columns only
# Step 1: Select numeric variables
df_numeric <- df_fin %>%
  select(where(is.numeric)) %>%
  drop_na()  # drop rows with NAs for clean correlation and VIF

# Step 2: Correlation matrix
cor_mat <- cor(df_numeric)

# Step 3: Visualize with corrplot
corrplot(cor_mat,
         method = "circle",
         type = "lower",
         tl.col = "black",
         tl.cex = 0.8,
         number.cex = 0.7,
         title = "Correlation Matrix (Field & Drone Metrics)",
         mar = c(0,0,1,0),
         addCoef.col = "black")


# is tehre a correlation betwen field and drone derived vertical strcuture?
df_model %>% 
  ggplot(aes(x = shannon_height,
             y = drone_cv_height_20)) + 
  geom_point() + 
  geom_smooth(method = "loess")

# Fit the linear model --------------------------------
df_model <- df_fin %>%
  dplyr::filter(drone_cv_height_20 < 10)  # example threshold


model <- lm(shannon_height ~ drone_cv_height_20, data = df_model)

# Summary output
summary(model)


# compare field vs drone metrics: ------------------------------------------------

# Define columns
field_metrics <- df_fin %>%
  select(stem_density, basal_area_ha_m2, shannon_species, shannon_height)

drone_metrics <- df_fin %>%
  select(drone_mean_height_20, drone_sd_height_20, drone_cv_height_20,
         drone_max_height_20, drone_min_height_20)

# Combine and track origin
combined_df <- bind_cols(field_metrics, drone_metrics)

# Label columns for plot clarity
colnames(combined_df) <- c(
  "StemDensity", "BA_ha", "Shannon_Species", "Shannon_Height",
  "Drone_Mean", "Drone_SD", "Drone_CV", "Drone_Max", "Drone_Min"
)

# Custom ggpairs plot: compare only Field vs Drone (upper part blank)
GGally::ggpairs(
  combined_df,
  columns = 1:9,
  lower = list(
    continuous = wrap("smooth", alpha = 0.6, size = 0.5)
  ),
  upper = list(continuous = "blank"),
  diag = list(continuous = "blankDiag"),
  mapping = aes(color = NULL)
) +
  ggtitle("Field vs Drone Metrics: Pairwise Comparison")

