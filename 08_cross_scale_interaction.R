
# Cross-scale interaction: field & drone variation in vertical structures -------------------

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


## Sumarize field observation per cluster  --------------------------------------

# guestimate dbh and ba per individual based on height distribution 
dat23_subplot_recode <- dat23_subplot %>% 
  mutate(
    # Create a numeric height estimate (keep original height class string)
    hgt_est = case_when(
      vegtype == "mature" & dbh == "10–20cm" ~ 10.0,
      vegtype == "mature" & dbh == "20–40cm" ~ 20.0,
      vegtype == "mature" & dbh == "40–60cm" ~ 30.0,
      hgt == "0.2–0.4"                       ~ 0.3,
      hgt == "0.4–0.6"                       ~ 0.5,
      hgt == "0.6–0.8"                       ~ 0.7,
      hgt == "0.8–1.0"                       ~ 0.9,
      hgt == "1.0–1.3"                       ~ 1.2,
      hgt == "1.3–2.0"                       ~ 1.7,
      hgt == "2–4"                           ~ 3.0,
      hgt == ">4"                            ~ 5.0,
      TRUE                                   ~ NA_real_
    ),
    
    # Estimate DBH (already numeric)
    dbh_est = case_when(
      vegtype == "mature" & dbh == "10–20cm" ~ 15.0,
      vegtype == "mature" & dbh == "20–40cm" ~ 30.0,
      vegtype == "mature" & dbh == "40–60cm" ~ 50.0,
      hgt == "0.2–0.4"                       ~ 0.3,
      hgt == "0.4–0.6"                       ~ 0.5,
      hgt == "0.6–0.8"                       ~ 0.7,
      hgt == "0.8–1.0"                       ~ 0.9,
      hgt == "1.0–1.3"                       ~ 1.2,
      hgt == "1.3–2.0"                       ~ 1.7,
      hgt == "2–4"                           ~ 3.0,
      hgt == ">4"                            ~ 5.0,
      TRUE                                   ~ NA_real_
    )
  ) %>% 
  # Calculate basal area
  mutate(
    basal_area_cm2 = pi * (dbh_est / 2)^2,
    ba_total_cm2   = basal_area_cm2 * n,
    ba_total_m2    = ba_total_cm2 / 10000,
    ba_ha_m2       = ba_total_m2 * scaling_factor
  )


### get stem density and shannon for field data --------------
df_cluster <- dat23_subplot_recode %>% 
  group_by(cluster, manag_intensity, salvage_intensity, protection_intensity) %>% 
  summarise(stem_density = sum(stem_density, na.rm = T),
            basal_area_ha_m2 = sum(ba_ha_m2, na.rm = TRUE))#,
           # mean_hgt  = mean)

# weighted average height --------------
#summarise(weighted_hgt = sum(hgt_est * basal_area_cm2) / sum(basal_area_cm2))


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





# Analysis ----------------------------------------------

## Calculate alpha, beta and gamm diversity, Jackard index --------------------------------------

library(vegan)

# Clean data: height class must be non-missing
dat_clean2 <- dat23_subplot_recode %>%
  filter(hgt_est != "" & !is.na(hgt_est)) %>% 
  mutate(hgt_est = as.factor(hgt_est))

# Create a presence-absence matrix: rows = plots, cols = height classes
height_pa_matrix <- dat_clean2 %>%
  mutate(present = 1) %>%
  distinct(ID, cluster, hgt_est, .keep_all = TRUE) %>%
  pivot_wider(
    id_cols = c(ID, cluster),
    names_from = hgt_est,
    values_from = present,
    values_fill = 0
  )

# Split into metadata and numeric presence-absence 
metadata <- height_pa_matrix %>%
  select(ID, cluster)

height_data <- height_pa_matrix %>%
  tibble::column_to_rownames("ID") %>%
  select(-cluster)

# Check that height_data is now numeric
stopifnot(all(sapply(height_data, is.numeric)))

### ---- Alpha diversity 
alpha_per_plot <- rowSums(height_data)
alpha_mean <- mean(alpha_per_plot)

# ---- Gamma diversity (total per cluster) 
gamma_per_cluster <- height_pa_matrix %>%
  pivot_longer(cols = -c(ID, cluster), names_to = "hgt_est", values_to = "present") %>%
  group_by(cluster, hgt_est) %>%
  summarise(present = max(present), .groups = "drop") %>%
  group_by(cluster) %>%
  summarise(gamma_diversity = sum(present), .groups = "drop")

### ---- Beta diversity 
beta_per_cluster <- metadata %>%
  mutate(alpha = alpha_per_plot) %>%
  left_join(gamma_per_cluster, by = "cluster") %>%
  group_by(cluster) %>%
  summarise(
    alpha_mean = mean(alpha),
    gamma = first(gamma_diversity),
    beta = gamma / alpha_mean
  )

###  Optional: Jaccard dissimilarity among plots
jaccard_dist <- vegdist(height_data, method = "jaccard", binary = TRUE)

# Convert to dataframe for summary if needed
jaccard_df <- as.data.frame(as.matrix(jaccard_dist))
jaccard_df$ID <- rownames(jaccard_df)


# Step 1: Convert distance matrix to long format
jaccard_df <- as.data.frame(as.matrix(jaccard_dist)) %>%
  tibble::rownames_to_column("plot1") %>%
  tidyr::pivot_longer(-plot1, names_to = "plot2", values_to = "dissimilarity")

# Remove self-comparisons
jaccard_df <- jaccard_df %>%
  filter(plot1 != plot2)

# Step 2: Add cluster info to both plots
cluster_lookup <- metadata %>% distinct(ID, cluster)

jaccard_df <- jaccard_df %>%
  left_join(cluster_lookup, by = c("plot1" = "ID")) %>%
  rename(cluster1 = cluster) %>%
  left_join(cluster_lookup, by = c("plot2" = "ID")) %>%
  rename(cluster2 = cluster)

# Step 3: Keep only within-cluster comparisons
jaccard_within_cluster <- jaccard_df %>%
  filter(cluster1 == cluster2) %>%
  group_by(cluster = cluster1) %>%
  summarise(
    mean_jaccard = mean(dissimilarity, na.rm = TRUE),
    n_comparisons = n()
  )

# View result
jaccard_within_cluster


#### merge dissimilarity between clusters --------------------------------------------------
df_similarity <- jaccard_within_cluster %>% 
  left_join(beta_per_cluster)


df_similarity %>% 
  ggplot(aes(x = beta,
             y = mean_jaccard)) +
  geom_point() + 
  geom_smooth()


# Cross scale variation ----------------------------------------------
##  variation of height classes across scales: plot, cluster, landscape -----------------------

# 

# Ensure height is treated as a factor, ignore missing/empty strings
dat_clean <- dat23_subplot_recode %>%
  filter(hgt_est != "" & !is.na(hgt_est))

# 1. Number of unique height classes per plot
height_div_plot <- dat_clean %>%
  group_by(ID, cluster) %>%
  summarise(n_height_classes = n_distinct(hgt_est)) %>%
  ungroup()

# 2. Number of unique height classes per cluster
height_div_cluster <- dat_clean %>%
  group_by(cluster) %>%
  summarise(n_height_classes = n_distinct(hgt_est)) %>%
  ungroup()

# 3. Number of unique height classes across landscape
n_height_classes_landscape <- dat_clean %>%
  summarise(n_height_classes = n_distinct(hgt_est)) %>%
  pull(n_height_classes)



# Compute medians and IQRs
summary_stats <- tibble(
  spat_scale = factor(c("Plot", "Cluster", "Landscape"),
                    level = c("Plot", "Cluster", "Landscape")),
  median = c(
    median(height_div_plot$n_height_classes),
    median(height_div_cluster$n_height_classes),
    n_height_classes_landscape
  ),
  IQR_low = c(
    quantile(height_div_plot$n_height_classes, 0.25),
    quantile(height_div_cluster$n_height_classes, 0.25),
    n_height_classes_landscape
  ),
  IQR_high = c(
    quantile(height_div_plot$n_height_classes, 0.75),
    quantile(height_div_cluster$n_height_classes, 0.75),
    n_height_classes_landscape
  )
)


ggplot(summary_stats, aes(y = spat_scale, x = median)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = IQR_low, xmax = IQR_high), 
                 height = 0) +
  labs(
    x = "Number of Height Classes",
    y = "Spatial Scale",
    title = "Height Class Diversity Across Spatial Scales"
  ) +
  theme_classic2()

## Plots vs cluster --------------------------------------------------------------
## HOw value in one plot affect clusre level mean? sd?
# split table in towo group: 
# select one random subplot
# subset remaining subplots - calculate means and sd

# One random row per cluster
height_div_subset <- height_div_plot %>%
  group_by(cluster) %>%
  slice_sample(n = 1) %>%
  ungroup() %>% 
  rename(plt_n_hgt = n_height_classes)

# Remaining rows
height_div_remaining <- anti_join(height_div_plot, height_div_subset, by = c("ID", "cluster"))

height_div_remaining_summary <- height_div_remaining %>% 
  group_by(cluster) %>% 
  summarize(cl_mean_hgt_n = mean(n_height_classes, na.rm = T),
            cl_median_hgt_n = median(n_height_classes, na.rm = T),
            cl_sd_hgt_n = sd(n_height_classes, na.rm = T),
            cl_cv_hgt_n = cl_sd_hgt_n/cl_mean_hgt_n) %>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  left_join(height_div_subset) #%>% # replace NA by 0s
  
# how does the plot value predict the cluster value? 
# predcit average value per plot
height_div_remaining_summary %>% 
  ggplot(aes(x = plt_n_hgt,
             y = cl_mean_hgt_n )) + 
  geom_point() +
  geom_smooth(method = 'loess') +
  theme_classic()

### predict coeffcient of variation -
height_div_remaining_summary %>% 
  ggplot(aes(x = plt_n_hgt,
             y = cl_cv_hgt_n )) + 
  geom_point() +
  geom_smooth(method = 'loess') +
  theme_classic()



### bootstrap analysis of how the random selection of one point represent the values per ----
# cluster

# Set number of bootstrap iterations
n_boot <- 1000

# Store results
bootstrap_results <- map_dfr(1:n_boot, function(i) {
  
  # 1. Random plot per cluster
  height_div_subset <- height_div_plot %>%
    group_by(cluster) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    rename(plot_ID = ID,
           plt_n_hgt = n_height_classes)
  
  # 2. Remaining plots
  height_div_remaining <- anti_join(
    height_div_plot,
    height_div_subset,
    by = c("ID" = "plot_ID", "cluster", "n_height_classes" = "plt_n_hgt")
  )
  
  # 3. Mean height richness of remaining plots
  cluster_remaining_summary <- height_div_remaining %>%
    group_by(cluster) %>%
    summarize(cl_mean_hgt_n = mean(n_height_classes, na.rm = TRUE)) %>%
    ungroup()
  
  # 4. Merge back to random plots
  paired <- height_div_subset %>%
    left_join(cluster_remaining_summary, by = "cluster") %>%
    mutate(iteration = i)
  
  return(paired)
})


### calculate average LOESS fit -----------------
# Define prediction grid (plot-level richness: 1, 2, 3)
pred_grid <- tibble(plt_n_hgt = seq(1, 3, by = 0.1))

# Apply LOESS fit to each iteration and predict over the grid
loess_fits <- bootstrap_results %>%
  filter(!is.na(cl_mean_hgt_n)) %>%
  group_by(iteration) %>%
  group_split() %>%
  map_dfr(function(df) {
    loess_model <- tryCatch(loess(cl_mean_hgt_n ~ plt_n_hgt, data = df, span = 0.75), error = function(e) NULL)
    if (!is.null(loess_model)) {
      pred <- predict(loess_model, newdata = pred_grid)
      tibble(plt_n_hgt = pred_grid$plt_n_hgt, cl_mean_pred = pred, iteration = unique(df$iteration))
    } else {
      NULL
    }
  })

# Aggregate: mean and 95% CI of LOESS fits
loess_summary <- loess_fits %>%
  group_by(plt_n_hgt) %>%
  summarize(
    mean_fit = mean(cl_mean_pred, na.rm = TRUE),
    lower_CI = quantile(cl_mean_pred, 0.025, na.rm = TRUE),
    upper_CI = quantile(cl_mean_pred, 0.975, na.rm = TRUE)
  )


ggplot(loess_summary, aes(x = plt_n_hgt, y = mean_fit)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray") +
  labs(
    x = "Height Class Richness in Random Plot",
    y = "Mean Height Class Richness in Remaining Plots",
    title = "Average LOESS Fit Across Bootstraps"
  ) +
  theme_classic2()



### calculate pseudo R2 ----------------------- 
# Function to compute pseudo-R² for each LOESS fit
compute_loess_r2 <- function(df) {
  loess_model <- tryCatch(loess(cl_mean_hgt_n ~ plt_n_hgt, data = df, span = 0.75), error = function(e) NULL)
  if (!is.null(loess_model)) {
    preds <- predict(loess_model, newdata = df)
    y <- df$cl_mean_hgt_n
    rss <- sum((y - preds)^2, na.rm = TRUE)
    tss <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- 1 - (rss / tss)
    return(tibble(iteration = unique(df$iteration), r2 = r2))
  } else {
    return(tibble(iteration = unique(df$iteration), r2 = NA_real_))
  }
}

# Apply to each bootstrap iteration
r2_results <- bootstrap_results %>%
  filter(!is.na(cl_mean_hgt_n)) %>%
  group_by(iteration) %>%
  group_split() %>%
  map_dfr(compute_loess_r2)

# Summarize R² distribution
r2_summary <- r2_results %>%
  summarize(
    mean_r2 = mean(r2, na.rm = TRUE),
    median_r2 = median(r2, na.rm = TRUE),
    r2_2.5 = quantile(r2, 0.025, na.rm = TRUE),
    r2_97.5 = quantile(r2, 0.975, na.rm = TRUE)
  )


ggplot(r2_results, aes(x = r2)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  labs(
    x = expression(R^2 ~ " of LOESS fits across bootstraps"),
    y = "Frequency",
    title = "Distribution of Pseudo-" ~ R^2 ~ " Values"
  ) +
  theme_minimal()





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

## merge field data with drone estimation ------------------------------

### Sensitivity analysis ----------

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

# convert to wide format as i have 4 different buffer sizes - need to do this for cluster level as well!
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
 # right_join(drone_cv_wide_renamed, by = "cluster") %>%
  right_join(df_similarity, by = "cluster") %>% 
  right_join(height_div_remaining_summary)
  

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


#s does vertical diversity between plots predict stem density? -------------

df_fin %>% 
  ggplot(aes(x = shannon_height     ,
             y = stem_density )) + 
  geom_point() + 
  geom_smooth(method = "lm")




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



# how does cluster variability affects foresrt resilience? ---------------
head(df_fin)


cor(df_fin$cl_cv_hgt_n , df_fin$stem_density)

hist(df_fin$stem_density)

library(mgcv)
# library(tweedie)
# library(statmod)

# Start with null model
m_null <- gam(stem_density ~ 1, data = df_fin, family = tw(link = "log"))

# Try candidate models with individual smooth terms
m1 <- gam(stem_density ~ s(pre_mean_cover_dens, k = 3), data = df_fin, family = tw(link = "log"))
m2 <- gam(stem_density ~ s(pre_share_coniferous, k = 3), data = df_fin, family = tw(link = "log"))
m3 <- gam(stem_density ~ s(pre_shannon, k = 3), data = df_fin, family = tw(link = "log"))
m4 <- gam(stem_density ~ s(mean_jaccard, k = 3), data = df_fin, family = tw(link = "log"))
m5 <- gam(stem_density ~ s(alpha_mean, k = 3), data = df_fin, family = tw(link = "log"))
m6 <- gam(stem_density ~ s(beta, k = 3), data = df_fin, family = tw(link = "log"))
m7 <- gam(stem_density ~ s(manag_intensity, k = 3), data = df_fin, family = tw(link = "log"))
m8 <- gam(stem_density ~ s(salvage_intensity, k = 3), data = df_fin, family = tw(link = "log"))
m9 <- gam(stem_density ~ s(protection_intensity, k = 3), data = df_fin, family = tw(link = "log"))

# Compare AICs
AIC(m_null, m1, m2, m3, m4, m5, m6, m7, m8, m9)

df_fin %>%
  ungroup(.) %>% 
  dplyr::select(beta, mean_jaccard, salvage_intensity, protection_intensity) %>%
  cor(use = "complete.obs")



best_model <- gam(
  stem_density ~ 
  #  s(alpha_mean, k = 3) +
    s(cl_cv_hgt_n , k = 3) +
    #s(mean_jaccard, k = 3) +
    #s(manag_intensity, k = 3) ,#+
    s(salvage_intensity, k = 3) +
    s(protection_intensity, k = 3) +
    s(pre_share_coniferous, k = 3) +
    s(pre_mean_cover_dens, k = 3),
    
  data = df_fin,
  family = tw(link = "log")
)
summary(best_model)
plot(best_model,page = 1)

library(gratia)
appraise(best_model)
