# test deltaes between 2023 and 2025

# needs to be complated!!!


# get estimation of growth  and mortality ---------------------------

intensity_vars <- c("planting_intensity","anti_browsing_intensity",
                    "grndwrk_intensity","clear_intensity")

response_vars  <- c("mean_hgt",
                    "dens_ha",
                    "cv_hgt",
                    "effective_numbers",
                    "sp_richness",
                    "IV",
                    "spruce_share",
                    "CWM_shade",
                    "CWM_drought"
)

# identify plots with both years present (NOW based on the clean table)
plots_both <- plot_df %>%
  distinct(plot_id, year) %>%
  count(plot_id, name = "n_years") %>%
  filter(n_years == 2) %>%
  pull(plot_id)


# get only one management per plot; select the higher one: collapse to one row per plot_id × year (resolve duplicates)
plot_df_max_intensities <- plot_df %>%
  group_by(plot_id) %>%
  summarise(
    across(all_of(intensity_vars), safe_max),               # higher management intensity
    time_snc23 = min(time_snc_full_disturbance, na.rm = TRUE), # lower time since disturbnace - in 2023
    .groups = "drop"
  )


# convert to wide widen + compute deltas
plot_wide_change <- 
  plot_df %>% 
  filter(plot_id %in% plots_both) %>% # filter only plots that have both years of recording
  mutate(
    year_short = case_when(
      year == 2023 ~ "23",
      year == 2025 ~ "25",
      TRUE ~ NA_character_
    )
  ) %>% 
  
  select(plot_id, year_short, all_of(response_vars)) %>% 
  mutate(  # change numeric variables into 0
    across(where(is.numeric), ~replace_na(.x, 0))
  ) %>% 
  pivot_wider(
    names_from  = year_short,
    values_from = c(all_of(response_vars)),
    names_sep   = ""
  ) %>%
  mutate(
    dt_mean_hgt          = mean_hgt25 - mean_hgt23,
    dt_dens_ha           = dens_ha25  - dens_ha23,
    dt_cv_hgt            = cv_hgt25   - cv_hgt23,
    dt_effective_numbers = effective_numbers25 - effective_numbers23,
    dt_sp_richness       = sp_richness25 - sp_richness23,
    dt_spruce_share      = spruce_share25 - spruce_share23,
    dt_CWM_shade         = CWM_shade25 - CWM_shade23,
    dt_CWM_drought       = CWM_drought25 - CWM_drought23,
    dt_IV                = IV25 - IV23
  ) %>% 
  left_join(plot_df_max_intensities)


# get all deltas
delta_vars <- names(plot_wide_change)[stringr::str_detect(names(plot_wide_change), "^dt_")]


# long table of deltas - quicjk ggplot f delta's distribution
delta_long <- 
  plot_wide_change %>%
  select(plot_id, all_of(delta_vars)) %>%
  pivot_longer(cols = all_of(delta_vars),
               names_to = "metric",
               values_to = "delta") %>% 
  left_join(plot_df_max_intensities)


delta_long %>%
  ggplot(aes(x = factor(time_snc23), y = delta)) +
  #geom_boxplot(outlier.alpha = 0.3, width = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_jitter(alpha = 0.3, size = 0.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    # x = NULL,
    y = "Change (2025 – 2023)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 9)
  )


# by management
delta_long %>%
  ggplot(aes(x = factor(clear_intensity), y = delta)) +
  #geom_boxplot(outlier.alpha = 0.3, width = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_jitter(alpha = 0.3, size = 0.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    #x = NULL,
    y = "Change (2025 – 2023)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 9)
  )


delta_long %>%
  ggplot(aes(x = factor(planting_intensity), y = delta)) +
  #geom_boxplot(outlier.alpha = 0.3, width = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_jitter(alpha = 0.3, size = 0.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    #x = NULL,
    y = "Change (2025 – 2023)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 9)
  )


delta_long %>%
  ggplot(aes(x = factor(anti_browsing_intensity), y = delta)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.4) +
  #geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_jitter(alpha = 0.3, size = 0.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    #x = NULL,
    y = "Change (2025 – 2023)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 9)
  )



delta_long %>%
  count(time_snc23)

# i have support only in 1:4 years
delta_long_main <- delta_long %>%
  filter(time_snc23 >= 1, time_snc23 <= 4)

# keep 0 year separated
delta_long_zero <- delta_long %>%
  filter(time_snc23 == 0)

delta_long_main %>%
  ggplot(aes(x = "", y = delta)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(alpha = 0.3, size = 1) +
  #  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3)


# get histogram for distributions

# compute summary statistics per metric
delta_stats <- delta_long_main %>%
  group_by(metric) %>%
  summarise(
    mean_delta   = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    .groups = "drop"
  )


delta_long_main %>%
  ggplot(aes(x = delta)) +
  geom_density(fill = "grey70", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ metric, scales = "free_x", ncol = 3) +
  theme_classic()

delta_long %>%
  ggplot(aes(x = delta)) +
  geom_histogram() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(data = delta_stats,
             aes(xintercept = mean_delta),
             linetype = "solid",
             color = 'red',
             linewidth = 0.6) +
  geom_vline(data = delta_stats,
             aes(xintercept = median_delta),
             linetype = "dotted",
             color = 'grey',
             linewidth = 0.6) +
  facet_wrap(~ metric, scales = "free_x", ncol = 3) +
  labs(
    
    y = "Change (2025 − 2023)"
  ) 

delta_long_main %>%
  filter(metric %in% c("dt_mean_hgt",
                       "dt_sp_richness")) %>%
  ggplot(aes(x = planting_intensity, y = delta)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  facet_wrap(~ metric, scales = "free_y")


delta_long_main %>%
  ggplot(aes(x = time_snc23, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    x = "Time since disturbance (years)",
    y = "Change (2025 − 2023)"
  ) +
  theme_classic(base_size = 10)


delta_long_main %>%
  ggplot(aes(x = planting_intensity, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    y = "Change (2025 − 2023)"
  ) +
  theme_classic(base_size = 10)

delta_long_main %>%
  ggplot(aes(x = anti_browsing_intensity, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    y = "Change (2025 − 2023)"
  ) +
  theme_classic(base_size = 10)

delta_long_main %>%
  ggplot(aes(x = grndwrk_intensity, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "glm", se = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(
    y = "Change (2025 − 2023)"
  ) +
  theme_classic(base_size = 10)



# Get models

m_hgt <- gam(
  delta ~ 
    s(time_snc23, k = 4) +
    planting_intensity*anti_browsing_intensity +
    grndwrk_intensity,
  data = delta_long %>% filter(metric == "dt_mean_hgt")
)

appraise(m_hgt)




