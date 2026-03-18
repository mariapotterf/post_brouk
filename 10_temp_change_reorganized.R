# ============================================================
# Forest regeneration: temporal change analysis 2023-2025
# ============================================================
# Sections:
#   0. Setup
#   1. Data ingestion & QC
#   2. Management data cleaning
#   3. Derived variables (IV, CWM, spruce, subplot/plot tables)
#   4. Species composition & beta diversity
#   5. Models (GAMs)
#   6. Figures
#   7. Tables & export
# ============================================================


# 0. Setup -------------------------------------------------------

gc()
set.seed(1234)

## Libraries
library(terra);      library(sf)
library(data.table); library(dplyr);   library(tidyr)
library(stringr);    library(purrr)
library(ggplot2);    library(ggpubr);  library(ggridges); library(gghalves)
library(scales);     library(forcats); library(RColorBrewer)
library(corrplot);   library(GGally)
library(vegan)
library(mgcv);       library(ggeffects); library(gratia)
library(broom);      library(emmeans)
library(tibble)
library(spdep)
library(sjPlot)
library(patchwork)

## Project variables (palettes, labels, species lists, species_class, earlyspecs_laura)
source('my_variables.R')

## ggplot theme
theme_set(
  theme_classic2(base_size = 10) +
    theme(axis.title = element_text(size = 10),
          axis.text  = element_text(size = 10))
)

## Constants
area_subplot_m2 <- 4
area_plot_m2    <- 5 * 4

## Shared axis label
x_lab_time_snc_full_dist <- "Time since stand\nreplacing disturbance (years)"

## Shared factor level labels
# intensity_levels <- c("0\u201319", "20\u201339", "40\u201359", "60\u201379", "80\u2013100")
# low_classes      <- c("0\u201319", "20\u201339")


intensity_levels <-  c("0–19", "20–39", "40–59", "60–79", "80–100")
low_classes <- c("0–19", "20–39")
applied_mng_intens_order <- c("0–19", "20–39","40–59", "60–79", "80–100") 



activity_intens_labels <- c(
  "clear_intensity"         = "Salvage\nlogging",
  "grndwrk_intensity"       = "Soil\npreparation",
  "planting_intensity"      = "Planting",
  "anti_browsing_intensity" = "Browsing\nprotection"
)

## Helper functions

safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

wfun <- function(x) ifelse(is.na(x) | x < 0, 0, x)

calc_iv_core <- function(data, size_var, ...) {
  data %>%
    group_by(..., species) %>%
    summarise(
      n_sp    = sum(n, na.rm = TRUE),
      size_sp = sum({{ size_var }}, na.rm = TRUE),
      .groups = "drop_last"
    ) %>%
    mutate(
      RA       = n_sp / sum(n_sp, na.rm = TRUE),
      size_tot = sum(size_sp, na.rm = TRUE),
      RS       = dplyr::if_else(size_tot > 0, size_sp / size_tot, NA_real_),
      IV       = (RA + RS) / 2
    ) %>%
    select(-size_tot) %>%
    ungroup()
}

calc_iv_subplot <- function(data, size_var) calc_iv_core(data, {{ size_var }}, plot, year, subplot)
calc_iv_plot    <- function(data, size_var) calc_iv_core(data, {{ size_var }}, plot, year)

## Plotting helpers for model predictions (used in §5-6)
pp <- function(model, terms, xlab = NULL, ylab = NULL,
               annot = NULL, scale_y = 1,
               annot_x = 3.5, annot_y = NULL) {
  pr <- ggpredict(model, terms = terms, exclude = "s(plot_id)") %>%
    as.data.frame() %>%
    mutate(across(c(predicted, conf.low, conf.high), ~ . * scale_y))

  p <- ggplot(pr, aes(x = x, y = predicted, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.8) +
    labs(x = xlab, y = ylab, colour = "Scale") +
    scale_color_manual(values = c("subplot" = "grey80", "plot" = "grey30")) +
    scale_fill_manual( values = c("subplot" = "grey80", "plot" = "grey30"), guide = "none") +
    theme_classic() +
    theme(axis.title = element_text(size = 7))

  if (!is.null(annot))
    p <- p + annotate("text", x = annot_x, y = annot_y,
                      label = annot, size = 2.5, hjust = 0.5)
  p
}

pp_inset_model <- function(model, scale_y = 1, p_lab = NULL,
                           annot_x = 1.5, annot_y = NULL) {
  pr <- ggpredict(model, terms = "level") %>%
    as.data.frame() %>%
    mutate(across(c(predicted, conf.low, conf.high), ~ . * scale_y))

  ggplot(pr, aes(x = x, y = predicted, colour = x)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .08) +
    geom_point(size = 2.4) +
    annotate("text", x = annot_x, y = annot_y, label = p_lab, size = 2.5) +
    scale_color_manual(values = c("subplot" = "grey80", "plot" = "grey30"), guide = "none") +
    theme_classic(base_size = 8) +
    theme(axis.title      = element_blank(),
          axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          axis.line.y     = element_blank(),
          axis.text.x     = element_text(size = 7),
          legend.position = "none",
          plot.margin     = margin(5, 5, 5, 0))
}


# 1. Data ingestion & QC -----------------------------------------

dat_overlap <- fread('outData/full_table_23_25.csv')

dat_master_subplot <- dat_overlap %>%
  dplyr::select(plot, subplot, year, status) %>%
  distinct()

table(dat_master_subplot$year)
table(dat_master_subplot$status)

n_plots_total    <- length(unique(dat_master_subplot$plot))
n_subplots_total <- length(unique(dat_master_subplot$subplot))

cat("Plots:", n_plots_total, "| Subplots:", n_subplots_total, "\n")

## Quick counts
total_trees_per_year <- dat_overlap %>%
  group_by(year) %>%
  summarise(total_trees = sum(n, na.rm = TRUE), .groups = "drop")

tree_summary <- dat_overlap %>%
  group_by(year, seral_stage) %>%
  summarise(n_trees_recovery = sum(n, na.rm = TRUE), .groups = "drop") %>%
  left_join(total_trees_per_year) %>%
  mutate(share = n_trees_recovery / total_trees * 100)


# 2. Management data cleaning ------------------------------------

management_vars <- c("clear", "grndwrk", "logging_trail", "planting", "anti_browsing")

management_vars_rmv <- c(
  "salvage_sum", "protection_sum", "clear_sum", "grndwrk_sum",
  "logging_trail_sum", "planting_sum", "anti_browsing_sum",
  "management_sum", "clear_intensity", "grndwrk_intensity",
  "logging_trail_intensity", "planting_intensity",
  "anti_browsing_intensity", "salvage_intensity",
  "protection_intensity", "management_intensity"
)

dat_overlap_mng <- dat_overlap %>%
  dplyr::select(-all_of(management_vars_rmv))

## Identify plots surveyed in both years
plots_with_both <- dat_overlap_mng %>%
  dplyr::filter(status == "both", year %in% c(2023, 2025)) %>%
  distinct(plot, year) %>%
  count(plot) %>%
  filter(n > 1) %>%
  pull(plot)

## 2023 management (reference)
mng_both23 <- dat_overlap_mng %>%
  filter(plot %in% plots_with_both, status == "both", year == 2023) %>%
  dplyr::select(plot, subplot, all_of(management_vars), year) %>%
  distinct() %>%
  group_by(plot, year) %>%
  arrange(plot, subplot) %>%
  mutate(row_id = row_number())

## 2025 management — resolve duplicates by taking max presence
mng_both25 <- dat_overlap_mng %>%
  filter(plot %in% plots_with_both, status == "both", year == 2025) %>%
  dplyr::select(plot, subplot, all_of(management_vars), year) %>%
  group_by(plot, subplot, year) %>%
  summarise(across(all_of(management_vars), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>%
  group_by(plot) %>%
  arrange(plot, subplot) %>%
  mutate(row_id = row_number())

## Combine by row order, take max per management variable
mng_both25_upd <- bind_rows(mng_both23, mng_both25) %>%
  group_by(plot, row_id) %>%
  summarise(
    subplot = subplot[year == 2025][1],
    across(all_of(management_vars), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  select(-row_id) %>%
  mutate(year = 2025)

## Rebuild full dataset with reconciled 2025 management
dat_overlap_mng_25 <- dat_overlap_mng %>%
  filter(status == "both" & year == 2025) %>%
  dplyr::select(-all_of(management_vars)) %>%
  left_join(mng_both25_upd) %>%
  dplyr::select(all_of(colnames(dat_overlap_mng)))

dat_overlap_mng_rest <- dat_overlap_mng %>%
  dplyr::filter(!(status == "both" & year == 2025))

dat_overlap_mng_upd <- bind_rows(dat_overlap_mng_rest, dat_overlap_mng_25)

## Management intensity per plot × year (proportion of subplots with each activity)
mng_subplot <- dat_overlap_mng_upd %>%
  dplyr::select(plot, subplot, year, n_subplots, all_of(management_vars)) %>%
  distinct()

management_intensity <- mng_subplot %>%
  group_by(plot, year) %>%
  summarise(
    n_subplots              = first(n_subplots),
    clear_intensity         = sum(clear         == 1, na.rm = TRUE) / n_subplots,
    grndwrk_intensity       = sum(grndwrk       == 1, na.rm = TRUE) / n_subplots,
    logging_trail_intensity = sum(logging_trail == 1, na.rm = TRUE) / n_subplots,
    planting_intensity      = sum(planting      == 1, na.rm = TRUE) / n_subplots,
    anti_browsing_intensity = sum(anti_browsing == 1, na.rm = TRUE) / n_subplots,
    .groups = "drop"
  )

dat_overlap_mng_upd2 <- dat_overlap_mng_upd %>%
  left_join(management_intensity, by = c("plot", "year", "n_subplots"))

## Subplot and plot management tables
df_mng_sub <- dat_overlap_mng_upd2 %>%
  distinct(plot, subplot, year,
           clear, grndwrk, logging_trail, planting, anti_browsing,
           clear_intensity, grndwrk_intensity, logging_trail_intensity,
           planting_intensity, anti_browsing_intensity)

df_mng_plot <- df_mng_sub %>%
  select(plot, year,
         clear_intensity, grndwrk_intensity, logging_trail_intensity,
         planting_intensity, anti_browsing_intensity) %>%
  distinct()

## Management intensity summary (for figures)
df_master_mng_intensity <- dat_overlap_mng_upd2 %>%
  distinct(plot, year,
           clear_intensity, grndwrk_intensity, logging_trail_intensity,
           planting_intensity, anti_browsing_intensity)

mng_intensity_props <- df_master_mng_intensity %>%
  pivot_longer(cols = c(clear_intensity, grndwrk_intensity, logging_trail_intensity,
                        planting_intensity, anti_browsing_intensity),
               names_to = "activity", values_to = "applied") %>%
  group_by(activity, applied, year) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(activity) %>%
  mutate(
    proportion = n / sum(n) * 100,
    intensity_class = cut(
      applied,
      breaks = c(-Inf, 0.19, 0.39, 0.59, 0.79, Inf),
      labels = intensity_levels,
      right  = TRUE
    ),
    intensity_binary = case_when(intensity_class == intensity_levels[1] ~ "no",
                                 TRUE ~ "yes")
  )

## Binary summary (collapsed across years)
mng_sub_conv <- mng_intensity_props %>%
  mutate(intensity_binary = factor(intensity_binary, levels = c("no", "yes")),
         intensity_class  = factor(intensity_class, levels = intensity_levels)) %>%
  group_by(activity, intensity_binary) %>%
  summarise(proportion = sum(proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(proportion_plot = if_else(intensity_binary == "no", -proportion, proportion))

## Intensity summary (shifted for diverging bar)
mng_shifted <- mng_intensity_props %>%
  mutate(intensity_binary = factor(intensity_binary, levels = c("no", "yes")),
         intensity_class  = factor(intensity_class, levels = intensity_levels)) %>%
  group_by(activity, intensity_class) %>%
  summarise(proportion = sum(proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    proportion_shifted   = if_else(intensity_class %in% low_classes, -proportion, proportion),
    intensity_class_plot = fct_relevel(intensity_class,
                                       c(rev(intensity_levels[3:5]),
                                         intensity_levels[1:2])),
    activity = factor(activity, levels = rev(c("clear_intensity", "grndwrk_intensity",
                                               "planting_intensity", "anti_browsing_intensity",
                                               "logging_trail_intensity")))
  )

## Intensity class summary table
intensity_class_summary_counts <- mng_intensity_props %>%
  group_by(activity, intensity_class) %>%
  summarise(n_plots = sum(n), .groups = "drop")

intensity_class_summary_percent <- intensity_class_summary_counts %>%
  group_by(activity) %>%
  mutate(share = round(100 * n_plots / sum(n_plots), 1)) %>%
  unite(col = "value", n_plots, share, sep = " (") %>%
  mutate(value = paste0(value, "%)")) %>%
  pivot_wider(names_from = intensity_class, values_from = value)

low_high_summary <- intensity_class_summary_counts %>%
  mutate(group = case_when(
    intensity_class %in% intensity_levels[1:2] ~ "Low",
    TRUE ~ "High"
  )) %>%
  group_by(activity, group) %>%
  summarise(n_plots = sum(n_plots), .groups = "drop") %>%
  group_by(activity) %>%
  mutate(share = round(100 * n_plots / sum(n_plots), 1)) %>%
  unite(col = "value", n_plots, share, sep = " (") %>%
  mutate(value = paste0(value, "%)")) %>%
  pivot_wider(names_from = group, values_from = value)

intensity_class_summary_final <- intensity_class_summary_percent %>%
  left_join(low_high_summary, by = "activity") %>%
  mutate(high_n = as.numeric(str_extract(High, "^[0-9]+"))) %>%
  arrange(desc(high_n)) %>%
  select(-high_n)

intensity_class_summary_final


# 3. Derived variables -------------------------------------------

## 3a. Importance values (IV)
dat_iv <- dat_overlap_mng_upd2 %>%
  mutate(size_height = hgt_est * n)

iv_sub  <- calc_iv_subplot(dat_iv, size_height)
iv_plot <- calc_iv_plot(dat_iv, size_height)

iv_max_sub <- iv_sub %>%
  group_by(plot, year, subplot) %>%
  slice_max(order_by = IV, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(IVmax_species = species) %>%
  select(plot, year, subplot, IV, IVmax_species)

iv_max_plot <- iv_plot %>%
  group_by(plot, year) %>%
  slice_max(order_by = IV, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(IVmax_species = species) %>%
  select(plot, year, IV, IVmax_species)

## IV by leaf type (coniferous vs deciduous)
iv_leaf_sub <- iv_sub %>%
  left_join(species_class, by = "species") %>%
  filter(!is.na(leaf_type)) %>%
  group_by(plot, year, subplot, leaf_type) %>%
  summarise(IV = sum(IV, na.rm = TRUE),
            RA = sum(RA, na.rm = TRUE),
            RS = sum(RS, na.rm = TRUE),
            .groups = "drop")

iv_leaf_plot <- iv_plot %>%
  left_join(species_class, by = "species") %>%
  filter(!is.na(leaf_type)) %>%
  group_by(plot, year, leaf_type) %>%
  summarise(IV = sum(IV, na.rm = TRUE),
            RA = sum(RA, na.rm = TRUE),
            RS = sum(RS, na.rm = TRUE),
            .groups = "drop")

iv_leaf_plot_wide <- iv_leaf_plot %>%
  dplyr::select(plot, year, leaf_type, IV) %>%
  tidyr::pivot_wider(names_from = leaf_type, values_from = IV, names_prefix = "IV_")

## 3b. Community-weighted means (CWM: shade & drought tolerance)
wvar <- "n"   # weight by stem counts; alternatives: "basal_area_cm2" or "hgt_est"

cwm_subplot <- dat_overlap_mng_upd2 %>%
  mutate(w = wfun(.data[[wvar]])) %>%
  filter(w > 0) %>%
  group_by(plot, subplot, year, time_snc_full_disturbance) %>%
  summarise(
    stems_with_traits = sum(w[!is.na(Shade_tolerance) & !is.na(Drought_tolerance)], na.rm = TRUE),
    stems_total       = sum(w, na.rm = TRUE),
    CWM_shade   = ifelse(stems_with_traits > 0,
                         sum(w * Shade_tolerance,  na.rm = TRUE) / stems_with_traits, NA_real_),
    CWM_drought = ifelse(stems_with_traits > 0,
                         sum(w * Drought_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    trait_coverage = stems_with_traits / pmax(stems_total, 1e-9),
    .groups = "drop"
  ) %>%
  select(-stems_total, -stems_with_traits, -trait_coverage)

cwm_plot <- dat_overlap_mng_upd2 %>%
  filter(species != "ots1") %>%
  mutate(w = wfun(.data[[wvar]])) %>%
  filter(w > 0) %>%
  group_by(plot, year, time_snc_full_disturbance) %>%
  summarise(
    stems_with_traits = sum(w[!is.na(Shade_tolerance) & !is.na(Drought_tolerance)], na.rm = TRUE),
    stems_total       = sum(w, na.rm = TRUE),
    CWM_shade   = ifelse(stems_with_traits > 0,
                         sum(w * Shade_tolerance,  na.rm = TRUE) / stems_with_traits, NA_real_),
    CWM_drought = ifelse(stems_with_traits > 0,
                         sum(w * Drought_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    trait_coverage = stems_with_traits / pmax(stems_total, 1e-9),
    .groups = "drop"
  ) %>%
  select(-stems_total, -stems_with_traits, -trait_coverage)

## 3c. Spruce share
spruce_share_plot <- dat_overlap_mng_upd2 %>%
  filter(!is.na(n)) %>%
  group_by(plot, year) %>%
  summarise(total_stems  = n(),
            spruce_stems = sum(species == "piab"),
            spruce_share = spruce_stems / total_stems,
            .groups = "drop") %>%
  select(plot, year, spruce_share)

spruce_share_sub <- dat_overlap_mng_upd2 %>%
  filter(!is.na(n)) %>%
  group_by(plot, subplot, year) %>%
  summarise(total_stems  = n(),
            spruce_stems = sum(species == "piab"),
            spruce_share = spruce_stems / total_stems,
            .groups = "drop") %>%
  select(plot, subplot, year, spruce_share)

## 3d. Subplot summary table
field_sub_summ <- dat_overlap_mng_upd2 %>%
  mutate(n = coalesce(n, 0L)) %>%
  group_by(plot, subplot, year,
           time_snc_full_disturbance, time_snc_part_disturbance,
           disturbance_year, forest_year, disturbance_length,
           clear, grndwrk, logging_trail, planting, anti_browsing) %>%
  summarise(
    stems_total       = sum(n, na.rm = TRUE),
    sp_richness       = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp        = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp       = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
    mean_hgt          = if (stems_total > 0) weighted.mean(hgt_est[n > 0], n[n > 0], na.rm = TRUE) else NA_real_,
    var_hgt = {
      if (stems_total > 1) {
        sel <- (n > 0) & is.finite(hgt_est)
        h <- hgt_est[sel]; ww <- n[sel]
        if (length(h) >= 1 && sum(ww) > 1) {
          mu <- weighted.mean(h, ww)
          v  <- sum(ww * (h - mu)^2) / (sum(ww) - 1)
          if (is.nan(v)) NA_real_ else v
        } else NA_real_
      } else NA_real_
    },
    cv_hgt = if (is.finite(mean_hgt) && mean_hgt > 0 && !is.na(var_hgt))
      sqrt(var_hgt) / mean_hgt else
        if (stems_total > 1 && is.finite(mean_hgt) && mean_hgt > 0) 0 else NA_real_,
    .groups = "drop"
  ) %>%
  left_join(cwm_subplot)

field_sub_summ_cleaned <- field_sub_summ %>%
  group_by(plot, subplot, year) %>%
  summarise(across(where(is.numeric), safe_max), .groups = "drop")

## 3e. Plot summary table (pooled across subplots)
plot_metrics_pooled <- dat_overlap_mng_upd2 %>%
  mutate(n = coalesce(n, 0L)) %>%
  group_by(plot, year) %>%
  summarise(
    stems_total       = sum(n, na.rm = TRUE),
    sp_richness       = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp        = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp       = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    mean_hgt          = if (stems_total > 0) weighted.mean(hgt_est, n, na.rm = TRUE) else NA_real_,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
    var_hgt = {
      if (stems_total > 1) {
        sel <- (n > 0) & is.finite(hgt_est)
        h <- hgt_est[sel]; ww <- n[sel]
        if (length(h) >= 1 && sum(ww) > 1) {
          mu <- weighted.mean(h, ww)
          v  <- sum(ww * (h - mu)^2) / (sum(ww) - 1)
          if (is.nan(v)) NA_real_ else v
        } else NA_real_
      } else NA_real_
    },
    cv_hgt = if (is.finite(mean_hgt) && mean_hgt > 0 && !is.na(var_hgt))
      sqrt(var_hgt) / mean_hgt else
        if (stems_total > 1 && is.finite(mean_hgt) && mean_hgt > 0) 0 else NA_real_,
    .groups = "drop"
  ) %>%
  left_join(cwm_plot,          by = join_by(plot, year)) %>%
  left_join(df_mng_plot,       by = join_by(plot, year)) %>%
  left_join(spruce_share_plot, by = join_by(plot, year)) %>%
  left_join(iv_max_plot,       by = join_by(plot, year))

## 3f. Plot context
df_plot_context <- dat_overlap_mng_upd2 %>%
  dplyr::select(plot, year,
                pre_dist_trees_n, area_m2, pre_dist_dens_ha,
                time_snc_full_disturbance,
                disturbance_year, forest_year, disturbance_length,
                x, y) %>%
  group_by(plot, year) %>%
  mutate(x = mean(x, na.rm = TRUE), y = mean(y, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct()

## 3g. Bind subplot + plot into unified analysis table
sub_df <- field_sub_summ_cleaned %>%
  transmute(
    ID      = subplot, plot_id = plot, year, level = "subplot",
    dens_m2 = stems_total / area_subplot_m2,
    cv_hgt, mean_hgt, sp_richness, shannon_sp, evenness_sp,
    effective_numbers, w = stems_total
  ) %>%
  mutate(dens_ha = dens_m2 * 10000) %>%
  left_join(df_plot_context,  by = c("plot_id" = "plot", "year")) %>%
  left_join(cwm_subplot,      by = c("plot_id" = "plot", "year",
                                     "ID" = "subplot",
                                     "time_snc_full_disturbance")) %>%
  left_join(df_mng_plot,      by = c("plot_id" = "plot", "year")) %>%
  left_join(spruce_share_sub, by = c("plot_id" = "plot", "year", "ID" = "subplot")) %>%
  left_join(iv_max_sub,       by = c("plot_id" = "plot", "year", "ID" = "subplot")) %>%
  distinct()

plot_df <- plot_metrics_pooled %>%
  transmute(
    ID      = plot, plot_id = plot, year, level = "plot",
    dens_m2 = stems_total / area_plot_m2,
    cv_hgt, mean_hgt, sp_richness, shannon_sp, evenness_sp,
    effective_numbers, w = stems_total
  ) %>%
  mutate(dens_ha  = dens_m2 * 10000,
         mean_hgt = replace_na(mean_hgt, 0)) %>%
  left_join(df_plot_context,   by = c("plot_id" = "plot", "year")) %>%
  left_join(cwm_plot,          by = c("plot_id" = "plot", "year",
                                      "time_snc_full_disturbance")) %>%
  left_join(df_mng_plot,       by = c("plot_id" = "plot", "year")) %>%
  left_join(spruce_share_plot, by = c("plot_id" = "plot", "year")) %>%
  left_join(iv_max_plot,       by = c("plot_id" = "plot", "year")) %>%
  distinct()

both_levels_re2 <- bind_rows(sub_df, plot_df) %>%
  mutate(
    time_snc_full_disturbance = pmin(time_snc_full_disturbance, 8),
    level   = factor(level, levels = c("subplot", "plot")),
    plot_id = factor(plot_id),
    w       = pmin(pmax(w, 1), 50),
    year_f  = factor(year)
  ) %>%
  mutate(across(all_of(c("cv_hgt", "mean_hgt", "sp_richness", "effective_numbers")),
                ~ replace_na(.x, 0))) %>%
  mutate(
    cv_hgt_pos        = ifelse(cv_hgt <= 0, 1e-4, cv_hgt),
    effective_numbers = ifelse(effective_numbers == 0, 1e-4, effective_numbers),
    cv_hgt_present    = as.integer(cv_hgt > 0),
    active_regeneration = (planting_intensity + anti_browsing_intensity) / 2
  )

table(both_levels_re2$cv_hgt_present)

df_plot_clean <- both_levels_re2 %>% filter(level == "plot")
df_sub_clean  <- both_levels_re2 %>% filter(level == "subplot")

## 3h. AEF export table (beta added after §4)
plot_df_AEF <- plot_df %>%
  select(-level, -ID, -w) %>%
  mutate(year    = as.integer(year),
         plot_id = as.character(plot_id)) %>%
  rename(plot = plot_id) %>%
  left_join(iv_leaf_plot_wide, by = join_by(plot, year)) %>%
  mutate(
    time_snc_full_disturbance = pmin(time_snc_full_disturbance, 8),
    plot_id = factor(plot),
    year_f  = factor(year)
  ) %>%
  mutate(across(all_of(c("cv_hgt", "mean_hgt", "sp_richness", "effective_numbers")),
                ~ replace_na(.x, 0)))


# 4. Species composition & beta diversity -------------------------

## 4a. Stem share per species per year
df <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>%
  filter(!is.na(species) & species != "")

year_totals <- df %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(total_trees = sum(n, na.rm = TRUE),
                   n_plots     = dplyr::n_distinct(plot),
                   n_subplots  = dplyr::n_distinct(subplot),
                   .groups = "drop")

species_stem_share_year <- df %>%
  dplyr::group_by(year, species) %>%
  dplyr::summarise(stems         = sum(n, na.rm = TRUE),
                   plots_present = dplyr::n_distinct(plot[n > 0]),
                   .groups = "drop") %>%
  dplyr::left_join(year_totals %>% dplyr::select(year, total_trees), by = "year") %>%
  dplyr::mutate(share = round(100 * stems / total_trees, 2)) %>%
  dplyr::arrange(dplyr::desc(share), species)

top10_by_year <- species_stem_share_year %>%
  dplyr::group_by(year) %>%
  dplyr::slice_max(order_by = share, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

v_top_species <- top10_by_year %>%
  dplyr::distinct(species) %>%
  dplyr::pull(species) %>%
  sort()

species_levels <- top10_by_year %>%
  group_by(species) %>%
  summarise(share_all = sum(share, na.rm = TRUE), .groups = "drop") %>%
  arrange(share_all) %>%
  pull(species)

## 4b. Single-species plots
plots_with_species_counts <- df %>%
  filter(n >= 1) %>%
  group_by(plot) %>%
  summarise(n_species = n_distinct(species), .groups = "drop")

single_species_plots <- plots_with_species_counts %>%
  dplyr::filter(n_species == 1)

species_per_single_species_plot <- df %>%
  filter(n >= 1) %>%
  semi_join(single_species_plots, by = "plot") %>%
  group_by(plot, species) %>%
  summarise(n_trees = sum(n), .groups = "drop")

summary_species <- species_per_single_species_plot %>%
  count(species, name = "single_plots") %>%
  arrange(desc(single_plots)) %>%
  mutate(share = single_plots / n_plots_total)

print(summary_species)

## 4c. Stem density and species occurrence
df_stem_dens_species <- df %>%
  group_by(plot, year, species, n_subplots) %>%
  summarize(sum_n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  mutate(scaling_factor = 10000 / (n_subplots * 4),
         stem_dens      = sum_n * scaling_factor)

df_stem_dens_species_sum <- df_stem_dens_species %>%
  group_by(species, year) %>%
  summarise(stem_dens = sum(stem_dens, na.rm = TRUE), .groups = "drop") %>%
  mutate(stem_dens_avg = stem_dens / n_plots_total)

df_stem_dens_species <- df_stem_dens_species %>%
  filter(sum_n > 0, species %in% v_top_species) %>%
  dplyr::group_by(species, year) %>%
  dplyr::mutate(median_stem_density = median(stem_dens, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  mutate(species = factor(species, levels = rev(v_top_species)))

species_levels <- factor(v_top_species)

total_plots <- dat_overlap %>% pull(plot) %>% n_distinct()

species_occurence <- df_stem_dens_species %>%
  dplyr::filter(sum_n > 0) %>%
  distinct(species, year, plot) %>%
  count(species, year, name = "n_plots") %>%
  mutate(share_of_plots = n_plots / total_plots * 100)

## 4d. Beta diversity (Jaccard within-plot turnover)
comm_sub <- dat_overlap_mng_upd2 %>%
  mutate(n = coalesce(n, 0L)) %>%
  filter(!is.na(species), species != "") %>%
  group_by(plot, year, subplot, species) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = species, values_from = n, values_fill = 0)

beta_within_plotyear <- function(df_plotyear) {
  meta <- df_plotyear %>% dplyr::select(plot, year, subplot)
  mat  <- df_plotyear %>% dplyr::select(-plot, -year, -subplot) %>% as.data.frame()
  keep <- rowSums(mat) > 0
  mat  <- mat[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  n_sub <- nrow(mat)

  if (n_sub < 2) {
    return(tibble::tibble(
      plot = unique(meta$plot), year = unique(meta$year),
      n_subplots_nonempty = n_sub,
      beta_jaccard_mean = NA_real_, beta_bray_mean = NA_real_,
      alpha_rich_mean   = if (n_sub == 1) sum(mat[1, ] > 0) else NA_real_,
      gamma_rich        = if (n_sub == 1) sum(colSums(mat) > 0) else NA_real_,
      beta_whittaker    = NA_real_
    ))
  }

  mat_pa <- (mat > 0) * 1
  d_jac  <- vegan::vegdist(mat_pa, method = "jaccard", binary = TRUE)
  d_bray <- vegan::vegdist(mat,    method = "bray")
  alpha  <- rowSums(mat_pa)
  gamma  <- sum(colSums(mat_pa) > 0)

  tibble::tibble(
    plot = unique(meta$plot), year = unique(meta$year),
    n_subplots_nonempty = n_sub,
    beta_jaccard_mean   = mean(as.numeric(d_jac),  na.rm = TRUE),
    beta_bray_mean      = mean(as.numeric(d_bray), na.rm = TRUE),
    alpha_rich_mean     = mean(alpha),
    gamma_rich          = gamma,
    beta_whittaker      = gamma / mean(alpha)
  )
}

beta_plotyear <- comm_sub %>%
  group_by(plot, year) %>%
  group_split() %>%
  purrr::map_dfr(beta_within_plotyear)

## Add beta to AEF table
eps <- 0.001
plot_df_AEF2 <- plot_df_AEF %>%
  left_join(beta_plotyear) %>%
  mutate(beta_jaccard_mean = pmin(pmax(beta_jaccard_mean, eps), 1 - eps))

## Spearman correlations (quick checks)
cor.test(plot_df_AEF2$time_snc_full_disturbance, plot_df_AEF2$beta_jaccard_mean, method = "spearman")
cor.test(plot_df_AEF2$planting_intensity,         plot_df_AEF2$spruce_share,      method = "spearman")
cor.test(both_levels_re2$planting_intensity, both_levels_re2$anti_browsing_intensity, method = "spearman")

both_levels_re2 %>%
  mutate(planted   = planting_intensity      > 0,
         protected = anti_browsing_intensity > 0) %>%
  count(planted, protected)


# 5. Models (GAMs) -----------------------------------------------
# All GAMs: REML, random effect s(plot_id, bs="re")
# Families: tw(log) continuous, nb(log) counts, binomial(logit) binary

## 5a. Mean height
gam_mean_hgt_both <- gam(
  mean_hgt ~
    planting_intensity * anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 4) +
    grndwrk_intensity + year_f + level +
    s(plot_id, bs = "re"),
  data   = both_levels_re2 %>% dplyr::filter(mean_hgt < 6),
  family = tw(link = "log"),
  method = "REML"
)

## 5b. CV height — binary (structural heterogeneity present/absent)
gam_cv_hgt_bin_both <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    planting_intensity * anti_browsing_intensity +
    grndwrk_intensity + level + year_f +
    s(plot_id, bs = "re"),
  data   = both_levels_re2,
  family = binomial(link = "logit"),
  method = "REML"
)

## 5c. CV height — positive part
gam_cv_hgt_pos_both <- gam(
  cv_hgt_pos ~
    planting_intensity * anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity + level + year_f +
    s(plot_id, bs = "re"),
  data   = both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)

## 5d. Effective species number
gam_eff_both <- gam(
  effective_numbers ~
    planting_intensity + anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 3) +
    grndwrk_intensity + level + year_f +
    s(plot_id, bs = "re"),
  data   = both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)

## 5e. Species richness
gam_rich_both <- gam(
  sp_richness ~
    planting_intensity * anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity + level + year_f +
    s(plot_id, bs = "re"),
  data   = both_levels_re2,
  family = nb(link = "log"),
  method = "REML"
)

## 5f. Beta diversity (Jaccard)
m_beta_add <- mgcv::gam(
  beta_jaccard_mean ~
    time_snc_full_disturbance +
    planting_intensity + anti_browsing_intensity +
    grndwrk_intensity + year_f +
    s(plot_id, bs = "re"),
  data   = plot_df_AEF2,
  method = "REML"
)

## Final model list
fin.models <- list(
  hgt   = gam_mean_hgt_both,
  cvpos = gam_cv_hgt_pos_both,
  eff   = gam_eff_both,
  rich  = gam_rich_both,
  beta  = m_beta_add
)

lapply(fin.models, summary)

## 5g. emmeans: scale effect (subplot vs plot)
fin.models.levels <- fin.models[names(fin.models) != "beta"]

emm_results <- map_dfr(
  fin.models.levels,
  ~ broom::tidy(emmeans(.x, ~ level, type = "response")),
  .id = "model"
)

emm_change <- emm_results %>%
  select(model, level, response) %>%
  pivot_wider(names_from = level, values_from = response) %>%
  mutate(diff           = plot - subplot,
         ratio          = plot / subplot,
         percent_change = (plot - subplot) / subplot * 100)

emm_change

## 5h. emmeans: management effect (intensity 0 -> 1)
response_labels <- c(
  hgt   = "Mean height [m]",
  cvpos = "CV [%]",
  eff   = "Effective species [#]",
  rich  = "Species richness [#]",
  beta  = "Turnover [dim.]"
)

term_labels <- c(
  planting_intensity      = "Planting",
  anti_browsing_intensity = "Browsing\nprotection",
  grndwrk_intensity       = "Soil\npreparation"
)

get_mng_emm <- function(model, model_name, focal_term) {
  at_list <- switch(
    focal_term,
    "planting_intensity"      = list(planting_intensity      = c(0, 1)),
    "anti_browsing_intensity" = list(anti_browsing_intensity = c(0, 1)),
    "grndwrk_intensity"       = list(grndwrk_intensity       = c(0, 1))
  )
  emm_call <- if (model_name == "beta") {
    emmeans::emmeans(model,
                     specs   = stats::as.formula(paste0("~ ", focal_term)),
                     at      = at_list,
                     type    = "response",
                     exclude = "level")
  } else {
    emmeans::emmeans(model,
                     specs = stats::as.formula(paste0("~ ", focal_term)),
                     at    = at_list,
                     type  = "response")
  }
  emm_call %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    dplyr::mutate(model = model_name, focal_term = focal_term)
}

emm_mng <- map_dfr(
  names(fin.models),
  function(nm) {
    bind_rows(
      get_mng_emm(fin.models[[nm]], nm, "planting_intensity"),
      get_mng_emm(fin.models[[nm]], nm, "anti_browsing_intensity"),
      get_mng_emm(fin.models[[nm]], nm, "grndwrk_intensity")
    )
  }
)

pvals_mng <- map_dfr(fin.models, ~ broom::tidy(.x, parametric = TRUE), .id = "model") %>%
  filter(term %in% c("planting_intensity", "anti_browsing_intensity", "grndwrk_intensity")) %>%
  transmute(model, focal_term = term, p.value,
            p_label = case_when(p.value < 0.001 ~ "<0.001",
                                TRUE ~ formatC(p.value, format = "f", digits = 3)))

emm_mng2 <- emm_mng %>%
  mutate(
    intensity = case_when(
      focal_term == "planting_intensity"      ~ planting_intensity,
      focal_term == "anti_browsing_intensity" ~ anti_browsing_intensity,
      focal_term == "grndwrk_intensity"       ~ grndwrk_intensity
    ),
    response_plot = dplyr::coalesce(response, emmean),
    response_lab  = recode(model, !!!response_labels),
    term          = recode(focal_term, !!!term_labels)
  )

model_intensity_all_df_pct_mng <- emm_mng2 %>%
  select(model, focal_term, term, response_lab, intensity, response_plot, lower.CL, upper.CL) %>%
  rename(response = response_plot) %>%
  pivot_wider(names_from  = intensity,
              values_from = c(response, lower.CL, upper.CL),
              names_sep   = "_") %>%
  left_join(pvals_mng, by = c("model", "focal_term")) %>%
  mutate(
    estimate_pct = 100 * (response_1 - response_0) / response_0,
    lower_pct    = 100 * (lower.CL_1 - upper.CL_0) / upper.CL_0,
    upper_pct    = 100 * (upper.CL_1 - lower.CL_0) / lower.CL_0,
    sig_col      = ifelse(p.value < 0.05, "sig", "n.s."),
    response     = factor(response_lab,
                          levels = c("Mean height [m]", "CV [%]",
                                     "Effective species [#]",
                                     "Species richness [#]",
                                     "Turnover [dim.]"))
  )

model_intensity_summary <- model_intensity_all_df_pct_mng %>%
  dplyr::select(response, term, estimate_pct, lower_pct, upper_pct, p_label) %>%
  dplyr::rename(effect_pct  = estimate_pct,
                ci_low_pct  = lower_pct,
                ci_high_pct = upper_pct)
model_intensity_summary


# 6. Figures -----------------------------------------------------

## 6a. Species composition (stem share + occurrence)
top10_by_year_clean <- top10_by_year %>%
  dplyr::select(year, species, share) %>%
  tidyr::complete(year, species = species_levels, fill = list(share = 0)) %>%
  mutate(species = factor(species, levels = species_levels),
         year    = factor(year))

p_bar <- top10_by_year_clean %>%
  ggplot(aes(x = share, y = species, fill = species, alpha = year)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           aes(colour = factor(year))) +
  scale_x_continuous(labels = label_number(accuracy = 1, suffix = ""),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_discrete(limits = species_levels, labels = species_labels, drop = FALSE) +
  scale_colour_manual(values = c("2023" = "grey50", "2025" = "black"), guide = "none") +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_alpha_manual(name = "Survey year",
                     values  = c("2023" = 0.45, "2025" = 1.00),
                     guide   = "none",
                     breaks  = c("2025", "2023")) +
  labs(x = "Stems share [%]", y = "") +
  theme_classic2(base_size = 10) +
  theme(axis.text.y   = element_text(face = "italic", size = 8),
        plot.margin   = margin(t = 12, r = 5, b = 5, l = 5))

p_occurence <- species_occurence %>%
  filter(species %in% v_top_species) %>%
  mutate(species = factor(species, levels = species_levels)) %>%
  ggplot(aes(x = share_of_plots, y = species,
             fill = species, alpha = factor(year))) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           aes(colour = factor(year))) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_alpha_manual(name = "Survey year",
                     values = c("2025" = 1, "2023" = 0.45),
                     breaks = c("2025", "2023")) +
  scale_colour_manual(name = "Survey year",
                      values = c("2025" = "black", "2023" = "grey50"),
                      breaks = c("2025", "2023")) +
  labs(x = "Species occurence [%]", y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y        = element_blank(),
        legend.position    = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.key         = element_rect(fill = NA),
        legend.title       = element_text(size = 9),
        legend.text        = element_text(size = 8),
        legend.key.width   = unit(1, "lines"),
        legend.key.height  = unit(0.4, "lines"),
        plot.margin        = margin(t = 12, r = 5, b = 5, l = 5))

p_species_composition <- ggarrange(
  p_bar, p_occurence,
  ncol = 2, common.legend = FALSE, align = "h",
  widths = c(1.5, 1),
  labels = c("[a]", "[b]"),
  font.label = list(size = 10, face = "plain"),
  label.x = 0.02, label.y = 1.01
)

## 6b. Management intensity plot (simpler stacked bar)
fill_colors <- brewer.pal(length(intensity_levels), "YlOrRd")

p_management_intensity_plot_simpler <- mng_shifted %>%
  filter(activity != "logging_trail_intensity") %>%
  droplevels() %>%
  ggplot(aes(x = proportion, y = activity, fill = intensity_class_plot)) +
  geom_col(width = 0.4, color = "black") +
  scale_fill_manual(values = fill_colors,
                    name   = "Management intensity\nclass",
                    breaks = intensity_levels) +
  scale_y_discrete(labels = activity_intens_labels) +
  scale_x_continuous(labels = abs, name = "Plots share [%]") +
  ylab("") +
  theme_classic2() +
  theme(legend.position = "right", axis.text.y = element_text(size = 10))

p_management_bin_plot <- mng_sub_conv %>%
  filter(activity != "logging_trail_intensity") %>%
  droplevels() %>%
  ggplot(aes(x = proportion_plot, y = activity, fill = intensity_binary)) +
  geom_col(width = 0.2, color = "black") +
  scale_fill_manual(values = c("no" = "forestgreen", "yes" = "red")) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 0.8, lty = "dashed") +
  scale_y_discrete(labels = activity_intens_labels) +
  scale_x_continuous(labels = abs, name = "Plots share [%]") +
  annotate("text", x = -60, y = 4.5, label = "Not Applied", hjust = 0, size = 2.5, fontface = "bold") +
  annotate("text", x =  80, y = 4.5, label = "Applied",     hjust = 1, size = 2.5, fontface = "bold") +
  ylab("") +
  theme_classic2() +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 8))

## 6c. Disturbance history
plot_context_chars <- dat_overlap_mng_upd2 %>%
  dplyr::select(plot, year, status, disturbance_year, forest_year,
                disturbance_length, time_snc_full_disturbance, time_snc_part_disturbance,
                planting_intensity, clear_intensity, grndwrk_intensity,
                logging_trail_intensity, anti_browsing_intensity) %>%
  distinct() %>%
  mutate(disturbance_year = case_when(disturbance_year < 2018 ~ 2018,
                                      disturbance_year > 2022 ~ 2022,
                                      TRUE ~ disturbance_year),
         time_since_disturbance = year - disturbance_year)

summary_disturbance_year <- plot_context_chars %>%
  count(disturbance_year, name = "n") %>%
  mutate(share = round(100 * n / sum(n), 1))

summary_time_since <- plot_context_chars %>%
  count(time_since_disturbance, name = "n") %>%
  mutate(share = round(100 * n / sum(n), 1))

p_hist_dist_year <- plot_context_chars %>%
  count(disturbance_year) %>%
  ggplot(aes(x = disturbance_year, y = n)) +
  geom_col(fill = "grey80", color = "black", width = 0.8) +
  scale_x_continuous(breaks = seq(2018, 2022, 2)) +
  labs(x = "Disturbance Year\n", y = "Number of Plots [#]") +
  theme_classic2() +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8))

p_hist_time_since_dist <- plot_context_chars %>%
  count(time_since_disturbance) %>%
  ggplot(aes(x = time_since_disturbance, y = n)) +
  geom_col(fill = "grey80", color = "black", width = 0.8) +
  labs(x = "Time since disturbance\n[years]", y = "Number of Plots [#]") +
  theme_classic2() +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8))

p_combined_disturb_fig <- ggarrange(
  p_hist_dist_year, p_hist_time_since_dist,
  labels = c("[a]", "[b]"),
  font.label = list(size = 10, face = "plain"),
  widths = c(0.9, 1.1), ncol = 2, nrow = 1
)

p_combined_management_intens <- ggarrange(
  p_combined_disturb_fig,
  p_management_intensity_plot_simpler,
  labels = c(" ", "[c]"),
  font.label = list(size = 10, face = "plain"),
  ncol = 1, nrow = 2,
  heights = c(1.2, 1.2)
)

## 6d. Cross-scale violin distributions
response_vars <- c("mean_hgt", "cv_hgt", "effective_numbers", "sp_richness")
discrete_vars_plot <- c("Species richness\n[#]", "Species diversity\n[Effective #]")

both_levels_long <- both_levels_re2 %>%
  pivot_longer(cols = all_of(response_vars), names_to = "variable", values_to = "value")

both_levels_long_capped <- both_levels_long %>%
  group_by(variable) %>%
  mutate(value_capped = ifelse(value > quantile(value, 0.95, na.rm = TRUE),
                               quantile(value, 0.95, na.rm = TRUE), value)) %>%
  ungroup() %>%
  mutate(value_capped = ifelse(variable == "cv_hgt", value_capped * 100, value_capped),
         variable = factor(variable, levels = c("mean_hgt", "cv_hgt",
                                                "effective_numbers", "sp_richness")),
         variable = recode(variable,
                           effective_numbers = "Species diversity\n[Effective #]",
                           sp_richness       = "Species richness\n[#]",
                           mean_hgt          = "Mean height\n[m]",
                           cv_hgt            = "Height variability\n[CV, %]"),
         is_discrete = variable %in% discrete_vars_plot)

facet_ymins <- both_levels_long_capped %>%
  distinct(variable) %>%
  mutate(ymin = c(-0.5, -20, -1.8, -0.9), ymax = 0)

wilcox_results <- both_levels_long_capped %>%
  group_by(variable) %>%
  summarise(test = list(wilcox.test(value_capped ~ level, data = cur_data())),
            .groups = "drop") %>%
  mutate(tidied = map(test, broom::tidy)) %>%
  unnest(tidied) %>%
  mutate(p.value   = round(p.value, 4),
         p.adj     = round(p.adjust(p.value, method = "BH"), 4),
         p.display = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))) %>%
  select(variable, p.value, p.adj, p.display)

label_positions <- both_levels_long_capped %>%
  group_by(variable) %>%
  summarise(y = max(value_capped, na.rm = TRUE) * 1.20)

wilcox_plot_labels <- wilcox_results %>%
  left_join(label_positions, by = "variable")

p_violin <- ggplot() +
  gghalves::geom_half_violin(
    data = both_levels_long_capped %>% filter(!is_discrete),
    aes(x = level, y = value_capped, fill = level),
    side = "l", color = NA, trim = FALSE, scale = "width"
  ) +
  gghalves::geom_half_violin(
    data = both_levels_long_capped %>% filter(is_discrete),
    aes(x = level, y = value_capped, fill = level),
    side = "l", color = NA, trim = FALSE, scale = "width", adjust = 2
  ) +
  geom_boxplot(data = both_levels_long_capped,
               aes(x = level, y = value_capped),
               width = 0.1, outlier.shape = NA, color = "black") +
  geom_rect(data = facet_ymins, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "white", color = NA) +
  stat_summary(data = both_levels_long_capped,
               aes(x = level, y = value_capped),
               fun = mean, geom = "point",
               shape = 21, size = 2, fill = "red", color = "black") +
  scale_fill_manual(values = c("subplot" = "grey50", "plot" = "grey80")) +
  geom_text(data = wilcox_plot_labels, inherit.aes = FALSE,
            aes(x = 1.5, y = y, label = p.display), size = 3) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4, nrow = 1) +
  theme_classic(base_size = 10) +
  labs(x = NULL, y = NULL, fill = "Level") +
  theme(legend.position = "none",
        strip.text  = element_text(size = 9, face = "plain"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

summary_stats_violin <- both_levels_long_capped %>%
  group_by(variable, level) %>%
  summarise(n      = n(),
            mean   = mean(value_capped, na.rm = TRUE),
            median = median(value_capped, na.rm = TRUE),
            Q1     = quantile(value_capped, 0.25, na.rm = TRUE),
            Q3     = quantile(value_capped, 0.75, na.rm = TRUE),
            IQR    = IQR(value_capped, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(variable, level)

## 6e. Time since disturbance x scale (panel)
ylim_hgt  <- c(0.2, 2.0)
ylim_cv   <- c(10, 85)
ylim_eff  <- c(1, 12)
ylim_rich <- c(0.8, 7)

y_time_hgt  <- ylim_hgt[2]  * 0.98
y_time_cv   <- ylim_cv[2]   * 0.98
y_time_eff  <- ylim_eff[2]  * 0.98
y_time_rich <- ylim_rich[2] * 0.98

p_hgt <- pp(gam_mean_hgt_both,
            terms   = c("time_snc_full_disturbance[0:7]", "level"),
            ylab    = "Mean height [m]",
            annot   = "Time:\np = 0.004",
            annot_y = y_time_hgt) + coord_cartesian(ylim = ylim_hgt)

p_cv_pos <- pp(gam_cv_hgt_pos_both,
               terms   = c("time_snc_full_disturbance[0:7]", "level"),
               ylab    = "CV [%]",
               scale_y = 100,
               annot   = "Time:\np = 0.668",
               annot_y = y_time_cv) + coord_cartesian(ylim = ylim_cv)

p_eff <- pp(gam_eff_both,
            terms   = c("time_snc_full_disturbance[0:7]", "level"),
            xlab    = "Time since disturbance\n(years)",
            ylab    = "Effective species [#]",
            annot   = "Time:\np = 0.452",
            annot_y = y_time_eff) + coord_cartesian(ylim = ylim_eff)

p_rich <- pp(gam_rich_both,
             terms   = c("time_snc_full_disturbance[0:7]", "level"),
             xlab    = "Time since disturbance\n(years)",
             ylab    = "Species richness [#]",
             annot   = "Time:\np = 0.896",
             annot_y = y_time_rich) + coord_cartesian(ylim = ylim_rich)

p_inset_hgt  <- pp_inset_model(gam_mean_hgt_both,   scale_y = 1,   p_lab = "Scale:\np < 0.001", annot_y = y_time_hgt)
p_inset_cv   <- pp_inset_model(gam_cv_hgt_pos_both, scale_y = 100, p_lab = "Scale:\np < 0.001", annot_y = y_time_cv)
p_inset_eff  <- pp_inset_model(gam_eff_both,        scale_y = 1,   p_lab = "Scale:\np < 0.001", annot_y = y_time_eff)
p_inset_rich <- pp_inset_model(gam_rich_both,       scale_y = 1,   p_lab = "Scale:\np < 0.001", annot_y = y_time_rich)

pair_hgt  <- (p_hgt  + p_inset_hgt)  + plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[a]")
pair_cv   <- (p_cv_pos + p_inset_cv) + plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[b]")
pair_eff  <- (p_eff  + p_inset_eff)  + plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[c]")
pair_rich <- (p_rich + p_inset_rich) + plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[d]")

p_final <- (pair_hgt | pair_cv) / (pair_eff | pair_rich) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

## 6f. Management effect (% change 0 -> 1)
p_model_response <- ggplot(
  model_intensity_all_df_pct_mng,
  aes(x = response, y = estimate_pct, fill = response)
) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  geom_col(width = 0.7, color = NA) +
  geom_errorbar(aes(ymin = lower_pct, ymax = upper_pct),
                width = 0.15, linewidth = 0.5, color = "grey30") +
  geom_text(aes(label    = p_label,
                y        = upper_pct + 12,
                fontface = ifelse(p.value < 0.05, "bold", "plain")),
            size = 2.5) +
  facet_wrap(~ term, ncol = 3) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "", y = "Effect on response [%]") +
  theme_classic(base_size = 8) +
  theme(legend.position    = "none",
        strip.text         = element_text(size = 8),
        strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.6),
        axis.text.x        = element_text(angle = 30, hjust = 1),
        panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3))

## 6g. Spruce share by management
p_spruce_shares_boxplot <- ggarrange(
  both_levels_re2 %>% filter(level == "plot") %>%
    ggplot(aes(x = factor(planting_intensity), y = spruce_share * 100,
               fill = factor(planting_intensity))) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) +
    geom_jitter(alpha = 0.5, size = 0.5) +
    stat_compare_means(method = "kruskal.test", label.x = 2,
                       label.y = 1.05 * max(both_levels_re2$spruce_share * 100, na.rm = TRUE)) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "red", color = "white") +
    labs(x = "Planting Intensity\n", y = "Spruce share [%]") +
    theme_classic2() + theme(legend.position = "none"),

  both_levels_re2 %>% filter(level == "plot") %>%
    ggplot(aes(x = factor(anti_browsing_intensity), y = spruce_share * 100,
               fill = factor(anti_browsing_intensity))) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) +
    geom_jitter(alpha = 0.5, size = 0.5) +
    stat_compare_means(method = "kruskal.test", label.x = 2,
                       label.y = 1.05 * max(both_levels_re2$spruce_share * 100, na.rm = TRUE)) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "red", color = "white") +
    labs(x = "Browsing protection\nintensity", y = "Spruce share [%]") +
    theme_classic2() + theme(legend.position = "none"),

  labels = c("[a]", "[b]"), ncol = 2, align = "v",
  font.label = list(face = "plain")
)

## 6h. Height by seral stage x management (heatmap)
dat_clean <- dat_overlap_mng_upd2 %>%
  as_tibble() %>%
  dplyr::filter(!is.na(hgt_est), !is.na(n), n > 0)

height_species_mng_yr <- dat_clean %>%
  dplyr::group_by(year, species, planting_intensity, anti_browsing_intensity) %>%
  dplyr::summarise(total_stems          = sum(n, na.rm = TRUE),
                   mean_height_weighted = weighted.mean(hgt_est, w = n, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(seral_group = dplyr::if_else(species %in% earlyspecs_laura, "early", "late"))

height_seral_mng_year <- height_species_mng_yr %>%
  dplyr::group_by(year, seral_group, planting_intensity, anti_browsing_intensity) %>%
  dplyr::summarise(mean_h = weighted.mean(mean_height_weighted, w = total_stems, na.rm = TRUE),
                   stems  = sum(total_stems, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::filter(stems >= 3)

p_height_seral_mng <- ggplot(height_seral_mng_year,
       aes(x = planting_intensity, y = anti_browsing_intensity, fill = mean_h)) +
  geom_tile(color = "grey80") +
  facet_grid(seral_group ~ year) +
  scale_fill_viridis_c(limits = c(0, 2.5), option = "E",
                       direction = -1, end = 0.95, name = "Mean height (m)") +
  geom_text(aes(label = round(mean_h, 2)), color = "white", size = 2.5, fontface = "bold") +
  coord_equal() +
  theme_classic() +
  labs(x = "Planting intensity", y = "Anti-browsing intensity") +
  theme(legend.position  = "bottom",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"))


# 7. Tables & export ---------------------------------------------

## 7a. Save figures
ggsave("outFigsTest/p_species_composition.png",
       p_species_composition, width = 7, height = 3, dpi = 300)

ggsave("outFigsTest/p_violin.png",
       p_violin, width = 7, height = 2.5, dpi = 300)
ggsave("outFigsTest/p_violin.svg",
       p_violin, width = 7, height = 2.5, dpi = 300)

ggsave("outFigsTest/p_time_scale.png",
       p_final, width = 6, height = 5, dpi = 300, bg = "white")

ggsave("outFigsTest/p_model_response.png",
       p_model_response, width = 7, height = 5, dpi = 300)

ggsave("outFigsTest/p_management_intensity_plot_simpler.png",# - this one cuts
       p_management_intensity_plot_simpler, width = 5, height = 2.1, dpi = 300)

ggsave("outFigsTest/mng_intensity_plot.png",
       p_management_intensity_plot_simpler, width = 5, height = 2.1, dpi = 300)

ggsave("outFigsTest/p_combined_disturb_fig.png",
       p_combined_disturb_fig, width = 5, height = 2.5, dpi = 300)

ggsave("outFigsTest/p_combined_management_intens.png",
       p_combined_management_intens, width = 6, height = 5, dpi = 300)

ggsave("outFigsTest/spruce_share_boxplots.png",
       p_spruce_shares_boxplot, width = 7, height = 3.2, dpi = 300)

ggsave("outFigsTest/density_plot.png",
       p_height_seral_mng, width = 6, height = 4, dpi = 300)

## 7b. Model summary table -> Word
clean_term  <- function(x) x %>%
  str_replace_all("[:()`,]", "_") %>%
  str_replace_all("__+", "_") %>%
  str_remove_all("^_|_$")

p_to_signif <- function(p) case_when(
  is.na(p)   ~ NA_character_,
  p <= 0.001 ~ "***",
  p <= 0.01  ~ "**",
  p <= 0.05  ~ "*",
  p <= 0.1   ~ ".",
  TRUE       ~ "n.s."
)

format_pval <- function(p) ifelse(is.na(p), NA, formatC(p, format = "f", digits = 3))

extract_gam_summary <- function(model, model_name) {
  s        <- summary(model)
  param    <- as.data.frame(s$p.table)
  pval_col <- grep("Pr\\(>.*\\)", colnames(param), value = TRUE)
  param    <- rownames_to_column(param, "term")

  intercept_val <- param %>%
    filter(term == "(Intercept)") %>%
    pull(Estimate) %>%
    exp()

  param_clean <- param %>%
    filter(term != "(Intercept)") %>%
    select(term, p.value = all_of(pval_col))

  smooth <- if (!is.null(s$s.table)) {
    as.data.frame(s$s.table) %>%
      rownames_to_column("term") %>%
      select(term, p.value = `p-value`)
  } else {
    tibble(term = character(), p.value = numeric())
  }

  all_terms <- bind_rows(param_clean, smooth) %>%
    mutate(term_clean = clean_term(term),
           signif     = p_to_signif(p.value),
           pval_fmt   = format_pval(p.value),
           model      = model_name)

  model_metrics <- tibble(
    model              = model_name,
    response_intercept = intercept_val,
    r_squared          = s$r.sq,
    deviance_explained = s$dev.expl,
    n_samples          = s$n
  )

  list(pvalues = all_terms, metrics = model_metrics)
}

results <- map2(fin.models, names(fin.models), extract_gam_summary)

pvals_fmt <- map_dfr(results, "pvalues") %>%
  select(model, term_clean, pval_fmt) %>%
  pivot_wider(names_from = term_clean, values_from = pval_fmt, names_prefix = "p_")

metrics_df    <- bind_rows(map(results, "metrics"))
final_results <- left_join(metrics_df, pvals_fmt, by = "model")

tab_df(final_results,
       title = "Model Summary Table",
       file  = "outTable/models_summary_final_results.doc")

## 7c. Spatial export for collaborators
plot_sf_AEF <- plot_df_AEF2 %>%
  sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE)

sf::st_write(
  plot_sf_AEF,
  dsn          = "outDataShare/Karim_AEF/regeneration_chars_3035.gpkg",
  layer        = "plot",
  delete_layer = TRUE
)
