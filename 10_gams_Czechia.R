
# ------------------------------------
# Temporal analysis: Czechia
# ------------------------------------


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



# Filter data from both level: lot and subplot from EU 

both_levels_cz <- data.table::fread("outData/outData/both_levels_re2_all_countries.csv.csv") %>%
  filter(country_name == "Czechia")

plot_df_cz <- data.table::fread("outData/plot_df_AEF2_all_countries.csv") %>%
  filter(country_name == "Czechia")

sub_df_cz <- data.table::fread("outData/sub_df_AEF_all_countries.csv") %>%
  filter(country_name == "Czechia")





# 1. Models (GAMs) -----------------------------------------------
# All GAMs: REML, random effect s(plot_id, bs="re")
# Families: tw(log) continuous, nb(log) counts, binomial(logit) binary

## 5a. Mean height
gam_mean_hgt_both <- gam(
  mean_hgt ~
    planting_intensity * anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 4) +
    grndwrk_intensity + year_f + level +
    s(plot_id, bs = "re"),
  data   = both_levels_cz %>% dplyr::filter(mean_hgt < 6),
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
  data   = both_levels_cz,
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
  data   = both_levels_cz,
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
  data   = both_levels_cz,
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
  data   = both_levels_cz,
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

both_levels_long <- both_levels_cz %>%
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
  both_levels_cz %>% filter(level == "plot") %>%
    ggplot(aes(x = factor(planting_intensity), y = spruce_share * 100,
               fill = factor(planting_intensity))) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) +
    geom_jitter(alpha = 0.5, size = 0.5) +
    stat_compare_means(method = "kruskal.test", label.x = 2,
                       label.y = 1.05 * max(both_levels_cz$spruce_share * 100, na.rm = TRUE)) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "red", color = "white") +
    labs(x = "Planting Intensity\n", y = "Spruce share [%]") +
    theme_classic2() + theme(legend.position = "none"),
  
  both_levels_cz %>% filter(level == "plot") %>%
    ggplot(aes(x = factor(anti_browsing_intensity), y = spruce_share * 100,
               fill = factor(anti_browsing_intensity))) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) +
    geom_jitter(alpha = 0.5, size = 0.5) +
    stat_compare_means(method = "kruskal.test", label.x = 2,
                       label.y = 1.05 * max(both_levels_cz$spruce_share * 100, na.rm = TRUE)) +
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


