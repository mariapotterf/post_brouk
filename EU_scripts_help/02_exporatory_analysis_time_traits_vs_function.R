# =============================================================================
# Post-disturbance functional community reorganization in Central European forests
# Analysis script
# =============================================================================
# Structure:
#   0. Setup & data loading
#   1. Data preparation
#   2. Exploratory data analysis (EDA)
#   3. Functional trait state space (CWM_shade × CWM_drought)
#   4. LMMs: drivers of functional community state (CWM)
#   5. Conifer/deciduous ratio trajectory
#   6. Diversity–environment relationships
#   7. Management × disturbance interactions
#   8. Export results & model summaries
# =============================================================================


# ── 0. Setup ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)      # p-values for lmer
library(MuMIn)         # R² for mixed models (r.squaredGLMM)
library(emmeans)       # marginal means & contrasts
library(ggplot2)
library(ggeffects)     # model-predicted effect plots
library(patchwork)     # combine ggplots
library(vegan)         # permanova / ordination if needed
library(performance)   # check_model diagnostics
library(broom.mixed)   # tidy model output


# ── 1. Data loading & preparation ─────────────────────────────────────────────

dat <- fread("EU_scripts_help/both_levels_EU_comb.csv") %>% as_tibble()

# split levels
sub  <- dat %>% filter(level == "subplot")
plot <- dat %>% filter(level == "plot")

# ── 1a. Derived & cleaned variables ───────────────────────────────────────────

sub <- sub %>%
  mutate(
    # log-transform density (heavy right skew)
    log_dens_ha      = log1p(dens_ha),
    
    # management composite: mean of all intensity scores
    # (higher = more intensively managed post-disturbance)
    mgmt_composite   = rowMeans(
      select(., clear_intensity, grndwrk_intensity,
             logging_trail_intensity, planting_intensity,
             anti_browsing_intensity),
      na.rm = TRUE),
    
    # planted vs. not: planting_intensity > 0
    planted          = if_else(planting_intensity > 0, "planted", "unplanted"),
    
    # functional group dominance
    dom_functional   = case_when(
      IV_coniferous > 0.6  ~ "conifer_dom",
      IV_deciduous  > 0.6  ~ "deciduous_dom",
      TRUE                  ~ "mixed"
    ),
    
    # scale continuous predictors (z-scores) for model comparability
    z_severity   = scale(disturbance_severity)[,1],
    z_elevation  = scale(elevation)[,1],
    z_time       = scale(time_since)[,1],
    z_edge       = scale(dist_to_edge)[,1],
    z_mgmt       = scale(mgmt_composite)[,1],
    
    # factors
    country_name = factor(country_name),
    region       = factor(region),
    plot_id      = factor(plot_id),
    year_f       = factor(year_f)
  )

plot <- plot %>%
  mutate(
    log_dens_ha    = log1p(dens_ha),
    mgmt_composite = rowMeans(
      select(., clear_intensity, grndwrk_intensity,
             logging_trail_intensity, planting_intensity,
             anti_browsing_intensity),
      na.rm = TRUE),
    planted        = if_else(planting_intensity > 0, "planted", "unplanted"),
    z_severity     = scale(disturbance_severity)[,1],
    z_elevation    = scale(elevation)[,1],
    z_time         = scale(time_since)[,1],
    z_edge         = scale(dist_to_edge)[,1],
    z_mgmt         = scale(mgmt_composite)[,1],
    country_name   = factor(country_name),
    region         = factor(region),
    plot_id        = factor(plot_id),
    year_f         = factor(year_f)
  )

# working dataset: subplots with CWM available
sub_cwm <- sub %>% filter(!is.na(CWM_shade), !is.na(CWM_drought))

cat("Subplots total:", nrow(sub), "\n")
cat("Subplots with CWM:", nrow(sub_cwm), "\n")
cat("Plots total:", nrow(plot), "\n")
cat("Countries:", n_distinct(sub$country_name), "\n")

plot %>% 
  ggplot(aes(x = CWM_shade,
             y = CWM_drought,
             color = country_name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(.~country_name)


plot %>% 
  ggplot(aes(x = time_since,
             y = CWM_drought,
             color = country_name)) + 
  geom_point(alpha = 0.1) + 
  facet_wrap(.~country_name) + 
  geom_smooth(method = 'lm')


plot %>% 
  ggplot(aes(x = time_since,
             y = CWM_shade,
             color = country_name)) + 
  geom_point(alpha = 0.1) + 
  facet_wrap(.~country_name) + 
  geom_smooth(method = 'lm')

