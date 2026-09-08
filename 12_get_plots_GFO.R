# ============================================================
# Final model plots — Czechia
# Lightweight script: loads exported final models + summaries only.
# No data prep, no model fitting, no diagnostics code needed here.
# ============================================================

gc()
rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(mgcv)
library(ggplot2)
library(svglite)
library(patchwork)
library(marginaleffects)   # for back-transforming coefficients to % change later

# ---- load exported model bundle ----

fin.models.all           <- readRDS("outModel/final_gam_models.rds")
models_with_cluster      <- readRDS("outModel/final_gam_models_with_cluster.rds")

final_coefs              <- read.csv("outModel/final_gam_coefficients.csv")
cluster_sensitivity_wide <- read.csv("outModel/cluster_sensitivity_wide.csv")

plot_coords              <- readRDS("outModel/plot_coords.rds")

# sanity check everything loaded as expected
names(fin.models.all)
names(models_with_cluster)
str(final_coefs)
str(cluster_sensitivity_wide)


# get raw data --------------------------------------------------------

ok, lets start again. show me, how to access raw data from teh model, paryicularly richness

Since fin.models.all$rich is the GAM fit on species richness, its own model frame has exactly what went into that fit — this is the same idea as before, just applied specifically to richness so you can look at it directly rather than trusting a plot of it.


# ---- pull the raw data straight from the fitted model ----
rich_data <- fin.models.all$rich$model

# what's actually in there
str(rich_data)
head(rich_data)

# how many rows per scale, and does that match what you expect from
# your original data prep (n_overlap_plots / n_overlap_subplots)?
table(rich_data$level)

# the shape of richness itself, split by scale -- this is the part
# that explains why the histogram/violin looked odd: a lot of subplot
# rows sitting at very low integer values
rich_data %>%
  dplyr::group_by(level) %>%
  dplyr::summarise(
    n        = dplyr::n(),
    min      = min(sp_richness),
    max      = max(sp_richness),
    mean     = round(mean(sp_richness), 2),
    median   = median(sp_richness),
    n_zero   = sum(sp_richness == 0),
    pct_zero = round(100 * mean(sp_richness == 0), 1)
  )

# the actual value-by-value counts per scale -- this is the most direct
# way to see the discreteness that a KDE was smoothing over
table(level = rich_data$level, richness = rich_data$sp_richness)

rich_data %>% 
  ggplot(aes(x = sp_richness, fill = level)) +
  geom_histogram() + 
  facet_grid(~level)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

# GAM time-trend prediction by scale (level), excluding the
# plot-level random effect (s(plot_id)), management held at 0
# (unmanaged baseline), year_f held at its most common level.
predict_time_trend <- function(model,
                               time_var  = "time_snc_full_disturbance",
                               level_var = "level",
                               time_seq  = NULL,
                               zero_covariates = c("planting_pred", "browsing_pred", "grndwrk_pred")) {
  
  mf <- model$model
  
  if (is.null(time_seq)) {
    time_seq <- seq(min(mf[[time_var]], na.rm = TRUE),
                    max(mf[[time_var]], na.rm = TRUE),
                    length.out = 50)
  }
  
  newdat <- expand.grid(
    time  = time_seq,
    level = levels(mf[[level_var]]),
    stringsAsFactors = FALSE
  )
  names(newdat) <- c(time_var, level_var)
  
  for (v in intersect(zero_covariates, names(mf))) newdat[[v]] <- 0
  
  if ("year_f" %in% names(mf)) {
    ref_year <- names(sort(table(mf$year_f), decreasing = TRUE))[1]
    newdat$year_f <- factor(ref_year, levels = levels(mf$year_f))
  }
  
  # predict.gam needs a valid factor level for plot_id even though
  # it's excluded below -- any observed level works, value doesn't matter
  if ("plot_id" %in% names(mf)) newdat$plot_id <- mf$plot_id[1]
  
  pr <- predict(model, newdata = newdat, type = "link", se.fit = TRUE,
                exclude = "s(plot_id)")
  
  linkinv <- family(model)$linkinv
  newdat %>%
    mutate(
      fit   = linkinv(pr$fit),
      lower = linkinv(pr$fit - 1.96 * pr$se.fit),
      upper = linkinv(pr$fit + 1.96 * pr$se.fit)
    )
}

# Scale marginal means, time held at its mean -- same logic, but a
# single value per level instead of a curve.
predict_scale_means <- function(model,
                                time_var  = "time_snc_full_disturbance",
                                level_var = "level",
                                zero_covariates = c("planting_pred", "browsing_pred", "grndwrk_pred")) {
  
  mf <- model$model
  newdat <- data.frame(level = levels(mf[[level_var]]), stringsAsFactors = FALSE)
  names(newdat) <- level_var
  newdat[[time_var]] <- mean(mf[[time_var]], na.rm = TRUE)
  
  for (v in intersect(zero_covariates, names(mf))) newdat[[v]] <- 0
  
  if ("year_f" %in% names(mf)) {
    ref_year <- names(sort(table(mf$year_f), decreasing = TRUE))[1]
    newdat$year_f <- factor(ref_year, levels = levels(mf$year_f))
  }
  if ("plot_id" %in% names(mf)) newdat$plot_id <- mf$plot_id[1]
  
  pr <- predict(model, newdata = newdat, type = "link", se.fit = TRUE,
                exclude = "s(plot_id)")
  
  linkinv <- family(model)$linkinv
  newdat %>%
    mutate(
      fit   = linkinv(pr$fit),
      lower = linkinv(pr$fit - 1.96 * pr$se.fit),
      upper = linkinv(pr$fit + 1.96 * pr$se.fit)
    )
}

# Raw response distribution by scale, pulled straight from the
# fitted model's own model frame (model$model) -- no re-fitting,
# no re-loading of both_levels_crossscale needed.
#
# type = "density"   -- smooth KDE per scale, mirrored/overlaid, response
#                        turned vertical (orientation = "y"); best for
#                        continuous responses (CV, effective species)
# type = "histogram" -- true binned counts, response vertical; better for
#                        coarse/discrete responses (e.g. species richness,
#                        small integer counts) where a KDE can look too smooth
# type = "violin"    -- the original violin + boxplot + jitter version, kept
#                        here in case you want to switch back for any panel
#
# NOTE: bounds = c(0, Inf) on the density needs ggplot2 >= 3.4.0 (keeps the
# KDE from drawing a tail below zero for these non-negative responses).
plot_scale_distribution <- function(model, response_var, ylab, level_var = "level",
                                    type = c("density", "histogram", "violin"),
                                    bins = 20) {
  type <- match.arg(type)
  dat  <- model$model
  
  if (type == "violin") {
    p <- dat %>%
      ggplot(aes(x = .data[[level_var]], y = .data[[response_var]], fill = .data[[level_var]])) +
      geom_violin(alpha = 0.4, trim = TRUE, color = NA) +
      geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
      geom_jitter(width = 0.05, alpha = 0.3, size = 0.8, color = "grey30") +
      scale_fill_manual(values = scale_colors, guide = "none") +
      labs(x = "Scale", y = ylab)
  } else {
    p <- dat %>%
      ggplot(aes(y = .data[[response_var]], fill = .data[[level_var]], color = .data[[level_var]]))
    
    if (type == "density") {
      p <- p +
        geom_density(orientation = "y", alpha = 0.4, linewidth = 0.6, bounds = c(0, Inf)) +
        labs(x = "Density")
    } else {
      p <- p +
        geom_histogram(orientation = "y", position = "identity",
                       bins = bins, alpha = 0.5, color = NA) +
        labs(x = "Count")
    }
    
    p <- p +
      scale_fill_manual(values = scale_colors, name = "Scale") +
      scale_color_manual(values = scale_colors, name = "Scale") +
      labs(y = ylab)
  }
  
  # legend suppressed here -- panel [a] already carries the one "Scale"
  # legend for the whole figure; color mapping is identical (scale_colors)
  p +
    theme_classic(base_size = 9) +
    theme(axis.title = element_text(size = 8), legend.position = "none")
}

scale_colors <- c(subplot = "grey60", plot = "black")

# ------------------------------------------------------------
# Panel A: height -- time trend by scale + scale inset
# ------------------------------------------------------------

# compute both predictions once, so a and b can share one y-axis window
hgt_time_pred  <- predict_time_trend(fin.models.all$hgt)
hgt_scale_pred <- predict_scale_means(fin.models.all$hgt)

y_range_hgt <- range(c(hgt_time_pred$lower,  hgt_time_pred$upper,
                       hgt_scale_pred$lower, hgt_scale_pred$upper))

p_time_hgt <- hgt_time_pred %>%
  ggplot(aes(x = time_snc_full_disturbance, y = fit, color = level, fill = level)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = scale_colors, name = "Scale") +
  scale_fill_manual(values = scale_colors, name = "Scale") +
  coord_cartesian(ylim = y_range_hgt) +
  labs(x = "Time since disturbance\n(years)", y = "Mean height [m]") +
  theme_classic(base_size = 9) +
  theme(
    legend.position      = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.title         = element_text(size = 8),
    legend.text          = element_text(size = 8),
    legend.key.size      = unit(0.7, "lines")
  )

p_scale_hgt <- hgt_scale_pred %>%
  ggplot(aes(x = level, y = fit, color = level)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), linewidth = 0.8, size = 0.6) +
  scale_color_manual(values = scale_colors, guide = "none") +
  coord_cartesian(ylim = y_range_hgt) +
  labs(x = "Scale", y = NULL) +
  theme_classic(base_size = 9)

p_hgt_panel <- p_time_hgt + p_scale_hgt +
  plot_layout(widths = c(2.4, 1))

# ------------------------------------------------------------
# Panels C-E: CV, effective species, richness -- built twice, so both
# alternatives are kept side by side rather than picking one upfront.
#
# "dens" set: vertical density (CV, effective species) / histogram
#             (richness -- small-integer count, KDE over-smooths it)
# "vio"  set: the original violin + boxplot + jitter version
# ------------------------------------------------------------

p_cv_dens   <- plot_scale_distribution(fin.models.all$cvpos, "cv_hgt_pos",        "CV [%]",
                                       type = "histogram")
p_eff_dens  <- plot_scale_distribution(fin.models.all$eff,   "effective_numbers", "Effective species [#]",
                                       type = "histogram")
p_rich_dens <- plot_scale_distribution(fin.models.all$rich,  "sp_richness",       "Species richness [#]",
                                       type = "histogram", bins = 12)

p_cv_vio   <- plot_scale_distribution(fin.models.all$cvpos, "cv_hgt_pos",        "CV [%]",
                                      type = "violin")
p_eff_vio  <- plot_scale_distribution(fin.models.all$eff,   "effective_numbers", "Effective species [#]",
                                      type = "violin")
p_rich_vio <- plot_scale_distribution(fin.models.all$rich,  "sp_richness",       "Species richness [#]",
                                      type = "violin")

# ------------------------------------------------------------
# Combine into two 2x2 figures -- same panel [a]/[b] (height) in both,
# only c-e differ
# ------------------------------------------------------------

# 5 tags, not 4: p_hgt_panel is itself a 2-plot composite (time trend +
# scale inset), so it consumes the first two tag slots ([a], [b]) before
# c/d/e take the rest
fig3_tags <- list(c("[a]", "[b]", "[c]", "[d]", "[e]"))

p_fig3_density <-
  (p_hgt_panel | p_cv_dens) /
  (p_eff_dens  | p_rich_dens) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = fig3_tags, title = "Density / histogram") &
  theme(plot.tag = element_text(size = 10, face = "plain"))

p_fig3_violin <-
  (p_hgt_panel | p_cv_vio) /
  (p_eff_vio   | p_rich_vio) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = fig3_tags, title = "Violin + boxplot") &
  theme(plot.tag = element_text(size = 10, face = "plain"))

p_fig3_density
p_fig3_violin

# stacked, for a direct side-by-side comparison in one window
p_fig3_density / p_fig3_violin

# ------------------------------------------------------------
# Export -- both versions
# ------------------------------------------------------------

ggsave("outFigsCZ/p_fig3_density.png", p_fig3_density,
       width = 7, height = 6, dpi = 300, bg = "white")
ggsave("outFigsCZ/p_fig3_density.svg", p_fig3_density,
       width = 7, height = 6, bg = "white")

ggsave("outFigsCZ/p_fig3_violin.png", p_fig3_violin,
       width = 7, height = 6, dpi = 300, bg = "white")
ggsave("outFigsCZ/p_fig3_violin.svg", p_fig3_violin,
       width = 7, height = 6, bg = "white")

# ------------------------------------------------------------
# NOTES / things to check on your machine:
# - Column names (mean_hgt, cv_hgt_pos, effective_numbers, sp_richness,
#   time_snc_full_disturbance, level, plot_id, year_f) are taken from
#   your original model formulas -- double check they match if any
#   model was refit with renamed variables.
# - Management covariates are zeroed to represent the unmanaged
#   baseline; drop that block if your models don't include them.
# - year_f is fixed at its most frequent observed level as the
#   reference year; change `ref_year` logic if you'd rather average
#   predictions across years or use a specific one.
# - exclude = "s(plot_id)" assumes the random effect was specified
#   as s(plot_id, bs = "re") -- matches your original model code.
# ------------------------------------------------------------


# ============================================================
# "Results - key figures" panel, rebuilt with real model output
# ============================================================
# Recreates the 4-panel poster/talk summary (height trend, diversity,
# compositional dissimilarity, spruce share) in the same bar-chart
# style as the mockup, but:
#   - pulls real numbers + CIs from fin.models.all (mockup had none)
#   - panel 2: shows actual predicted values instead of an index, so
#     the axis label and the numbers on the bars finally agree
#     (mockup had y-axis "Change (%)" 0-80 next to bars labeled
#     1.00 / 1.59, which are an index, not a percent)
#   - panel 3: relabeled Jaccard dissimilarity, not Bray-Curtis --
#     matches beta_jaccard_mean / your Methods section
#
# Per the project notes, this is the quick "same bar style, real
# data" fix -- not the fuller redesign (continuous planting_intensity
# + GAM ribbon, paradox panels on a shared x-axis) that was also
# scoped out for the GfO talk deck, if you want to revisit that later.
# ============================================================

library(dplyr)
library(mgcv)
library(emmeans)
library(ggplot2)
library(patchwork)

# ---- load exported model bundle + coefficient table ----
fin.models.all <- readRDS("outModel/final_gam_models.rds")
final_coefs     <- read.csv("outModel/final_gam_coefficients.csv")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

# Unplanted (0) vs planted (1) marginal means + 95% CI for one model.
# term = "planting_pred" for the cross-scale models (hgt/cvpos/eff/rich),
# term = "planting_intensity" for the plot-only models (beta/adapt/spruce).
# Marginalizes over any other factor in the model (e.g. level, year_f)
# with equal weights, same as your original get_mng_emm() pipeline.
#
# NOTE: emmeans names the backtransformed estimate column differently
# depending on the model's link -- "response" for a log/logit-link
# family (tw/nb/betar: your eff/rich/spruce models), but "emmean" for
# an identity-link family (gaussian: your beta dissimilarity model,
# since there's nothing to backtransform). Only ONE of the two columns
# actually exists in any given result, so coalesce(response, emmean)
# errors the moment it references the missing one -- detect which
# column is present instead of assuming both are there.
get_planting_effect <- function(model, term = c("planting_pred", "planting_intensity")) {
  term <- match.arg(term)
  at_list <- setNames(list(c(0, 1)), term)
  
  df <- emmeans(model, specs = as.formula(paste("~", term)), at = at_list, type = "response") %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    rename(planting = all_of(term))
  
  est_col <- intersect(c("response", "emmean"), names(df))[1]
  
  df %>%
    mutate(
      status   = factor(ifelse(planting == 0, "Unplanted", "Planted"),
                        levels = c("Unplanted", "Planted")),
      estimate = .data[[est_col]]
    )
}

# p-value for the focal planting term, pulled from your already-exported
# parametric coefficient table (broom::tidy(model, parametric = TRUE))
get_planting_pvalue <- function(model_name, term = c("planting_pred", "planting_intensity")) {
  term <- match.arg(term)
  p <- final_coefs %>%
    filter(model == model_name, term == !!term) %>%
    pull(p.value)
  if (length(p) == 0) return(NA_real_)
  p
}

p_label <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) "p < 0.001" else paste0("p = ", formatC(p, format = "f", digits = 3))
}

bar_colors <- c(Unplanted = "grey70", Planted = "#2e7d32")

# One unplanted-vs-planted bar panel with CI error bars + value labels.
# `title` is the descriptive header text (poster-style, e.g. "Planting
# reduces compositional dissimilarity") -- printed above the p-value,
# in place of the placeholder numeric tags.
plot_planting_bar <- function(model, model_name, term, ylab, title = NULL,
                              digits = 2, fill_override = NULL, pct_scale = FALSE) {
  
  dat <- get_planting_effect(model, term)
  if (pct_scale) dat <- dat %>% mutate(estimate = estimate * 100, lower.CL = lower.CL * 100, upper.CL = upper.CL * 100)
  
  cols <- if (is.null(fill_override)) bar_colors else fill_override
  p_val <- get_planting_pvalue(model_name, term)
  
  ggplot(dat, aes(x = status, y = estimate, fill = status)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.12, linewidth = 0.5) +
    geom_text(aes(label = round(estimate, digits), y = upper.CL),
              vjust = -0.6, fontface = "bold", size = 3.2) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = "Management", y = ylab, title = title,
         subtitle = p_label(p_val)) +
    theme_classic(base_size = 9) +
    theme(
      plot.title    = element_text(size = 9, face = "bold", hjust = 0.5,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(hjust = 0.5, face = "bold",
                                   color = if (is.null(fill_override)) "grey20" else fill_override[["Planted"]]),
      axis.title.y  = element_text(size = 8)
    )
}

# ------------------------------------------------------------
# Panel 1: height vs time (reference: plot scale, unplanted baseline --
# change level_ref / planting_ref below if you'd rather show a
# different reference, e.g. averaged across scale)
# ------------------------------------------------------------

level_ref    <- "plot"
planting_ref <- 0

hgt_model <- fin.models.all$hgt
mf_hgt    <- hgt_model$model

time_seq <- seq(min(mf_hgt$time_snc_full_disturbance, na.rm = TRUE),
                max(mf_hgt$time_snc_full_disturbance, na.rm = TRUE),
                length.out = 50)

newdat_hgt <- data.frame(
  time_snc_full_disturbance = time_seq,
  level         = factor(level_ref, levels = levels(mf_hgt$level)),
  planting_pred = planting_ref,
  browsing_pred = 0,
  grndwrk_pred  = 0,
  year_f        = factor(names(sort(table(mf_hgt$year_f), decreasing = TRUE))[1],
                         levels = levels(mf_hgt$year_f)),
  plot_id       = mf_hgt$plot_id[1]
)

pr_hgt <- predict(hgt_model, newdata = newdat_hgt, type = "link", se.fit = TRUE,
                  exclude = "s(plot_id)")
linkinv_hgt <- family(hgt_model)$linkinv
newdat_hgt <- newdat_hgt %>%
  mutate(
    fit   = linkinv_hgt(pr_hgt$fit),
    lower = linkinv_hgt(pr_hgt$fit - 1.96 * pr_hgt$se.fit),
    upper = linkinv_hgt(pr_hgt$fit + 1.96 * pr_hgt$se.fit)
  )

# smooth-term p-value (not in final_coefs, which is parametric terms only)
p_time_hgt_val <- summary(hgt_model)$s.table["s(time_snc_full_disturbance)", "p-value"]

p_panel1 <- ggplot() +
  geom_point(data = mf_hgt %>% filter(level == level_ref),
             aes(x = time_snc_full_disturbance, y = mean_hgt),
             alpha = 0.25, size = 1, color = "grey40") +
  geom_ribbon(data = newdat_hgt, aes(x = time_snc_full_disturbance, ymin = lower, ymax = upper),
              fill = "#2e7d32", alpha = 0.2) +
  geom_line(data = newdat_hgt, aes(x = time_snc_full_disturbance, y = fit),
            color = "#2e7d32", linewidth = 1) +
  labs(x = "Time since disturbance (years)", y = "Mean tree height (m)",
       title = "Mean tree height",
       subtitle = p_label(p_time_hgt_val)) +
  theme_classic(base_size = 9) +
  theme(
    plot.title    = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    plot.subtitle = element_text(hjust = 0.5, face = "bold", color = "#2e7d32")
  )

# ------------------------------------------------------------
# Panel 2: diversity -- species richness only (effective species
# dropped per your request; swap fin.models.all$eff back in with
# its own plot_planting_bar() call if you want it back)
# ------------------------------------------------------------

p_panel2 <- plot_planting_bar(fin.models.all$rich, "rich", "planting_pred",
                              ylab  = "Species richness [#]",
                              title = "Species richness")

# ------------------------------------------------------------
# Panel 3: compositional dissimilarity -- Jaccard, not Bray-Curtis
# ------------------------------------------------------------

p_panel3 <- plot_planting_bar(fin.models.all$beta, "beta", "planting_intensity",
                              ylab  = "Jaccard dissimilarity",
                              title = "Jaccard dissimilarity",
                              digits = 2)

# ------------------------------------------------------------
# Panel 4: Norway spruce share -- the "paradox" panel, kept in the
# orange/red accent to flag it as the counterpoint to panel 2
# ------------------------------------------------------------

paradox_colors <- c(Unplanted = "grey70", Planted = "#d84315")

p_panel4 <- plot_planting_bar(fin.models.all$spruce, "spruce", "planting_intensity",
                              ylab  = "Norway spruce share [%]",
                              title = "Norway spruce share",
                              digits = 1, fill_override = paradox_colors, pct_scale = TRUE)

# ------------------------------------------------------------
# Combine into one row, matching the poster layout
# ------------------------------------------------------------

p_key_figures <-
  p_panel1 + p_panel2 + p_panel3 + p_panel4 +
  plot_layout(nrow = 1, widths = c(1.3, 1, 1, 1)) #+
  #plot_annotation(title = "")

p_key_figures

ggsave("outFigsCZ/p_key_figures_real.png", p_key_figures,
       width = 10, height = 3.2, dpi = 300, bg = "white")
ggsave("outFigsCZ/p_key_figures_real.svg", p_key_figures,
       width = 10, height = 3.2, bg = "white")

# editable vector PDF for Corel -- cairo_pdf keeps text as real,
# selectable/editable text objects (and paths as vector paths) rather
# than rasterizing, same device used elsewhere in your pipeline
# (e.g. p_combined_function.pdf)
ggsave("outFigsCZ/p_key_figures_real.pdf", p_key_figures,
       width = 10, height = 3.2, device = cairo_pdf)

# ------------------------------------------------------------
# NOTES:
# - Panel 1 is drawn at level = "plot", planting = 0 (unplanted) as
#   the reference. If your MS reports this trend marginalized across
#   scale instead, that's a bigger change (predict at both levels and
#   average the link-scale predictions before back-transforming --
#   not just averaging the two response-scale curves) -- flag if you
#   want that version instead.
# - get_planting_effect() marginalizes over every other factor (level,
#   year_f) with equal weights via emmeans' default reference grid --
#   same behavior as your original get_mng_emm() / emm_table pipeline,
#   so these numbers should match what's already in
#   outTable/management_effects_emmeans.doc if you want to cross-check.
# - This script only rebuilds the 4 charts, not the poster's dark
#   green header banner / caption text styling -- that's a slide/poster
#   layout layer, better done in your actual pptx build
#   (build_deck.js) or poster software than reproduced in ggplot.
# - Values plotted are the true unit (species count / Jaccard
#   dissimilarity / % spruce), not a fold-change index -- this was a
#   deliberate fix for the axis-label mismatch flagged in the poster
#   critique (y-axis said "Change (%)" 0-80 next to bars labeled as
#   an index, 1.00 / 1.59). Say the word if you'd rather keep the
#   "% change from unplanted" framing and I'll rebuild panel 2 that way
#   with a corrected, consistent axis label instead.
# ------------------------------------------------------------