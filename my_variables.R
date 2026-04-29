

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
               annot_x = 3.5, annot_y = NULL,
               ylim = NULL) {
  pr <- ggpredict(model, terms = terms, exclude = "s(plot_id)") %>%
    as.data.frame() %>%
    mutate(across(c(predicted, conf.low, conf.high), ~ . * scale_y))
  
  p <- ggplot(pr, aes(x = x, y = predicted, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.8) +
    labs(x = xlab, y = ylab, colour = "Scale") +
    scale_color_manual(values = c("subplot" = "grey80", "plot" = "grey30")) +
    scale_fill_manual(values = c("subplot" = "grey80", "plot" = "grey30"), guide = "none") +
    coord_cartesian(ylim = ylim, clip = "off") +   # clip off here, ylim passed in
    theme_classic() +
    theme(axis.title  = element_text(size = 7),
          plot.margin = margin(t = 25, r = 5, b = 5, l = 5))
  
  if (!is.null(annot))
    p <- p + annotate("text",
                      x     = annot_x,
                      y     = Inf,
                      label = annot,
                      size  = 2.5,
                      hjust = 0.5,
                      vjust = 1.1)
  p
}


# ── pp_inset_model: also add ylim parameter ───────────────────────────────────
pp_inset_model <- function(model, scale_y = 1, p_lab = NULL,
                           annot_x = 1.5, annot_y = NULL,
                           ylim = NULL) {
  pr <- ggpredict(model, terms = "level") %>%
    as.data.frame() %>%
    mutate(across(c(predicted, conf.low, conf.high), ~ . * scale_y))
  
  ggplot(pr, aes(x = x, y = predicted, colour = x)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .08) +
    geom_point(size = 2.4) +
    annotate("text",
             x     = annot_x,
             y     = Inf,
             label = p_lab,
             size  = 2.5,
             hjust = 0.5,
             vjust = 1.1) +
    scale_color_manual(values = c("subplot" = "grey80", "plot" = "grey30"),
                       guide = "none") +
    coord_cartesian(ylim = ylim, clip = "off") +
    theme_classic(base_size = 8) +
    theme(axis.title      = element_blank(),
          axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          axis.line.y     = element_blank(),
          axis.text.x     = element_text(size = 7),
          legend.position = "none",
          plot.margin     = margin(t = 25, r = 5, b = 5, l = 5))
}






# hold on all variables regarding species classifications
# update laurs's code for my species
# earlyspecs_laura <- c("lade","acps","frex","cabe","bepe","alin",
#                       "algl","acca","acpl","soau","soar",
#                       "coav","alvi","potr","poni","ulgl",
#                       "saca","rops")

# change the coding to match my data
earlyspecs_laura <- c("lade",
                      "acps",
                      "frex",
                      "cabe",
                      "besp", # from bepe
                      "alin",
                      "algl",
                      "acca",
                      "acpl",
                      "soau",
                      "soar",
                      "coav", # not present
                      "alvi",
                      "potr",
                      "poni", # not present
                      "ulsp",  # from ulgl
                      "saca",
                      "rops")


# from Claude: 
earlyspecs_cz <- c(
  # Pioneer conifers - light-demanding, open-ground colonizers
  "lade",   # Larix decidua
  "pisy",   # Pinus sylvestris  ← important addition for CZ bark beetle sites
  "psme",   # Pseudotsuga menziesii ← planted on salvage-logged sites
  
  # Classic pioneer broadleaves
  "besp",   # Betula spp.
  "potr",   # Populus tremula
  "posp",   # Populus spp.
  "saca",   # Salix caprea
  "sasp",   # Salix spp.
  "alin",   # Alnus incana
  "algl",   # Alnus glutinosa
  "alvi",   # Alnus viridis
  "soau",   # Sorbus aucuparia
  "soar",   # Sorbus aria
  
  # Secondary pioneers - fast on disturbed ground, moderate shade tolerance
  "acps",   # Acer pseudoplatanus
  "acpl",   # Acer platanoides
  "frex",   # Fraxinus excelsior
  "ulsp",   # Ulmus spp. - vigorous post-disturbance sprouting
  "rops",   # Robinia pseudoacacia - invasive but strongly early seral
  "aial"    # Ailanthus altissima - invasive, aggressive disturbed-ground colonizer
)



# > unique(species_class$species)
# [1] "absp" "piab" "pisy" "psme" "lade" "juni" "taba" "acpl" "acps" "algl" "alin"
# [12] "besp" "cabe" "fasy" "potr" "prav" "qusp" "saca" "sasp" "soar" "soau" "tisp"
# [23] "acca" "aehi" "aial" "alvi" "casa" "frex" "fror" "jure" "osca" "posp" "rops"
# [34] "soto" "ulsp"




# update species labels
species_labels <- c(
  piab = "Picea abies",
  besp = "Betula sp.",
  pisy = "Pinus sylvestris",
  qusp = "Quercus sp.",
  fasy = "Fagus sylvatica",
  lade = "Larix decidua",
  saca = "Salix caprea",
  soau = "Sorbus aucuparia",
  acps = "Acer pseudoplatanus",
  potr = "Populus tremula",
  absp = "Abies sp."#,
  #sasp = "Salix sp."
)

species_labels_other <- c(
  species_labels,          # your existing lookup from my_variables.R
  "other" = "Other species"
)

# add other to species_colors
species_colors_other <- c(
  species_colors,
  "other" = "#d9d9d9"
)

#species_levels <- rev(names(species_labels))  # Custom order, matching color palette and labels

# species_colors
# piab      besp      pisy      qusp      fasy      lade      saca      soau      acps      potr      absp      sasp 
#"#006837" "#17934D" "#58B65F" "#94D168" "#C6E77F" "#EDF7A7" "#FEF0A7" "#FDCD7B" "#FA9C58" "#EE613D" "#D22B26" "#A50026" 

species_colors <- c(
  piab = "#006837",
  besp = "#17934D",
  pisy = "#58B65F",
  qusp = "#94D168",
  fasy = "#C6E77F",
  lade = "#EDF7A7",
  saca = "#FEF0A7",
  soau = "#FDCD7B",
  acps = "#FA9C58",
  potr = "#EE613D",
  absp = "#D22B26",
  sasp = "#A50026"
)

# add other to species_colors
species_colors_other <- c(
  species_colors,
  "other" = "#d9d9d9"
)

# divid species on coniferous vs deciduousl

species_class <- tibble::tribble(
  ~species, ~leaf_type,
  "absp", "coniferous",
  "piab", "coniferous",
  "pisy", "coniferous",
  "psme", "coniferous",
  "lade", "coniferous",
  "juni", "coniferous",
  "taba", "coniferous",
  
  "acpl", "deciduous",
  "acps", "deciduous",
  "algl", "deciduous",
  "alin", "deciduous",
  "besp", "deciduous",
  "cabe", "deciduous",
  "fasy", "deciduous",
  "potr", "deciduous",
  "prav", "deciduous",
  "qusp", "deciduous",
  "saca", "deciduous",
  "sasp", "deciduous",
  "soar", "deciduous",
  "soau", "deciduous",
  "tisp", "deciduous",
  "acca", "deciduous",
  "aehi", "deciduous",
  "aial", "deciduous",
  "alvi", "deciduous",
  "casa", "deciduous",
  "frex", "deciduous",
  "fror", "deciduous",
  "jure", "deciduous",
  "osca", "deciduous",
  "posp", "deciduous",
  "rops", "deciduous",
  "soto", "deciduous",
  "ulsp", "deciduous"
)



# divide on functionnal groups
species_functional <- tribble(
  ~species, ~func_group,              ~rationale,
  # ── Spontaneous pioneers: low shade tolerance, rarely planted ──────────────
  "besp",   "spontaneous_pioneer",    "light-demanding, wind-dispersed, rarely planted",
  "potr",   "spontaneous_pioneer",    "clonal sprouter, not planted",
  "saca",   "spontaneous_pioneer",    "wind-dispersed, not planted",
  "sasp",   "spontaneous_pioneer",    "wind-dispersed, not planted",
  "soau",   "spontaneous_pioneer",    "bird-dispersed, not planted",
  "alin",   "spontaneous_pioneer",    "N-fixer, riparian, not planted",
  "algl",   "spontaneous_pioneer",    "N-fixer, not planted",
  "alvi",   "spontaneous_pioneer",    "N-fixer, not planted",
  "rops",   "spontaneous_pioneer",    "invasive, not planted in CZ policy",
  "aial",   "spontaneous_pioneer",    "invasive, not planted",
  
  # ── Economically valuable early-seral: planted AND naturally recruiting ────
  "pisy",   "economic_early",         "timber value, drought tolerant, planted + natural",
  "lade",   "economic_early",         "timber value, planted on salvage sites",
  "psme",   "economic_early",         "timber value, planted on salvage sites",
  "acps",   "economic_early",         "mixed forest target, moderate shade tolerance",
  
  # ── Climate-target late-seral: planted as policy species ──────────────────
  "fasy",   "planted_late_seral",     "shade tolerant, CZ subsidy target, low drought tol.",
  "qusp",   "planted_late_seral",     "drought tolerant, CZ subsidy target",
  "absp",   "planted_late_seral",     "shade tolerant, mixed forest target",
  "acpl",   "planted_late_seral",     "mixed forest target",
  "frex",   "planted_late_seral",     "mixed forest target, declining due to ash dieback",
  "ulsp",   "planted_late_seral",     "mixed forest target",
  
  # ── Legacy conifer: the pre-disturbance dominant ───────────────────────────
  "piab",   "legacy_conifer",         "pre-disturbance dominant, drought sensitive"
)



seral_colors <- c(
  "Early-dominated"  = "#4dac26",
  "Mixed"            = "#b8b8b8",
  "Late-dominated"   = "#d01c8b",
  # "Late seral only"  = "#e8c8e8",   # light purple — late but present
  "No regeneration"  = "#d9d9d9"    # grey — truly empty
)



# ── Shared theme ──────────────────────────────────────────────────────────────
theme_paper <- function(...) {
  theme_classic(base_size = 10) +
    theme(
      # axis
      axis.text        = element_text(size = 8),
      axis.title       = element_text(size = 8),
      # legend
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.key.width  = unit(1, "lines"),
      legend.key.height = unit(0.4, "lines"),
      # panel
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
      # default margin
      plot.margin      = margin(t = 5, r = 5, b = 5, l = 5),
      # pass any overrides
      ...
    )
}

