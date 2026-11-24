# ============================================================
# Forest regeneration: temporal change analysis 2023-2025
# ============================================================


# includes data from EU countries (2023) and for Czechia (2023&2025)

# Sections:
#   0. Setup
#   1. Data ingestion & QC
#   2. Management data cleaning, remove matue trees from CZ
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


# 1. Data ingestion & QC -----------------------------------------
# data is on subplot level, on individual species and 3 vertical layers
dat_overlap <- fread('outData/full_table_23_25.csv') # contacis data acros 10 Eu countries, CZ is collected 2x times

# mature trees were already removed from2023 estimation - can be useful for some other metrics? look them up again just to keep the databases
dat_counts_mature_species23 <- fread("C:/Users/potterf/OneDrive - CZU v Praze/Dokumenty/2023_PanEuropean/r_paneurop/outData/sub_plots_counts_DBH.csv") |> 
  dplyr::filter(VegType == "Mature" & count > 0)

# and export from 2025 data 
dat_counts_mature_species25 <- dat_overlap |> 
  dplyr::filter(vegtype == "mature" & n > 0)


# get a mastertable about plots and subplots for 23 data - keep all teh empty plots and subplots
dat_master23 <- fread("C:/Users/potterf/OneDrive - CZU v Praze/Dokumenty/2023_PanEuropean/r_paneurop/outData/sub_plots_counts_DBH.csv") |> 
  dplyr::select(ID, cluster, country) |> 
  group_by(cluster) |> 
  mutate(n_subplots = n_distinct(ID)) |> 
  distinct()

# 180/4221 - 4.26% of subplots have mature trees present

hist(dat_overlap$x)

dat_master_subplot <- dat_overlap %>%
  dplyr::select(plot, subplot, year, status, x, y) %>%
  distinct()

table(dat_master_subplot$year)
table(dat_master_subplot$status)

n_plots_total    <- length(unique(dat_master_subplot$plot))
n_subplots_total <- length(unique(dat_master_subplot$subplot))


n_plots_total
n_subplots_total


cat("Plots:", n_plots_total, "| Subplots:", n_subplots_total, "\n")


# Get subplot& plot xy
dat_sub_xy <- dat_master_subplot %>%
  dplyr::select(subplot, year, x, y) %>%
  dplyr::distinct()


dat_plot_xy <- dat_master_subplot %>%
  dplyr::group_by(plot, year) %>%
  dplyr::summarise(
    x = mean(x, na.rm = TRUE),
    y = mean(y, na.rm = TRUE),
    .groups = "drop"
  )


# 2. Management data cleaning ------------------------------------
# define variables for cumulative management
# idea - if teh site jhas salvage, site preparation in 2023, 
# needs to have it also in 2025. if planting was not present in 2025,
# it can occur in 2025
mgmt_vars_cum <- c("clear", "grndwrk", "planting", "anti_browsing")

# define variables on subplot level
management_vars <- c("clear", "grndwrk", "logging_trail", "planting", "anti_browsing")

# define cols to remove from overall regeneration table
management_vars_rmv <- c(
  "salvage_sum", "protection_sum", "clear_sum", "grndwrk_sum",
  "logging_trail_sum", "planting_sum", "anti_browsing_sum",
  "management_sum", "clear_intensity", "grndwrk_intensity",
  "logging_trail_intensity", "planting_intensity",
  "anti_browsing_intensity", "salvage_intensity",
  "protection_intensity", "management_intensity"
)

dat_overlap_mng <- dat_overlap %>%
  dplyr::select(-all_of(management_vars_rmv)) %>% 
  # if managemnet is NA -> convert to 0
  mutate(across(all_of(management_vars), ~ replace_na(.x, 0)))

## Identify plots surveyed in both years - needs to have partly consistent management; 
# ON SUBPLOT LEVEL
# if present in 2023 -> needs to be present in 2025;
# if not present in 2023 - > it can be recorded in 2025
# hence, I am modifying management in 2025 to be consistent with 2023
# if already recorded in 2023; some management can be indipendent (browing, planting) and can also occur later
# another problem - management should be consistent on plot level, 
# but i do not know how it fits on subplot level (which subplot belongs to which)
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
  arrange(plot, 
          desc(clear), 
          desc(grndwrk), 
          desc(planting),    # then fixed priority per activity
          desc(anti_browsing), 
          desc(logging_trail),
          subplot) %>%
  mutate(row_id = row_number())

## 2025 management — resolve duplicates by taking max presence
mng_both25 <- dat_overlap_mng %>%
  filter(plot %in% plots_with_both, status == "both", year == 2025) %>%
  dplyr::select(plot, subplot, all_of(management_vars), year) %>%
  group_by(plot, subplot, year) %>%
  summarise(across(all_of(management_vars), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>%
  group_by(plot) %>%
  arrange(plot, 
          desc(clear), 
          desc(grndwrk), 
          desc(planting),    # then fixed priority per activity
          desc(anti_browsing), 
          desc(logging_trail),
          subplot) %>%
  mutate(row_id = row_number())


### Get management separately: per PLOT per year — long format ────────────────────────────────
# to keep 
mng_by_year_plot <- dat_overlap %>%
  # if there is NA in management - > set to 0
  mutate(across(all_of(management_vars), ~ replace_na(.x, 0))) %>% 
  filter(status == "both", year %in% c(2023, 2025)) %>%
  distinct(plot, subplot, year,
           clear, grndwrk, logging_trail, planting, anti_browsing) %>%
  group_by(plot, subplot, year) %>%
  summarise(across(c(clear, grndwrk, logging_trail, planting, anti_browsing),
                   ~ max(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  # plot-level intensity per year
  group_by(plot, year) %>%
  summarise(
    clear_intensity          = mean(clear,          na.rm = TRUE),
    grndwrk_intensity        = mean(grndwrk,        na.rm = TRUE),
    logging_trail_intensity  = mean(logging_trail,  na.rm = TRUE),
    planting_intensity       = mean(planting,       na.rm = TRUE),
    anti_browsing_intensity  = mean(anti_browsing,  na.rm = TRUE),
    .groups = "drop"
  )

## ── Wide format: one row per plot, columns per year ───────────────────────────
mng_wide <- mng_by_year_plot %>%
  pivot_wider(
    names_from  = year,
    values_from = c(clear_intensity, grndwrk_intensity, logging_trail_intensity,
                    planting_intensity, anti_browsing_intensity),
    names_sep   = "_"
  ) %>%
  # delta: did management intensity change?
  mutate(
    delta_planting      = planting_intensity_2025       - planting_intensity_2023,
    delta_browsing      = anti_browsing_intensity_2025  - anti_browsing_intensity_2023,
    delta_grndwrk       = grndwrk_intensity_2025        - grndwrk_intensity_2023,
    delta_clear         = clear_intensity_2025          - clear_intensity_2023,
    delta_logging_trail = logging_trail_intensity_2025  - logging_trail_intensity_2023
  )

hist(mng_wide$delta_planting)
hist(mng_wide$delta_browsing)
hist(mng_wide$delta_grndwrk)
hist(mng_wide$delta_clear)
hist(mng_wide$delta_logging_trail)


## Combine by row order, take max per management variable - this is because I do not know which of the subplots 
# actually overlaps; i know only plots
# subplots are now ordered by management applied
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
 # dplyr::filter(vegtype != "mature") |>    # remove again mature tree layer to keep calculations on regeneration layer
  left_join(management_intensity, by = c("plot", "year", "n_subplots"))


# add indication of country and region
dat_overlap_mng_upd_cntries <- dat_overlap_mng_upd2 %>% 
  mutate(
  region = str_extract(plot, "^\\d+(?=_)"),
  country = case_when(
    region %in% c("11", "12", "14", "18", "19", "20", "25") ~ 11L,
    region == "17"                                           ~ 12L,
    region %in% c("15", "26")                               ~ 13L,
    region == "13"                                           ~ 14L,
    region == "16"                                           ~ 15L,
    region == "23"                                           ~ 16L,
    region == "21"                                           ~ 17L,
    region == "22"                                           ~ 18L,
    region %in% c("24", "27")                               ~ 19L,
    is.na(region)                                           ~ 13L,  # bare IDs -> Czechia
    TRUE                                                     ~ NA_integer_
  ),
  country_name = case_when(
    country == 11 ~ "Germany",
    country == 12 ~ "Poland",
    country == 13 ~ "Czechia",
    country == 14 ~ "Austria",
    country == 15 ~ "Slovakia",
    country == 16 ~ "Slovenia",
    country == 17 ~ "Italy",
    country == 18 ~ "Switzerland",
    country == 19 ~ "France",
    TRUE          ~ NA_character_
  )
) %>% 
  mutate(country_name   = factor(country_name),
         region = factor(region),
         country_name = factor(country_name))


dat_overlap_mng_upd2

# filter only regeneration data
dat_regen <- dat_overlap_mng_upd_cntries |> 
 dplyr::filter(vegtype != "mature")   # remove again mature tree layer to keep calculations on regeneration layer


### Subplot and plot management tables ----------------------------------------------
# subplot level - keep both rwa values and plot level merges
df_mng_sub <- dat_regen %>%
  distinct(plot, subplot, year,
           clear, 
           grndwrk, 
           logging_trail, 
           planting, 
           anti_browsing,
           clear_intensity, 
           grndwrk_intensity, 
           logging_trail_intensity,
           planting_intensity, 
           anti_browsing_intensity) %>% 
  mutate(across(all_of(management_vars), ~ replace_na(.x, 0)))

# summarize on plot level
df_mng_plot <- df_mng_sub %>%
  select(plot, year,
         clear_intensity, 
         grndwrk_intensity, 
         logging_trail_intensity,
         planting_intensity, 
         anti_browsing_intensity) %>%
  distinct() 


#####  how many stems i have?  ---------------------------------------------

# stems per year in Czechia (small + advanced regeneration; mature layer already removed)
n_stems_regen_czechia <- dat_regen %>%
  dplyr::filter(country_name == "Czechia") %>%
  group_by(country_name) %>%
  summarise(n_stems = sum(n, na.rm = TRUE), .groups = "drop")
n_stems_regen_czechia 
# 5074 -> 5090 after including 26_134


# including mature trees
n_stems_all_layers_czechia <- dat_overlap_mng_upd_cntries %>%
  dplyr::filter(country_name == "Czechia") %>%
  group_by(country_name) %>%
  summarise(n_stems = sum(n, na.rm = TRUE), .groups = "drop")
n_stems_all_layers_czechia 
# 5088 - difference is 29 (!) trees 


### identify cumulative management type for Karim - to correspond to Johanness  --------------------------------------------
# Subplot-level cumulative management classification — three schemes
# for EU level estimation !!!!
# 
##### Applied directly at the SUBPLOT level (no plot aggregation). -------------------
#
# Input : df_mng_sub  — one row per subplot, with 0/1 columns:
#           clear, grndwrk, logging_trail, planting, anti_browsing
#         (planting/anti_browsing already recoded to 0/1; no `== 2` artefact)
#
# Output: df_mng_sub_cum — same rows + three parallel classifications:
#
#   A) mgmt_type        — "ours": strict cumulative LADDER, gated on `clear`.
#                         passive < salvaged < planted < protected; off-ladder
#                         combos -> "uncommon" (9).
#   B) mgmt_type_joh    — Johannes: PRECEDENCE, highest-activity-wins.
#                         sitePrep = clear|grndwrk|logging_trail; then
#                         pmax(sitePrep, planting*2, anti_browsing*3).
#                         Never yields an "uncommon" bin.
#   C) establishment + protect — TWO-AXIS: an ordinal establishment pathway
#                         (passive < sitePrep < planted) PLUS an independent
#                         binary protection flag. Lets protection be estimated
#                         separately from planting; empties the "uncommon" bin.
#
# NA rule (uniform across all three schemes, for a fair comparison): if ANY of
# the five management variables is NA for a subplot, all three classifications
# are set to NA for that subplot. (You reported no NAs, so this is a no-op here,
# but keeps the schemes aligned if that ever changes.)
# ==============================================================================

mng_vars <- c("clear", "grndwrk", "logging_trail", "planting", "anti_browsing")

df_mng_sub_cum <- df_mng_sub %>%
  mutate(
    # convenience aliases (match earlier naming)
    plant  = planting,
    browse = anti_browsing,
    
    # shared flags
    any_na    = if_any(all_of(mng_vars), is.na),
    sitePrep  = clear == 1 | grndwrk == 1 | logging_trail == 1,  # Johannes / two-axis
    planted   = plant  == 1,
    protected = browse == 1,
    
    # ---- A) OURS: strict cumulative ladder (gated on `clear`) ----------------
    mgmt_type = case_when(
      any_na                                  ~ NA_integer_,
      clear == 0 & plant == 0 & browse == 0   ~ 0L,   # passive
      clear == 1 & plant == 0 & browse == 0   ~ 1L,   # salvaged
      clear == 1 & plant == 1 & browse == 0   ~ 2L,   # planted
      clear == 1 & plant == 1 & browse == 1   ~ 3L,   # protected
      TRUE                                    ~ 9L    # uncommon combos
    ),
    mgmt_type_label = factor(mgmt_type,
                             levels = c(0:3, 9),
                             labels = c("passive", "salvaged", "planted", "protected", "uncommon")),
    
    # ---- B) JOHANNES: precedence / highest-activity-wins ---------------------
    mgmt_type_joh = if_else(
      any_na,
      NA_integer_,
      pmax(as.integer(sitePrep),
           as.integer(planted)   * 2L,
           as.integer(protected) * 3L)
    ),
    mgmt_type_joh_label = factor(mgmt_type_joh,
                                 levels = 0:3,
                                 # level 1 spans clear/grndwrk/logging_trail, hence "sitePrep" not "salvaged"
                                 labels = c("passive_joh", 
                                            "salvaged_joh", 
                                            "planted_joh", 
                                            "protected_joh")),
    
    # ---- C) TWO-AXIS: establishment pathway + independent protection ---------
    # Axis 1 (ordinal): passive < sitePrep < planted; planting takes priority
    establishment = case_when(
      any_na   ~ NA_integer_,
      planted  ~ 2L,
      sitePrep ~ 1L,
      TRUE     ~ 0L
    ),
    establishment = factor(establishment,
                           levels = 0:2, labels = c("passive", "salvaged", "planted"),
                           ordered = TRUE),
    
    # Axis 2 (binary): browsing protection, independent of the pathway
    protect = if_else(any_na, NA_integer_, as.integer(protected)),
    
    # combined descriptive label (base x protection) — for tabulation/plots only
    mgmt_2axis = factor(
      if_else(is.na(establishment),
              NA_character_,
              paste0(as.character(establishment),
                     if_else(protect == 1L, "_protected", ""))),
      levels = c("passive", 
                 "passive_protected",
                 "salvaged", 
                 "salvaged_protected",
                 "planted",  
                 "planted_protected"))
  ) %>%
  select(-any_na)  # drop helper

# distribution of each scheme 
dist_ladder <- df_mng_sub_cum %>%
  count(mgmt_type, mgmt_type_label, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1))

dist_joh <- df_mng_sub_cum %>%
  count(mgmt_type_joh, mgmt_type_joh_label, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1))

dist_2axis <- df_mng_sub_cum %>%
  count(establishment, protect, mgmt_2axis, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1))

print(dist_ladder)
print(dist_joh)
print(dist_2axis)

#### Crosswalk: how the three schemes agree / diverge ----------------------------
# (off-diagonals here ARE your methods-section justification)
# ours vs Johannes
crosswalk_ladder_joh <- df_mng_sub_cum %>%
  count(mgmt_type_label, mgmt_type_joh_label) %>%
  pivot_wider(names_from = mgmt_type_joh_label, values_from = n, values_fill = 0)

# full three-way (where do the "uncommon" subplots land under the other schemes?)
crosswalk_all3 <- df_mng_sub_cum %>%
  count(mgmt_type_label, mgmt_type_joh_label, mgmt_2axis) %>%
  arrange(desc(n))

print(crosswalk_ladder_joh)
print(crosswalk_all3)

## export management types
df_mng_sub_cum_out <- df_mng_sub_cum %>%
  select(-ends_with("_intensity"))




### Management intensity summary (for figures, EU + CZK level) ------------------------------
df_master_mng_intensity <- dat_regen %>%
  distinct(plot, country_name,
           year,
           clear_intensity, 
           grndwrk_intensity, 
           logging_trail_intensity,
           planting_intensity, 
           anti_browsing_intensity)

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

## Intensity class summary table: for all data: EU from 2023 and czechia 2025!!
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

intensity_class_summary_EU_CZ <- intensity_class_summary_percent %>%
  left_join(low_high_summary, by = "activity") %>%
  mutate(high_n = as.numeric(str_extract(High, "^[0-9]+"))) %>%
  arrange(desc(high_n)) %>%
  select(-high_n)

intensity_class_summary_EU_CZ

#### management intensity plot - EU + CZ make a plot to illustrate management for all samples -------------------------

# yellow -> red, low -> high intensity (ColorBrewer YlOrRd, 5 classes)
col_intensity <- setNames(
  c("#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026"),
  intensity_levels            # "0–19", "20–39", "40–59", "60–79", "80–100"
)

activity_labels <- c(
  clear_intensity          = "Clearing",
  grndwrk_intensity        = "Ground work",
  planting_intensity       = "Planting",
  anti_browsing_intensity  = "Anti-browsing",
  logging_trail_intensity  = "Logging trails"
)

# ---- 1) overall: all plots pooled (EU 2023 + Czechia 2025) ----------------
plot_manag_intensity <- intensity_class_summary_counts %>%
  ungroup() %>%
  group_by(activity) %>%
  mutate(share = n_plots / sum(n_plots) * 100) %>%
  ungroup() %>%
  mutate(intensity_class = factor(intensity_class, levels = intensity_levels),
         # order activities by share of high classes (60-100)
         activity = fct_reorder(activity,
                                if_else(intensity_class %in% intensity_levels[4:5], share, 0),
                                .fun = sum))

p_manag_intensity_EU_CZ <- ggplot(plot_manag_intensity,
                aes(x = share, y = activity, fill = fct_rev(intensity_class))) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = if_else(share >= 4, sprintf("%.0f%%", share), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = col_intensity, breaks = intensity_levels,
                    name = "Intensity (%)") +
  scale_y_discrete(labels = activity_labels) +
  scale_x_continuous(labels = label_percent(scale = 1), expand = c(0, 0)) +
  labs(x = "Share of plots", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom")

p_manag_intensity_EU_CZ

#ggsave("outFigs/mng_intensity_all_countries.png", p_manag_intensity_EU_CZ,
#       width = 7, height = 5, dpi = 300, bg = 'white')


# management per country 

## ---- 1) intensity classes per country (and year for Czechia) 
mng_by_country <- df_master_mng_intensity %>%
  mutate(
    # Czechia has two surveys -> separate bars; other countries have one record
    country_lab = if_else(country_name == "Czechia",
                          paste0("Czechia ", year),
                          as.character(country_name))
  ) %>%
  pivot_longer(cols = c(clear_intensity, grndwrk_intensity, logging_trail_intensity,
                        planting_intensity, anti_browsing_intensity),
               names_to = "activity", values_to = "applied") %>%
  filter(!is.na(applied)) %>%
  mutate(intensity_class = cut(applied,
                               breaks = c(-Inf, 0.19, 0.39, 0.59, 0.79, Inf),
                               labels = intensity_levels, right = TRUE)) %>%
  count(country_lab, activity, intensity_class, name = "n_plots") %>%
  group_by(country_lab, activity) %>%
  mutate(share = n_plots / sum(n_plots) * 100) %>%
  ungroup() %>%
  mutate(intensity_class = factor(intensity_class, levels = intensity_levels),
         activity = factor(activity, levels = names(activity_labels)))

## sample size per bar (same for all activities within a country)
n_per_country <- df_master_mng_intensity %>%
  mutate(country_lab = if_else(country_name == "Czechia",
                               paste0("Czechia ", year), as.character(country_name))) %>%
  count(country_lab, name = "n_total")

mng_by_country <- mng_by_country %>%
  left_join(n_per_country, by = "country_lab") %>%
  mutate(country_n = paste0(country_lab, " (n=", n_total, ")"))

## order countries by overall share of high-intensity management (60-100 %)
country_order <- mng_by_country %>%
  filter(intensity_class %in% intensity_levels[4:5]) %>%
  group_by(country_n) %>%
  summarise(high = sum(share) / n_distinct(activity), .groups = "drop") %>%
  arrange(high) %>%
  pull(country_n)

mng_by_country <- mng_by_country %>%
  mutate(country_n = factor(country_n, levels = country_order))

#### ---- 2) plot: one panel per activity, one bar per country ------------------
p_manag_intensity_country <- ggplot(
  mng_by_country,
  aes(x = share, y = country_n, fill = fct_rev(intensity_class))) +
  geom_col(width = 0.75, colour = "grey30", linewidth = 0.15) +
  facet_wrap(~ activity, nrow = 1, labeller = labeller(activity = activity_labels)) +
  scale_fill_manual(values = col_intensity, breaks = intensity_levels,
                    name = "Intensity (%)") +
  scale_x_continuous(labels = label_percent(scale = 1, suffix = ""),
                     breaks = c(0, 50, 100), expand = c(0, 0)) +
  labs(x = "Share of plots (%)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.spacing.x   = unit(1, "lines"),
        strip.text        = element_text(face = "bold"),
        legend.position   = "bottom")

p_manag_intensity_country

#ggsave("outFigs/mng_intensity_by_country_EU_CZ.png", p_manag_intensity_country,
#        width = 13, height = 6, dpi = 300)


# 3. Derived variables -------------------------------------------

## 3a. Importance values (IV)
dat_iv <- dat_regen %>%
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


iv_leaf_sub_wide <- iv_leaf_sub %>%
  dplyr::select(subplot, year, leaf_type, IV) %>%
  tidyr::pivot_wider(names_from = leaf_type, values_from = IV, names_prefix = "IV_")




## 3b. Community-weighted means (CWM: shade & drought tolerance)
wvar <- "n"   # weight by stem counts; alternatives: "basal_area_cm2" or "hgt_est"

cwm_subplot <- dat_regen %>%
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

cwm_plot <- dat_regen %>%
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
spruce_share_plot <- dat_regen %>%
  filter(!is.na(n)) %>%
  group_by(plot, year) %>%
  summarise(total_stems  = n(),
            spruce_stems = sum(species == "piab"),
            spruce_share = spruce_stems / total_stems,
            .groups = "drop") %>%
  select(plot, year, spruce_share)

spruce_share_sub <- dat_regen %>%
  filter(!is.na(n)) %>%
  group_by(plot, subplot, year) %>%
  summarise(total_stems  = n(),
            spruce_stems = sum(species == "piab"),
            spruce_share = spruce_stems / total_stems,
            .groups = "drop") %>%
  select(plot, subplot, year, spruce_share)

## 3d. Subplot summary table
field_sub_summ <- dat_regen %>%
  mutate(n = coalesce(n, 0L)) %>%
  group_by(plot, subplot, year,
           time_snc_full_disturbance, time_snc_part_disturbance,
           disturbance_year, forest_year, disturbance_length,
           clear, grndwrk, logging_trail, planting, anti_browsing) %>%
  summarise(
    stems_total       = sum(n, na.rm = TRUE),
    sp_richness       = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp        = if (stems_total > 0) vegan::diversity(tapply(n, species, sum), index = "shannon") else 0,  # <- changed
    evenness_sp       = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
    #sp_richness       = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
   # shannon_sp        = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
   # evenness_sp       = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
   # effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
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

# is this step needed?
nrow(field_sub_summ) # 5460 - > 5480

#field_sub_summ_cleaned <- field_sub_summ %>%
#  group_by(plot, subplot, year) %>%
#  summarise(across(where(is.numeric), safe_max), .groups = "drop")
#nrow(field_sub_summ_cleaned) #5460

## 3e. Plot summary table (pooled across subplots)
plot_metrics_pooled <- dat_regen %>%
  mutate(n = coalesce(n, 0L)) %>%
  group_by(plot, year) %>%
  summarise(
    n_subplots_plot = n_distinct(subplot),     # account that i can have 4 or 5 subplots per plot ()
    stems_total       = sum(n, na.rm = TRUE),
    sp_richness       = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp        = if (stems_total > 0) vegan::diversity(tapply(n, species, sum), index = "shannon") else 0,  # <- changed
    evenness_sp       = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
   # sp_richness       = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
   # shannon_sp        = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
   # evenness_sp       = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    mean_hgt          = if (stems_total > 0) weighted.mean(hgt_est, n, na.rm = TRUE) else NA_real_,
  #  effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
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


# is not effcetive numbers correctly calculated?
plot(plot_metrics_pooled$sp_richness, plot_metrics_pooled$effective_numbers)
abline(0, 1, lty = 2)   # 1:1 line, no point should lie above it

## 3f. Plot context
df_plot_context <- dat_regen %>%
  dplyr::select(plot, year,
                pre_dist_trees_n, area_m2, pre_dist_dens_ha,
                time_snc_full_disturbance,
                disturbance_year, forest_year, disturbance_length,
                x, y) %>%
  group_by(plot, year) %>%
  mutate(x = mean(x, na.rm = TRUE), y = mean(y, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct()

nrow(df_plot_context)  # 1092  # EU + CZ

## 3g. Bind subplot + plot into unified analysis table
sub_df <- #field_sub_summ_cleaned %>%
  field_sub_summ |> 
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
    ID      = plot, 
    plot_id = plot, year, 
    level = "plot",
    dens_m2 = stems_total / (n_subplots_plot * 4),   # number of subplots per plot + 4 m2
    cv_hgt, 
    mean_hgt, 
    sp_richness, 
    shannon_sp, 
    evenness_sp,
    effective_numbers, 
    w = stems_total
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


# add country indicator
both_levels_re2 <- both_levels_re2 %>% 
  mutate(
  region = str_extract(plot_id, "^\\d+(?=_)"),
  country = case_when(
    region %in% c("11", "12", "14", "18", "19", "20", "25") ~ 11L,
    region == "17"                                           ~ 12L,
    region %in% c("15", "26")                               ~ 13L,
    region == "13"                                           ~ 14L,
    region == "16"                                           ~ 15L,
    region == "23"                                           ~ 16L,
    region == "21"                                           ~ 17L,
    region == "22"                                           ~ 18L,
    region %in% c("24", "27")                               ~ 19L,
    is.na(region)                                           ~ 13L,  # bare IDs -> Czechia
    TRUE                                                     ~ NA_integer_
  ),
  country_name = case_when(
    country == 11 ~ "Germany",
    country == 12 ~ "Poland",
    country == 13 ~ "Czechia",
    country == 14 ~ "Austria",
    country == 15 ~ "Slovakia",
    country == 16 ~ "Slovenia",
    country == 17 ~ "Italy",
    country == 18 ~ "Switzerland",
    country == 19 ~ "France",
    TRUE          ~ NA_character_
  )
) %>% 
  mutate(country_name   = factor(country_name),
         region = factor(region),
         country_name = factor(country_name))

cz_counts <- both_levels_re2 %>%
   dplyr::filter(country_name == "Czechia") %>%
  group_by(year) %>%
  dplyr::summarise(n_plots     = n_distinct(plot_id[level == "plot"]),
            n_subplots  = sum(level == "subplot"),
            .groups = "drop")
cz_counts

# plots surveyed in both years
both_levels_re2 %>%
  filter(level == "plot") %>%
  filter(country_name == "Czechia") %>%
  distinct(plot_id, year) %>%
  count(plot_id) %>%
  count(n, name = "n_plots")

# A tibble: 2 × 2
# n n_plots
# <int>   <int>
#   1     1      83
# 2     2     126


# plots and subplots per year in Czechia
both_levels_re2 %>%
  filter(country_name == "Czechia") %>%
  group_by(year) %>%
  summarise(n_plots    = n_distinct(plot_id[level == "plot"]),
            n_subplots = n_distinct(ID[level == "subplot"]),
            .groups = "drop")

# year n_plots n_subplots
# <dbl>   <int>      <int>
#   1  2023     134        670
# 2  2025     199        995



table(both_levels_re2$cv_hgt_present)

df_plot_clean <- both_levels_re2 %>% filter(level == "plot")
df_sub_clean  <- both_levels_re2 %>% filter(level == "subplot")

## 3h. AEF export table (beta added after §4)

sub_df_AEF <- df_sub_clean %>%
  select(-level, -w, 
         -x, -y) %>% # remove wrong coordinates
  mutate(year    = as.integer(year),
         plot_id = as.character(plot_id)) %>%
  rename(plot = plot_id,
         subplot = ID) %>%
  left_join(iv_leaf_sub_wide, by = join_by(subplot, year)) %>%
  dplyr::left_join(dat_sub_xy, by = dplyr::join_by(subplot, year)) %>%   # add new averaged xy
  mutate(
    time_snc_full_disturbance = pmin(time_snc_full_disturbance, 8),
    plot_id = factor(plot),
    year_f  = factor(year),
  ) %>%
  mutate(across(all_of(c("cv_hgt", "mean_hgt", "sp_richness", "effective_numbers")),
                ~ replace_na(.x, 0)))



plot_df_AEF <- df_plot_clean %>%
  select(-level, -ID, -w,
         -x, -y) %>%
  mutate(year    = as.integer(year),
         plot_id = as.character(plot_id)) %>%
  rename(plot = plot_id) %>%
  left_join(iv_leaf_plot_wide, by = join_by(plot, year)) %>%
  dplyr::left_join(dat_plot_xy, by = dplyr::join_by(plot, year)) %>%   # add new averaged xy
  mutate(
    time_snc_full_disturbance = pmin(time_snc_full_disturbance, 8),
    plot_id = factor(plot),
    year_f  = factor(year)
  ) %>%
  mutate(across(all_of(c("cv_hgt", "mean_hgt", "sp_richness", "effective_numbers")),
                ~ replace_na(.x, 0)))


# 4. Species composition & beta diversity -------------------------
# remove mature treees!!
## 4a. Stem share per species per year
df <- dat_overlap %>%
  dplyr::filter(vegtype != 'mature') |> 
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
comm_sub <- dat_regen %>%
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


# 5. Export EU level data ----------------------------------

# ── Save objects for EU& CZ analysis input data 
data.table::fwrite(both_levels_re2, "outData/both_levels_re2_all_countries.csv")  # mature stand is already removed
data.table::fwrite(plot_df_AEF2,    "outData/plot_df_AEF2_all_countries.csv")
data.table::fwrite(sub_df_AEF,      "outData/sub_df_AEF_all_countries.csv")
data.table::fwrite(beta_plotyear,   "outData/beta_plotyear_all_countries.csv")
data.table::fwrite(dat_regen,       "outData/dat_regen_species_all_countries.csv")

# mature trees on EU level - 2023
data.table::fwrite(dat_counts_mature_species23, 'outData/mature_trees_count_EU_23.csv') # export Mature trees position at 2023 data
data.table::fwrite(dat_counts_mature_species25, 'outData/mature_trees_count_EU_25.csv') # export Mature trees position at 2025 data (only czechia)

# export for Karim - types of management
data.table::fwrite(df_mng_sub_cum_out, "outDataShare/Karim_AEF/cleaned/subplot_cumul_management_20260727.csv")

#mng_wide
data.table::fwrite(mng_wide,   "outData/mng_intensity_year_CZ.csv")

## 6. Spatial export for collaborators -------------------------------------
plot_sf_AEF <- plot_df_AEF %>%
  sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE)

sf::st_write(
  plot_sf_AEF,
  dsn          = "outDataShare/Karim_AEF/regeneration_chars_plot_3035.gpkg",
  layer        = "plot",
  delete_layer = TRUE
)

sub_sf_AEF <- sub_df_AEF %>%
  sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE)

sf::st_write(
  sub_sf_AEF,
  dsn          = "outDataShare/Karim_AEF/regeneration_chars_subplot_3035.gpkg",
  layer        = "subplot",
  delete_layer = TRUE
)

# check for coordinate systems
summary(plot_df_AEF$x)
summary(plot_df_AEF$y)
summary(sub_df_AEF$x)
summary(sub_df_AEF$y)



# Check CRS
sf::st_crs(plot_sf_AEF)
sf::st_crs(sub_sf_AEF)

# Plot together
plot(sf::st_geometry(plot_sf_AEF), col = "blue")
plot(sf::st_geometry(sub_sf_AEF), col = "red", add = TRUE)



# Sanity check 2 - on map

library(rnaturalearth)

europe <- ne_countries(continent = "Europe", returnclass = "sf") %>%
  st_transform(3035)

ggplot() +
  geom_sf(data = europe, fill = "grey90", colour = "white") +
  geom_sf(data = plot_sf_AEF, colour = "steelblue", size = 0.8) +
  coord_sf(xlim = c(3800000, 5400000), ylim = c(2400000, 3500000))

