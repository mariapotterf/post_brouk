
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

## Project variables (palettes, labels, species lists, species_class, earlyspecs_cz)
source('my_variables.R')



# Filter data from both level: plot and subplot from EU  ----------------------
# these are summarized on plot and subplot level 
both_levels_cz <- fread("outData/both_levels_re2_all_countries.csv") %>%
  filter(country_name == "Czechia") %>%
  mutate(
    level        = factor(level, levels = c("subplot", "plot")),
    plot_id      = factor(plot_id),
    year_f       = factor(year),
    country_name = factor(country_name),
    region       = factor(region)
  )

plot_df_cz <- fread("outData/plot_df_AEF2_all_countries.csv") %>%
  filter(country_name == "Czechia") %>%
  mutate(
    plot_id      = factor(plot),
    year_f       = factor(year),
    country_name = factor(country_name),
    region       = factor(region)
  )

sub_df_cz <- fread("outData/sub_df_AEF_all_countries.csv") %>%
  filter(country_name == "Czechia") %>%
  mutate(
    plot_id      = factor(plot),
    year_f       = factor(year),
    country_name = factor(country_name),
    region       = factor(region)
  )

# data with individual species per subplot and plots each species
dat_overlap <- fread("outData/dat_full_species_all_countries.csv") %>%
  filter(country_name == "Czechia") %>%
  mutate(
    plot      = factor(plot)
  )


# ── Sanity check ───────────────────────────────────────────────
cat("Plots CZ:   ", n_distinct(both_levels_cz$plot_id), "\n")
cat("Years:      ", sort(unique(both_levels_cz$year)), "\n")
cat("Levels:     ", levels(both_levels_cz$level), "\n")
cat("plot_df rows:", nrow(plot_df_cz), "\n")
cat("sub_df rows: ", nrow(sub_df_cz), "\n")


# Species composition -------------------------------------------------------------------


# Summary --------------------------------------
## Analyze data: first check up ----------------------
# get master table, having all unique plots and subplots - even teh empty ones
dat_master_subplot <- dat_overlap %>% 
  dplyr::select(plot, subplot, year, status) %>% 
  distinct()

table(dat_master_subplot$year)  
table(dat_master_subplot$status)  

n_plots_total <- length(unique(dat_master_subplot$plot))     # 208, 333
n_plots_total # 208,. 333 over 2 years, 125 recoccurs

n_subplots_total <-length(unique(dat_master_subplot$subplot))  # 1665
n_subplots_total # 1665


# filter only data with stem values 
dat_overlap_populated <- dat_overlap %>% 
  dplyr::filter(n>0)

# make a dataset where if n is NA, replace by 0
dat_overlap_n0 <- dat_overlap %>% 
  mutate(n = ifelse(is.na(n), 0, n))  


## Species composition --------------------------------------------------


### Get summary across all trees and study sites --------------------------
# 0) Safe counts (treat NA counts as 0)
df <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>%
  filter(!is.na(species) & species != "")


# find species with the highest share of stems per year
# Year-specific totals (denominator per year)
year_totals <- df %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    total_trees = sum(n, na.rm = TRUE),
    n_plots = dplyr::n_distinct(plot),
    n_subplots = dplyr::n_distinct(subplot),
    .groups = "drop"
  )

# 2) Species x year summary + year-specific share
species_stem_share_year <- df %>%
  dplyr::group_by(year, species) %>%
  dplyr::summarise(
    stems = sum(n, na.rm = TRUE),
    plots_present = dplyr::n_distinct(plot[n > 0]),
    .groups = "drop"
  ) %>%
  dplyr::left_join(year_totals %>% dplyr::select(year, total_trees), by = "year") %>%
  dplyr::mutate(
    share = round(100 * stems / total_trees, 2)
  ) %>%
  dplyr::arrange(dplyr::desc(share), species)

species_stem_share_year

# 3) Top species: option 1 = top 10 PER YEAR
top10_by_year <- species_stem_share_year %>%
  dplyr::group_by(year) %>%
  dplyr::slice_max(order_by = share, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

# if i have two years, thise leads to 11 species: "absp" "acps" "besp" "fasy" "lade" "piab" "pisy" "potr" "qusp" "saca" "soau"
# 'absp' has the lowest share: 1.95%, only on 10 plots - can i remove it?

# 4) Palette for plotting (consistent colors across ALL top species across years)
# order by total share across both years (ascending = small top, big bottom for coord flip)
species_levels <- top10_by_year %>%
  group_by(species) %>%
  summarise(share_all = sum(share, na.rm = TRUE), .groups = "drop") %>%
  arrange(share_all) %>%
  pull(species)

# v_top_species follows the same order (not alphabetical)
v_top_species <- species_levels

species_levels



# Between years compositional change at PLOT level -----------------
plots_both_years <- dat_overlap %>%
  filter(status == "both") %>%
  distinct(plot, year)

# Species presence per plot × year (at least 1 stem)
sp_presence <- dat_overlap %>%
  filter(status == "both", n > 0) %>%
  distinct(plot, year, species)

# Wide: one row per plot, columns = species, values = year presence
sp_wide <- sp_presence %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = year, 
              values_from = present,
              values_fill = 0, names_prefix = "yr_") %>%
  rename(presence_2023 = yr_2023, 
         presence_2025 = yr_2025)

# Classify each species × plot as: persistent, new arrival, local loss
sp_change <- sp_wide %>%
  mutate(change = case_when(
    presence_2023 == 1 & presence_2025 == 1 ~ "persistent",
    presence_2023 == 0 & presence_2025 == 1 ~ "new_arrival",
    presence_2023 == 1 & presence_2025 == 0 ~ "local_loss",
    TRUE ~ NA_character_
  ),
  seral_stage = case_when(
    species %in% earlyspecs_cz ~ "early",
    TRUE ~ "late"          # includes spruce, fir, beech, oak etc.
  )
  ) %>%
  left_join(species_class, by = "species")  # brings in seral_stage / leaf_type

# Per-plot summary: n gains, n losses, net change
plot_turnover <- sp_change %>%
  group_by(plot) %>%
  summarise(
    n_persistent  = sum(change == "persistent"),
    n_new_arrival = sum(change == "new_arrival"),
    n_local_loss  = sum(change == "local_loss"),
    net_change    = n_new_arrival - n_local_loss,
    .groups = "drop"
  )

# Seral-stage breakdown of gains/losses
seral_change <- sp_change %>%
  filter(!is.na(seral_stage)) %>%
  group_by(plot, seral_stage, change) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = change, values_from = n, values_fill = 0)


# Seral stage stem share per plot × year
seral_share_plot <- dat_overlap %>%
  filter(status == "both", !is.na(n), n > 0) %>%
  mutate(seral_stage = if_else(species %in% earlyspecs_cz, "early", "late")) %>%
  group_by(plot, year) %>%
  mutate(total_stems = sum(n)) %>%
  group_by(plot, year, seral_stage) %>%
  summarise(
    share = sum(n) / first(total_stems),
    .groups = "drop"
  )
# Did the early:late ratio change between years?
seral_ratio <- seral_share_plot %>%
  pivot_wider(
    names_from  = c(seral_stage, year),
    values_from = share,
    values_fill = 0
  ) %>%
  # now one row per plot, columns: early_2023, late_2023, early_2025, late_2025
  mutate(early_share_2023 = early_2023,
         early_share_2025 = early_2025)

# check - should be ~125 rows (plots surveyed both years), 6 columns
glimpse(seral_ratio)

# only keep plots that have data in BOTH years
seral_ratio_both <- seral_ratio %>%
  filter(!is.na(early_2023) | late_2023 > 0,
         !is.na(early_2025) | late_2025 > 0)

wilcox.test(seral_ratio_both$early_share_2023,
            seral_ratio_both$early_share_2025,
            paired = TRUE)



# What are the actual median shares?
seral_ratio_both %>%
  summarise(
    median_early_2023 = median(early_share_2023, na.rm = TRUE),
    median_early_2025 = median(early_share_2025, na.rm = TRUE),
    mean_early_2023   = mean(early_share_2023, na.rm = TRUE),
    mean_early_2025   = mean(early_share_2025, na.rm = TRUE)
  )

# And how many plots are early-dominated vs late-dominated in each year?
seral_ratio_both %>%
  summarise(
    early_dom_2023 = sum(early_share_2023 > 0.5, na.rm = TRUE),
    late_dom_2023  = sum(early_share_2023 < 0.5, na.rm = TRUE),
    early_dom_2025 = sum(early_share_2025 > 0.5, na.rm = TRUE),
    late_dom_2025  = sum(early_share_2025 < 0.5, na.rm = TRUE)
  )

# And look at the distribution shape
par(mfrow = c(1,2))
hist(seral_ratio_both$early_share_2023, breaks = 20, 
     main = "2023", xlab = "Early seral share")
hist(seral_ratio_both$early_share_2025, breaks = 20,
     main = "2025", xlab = "Early seral share")


# Quantify the bimodality more explicitly
library(ggalluvial)

# Build the data
traj_alluvial <- seral_ratio_both %>%
  mutate(
    traj_2023 = case_when(
      early_share_2023 < 0.2 ~ "Late-dominated",
      early_share_2023 > 0.8 ~ "Early-dominated",
      TRUE                   ~ "Mixed"
    ),
    traj_2025 = case_when(
      early_share_2025 < 0.2 ~ "Late-dominated",
      early_share_2025 > 0.8 ~ "Early-dominated",
      TRUE                   ~ "Mixed"
    ),
    traj_2023 = factor(traj_2023, levels = c("Early-dominated", "Mixed", "Late-dominated")),
    traj_2025 = factor(traj_2025, levels = c("Early-dominated", "Mixed", "Late-dominated"))
  ) %>%
  count(traj_2023, traj_2025, name = "n")

# Alluvial plot
p_alluvial <- ggplot(traj_alluvial,
                     aes(axis1 = traj_2023, axis2 = traj_2025, y = n)) +
  geom_alluvium(aes(fill = traj_2023), width = 1/4,
                alpha = 0.6, knot.pos = 0.4) +
  geom_stratum(aes(fill = after_stat(stratum)),  # ← color bars by class
               width = 1/4, color = "white", linewidth = 0.5) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3.2) +
  scale_x_discrete(limits = c("2023", "2025"),
                   expand = expansion(mult = c(0.12, 0.12))) +
  scale_fill_manual(
    values = c(
      "Early-dominated" = "#4dac26",
      "Mixed"           = "#b8b8b8",
      "Late-dominated"  = "#d01c8b"
    ),
    guide = "none"
  ) +
  labs(x = NULL, y = "Number of plots") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(size = 11, face = "bold"),
    axis.text.y  = element_text(size = 9),
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

p_alluvial


# does teh management pushes teh change? 

traj_with_mng <- seral_ratio_both %>%
  mutate(
    traj_2023 = case_when(
      early_share_2023 < 0.2 ~ "late_dominated",
      early_share_2023 > 0.8 ~ "early_dominated",
      TRUE                   ~ "mixed"
    ),
    traj_2025 = case_when(
      early_share_2025 < 0.2 ~ "late_dominated",
      early_share_2025 > 0.8 ~ "early_dominated",
      TRUE                   ~ "mixed"
    ),
    # Did plot stay stable or change?
    stable = traj_2023 == traj_2025,
    # Did it move toward late or early seral?
    direction = case_when(
      traj_2023 == "early_dominated" & traj_2025 == "mixed"          ~ "toward_late",
      traj_2023 == "mixed"           & traj_2025 == "late_dominated"  ~ "toward_late",
      traj_2023 == "mixed"           & traj_2025 == "early_dominated" ~ "toward_early",
      traj_2023 == "late_dominated"  & traj_2025 == "mixed"           ~ "toward_early",
      stable                                                           ~ "stable",
      TRUE                                                             ~ "other"
    )
  ) %>%
  left_join(
    plot_df_cz %>%
      filter(year == 2023) %>%
      select(plot, planting_intensity, time_snc_full_disturbance),
    by = "plot"
  )

# Is planting intensity different between stable and changing plots?
kruskal.test(planting_intensity ~ direction, data = traj_with_mng)

# Pairwise if significant
pairwise.wilcox.test(traj_with_mng$planting_intensity,
                     traj_with_mng$direction,
                     p.adjust.method = "BH")

# Visualize
ggplot(traj_with_mng, aes(x = direction, y = planting_intensity, fill = direction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c(
    "toward_early" = "#4dac26",
    "toward_late"  = "#d01c8b",
    "stable"       = "#b8b8b8",
    "other"        = "#dddddd"
  )) +
  labs(x = NULL, y = "Planting intensity") +
  theme_classic() +
  theme(legend.position = "none")



# Delta early share per plot (continuous, not categorical)
delta_df <- seral_ratio_both %>%
  mutate(delta_early = early_share_2025 - early_share_2023) %>%
  left_join(
    plot_df_cz %>%
      filter(year == 2023) %>%
      select(plot, planting_intensity, time_snc_full_disturbance),
    by = "plot"
  )

# Does planting predict direction of change?
cor.test(delta_df$planting_intensity, delta_df$delta_early, method = "spearman")

# Simple model
lm_delta <- lm(delta_early ~ planting_intensity + time_snc_full_disturbance,
               data = delta_df)
summary(lm_delta)

# Plot
ggplot(delta_df, aes(x = planting_intensity, y = delta_early)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", color = "#d01c8b", se = TRUE) +
  labs(x = "Planting intensity (2023)",
       y = "Change in early-seral share\n(2025 − 2023)") +
  theme_classic()












# connect transition between early and seral to  management
seral_ratio_both %>%
  mutate(
    traj_2023 = case_when(
      early_share_2023 < 0.2 ~ "late_dominated",
      early_share_2023 > 0.8 ~ "early_dominated",
      TRUE                   ~ "mixed"
    ),
    traj_2025 = case_when(
      early_share_2025 < 0.2 ~ "late_dominated",
      early_share_2025 > 0.8 ~ "early_dominated",
      TRUE                   ~ "mixed"
    ),
    transition = paste(traj_2023, traj_2025, sep = " → ")
  ) %>%
  left_join(
    plot_df_cz %>% 
      filter(year == 2023) %>%
      select(plot, planting_intensity, anti_browsing_intensity, 
             clear_intensity, time_snc_full_disturbance),
    by = "plot"
  ) %>%
  group_by(traj_2023, traj_2025) %>%
  summarise(
    n                  = n(),
    mean_planting      = mean(planting_intensity, na.rm = TRUE),
    mean_antibrowsing  = mean(anti_browsing_intensity, na.rm = TRUE),
    mean_time          = mean(time_snc_full_disturbance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(traj_2023, traj_2025)


# Quick visualization of the key contrast
traj_planting <- seral_ratio_both %>%
  mutate(
    traj_2023 = case_when(
      early_share_2023 < 0.1 ~ "late-dominated",
      early_share_2023 > 0.9 ~ "early-dominated",
      TRUE                   ~ "mixed"
    ),
    traj_2025 = case_when(
      early_share_2025 < 0.1 ~ "late-dominated",
      early_share_2025 > 0.9 ~ "early-dominated",
      TRUE                   ~ "mixed"
    )
  ) %>%
  left_join(
    plot_df_cz %>%
      filter(year == 2023) %>%
      select(plot, planting_intensity),
    by = "plot"
  ) %>%
  group_by(traj_2023, traj_2025) %>%
  summarise(
    n             = n(),
    mean_planting = mean(planting_intensity, na.rm = TRUE),
    sd_planting   = sd(planting_intensity,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n > 1) %>%   # drop the single early→late plot
  mutate(
    traj_2023 = factor(traj_2023, levels = c("early-dominated", "mixed", "late-dominated")),
    traj_2025 = factor(traj_2025, levels = c("early-dominated", "mixed", "late-dominated")),
    label = paste0("n=", n)
  )

ggplot(traj_planting,
       aes(x = traj_2023, y = traj_2025, fill = mean_planting)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = paste0(round(mean_planting, 2), "\n", label)),
            size = 3.5, lineheight = 1.3) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Mean planting\nintensity",
                       limits = c(0, 1)) +
  labs(x = "Trajectory 2023", y = "Trajectory 2025",
       title = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(size = 10),
        legend.position = "right")



# spruce change: legacy vs planted???
# Does spruce occurrence increase more on planted plots? ------------------------
spruce_occ_change <- sp_change %>%
  filter(species == "piab") %>%
  left_join(plot_df_cz %>% 
              filter(year == 2023) %>% 
              select(plot, planting_intensity), 
            by = "plot") %>%
  mutate(planted = planting_intensity > 0)

table(spruce_occ_change$change, spruce_occ_change$planted)
# if new_arrival concentrates on planted = TRUE, that confirms the planting signal

table(spruce_occ_change$change, spruce_occ_change$planting_intensity)




#### Single species plots chars ----------------------

# Summarize number of species per plot, only for n > 1 trees - why is that??
plots_with_species_counts <- df %>%
  filter(n >= 1) %>% 
  group_by(plot) %>%
  summarise(
    n_species = n_distinct(species),
    .groups = "drop"
  )

# Filter to plots with exactly 1 species
single_species_plots <- plots_with_species_counts %>%
  dplyr::filter(n_species == 1)

# Step 4: Join back to original to get the species name
species_per_single_species_plot <- df %>%
  filter(n >= 1) %>%
  semi_join(single_species_plots, by = "plot") %>%
  group_by(plot, species) %>%
  summarise(n_trees = sum(n), .groups = "drop")  # just to show tree count if needed

# Count how many plots per species
summary_species <- species_per_single_species_plot %>%
  count(species, name = "single_plots") %>%
  arrange(desc(single_plots)) %>% 
  mutate(share = single_plots/333)

# Print result
print(summary_species)


#### Stem density per species per top X species --------------------------------
df_stem_dens_species <- df %>% 
  group_by(plot, year, species, n_subplots ) %>%
  summarize(sum_n = sum(n, na.rm =T)) %>% 
  mutate(scaling_factor = 10000/(n_subplots * 4),
         stem_dens = sum_n*scaling_factor) #%>% 
# mutate(log_sum_stem_density = log10(stem_dens + 1)) #%>%  # Adding 1 to avoid log(0)
#ungroup()

# get total sum and calculate as average value over all sites 
df_stem_dens_species_sum <- 
  df_stem_dens_species %>% 
  group_by(species, year) %>% #, year
  summarise(stem_dens = sum(stem_dens, na.rm = T)) %>% #,
  mutate(stem_dens_avg = stem_dens/n_plots_total) #,


df_stem_dens_species <- 
  df_stem_dens_species %>% 
  ungroup(.) %>% 
  filter(sum_n >0) %>% 
  filter(species %in% v_top_species) %>% 
  dplyr::group_by(species, year) %>%
  dplyr::mutate(median_stem_density = median(stem_dens, na.rm = TRUE)) %>% 
  dplyr::ungroup(.) %>%
  mutate(species = factor(species, levels = rev(v_top_species))) # Set custom order



# make a barplot of stem occurence ---------------------------

#Make sure every species appears in both years (missing gets 0)
top10_by_year_clean <- top10_by_year %>%
  dplyr::select(year, species, share) %>%
  tidyr::complete(year, species = species_levels, fill = list(share = 0)) %>%
  mutate(
    species = factor(species, levels = species_levels),
    year = factor(year)  # for alpha legend control
  )

# 3) Plot: same hue per species, lighter/darker by year
p_bar <- top10_by_year_clean %>%
  ggplot(aes(x = share, y = species, fill = species, alpha = year)) +
  geom_col(position = position_dodge(width = 0.7), 
           width = 0.6, 
           aes(colour = factor(year))) + # "grey50"
  scale_x_continuous(labels = label_number(accuracy = 1, suffix = ""),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_discrete(
    limits = species_levels,
    labels = species_labels,   # keep your mapping
    drop = FALSE
  ) +
  scale_colour_manual(
    values = c("2023" = "grey50",
               "2025" = "black"),
    guide = "none"   # prevents second legend
  )+
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_alpha_manual(
    name = "Survey year",
    values = c("2023" = 0.45, "2025" = 1.00),  # lighter vs darker
    guide = "none",                              # set to "legend" if you want it
    breaks = c("2025", "2023")  # controls legend order
  ) +
  
  labs(x = "Stems share [%]", y = "") +
  theme_classic2(base_size = 10) +
  theme(
    axis.text.y = element_text(
      face = "italic",
      size = 8
    )
  )
p_bar

####  Get species occurence from total number of plots -------------
# Total number of unique plots
total_plots <- dat_overlap %>%
  pull(plot, year) %>%
  n_distinct()

# Share of plots per species (where species are present = non-zero stems)
species_occurence <- 
  df_stem_dens_species %>%
  ungroup(.) %>% 
  dplyr::filter(sum_n > 0) %>%                 # Only where species occurred
  distinct(species, year, plot) %>%    #year,       # Unique species × plot combos
  count(species, year, name = "n_plots")  %>% #year, %>%  # Count number of plots per species
  mutate(share_of_plots = n_plots / total_plots*100) %>% 
  arrange()

species_occurence

# Plot
p_occurence <- 
  species_occurence %>% 
  filter(species %in% v_top_species ) %>% 
  mutate(species = factor(species, levels = species_levels)) %>% 
  ggplot(aes(x = share_of_plots, 
             y = species, 
             fill = species,
             alpha = factor(year))) +
  geom_col(position = position_dodge(width = 0.7), 
           width = 0.6,
           aes(colour = factor(year))) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), 
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_manual(values = species_colors, 
                    guide = "none") +
  scale_alpha_manual(
    name = "Survey year",
    values = c("2025" = 1,
               "2023" = 0.45),
    breaks = c("2025", "2023")
  ) +
  
  scale_colour_manual(
    name = "Survey year",
    values = c("2025" = "black",
               "2023" = "grey50"),
    breaks = c("2025", "2023")
  )+
  labs(
    x = "Species occurence [%]",
    y = NULL#,
    # fill = "Year"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    
    # legend inside bottom-right
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    
    # make legend box readable inside plot
    #legend.background = element_rect(fill = "white", colour = "grey70"),
    legend.key = element_rect(fill = NA),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  theme(
    legend.key.width  = unit(1, "lines"),   # wider
    legend.key.height = unit(0.4, "lines")    # shorter → not square
  )


p_occurence

p_bar <- p_bar + 
  theme(plot.margin = margin(t = 12, 
                             r = 5, 
                             b = 5, 
                             l = 5))
p_occurence <- p_occurence + 
  theme(plot.margin = margin(t = 12, 
                             r = 5, 
                             b = 5, 
                             l = 5))


p_species_composition <- ggarrange(p_bar, p_occurence,# p_density,  
                                   ncol = 2, common.legend = F, 
                                   align = "h",        # <-- critical
                                   widths = c(1.5, 1),
                                   labels = c("[a]", "[b]"),
                                   font.label = list(size = 10, face = "plain"),
                                   label.x = 0.02,    # near left edge
                                   label.y = 1.01     # just above plot
)

p_species_composition

# Export to PNG
# ggsave("outFigs/p_species_composition.png", 
#        plot = p_species_composition,
#        width = 7, height = 3, units = "in", dpi = 300)



### Total trees per year --------------------------------------------------------
total_trees_per_year <- dat_overlap %>%
  group_by(year) %>%
  summarise(
    total_trees = sum(n, na.rm = TRUE),
    .groups = "drop"
  )

tree_summary <- dat_overlap %>%
  group_by(year, seral_stage ) %>%
  summarise(
    n_trees_recovery = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(total_trees_per_year) %>% 
  mutate(share = n_trees_recovery /total_trees * 100)


summary(tree_summary)





## 6b. Management intensity plot (simpler stacked bar)
fill_colors <- brewer.pal(length(intensity_levels), "YlOrRd")


# START
## Management Subplot level -----------------------------------------------------------

df_mng_sub <- dat_overlap %>% 
  #dplyr::filter(year == "2023") %>% # keep management oionly frm 2023 for consistency
  distinct(plot, subplot,year,
           clear,
           grndwrk,
           logging_trail,
           planting,
           anti_browsing,
           clear_intensity,           
           grndwrk_intensity,         
           logging_trail_intensity,  
           planting_intensity ,       
           anti_browsing_intensity#,   
          
  )


df_mng_plot <- df_mng_sub %>% 
  select(plot, year, 
         clear_intensity,           
         grndwrk_intensity,         
         logging_trail_intensity,  
         planting_intensity ,       
         anti_browsing_intensity) %>% 
  distinct()


mng_sub_props <- 
  df_mng_sub %>%
  # group by year???
  pivot_longer(cols = c(clear, grndwrk, logging_trail, planting, anti_browsing),
               names_to = "activity",
               values_to = "applied") %>%
  mutate(applied = replace_na(applied, 0)) %>% 
  group_by(activity, applied, year) %>%
  # View()
  summarise(n = n(), .groups = "drop") %>%
  group_by(activity) %>%
  mutate(proportion = n / sum(n)*100)

mng_sub_props


# First extract proportion for applied == 1 per activity
mng_sub_applied_order <- mng_sub_props %>%
  filter(applied == 1) %>%
  arrange(desc(proportion)) %>%
  pull(activity) %>% 
  unique()


# Prepare data for diverging bar plot
mng_sub_conv <-  mng_sub_props %>%
  mutate(proportion = ifelse(applied == 0, -proportion, proportion),
         applied = ifelse(applied == 1, "Presence", "Absence")) %>% 
  arrange(desc(proportion)) %>% 
  mutate(activity = factor(activity, levels = rev(mng_sub_applied_order))) #%>%

#library(forcats)

activity_labels <- c(
  "clear" = "Salvage logging",
  "grndwrk" = "Soil\npreparation", # "Groundwork", #
  "planting" = "Planting",
  "anti_browsing" = "Browsing\nprotection",
  "logging_trail" = "Logging trail"
)

# Create diverging bar plot - need to check for divide where site is salvage logging = 1 in 2023 but 0 in 2025
ggplot(mng_sub_conv, aes(x = proportion, y = activity, fill = applied)) +
  geom_col(width = 0.2, col = 'black') +
  scale_x_continuous(#labels = scales::percent_format(accuracy = 1),
    breaks = seq(-100, 100, 50),
    limits = c(-100, 100),
    labels = function(x) paste0(abs(x), "%")) +
  scale_fill_manual(values = c("Presence" = "red", 
                               "Absence" = "darkgreen")) +
  scale_y_discrete(labels = activity_labels) +
  geom_vline(xintercept = 0, color = 'darkgrey', lty = 'dashed')+
  labs(x = "Share of subplots [%]", 
       y = "",
       title = "",
       fill = "") +
  theme_classic2() +
  facet_grid(.~year) +
  annotate("text", x = -80, y = 5.5, label = "Absence", hjust = 0, size = 2.5, fontface = "bold") +
  annotate("text", x =  80, y = 5.5, label = "Presence",     hjust = 1, size = 2.5, fontface = "bold") +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 10, face = "italic"),
        legend.position = "none",
        text = element_text(size = 10))


# 
### Plot level: binary  ------------------------------------
df_master_mng_intensity <- dat_overlap %>% 
  #dplyr::filter(year == "2023") %>% # keep management oionly frm 2023 for consistency
  distinct(plot, year,
           # intensities
           clear_intensity,           
           grndwrk_intensity,         
           logging_trail_intensity,  
           planting_intensity ,      
           anti_browsing_intensity
           
  )
nrow(df_master_mng_intensity)

mng_intensity_props <- 
  df_master_mng_intensity %>%
  pivot_longer(cols = c(clear_intensity,           
                        grndwrk_intensity,         
                        logging_trail_intensity,  
                        planting_intensity ,       
                        anti_browsing_intensity
  ),
  names_to = "activity",
  values_to = "applied") %>%
  group_by(activity, applied, year) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(activity) %>%
  mutate(proportion = n / sum(n)*100)

mng_intensity_props

mng_intensity_props <- mng_intensity_props %>%
  mutate(
    intensity_class = cut(
      applied,
      breaks = c(-Inf, 0.19, 0.39, 0.59, 0.79, Inf),
      labels = c("0–19", "20–39", "40–59", "60–79", "80–100"),
      right = TRUE
    ),
    intensity_binary = case_when(
      intensity_class == "0–19" ~ "no",
      TRUE ~ "yes"
    )
  )


# prepare for binary classifucation  - need to group first, as I have two years
mng_sub_conv <- mng_intensity_props %>%
  mutate(
    intensity_binary = factor(intensity_binary, levels = c('no', 'yes')),
    intensity_class = factor(intensity_class, levels = intensity_levels)
  ) %>% 
  group_by(activity, intensity_binary) %>% 
  summarise(proportion = sum(proportion, na.rm = T))

# prepare for intensity classification
mng_sub_intensity <- mng_intensity_props %>%
  mutate(
    intensity_binary = factor(intensity_binary, levels = c('no', 'yes')),
    intensity_class = factor(intensity_class, levels = intensity_levels)#,
    #activity = factor(activity, levels = applied_mng_intens_order)
  ) %>% 
  group_by(activity, intensity_class) %>% 
  summarise(proportion = sum(proportion, na.rm = T))


# define colors 
fill_colors <- brewer.pal(length(intensity_levels), "YlOrRd")

# for binary: Flip "no" values to negative
mng_sub_conv_plot <- mng_sub_conv %>%
  mutate(proportion_plot = if_else(intensity_binary == "no", -proportion, proportion))

# for intensity - flip low intensity to left side of the plot
mng_shifted <- mng_sub_intensity %>%
  mutate(
    intensity_class = factor(intensity_class, levels = intensity_levels),
    proportion_shifted = if_else(intensity_class %in% low_classes, -proportion, proportion),
    # reverse class order for right-hand side only
    intensity_class_plot = fct_relevel(
      intensity_class,
      c("80–100","60–79", "40–59","0–19","20–39" )
    ),
    activity = factor(activity, levels = rev(c(
      "clear_intensity",
      "grndwrk_intensity",      
      "planting_intensity",     
      "anti_browsing_intensity",
      "logging_trail_intensity"
    )))
  )


activity_intens_labels <- c(
  "clear_intensity" = "Salvage\nlogging",
  "grndwrk_intensity" =  "Soil\npreparation", #"Groundwork", 
  "planting_intensity" = "Planting",
  "anti_browsing_intensity" = "Browsing\nprotection"#,
  #"logging_trail_intensity" = "Logging trail"
)

p_management_bin_plot <- mng_sub_conv_plot %>% 
  filter(activity != "logging_trail_intensity") %>% 
  droplevels(.) %>% 
  ggplot(aes(x = proportion_plot, y = activity,
             fill = intensity_binary )) +
  geom_col(width = 0.2, color = "black") +
  scale_fill_manual(values =  c("no" = "forestgreen", "yes" = "red")#,
                    #name = "Intensity class",
                    #breaks = intensity_levels
  ) +
  geom_vline(xintercept = 0, color = "grey", 
             linewidth = 0.8, lty = 'dashed') +
  ylab('') +
  scale_y_discrete(labels = activity_intens_labels) +   # 👈 this does the relabeling
  scale_x_continuous(labels = abs, name = "Plots share [%]") +
  annotate("text", x = -60, y = 4.5, label = "Not Applied", hjust = 0, size = 2.5, fontface = "bold") +
  annotate("text", x =  80, y = 4.5, label = "Applied",     hjust = 1, size = 2.5, fontface = "bold") +
  theme_classic2() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.grid.major.y = element_blank()
  )

p_management_bin_plot


p_management_intensity_plot <- mng_shifted %>% 
  filter(activity != "logging_trail_intensity") %>% 
  droplevels(.) %>% 
  ggplot(aes(x = proportion_shifted, 
             y = activity,
             fill = intensity_class_plot)) +
  geom_vline(xintercept = 0, color = "grey", 
             linewidth = 0.5, lty = 'dashed') +
  geom_col(width = 0.4, color = "black") +
  scale_fill_manual(values = fill_colors, name = "Intensity class",
                    breaks = intensity_levels) +
  ylab('') +
  scale_y_discrete(labels = activity_intens_labels) +   # 👈 this does the relabeling
  scale_x_continuous(labels = abs, name = "Plots share [%]") +
  annotate("text", x = -80, y = 4.5, label = "Low intensity", hjust = 0, size = 2.5, fontface = "bold") +
  annotate("text", x =  80, y = 4.5, label = "High intensity",     hjust = 1, size = 2.5, fontface = "bold") +
  theme_classic2() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )

p_management_intensity_plot

# # Save as PNG
ggsave("outFigs/mng_intensity_plot.png", plot = p_management_intensity_plot,
       width = 5, height = 2.1, units = "in", dpi = 300)



# make management intensity plot simpler
p_management_intensity_plot_simpler <- 
  mng_shifted %>% 
  filter(activity != "logging_trail_intensity") %>% 
  droplevels(.) %>% 
  ggplot(aes(x = proportion , 
             y = activity,
             fill = intensity_class_plot)) +
  geom_col(width = 0.4, color = "black") +
  scale_fill_manual(values = fill_colors, 
                    name = "Management intensity\nclass",
                    breaks = intensity_levels) +
  ylab('') +
  scale_y_discrete(labels = activity_intens_labels) +   # 👈 this does the relabeling
  scale_x_continuous(labels = abs, 
                     name = "Plots share [%]") +
  theme_classic2() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )

p_management_intensity_plot_simpler

# # Save as PNG
ggsave("outFigs/p_management_intensity_plot_simpler.png", 
       plot = p_management_intensity_plot_simpler,
       width = 5, height = 1.9, units = "in", dpi = 300)


# Summarize counts per class
intensity_class_summary_counts <- mng_intensity_props %>%
  group_by(activity, intensity_class) %>%
  summarise(n_plots = sum(n), .groups = "drop")

#Full class breakdown with % and formatted "n (%)"
intensity_class_summary_percent <- intensity_class_summary_counts %>%
  group_by(activity) %>%
  mutate(share = round(100 * n_plots / sum(n_plots), 1)) %>%
  unite(col = "value", n_plots, share, sep = " (") %>%
  mutate(value = paste0(value, "%)")) %>%
  pivot_wider(names_from = intensity_class, values_from = value)

# Add low/high total columns (counts + %)
low_high_summary <- intensity_class_summary_counts %>%
  mutate(group = case_when(
    intensity_class %in% c("0–19", "20–39") ~ "Low",
    intensity_class %in% c("40–59", "60–79", "80–100") ~ "High"
  )) %>%
  group_by(activity, group) %>%
  summarise(n_plots = sum(n_plots), .groups = "drop") %>%
  group_by(activity) %>%
  mutate(share = round(100 * n_plots / sum(n_plots), 1)) %>%
  unite(col = "value", n_plots, share, sep = " (") %>%
  mutate(value = paste0(value, "%)")) %>%
  pivot_wider(names_from = group, values_from = value)

# Combine and sort by High intensity count
intensity_class_summary_final <- intensity_class_summary_percent %>%
  left_join(low_high_summary, by = "activity") %>%
  mutate(
    high_n = as.numeric(str_extract(High, "^[0-9]+"))  # extract count before "("
  ) %>%
  arrange(desc(high_n)) %>%
  select(-high_n)  # optional: remove helper column

# View final table
intensity_class_summary_final 






## 6c. Disturbance history
plot_context_chars <- df %>%
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
dat_clean <- df %>%
  as_tibble() %>%
  dplyr::filter(!is.na(hgt_est), !is.na(n), n > 0)

height_species_mng_yr <- dat_clean %>%
  dplyr::group_by(year, species, planting_intensity, anti_browsing_intensity) %>%
  dplyr::summarise(total_stems          = sum(n, na.rm = TRUE),
                   mean_height_weighted = weighted.mean(hgt_est, w = n, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(seral_group = dplyr::if_else(species %in% earlyspecs_cz, "early", "late"))

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
ggsave("outFigsCZ/p_species_composition.png",
       p_species_composition, width = 7, height = 3, dpi = 300)

ggsave("outFigsCZ/p_violin.png",
       p_violin, width = 7, height = 2.5, dpi = 300)
ggsave("outFigsCZ/p_violin.svg",
       p_violin, width = 7, height = 2.5, dpi = 300)

ggsave("outFigsCZ/p_time_scale.png",
       p_final, width = 6, height = 5, dpi = 300, bg = "white")

ggsave("outFigsCZ/p_model_response.png",
       p_model_response, width = 7, height = 5, dpi = 300)

ggsave("outFigsCZ/p_management_intensity_plot_simpler.png",# - this one cuts
       p_management_intensity_plot_simpler, width = 5, height = 2.1, dpi = 300)

ggsave("outFigsCZ/p_combined_disturb_fig.png",
       p_combined_disturb_fig, width = 5, height = 2.5, dpi = 300)

ggsave("outFigsCZ/spruce_share_boxplots.png",
       p_spruce_shares_boxplot, width = 7, height = 3.2, dpi = 300)

ggsave("outFigsCZ/density_plot.png",
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


