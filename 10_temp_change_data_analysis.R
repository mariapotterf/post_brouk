
# Analyze data temporal change
#   

# read overlapping data: from 223 and 2025
# investigate how they change with time since disturbnace
# structure (height, cv, stem density), composition (shannon, richness, eveness)
# run models:


gc()

set.seed(1234) # assure reproducibility

library(terra)
library(sf)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggpubr)
library(corrplot)
library(GGally)
library(vegan) # for diversity indices

library(mgcv)
library(ggeffects)
library(gratia)

library(RColorBrewer)

library(ggridges)
library(scales)
library(forcats) # order factors
library(broom)


source('my_variables.R')


theme_set(theme_classic2(base_size = 10) +
            theme(axis.title = element_text(size = 10),
                  axis.text  = element_text(size = 10)))


# Read data -----------------------------
# dat_overlap_mng_upd2 <- fread('outData/full_table_overlap_23_25.csv')

# test 2025/11/03 -> then needs to rename the layer to 'dat'!
dat_overlap  <- fread('outData/full_table_23_25.csv')  # accound for all data points, not just the ovelapping ones
# this can help me to benefit from all disturbnace history sites



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
v_top_species <- top10_by_year %>%
  dplyr::distinct(species) %>%
  dplyr::pull(species) %>%
  sort()

# order species across both years 
species_levels <- top10_by_year %>%
  group_by(species) %>%
  summarise(share_all = sum(share, na.rm = TRUE), .groups = "drop") %>%
  arrange(share_all) %>%                 # small at top, big at bottom (since we reverse later)
  pull(species)

species_levels
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
            #log_sum_stem_density = sum(log_sum_stem_density, na.rm = T)) %>%
  mutate(stem_dens_avg = stem_dens/n_plots_total) #,
      #   log_sum_stem_density_avg = log_sum_stem_density/n_plots_total)


df_stem_dens_species <- 
  df_stem_dens_species %>% 
  ungroup(.) %>% 
  filter(sum_n >0) %>% 
  filter(species %in% v_top_species) %>% 
  dplyr::group_by(species, year) %>%
  dplyr::mutate(median_stem_density = median(stem_dens, na.rm = TRUE)) %>% 
  dplyr::ungroup(.) %>%
  mutate(species = factor(species, levels = rev(v_top_species))) # Set custom order

species_levels = factor(v_top_species)

#p_density<-
  df_stem_dens_species_sum %>% #df_stem_dens_species_year2 %>% 
  mutate(species = factor(species, levels = species_levels)) %>% 
  filter(!is.na(species)) %>% 
  ggplot(aes(x = stem_dens ,
             y = species,
             fill = factor(year))) +
  geom_boxplot(
    aes(group = interaction(species, year), 
      alpha = factor(year)),
    position = position_dodge(width = 0.6),
    outlier.shape = NA,
    width = 0.45#,
    # color = "black"
  ) +
  
  # coord_flip() +
  labs(
    x = "sum stem density",
    y = "",  
    fill = "Year"
  ) +
  scale_fill_manual(values= species_colors) +
  theme_classic(base_size = 10) +
  scale_y_discrete(labels = species_labels) +
  theme(axis.text.y = element_text(face = "italic", size = 8),
        legend.position = 'none')


#p_density

# make a barplot of stem occurence
  
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

# Optional: order species by max share across years
# species_order <- species_occurence %>%
#   group_by(species) %>%
#   summarise(max_share = max(share_of_plots)) %>%
#   arrange(desc(max_share)) %>%
#   pull(species)
# 
# species_plot_share <- species_occurence %>%
#   group_by(species) %>% 
#   summarize(share_of_plots_avg = mean(share_of_plots)) %>% 
#   mutate(species = factor(species, levels = rev(species_order)))
# 
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

# p_density with no y labels
# p_density <- p_density +
#   labs(y = NULL) +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())
# 

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

## Get species importance values (IV) - per species --------------------------------
# based on relative counts and relative basal areas/height
# 
# IV = mean(relative abundance, relative size: can be height- for regeneration or BA - for mature)
# Relative abundance = species stem count / total stem count (per plot or subplot)
# Relative size = (n × max height per species) / sum(n × max height) within unit
# IV ranges from 0–1 and is calculated separately at subplot and plot level

# double check how species dominance changes if using guestimated DBH for regeneration of height values:
# spearman correlation of 0.99, i can do either, or: teh most correct is using height fr or small regenerating trees 
dat_iv <- dat_overlap %>%
  mutate(
    size_height   = hgt_est * n#,
   # size_ba_guess = basal_area_cm2 * n
  )

calc_iv_core <- function(data, size_var, ...) {
  data %>%
    group_by(..., species) %>%
    summarise(
      n_sp    = sum(n, na.rm = TRUE),
      size_sp = sum({{ size_var }}, na.rm = TRUE),
      .groups = "drop_last"
    ) %>%
    mutate(
      RA = n_sp / sum(n_sp, na.rm = TRUE),
      size_tot = sum(size_sp, na.rm = TRUE),
      RS = dplyr::if_else(size_tot > 0, size_sp / size_tot, NA_real_),
      IV = (RA + RS) / 2
    ) %>%
    select(-size_tot) %>%
    ungroup()
}

calc_iv_subplot <- function(data, size_var) {
  calc_iv_core(data, {{ size_var }}, plot, year, subplot)
}
calc_iv_plot <- function(data, size_var) {
  calc_iv_core(data, {{ size_var }}, plot, year)
}

# calculate rIV per subplot, plot
iv_sub  <- calc_iv_subplot(dat_iv, size_height)
iv_plot <- calc_iv_plot(dat_iv, size_height)


# identify top species (max rIVI) per plot, subplot - if it is equal, choose a random one
iv_max_sub <- iv_sub %>%
  group_by(plot, year, subplot) %>%
  slice_max(order_by = IV, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  rename(IVmax_species = species) %>% 
  select(plot,   year, subplot, IV,   IVmax_species)


iv_max_plot <- iv_plot %>%
  group_by(plot, year) %>%
  slice_max(order_by = IV, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  rename(IVmax_species = species) %>% 
  select(plot,   year, IV, IVmax_species)

iv_max_sub
iv_max_plot



#### rIVI: coniferous vs deciduous (subplot)
iv_leaf_sub <- iv_sub %>%
  left_join(species_class, by = "species") %>%
  filter(!is.na(leaf_type)) %>%
  group_by(plot, year, subplot, leaf_type) %>%
  summarise(
    IV = sum(IV, na.rm = TRUE),
    RA = sum(RA, na.rm = TRUE),
    RS = sum(RS, na.rm = TRUE),
    .groups = "drop"
  )

# coniferous vs deciduous (plot)  <-- fixed join order
iv_leaf_plot <- iv_plot %>%
  left_join(species_class, by = "species") %>%
  filter(!is.na(leaf_type)) %>%
  group_by(plot, year, leaf_type) %>%
  summarise(
    IV = sum(IV, na.rm = TRUE),
    RA = sum(RA, na.rm = TRUE),
    RS = sum(RS, na.rm = TRUE),
    .groups = "drop"
  )

iv_leaf_plot

iv_leaf_plot %>% 
  ggplot(aes(x = leaf_type, y = IV)) + 
  geom_boxplot()

# convert to wide format
iv_leaf_plot_wide <- iv_leaf_plot %>%
  dplyr::select(plot, year, leaf_type, IV) %>%
  tidyr::pivot_wider(
    names_from  = leaf_type,
    values_from = IV,
    names_prefix = "IV_"
  )

iv_leaf_plot_wide

iv_leaf_plot_wide %>%
  ggplot(aes(x = IV_coniferous,
             y = IV_deciduous)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 1, slope = -1, linetype = "dashed") +
  coord_equal() +
  labs(
    x = "Importance value (Coniferous)",
    y = "Importance value (Deciduous)",
    title = "Leaf-type dominance at plot level"
  ) +
  theme_minimal()


### Tree heights by species -------------------------------------------------------
# get heights by species
dat_overlap %>% 
  dplyr::filter(hgt_est>0) %>% 
  mutate(species = fct_reorder(species, hgt_est, .fun = median, .desc = TRUE)) %>%
  ggplot(aes(x = species,
             y = hgt_est,
             fill = factor(year))) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,7.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# compare change in feights pioneer vs late
dat_overlap %>% 
  dplyr::filter(seral_stage != "other") %>%
  dplyr::filter(hgt_est>0) %>% 
  #mutate(species = fct_reorder(species, hgt_est, .fun = median, .desc = TRUE)) %>%
  ggplot(aes(x = seral_stage,
             y = hgt_est,
             fill = factor(year))) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,6)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##  Update management for overlapping sites: combine management ---------------------------------
#library(dplyr)

# as subplots names do not fit, link subplots by row_id - almost randomly 
# (as randomness was done in terrain - everyobe could choose where to start with first dot )

#  Define management variables
management_vars <- c(
  "clear", "grndwrk", "logging_trail", "planting", "anti_browsing"
)

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

#  Identify plots with both 2023 and 2025 entries for status == "both"
plots_with_both <- dat_overlap_mng %>%
  dplyr::filter(status == "both", year %in% c(2023, 2025)) %>%
  distinct(plot, year) %>%
  count(plot) %>%
  filter(n > 1) %>%
  pull(plot)

# Get 2023 & 2025 management data
mng_both23 <- dat_overlap_mng %>%
  filter(plot %in% plots_with_both, status == "both", year == 2023) %>%
  dplyr::select(plot, subplot, all_of(management_vars), year) %>%
  distinct() %>% 
  group_by(plot, year) %>%
  arrange(plot, subplot) %>%
  mutate(row_id = row_number()) 

head(mng_both23)

# somehow i have duplicated values - keep only the higher value 
# (presence, if both presence and absence are recorded on 
# same plot. subplots are linked by row_id)
mng_both25 <- dat_overlap_mng %>%
  filter(plot %in% plots_with_both, status == "both", year == 2025) %>%
  dplyr::select(plot, subplot, all_of(management_vars), year) %>%
 # distinct() %>% 
  group_by(plot, subplot, year) %>%
  summarise(across(all_of(management_vars), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>% 
  group_by(plot) %>%
  arrange(plot, subplot) %>%
  mutate(row_id = row_number()) 

nrow(mng_both25)
head(mng_both25)

# Step 4: Combine by row order and take max per variable
mng_both25_upd <- bind_rows(mng_both23, mng_both25) %>%
  group_by(plot, row_id) %>%
  summarise(
    subplot = subplot[year == 2025][1],  # use 2025 subplot name
    across(all_of(management_vars), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>% 
  select(-row_id) %>% 
  mutate(year = 2025)


# remove rows from 2025 from overall table, update management ---
dat_overlap_mng_25 <-
  dat_overlap_mng %>% 
  filter(status == 'both' & year == 2025) %>% 
  dplyr::select(-all_of(management_vars)) %>%
#  nrow()
  left_join(mng_both25_upd) %>% 
  dplyr::select(all_of(colnames(dat_overlap_mng))) # reorder columns in proper order
  
# Keep all other rows (i.e., everything except 2025 "both")
dat_overlap_mng_rest <- dat_overlap_mng %>%
  dplyr::filter(!(status == "both" & year == 2025))

# Bind updated 2025 rows back to the full dataset
dat_overlap_mng_upd <- bind_rows(dat_overlap_mng_rest, 
                                 dat_overlap_mng_25)


# Calculate intensity per plot & year 
# Collapse to one row per subplot (presence-based)
mng_subplot <- dat_overlap_mng_upd %>%
  dplyr::select(plot, subplot, year, n_subplots, all_of(management_vars)) %>%
  distinct() 

# plot based management intensity with updates 2025 values 
management_intensity <- mng_subplot %>%
  group_by(plot, year) %>%
  summarise(
    n_subplots = first(n_subplots),
    clear_intensity = sum(clear == 1, na.rm = TRUE) / n_subplots,
    grndwrk_intensity = sum(grndwrk == 1, na.rm = TRUE) / n_subplots,
    logging_trail_intensity = sum(logging_trail == 1, na.rm = TRUE) / n_subplots,
    planting_intensity = sum(planting == 1, na.rm = TRUE) / n_subplots,
    anti_browsing_intensity = sum(anti_browsing == 1, na.rm = TRUE) / n_subplots,
    .groups = "drop"
  )

# Optional: Join back to main table (if needed per-row)
dat_overlap_mng_upd2 <- dat_overlap_mng_upd %>%
  left_join(management_intensity, by = c("plot", "year", "n_subplots"))

## Management Subplot level -----------------------------------------------------------

df_mng_sub <- dat_overlap_mng_upd2 %>% 
  #dplyr::filter(year == "2023") %>% # keep management oionly frm 2023 for consistency
  distinct(plot, subplot,year,
           clear,
           grndwrk,
           logging_trail,
           planting,
           anti_browsing,
           
           # intensities = % of subplots with that activity
           clear_intensity,           
           grndwrk_intensity,         
           logging_trail_intensity,  
           planting_intensity ,       
           anti_browsing_intensity#,   
           # salvage_intensity,         
           # protection_intensity,      
           # management_intensity
           # 
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
df_master_mng_intensity <- dat_overlap_mng_upd2 %>% 
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



# Define the factor levels in correct order
intensity_levels <-  c("0–19", "20–39", "40–59", "60–79", "80–100")
low_classes <- c("0–19", "20–39")
applied_mng_intens_order <- c("0–19", "20–39","40–59", "60–79", "80–100") 


# prepare for binary classifucation  - need to group first, as I have two years
mng_sub_conv <- mng_intensity_props %>%
  mutate(
    intensity_binary = factor(intensity_binary, levels = c('no', 'yes')),
    intensity_class = factor(intensity_class, levels = intensity_levels)#,
    #activity = factor(activity, levels = applied_mng_intens_order)
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
 
 



# Generate new variables -----------------------------------------------------------
### Early vs late ( plot, subplot) ---------------
# 
# # Calculate stem counts by recovery type at the plot level
# share_early_vs_late <- 
#   dat_overlap_mng_upd2%>%
#   filter(!is.na(n)) %>%  # Optional: remove NAs if present
#   group_by(plot, seral_stage, time_snc_full_disturbance) %>%
#   summarise(n_stems = n(), .groups = "drop") %>%
#   pivot_wider(names_from = seral_stage,
#               values_from = n_stems,
#               values_fill = 0) %>% 
#   mutate(total = early + late,
#          share_early = early/total*100,
#          share_late = late/total*100) %>% 
#   select(plot, time_snc_full_disturbance, share_early, share_late) %>%
#   pivot_longer(cols = starts_with("share_"),
#                names_to = "seral_stage",
#                values_to = "share")# %>%
# 
# ### Early vs late: plot level  
# df_plot_share_early <-  
#   dat_overlap_mng_upd2%>%
#   filter(!is.na(n)) %>%  # Optional: remove NAs if present
#   group_by(plot,year, seral_stage, time_snc_full_disturbance) %>%
#   summarise(n_stems = n(), .groups = "drop") %>%
#   pivot_wider(names_from = seral_stage,
#               values_from = n_stems,
#               values_fill = 0) %>% 
#   mutate(total = early + late,
#          share_early = early/total*100,
#          share_late = late/total*100) %>% 
#   select(plot, year, share_early, share_late,time_snc_full_disturbance) 
# 
# ### early vs late : subplot level 
# df_sub_share_early <-  dat_overlap_mng_upd2%>%
#   filter(!is.na(n)) %>%  # Optional: remove NAs if present
#   group_by(subplot, plot,year, seral_stage, time_snc_full_disturbance) %>%
#   summarise(n_stems = n(), .groups = "drop") %>%
#   pivot_wider(names_from = seral_stage,
#               values_from = n_stems,
#               values_fill = 0) %>% 
#   mutate(total = early + late,
#          share_early = early/total*100,
#          share_late = late/total*100) %>% 
#   select(subplot, plot, year, share_early, share_late,time_snc_full_disturbance) 
# 
# 
# # Plot
# ggplot(share_early_vs_late, 
#        aes(x = factor(time_snc_full_disturbance),
#            y = share,
#            fill = seral_stage)) +
#   geom_boxplot() +
#   #geom_jitter() +
#   #geom_bar(stat = "identity", position = "stack") +
#   labs(x = "Time since stand replacing\ndisturbance (years)",
#        y = "Share of stems (%)",
#        fill = "Seral stage") 
# 

## Spruce share ( plot, subplot) ----------------------------------------------------------------------


# Calculate spruce share at the plot level
spruce_share_plot <- dat_overlap_mng_upd2 %>%
  filter(!is.na(n)) %>%  # Optional: remove NAs if present
  group_by(plot, year) %>%
  summarise(
    total_stems = n(),
    spruce_stems = sum(species == "piab"),
    spruce_share = spruce_stems / total_stems
  ) %>%
  select(plot, year, spruce_share)

hist(spruce_share_plot$spruce_share, breaks = 50)

# Calculate spruce share at the subplot level
spruce_share_sub <- dat_overlap_mng_upd2 %>%
  filter(!is.na(n)) %>%  # Optional: remove NAs if present
  group_by(plot, subplot, year) %>%
  summarise(
    total_stems = n(),
    spruce_stems = sum(species == "piab"),
    spruce_share = spruce_stems / total_stems,
    .groups = "drop"
  ) %>%
  select(plot, subplot, year, spruce_share)

hist(spruce_share_sub$spruce_share, breaks = 50)


## Traits: Weighted community mean (subplot & plot) ------------------

# choose  weighting : by stems (for early communities), by structure?
# Option A (default): weight by stem counts
wvar <- "n"

# Option B: weight by structure (uncomment ONE)
# wvar <- "basal_area_cm2"
# wvar <- "hgt_est"

# helper to pull a numeric weight safely
wfun <- function(x) ifelse(is.na(x) | x < 0, 0, x)

# Continue from here to replace dat_overlap_mng_upd2to dat_overlap_cleaned!
# Subplot × year CWMs 
cwm_subplot <- dat_overlap_mng_upd2 %>%
  #filter(species != 'ots1') %>% 
  mutate(w = wfun(.data[[wvar]])) %>%
  # keep only rows that contribute weight and have trait scores
  filter(w > 0) %>%
  group_by(plot, subplot, year, time_snc_full_disturbance) %>%
  summarise(
    stems_with_traits = sum(w[!is.na(Shade_tolerance) & !is.na(Drought_tolerance)], na.rm = TRUE),
    stems_total       = sum(w, na.rm = TRUE),
    CWM_shade   = ifelse(stems_with_traits > 0,
                         sum(w * Shade_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    CWM_drought = ifelse(stems_with_traits > 0,
                         sum(w * Drought_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    trait_coverage = stems_with_traits / pmax(stems_total, 1e-9),
    .groups = "drop"
  ) %>% 
  select(-stems_total, - stems_with_traits, -trait_coverage)

# 2) Plot × year CWMs
cwm_plot <- dat_overlap_mng_upd2%>%
  filter(species != 'ots1') %>% 
  mutate(w = wfun(.data[[wvar]])) %>%
  filter(w > 0) %>%
  group_by(plot, year, time_snc_full_disturbance) %>%
  summarise(
    stems_with_traits = sum(w[!is.na(Shade_tolerance) & !is.na(Drought_tolerance)], na.rm = TRUE),
    stems_total       = sum(w, na.rm = TRUE),
    CWM_shade   = ifelse(stems_with_traits > 0,
                         sum(w * Shade_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    CWM_drought = ifelse(stems_with_traits > 0,
                         sum(w * Drought_tolerance, na.rm = TRUE) / stems_with_traits, NA_real_),
    trait_coverage = stems_with_traits / pmax(stems_total, 1e-9),
    .groups = "drop"
  ) %>% 
  select(-stems_total, - stems_with_traits, -trait_coverage)

### quick plots on subplot and plot level 
ord_year <- function(x) factor(as.character(x), levels = c("2023","2025"))


# Long format for subplot
cwm_subplot_long <- cwm_subplot %>%
  transmute(level = "subplot",
            plot, subplot,
            time_snc_full_disturbance,
            year = ord_year(year),
            CWM_shade, CWM_drought) %>%
  pivot_longer(c(CWM_shade, CWM_drought),
               names_to = "trait", values_to = "CWM")

# Long format for plot
cwm_plot_long <- cwm_plot %>%
  transmute(level = "plot",
            plot,
            time_snc_full_disturbance,
            year = ord_year(year),
            CWM_shade, CWM_drought) %>%
  pivot_longer(c(CWM_shade, CWM_drought),
               names_to = "trait", values_to = "CWM")

# Combine
cwm_all_long <- bind_rows(cwm_subplot_long, cwm_plot_long) %>%
  filter(!is.na(CWM))

# Nice facet labels
trait_labs <- c(CWM_shade = "Shade tolerance", 
                CWM_drought = "Drought tolerance")

# ggplot(cwm_all_long, aes(x = year, y = CWM, fill = year)) +
#   geom_boxplot(width = 0.6, outlier.alpha = 0.15) +
#   geom_jitter(width = 0.1, alpha = 0.25, size = 0.7) +
#   facet_grid(level ~ trait, labeller = labeller(trait = trait_labs)) +
#   labs(x = "Year", y = "Community-weighted mean (CWM)",
#        title = "Trait CWMs by level and year") +
#   theme_classic2(base_size = 10)
# 
# cwm_all_long %>% 
#   ggplot(aes(x = time_snc_full_disturbance, 
#              y = CWM, 
#              fill = factor(time_snc_full_disturbance))) +
#   geom_boxplot(width = 0.6, outlier.shape = NA) +
#  # geom_jitter(width = 0.1, alpha = 0.25, size = 0.7) +
#   facet_grid(level ~ trait, 
#              labeller = labeller(trait = trait_labs),
#              scales = 'free_y') +
#   labs(x = "Year", y = "Community-weighted mean (CWM)",
#        title = "Trait CWMs",
#        subtitle = "by level and time since disturbance") +
#   theme_bw(base_size = 10) + 
#   theme(legend.position = "none")


# Create table on subplot and plot level ----------------------------------------
## Field data summary: subplot metrics ----------------------------------------------------------
field_sub_summ <- dat_overlap_mng_upd2%>%
  mutate(n = coalesce(n, 0L)) %>% # change NAs to 0 (0L = literrary 0 = integer (not double))
  # mutate(year = as.factor(year)) %>% 
  group_by(plot, subplot, year, 
           time_snc_full_disturbance, 
           time_snc_part_disturbance,
           disturbance_year, 
           forest_year, 
           disturbance_length,
           clear,
           grndwrk,
           logging_trail,
           planting,
           anti_browsing) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE),
    sp_richness = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp  = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
    # weighted mean height using only present stems
    mean_hgt    = if (stems_total > 0) weighted.mean(hgt_est[n > 0], n[n > 0], na.rm = T) else NA_real_,
    
    # compute weighted variance whenever we have ≥2 stems and any finite heights
    var_hgt = {
      if (stems_total > 1) {
        sel <- (n > 0) & is.finite(hgt_est)
        h   <- hgt_est[sel]; ww <- n[sel]
        if (length(h) >= 1 && sum(ww) > 1) {
          mu <- weighted.mean(h, ww)
          v  <- sum(ww * (h - mu)^2) / (sum(ww) - 1)  # freq-weighted, Bessel corrected
          if (is.nan(v)) NA_real_ else v              # will be 0 if all h are equal
        } else NA_real_
      } else NA_real_
    },
    
    cv_hgt = if (is.finite(mean_hgt) && mean_hgt > 0 && !is.na(var_hgt))
      sqrt(var_hgt) / mean_hgt else
        if (stems_total > 1 && is.finite(mean_hgt) && mean_hgt > 0) 0 else NA_real_,
    .groups = "drop"  #,
    #range_hgt = if (sum(n > 0) >= 2) diff(range(hgt_est[n > 0], na.rm = TRUE)) else NA_real_,
    
  )  %>% 
  left_join(cwm_subplot)

  
# Problem: issue that some sites have different management from 2023 - need to keep management 
  # from 2023

field_sub_summ %>% 
  filter(subplot == "628_T2_AH_20250827") %>% 
  str()

safe_max <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    max(x, na.rm = TRUE)
  }
}

field_sub_summ_cleaned <- field_sub_summ %>%
  group_by(plot, subplot, year) %>%
  summarise(
    across(
      where(is.numeric),
      safe_max
    ),
    .groups = "drop"
  )

field_sub_summ_cleaned %>% 
  filter(subplot == "628_T2_AH_20250827") %>% 
  str()



# is teh change in community shading/drought tolerance driven by planting????
# Plot: Shade ~ Time since disturbance, by planting
x_lab_time_snc_full_dist = "Time since stand\nreplacing disturbance (years)"
# 
# #### Effect of planting? 
# p_shade_planting <- field_sub_summ_cleaned%>% 
#   ggplot(aes(x = as.factor(time_snc_full_disturbance),
#              y = CWM_shade,
#              fill = factor(planting))) +
#  #
#   geom_boxplot() +
#   labs(x = x_lab_time_snc_full_dist,
#        y = "CWM shade",
#        fill = "Planting") +
#   theme_classic2() +
#   theme(text  = element_text(size = 10))
# 
# # Plot: Drought ~ Time since disturbance, by planting
# p_shade_drought <- field_sub_summ_cleaned%>% 
#   ggplot(aes(x = as.factor(time_snc_full_disturbance),
#              y = CWM_drought,
#              fill = factor(planting))) +
#   geom_boxplot(outlier.shape = NA) +
#   labs(x = x_lab_time_snc_full_dist,
#        y = "CWM drought",
#        fill = "Planting") +
#   theme_classic2()+
#   theme(text  = element_text(size = 10))
# 
# # Plot: Shade ~ Planting
# p_shade_total <- field_sub_summ_cleaned%>% 
#   ggplot(aes(x = factor(planting),
#              y = CWM_shade,
#              fill = factor(planting))) +
#   #
#   geom_boxplot(outlier.shape = NA) +
#   stat_compare_means(method = "wilcox.test", label = "p.format", 
#                      label.y = 4, size = 3) +
#   labs(x = "Planting",
#        y = "CWM shade",
#        fill = "Planting") +
#   theme_classic2() +
#   theme(text  = element_text(size = 10))
# 
# # Plot: Drought ~ Planting
# p_drought_total <- field_sub_summ_cleaned%>% 
#   ggplot(aes(x = factor(planting),
#              y = CWM_drought,
#              fill = factor(planting))) +
#   geom_boxplot() +
#   stat_compare_means(method = "wilcox.test", label = "p.format", 
#                      label.y = 4, size = 3) +
#   labs(x = "Planting",
#        y = "CWM drought",
#        fill = "Planting") +
#   theme_classic2() +
#   theme(text  = element_text(size = 10))
# 
# # Arrange all plots
# annotate_figure(
#   ggarrange(p_shade_planting, p_shade_total,
#             p_shade_drought, p_drought_total,
#             ncol = 2, nrow = 2,
#             common.legend = TRUE, legend = "bottom"),
#   top = text_grob("Subplot level", 
#                   face = "bold", size = 12)
# )
# 


#### Subplot quick plotting: all vars ---------------------
df_sub_long <- field_sub_summ_cleaned%>%
  dplyr::filter(stems_total > 0) %>% 
  filter(cv_hgt >0) %>% 
  #ungroup() %>%
  select(year,
         clear,
         grndwrk,
         logging_trail,
         planting,
         anti_browsing,
         time_snc_full_disturbance, time_snc_part_disturbance, 
         stems_total,
         mean_hgt, cv_hgt, shannon_sp, sp_richness,
         CWM_shade ,
         CWM_drought ) %>%
  pivot_longer(-c(year,
                  time_snc_full_disturbance,
                  time_snc_part_disturbance,
                  #CWM_shade ,
                  #CWM_drought,
                  clear,
                  grndwrk,
                  logging_trail,
                  planting,
                  anti_browsing),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value))# %>%



# Drivers -----------------------------------------------------------------

## Summary stats on plot level stem density per m² --------------------
area_subplot_m2 <- 4      # 4 m²
area_plot_m2    <- 5*4    # 20 m²



## --- Plot-level metrics (aggregate over subplots) ---
plot_metrics_mean <- field_sub_summ_cleaned%>%
  group_by(plot, year) %>%
  summarise(
    mean_sp_richness = mean(sp_richness, na.rm = TRUE),
    #var_sp_richness  = var(sp_richness,  na.rm = TRUE),
    mean_shannon_sp  = mean(shannon_sp,  na.rm = TRUE),
    mean_evenness_sp = mean(evenness_sp, na.rm = TRUE),
    mean_mean_hgt    = mean(mean_hgt,    na.rm = TRUE),
    mean_cv_hgt      = mean(cv_hgt,      na.rm = TRUE),
    mean_eff_numb    = mean(effective_numbers, na.rm = TRUE),    
    #var_cv_hgt       = var(cv_hgt,       na.rm = TRUE),
    #mean_range_hgt   = mean(range_hgt,   na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(cwm_plot, by = join_by(plot, year)) %>% 
  left_join(df_mng_plot,by = join_by(plot, year)) %>% 
  left_join(spruce_share_plot, by = join_by(plot, year)) %>%
  left_join(iv_max_plot, by = join_by(plot, year)) #%>% 


# --- pooled CV directly from dat23_subplot_recode ---
plot_metrics_pooled  <- dat_overlap_mng_upd2 %>%
  mutate(n = coalesce(n, 0L)) %>%                # no NA counts
  group_by(plot, year) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE),
    sp_richness = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp  = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    mean_hgt    = if (stems_total > 0) weighted.mean(hgt_est, n, na.rm = TRUE) else NA_real_,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
    # compute weighted variance whenever we have ≥2 stems and any finite heights
    var_hgt = {
      if (stems_total > 1) {
        sel <- (n > 0) & is.finite(hgt_est)
        h   <- hgt_est[sel]; ww <- n[sel]
        if (length(h) >= 1 && sum(ww) > 1) {
          mu <- weighted.mean(h, ww)
          v  <- sum(ww * (h - mu)^2) / (sum(ww) - 1)  # freq-weighted, Bessel corrected
          if (is.nan(v)) NA_real_ else v              # will be 0 if all h are equal
        } else NA_real_
      } else NA_real_
    },
    
    cv_hgt = if (is.finite(mean_hgt) && mean_hgt > 0 && !is.na(var_hgt))
      sqrt(var_hgt) / mean_hgt else
        if (stems_total > 1 && is.finite(mean_hgt) && mean_hgt > 0) 0 else NA_real_,
    .groups = "drop" 
  ) %>% 
  left_join(cwm_plot, by = join_by(plot, year)) %>% 
  left_join(df_mng_plot,by = join_by(plot, year)) %>% 
  left_join(spruce_share_plot, by = join_by(plot, year)) %>%
  left_join(iv_max_plot, by = join_by(plot, year)) #%>%  

#mutate(cv_hgt = ifelse(is.na(cv_hgt), 0L, cv_hgt)) # replace NA by 0 if stems are missing

# get only context and disturbance information on plot level
df_plot_context <- dat_overlap_mng_upd2%>%
  dplyr::select(plot,year,  
                pre_dist_trees_n,  
                area_m2, 
                pre_dist_dens_ha, 
                time_snc_full_disturbance,
               # time_snc_part_disturbance,
                disturbance_year, 
                forest_year, 
                disturbance_length,
               x, y) %>% 
  group_by(plot, year) %>% 
  mutate(x = mean(x, na.rm =T),
            y = mean(y, na.rm =T)) %>% 
  ungroup(.) %>% 
  distinct() 



## create final table for both levels -----------------------------------------------------
# 1) Subplot table (has subplot mean_hgt and stems_total as weights)
sub_df <- field_sub_summ_cleaned %>%
 # filter(stems_total > 0) %>% #, cv_hgt > 0
  transmute(
    ID       = subplot,
    plot_id  = plot, #str_replace(subplot, "^[^_]+_([^_]+_[^_]+)_.*$", "\\1"),
    year     = year,
    level    = "subplot",
    dens_m2  = stems_total / area_subplot_m2,     # 4 m² subplot
    cv_hgt   = cv_hgt,
    mean_hgt = mean_hgt,
    sp_richness = sp_richness,
    shannon_sp  = shannon_sp ,
    evenness_sp = evenness_sp ,
    effective_numbers = effective_numbers,
    w        = stems_total
  ) %>% 
  mutate(dens_ha = dens_m2*10000) %>% 
  left_join(df_plot_context, by = c("plot_id" = "plot",
                                    "year" = "year")) %>% 
  left_join(cwm_subplot, by = c("plot_id" = "plot",
                                "year" = "year",
                                "ID" = 'subplot',
                                "time_snc_full_disturbance" = "time_snc_full_disturbance") ) %>% 
  left_join(df_mng_plot, by = c("plot_id" = "plot",
                                       "year" = "year")) %>%  # need to separate here by year?
  left_join(spruce_share_sub, by = c("plot_id" = "plot",
                                     "year" = "year",
                                     "ID" = 'subplot')) %>%
  left_join(iv_max_sub, by = c("plot_id" = "plot",
                                     "year" = "year",
                                     "ID" = 'subplot')) %>%  
  distinct()

nrow(sub_df)

names(sub_df)
# 2) Plot table (pooled metrics already computed)
plot_df <- plot_metrics_pooled %>%
  transmute(
    ID       = plot,
    plot_id  = plot,
    year     = year,
    level    = "plot",
    dens_m2  = stems_total / area_plot_m2,    # 5×4 m² = 20 m²
    cv_hgt   = cv_hgt,
    mean_hgt = mean_hgt,
    sp_richness = sp_richness,
    shannon_sp  = shannon_sp ,
    evenness_sp = evenness_sp ,
    effective_numbers = effective_numbers,
    w        = stems_total
  ) %>%
  #filter(!is.na(cv_hgt), cv_hgt > 0) %>% 
  mutate(dens_ha = dens_m2*10000,
         mean_hgt = replace_na(mean_hgt, 0)) %>% # Replace NA with 0)
  left_join(df_plot_context, by = c("plot_id" = "plot",
                                    "year" = "year")) %>% 
  left_join(cwm_plot, by = c("plot_id" = "plot",
                             "year" = "year",
                             "time_snc_full_disturbance" = "time_snc_full_disturbance") ) %>% 
  left_join(df_mng_plot, by = c("plot_id" = "plot",
                                       "year" = "year")) %>% 
  left_join(spruce_share_plot, by = c("plot_id" = "plot",
                                       "year" = "year")) %>% 
  left_join(iv_max_plot, by = c("plot_id" = "plot",
                               "year" = "year")) %>%  
  distinct() # keep only unique rows


names(plot_df)
nrow(plot_df)

##  Export data for Karim for AEF testing, 20260212 -------------------------

# add IV coniferous vs deciduous
plot_df_AEF <- plot_df %>% 
  dplyr::select(-level, -ID, -w) %>% 
  mutate(
    year = as.integer(year),           # from dbl -> int
    plot_id = as.character(plot_id)    # keep key as character
  ) %>%
  rename(plot = plot_id) %>% 
  left_join(iv_leaf_plot_wide, by = join_by(plot, year)) %>% 
  mutate(
    # cap time since disturbance
    time_snc_full_disturbance = pmin(time_snc_full_disturbance, 8),
    plot_id = factor(plot),
 #   w       = pmin(pmax(w, 1), 50),   # cap weights so a few dense plots don't dominate
    year_f = factor(year)
  ) %>% 
  mutate(
    across(all_of(c("cv_hgt", "mean_hgt", "sp_richness", "effective_numbers")), ~ replace_na(.x, 0))
  ) #%>%


plot_sf_AEF <- plot_df_AEF %>%
  sf::st_as_sf(coords = c("x", "y"), crs = 3035, remove = FALSE)

# sf::st_write(
#   plot_sf_AEF,
#   dsn = 'outDataShare/Karim_AEF/regeneration_chars_3035.gpkg',
#   layer = "plot",
#   delete_layer = TRUE  # overwrite layer if it already exists
# )



#### Bind & clean final table with both levels ----------
both_levels_re2 <- bind_rows(sub_df, plot_df) %>%
    mutate(
      # cap time since disturbance
      time_snc_full_disturbance = pmin(time_snc_full_disturbance, 8),
      level   = factor(level, levels = c("subplot","plot")),
      plot_id = factor(plot_id),
      w       = pmin(pmax(w, 1), 50),   # cap weights so a few dense plots don't dominate
      year_f = factor(year)
  ) %>% 
  mutate(
    across(all_of(c("cv_hgt", "mean_hgt", "sp_richness", "effective_numbers")), ~ replace_na(.x, 0))
  ) %>%
 # add small value to CV = 0
  mutate(
    cv_hgt_pos = ifelse(cv_hgt <= 0, 1e-4, cv_hgt),
    effective_numbers = ifelse(effective_numbers == 0, 1e-4, effective_numbers)#is.na(cv_hgt) | 
  ) %>%
  mutate(cv_hgt_present = as.integer(cv_hgt > 0)) %>% 
  mutate(active_regeneration =
           (planting_intensity + anti_browsing_intensity)/2) #%>% 
  
  #filter(!is.na(mean_hgt), !is.na(w), !is.na(cv_hgt)) #, !is.na(dens_m2)


table(both_levels_re2$cv_hgt_present)


# get data for plot and subplot 
df_plot_clean <- both_levels_re2 %>% 
  filter( level == 'plot')

df_sub_clean <- both_levels_re2 %>% 
  filter( level == 'subplot')






#### Test correlation between planting and anti-browsing intensity ---------------

# does planting and antibrowsing cooccurs?
both_levels_re2 %>% 
  dplyr::filter(level == 'plot') %>% 
  ggplot(aes(x = planting_intensity, 
             y = anti_browsing_intensity,
             color = factor(year),
             fill = factor(year))) +
  geom_jitter() +
  geom_smooth(method = "lm", se = TRUE) 


##### how does spruce share behaves? ------------------------------------------------

# 1. Summarise spruce cover by intensity grid
spruce_summary <- both_levels_re2 %>%
  filter(level == "plot") %>%
  group_by(planting_intensity, anti_browsing_intensity) %>%
  summarise(spruce_share = mean(spruce_share, na.rm = TRUE),
            n = n(),
            .groups = "drop")

ggplot(spruce_summary, aes(x = planting_intensity, y = anti_browsing_intensity)) +
  geom_tile(aes(fill = spruce_share)) +
 # geom_point(aes(size = n), shape = 21, fill = "black", color = "white", stroke = 0.3) +
  scale_fill_viridis_c(option = "magma", 
                       begin = 0.2,
                       direction = -1,
                       name = "Spruce share", limits = c(0, 1)) +
  scale_size_continuous(range = c(2, 10), name = "Plots [#]") +
  labs(
    x = "Planting intensity",
    y = "Anti-browsing intensity"
  ) +
  theme_classic(base_size = 10) +
  guides(
    fill = guide_colourbar(barwidth = 0.5, barheight = 8),
    size = guide_legend(ncol = 1)
  ) +
  theme(
    legend.position = "right",
    legend.box = "horizontal"
  )







### Spearman correlation between management intensities using base R ----------------------------------
cor.test(both_levels_re2$planting_intensity,
         both_levels_re2$anti_browsing_intensity,
         method = "spearman")

# Count combinations: any planting without protection or vice versa?
both_levels_re2 %>%
  mutate(
    planted = planting_intensity > 0,
    protected = anti_browsing_intensity > 0
  ) %>%
  count(planted, protected)


# planted protected     n
# <lgcl>    <lgcl> <int>
#   1:   FALSE     FALSE    70
# 2:   FALSE      TRUE     7
# 3:    TRUE     FALSE    70
# 4:    TRUE      TRUE   186

# # A tibble: 4 × 3
# planted protected     n
# <lgl>   <lgl>     <int>
#   1 FALSE   FALSE       545
# 2 FALSE   TRUE         73
# 3 TRUE    FALSE       781
# 4 TRUE    TRUE       1978




##### Boxplot: which management is most important?  --------------------
# Create long table for plotting across factors
both_levels_long <- both_levels_re2 %>%
  filter(level == 'plot') %>% 
  select(plot_id, level,
         cv_hgt, mean_hgt, sp_richness, effective_numbers,
         planting_intensity,
         clear_intensity,
         anti_browsing_intensity,
         grndwrk_intensity
         ) %>%
  pivot_longer(
    cols = c( mean_hgt, cv_hgt, effective_numbers, sp_richness),
    names_to = "variable",
    values_to = "value"
  )

# check which one of management activities is more important for species richness 
# and vertical structure?
p1<-ggplot(both_levels_long %>% filter(level == "plot"),
       aes(x = factor(planting_intensity) , y = value, fill = factor(planting_intensity) )) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),  # nicely formatted p-values
    method = "kruskal.test",                     # or "t.test" depending on data
    label.y.npc = 0.95,                         # position near top of panel
    size = 3
  ) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  theme_bw() + 
  theme(legend.position = 'none')

p2<-ggplot(both_levels_long %>% filter(level == "plot"),
       aes(x = factor(anti_browsing_intensity) , y = value, fill = factor(anti_browsing_intensity) )) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),  # nicely formatted p-values
    method = "kruskal.test",                     # or "t.test" depending on data
    label.y.npc = 0.95,                         # position near top of panel
    size = 3
  ) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +  
  theme_bw() +
  theme(legend.position = 'none')

p3<-ggplot(both_levels_long %>% filter(level == "plot"),
           aes(x = factor(grndwrk_intensity ) , y = value, fill = factor(grndwrk_intensity ) )) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),  # nicely formatted p-values
    method = "kruskal.test",                     # or "t.test" depending on data
    label.y.npc = 0.95,                         # position near top of panel
    size = 3
  ) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(legend.position = 'none')


ggarrange(p1,p2,p3, ncol = 1, nrow = 3)


#### GAM management interaction + time since disturbance smooths --------------



##### Mean height -----------------------------------------------------

# test linear model
gam_mean_hgt_plot <- gam(
  mean_hgt ~ 
    planting_intensity*anti_browsing_intensity +
     s(time_snc_full_disturbance, k = 4) +
    grndwrk_intensity+
    year_f +
    # level +
    s(plot_id, bs = "re"),
  data = df_plot_clean,
  family = tw(link = "log"),
  method = "REML"
)


# run on botj plot and sublot level
# test linear model
gam_mean_hgt_both <- gam(
  mean_hgt ~ 
    planting_intensity*anti_browsing_intensity +
      s(time_snc_full_disturbance, k = 4) +
    grndwrk_intensity+
    year_f +
     level +
    s(plot_id, bs = "re"),
  data = both_levels_re2 %>% dplyr::filter(mean_hgt < 6), #both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)


pred_hgt_level <- ggpredict(
  gam_mean_hgt_both,
  terms = c("time_snc_full_disturbance[0:8]", "level"),   # both levels
  #terms = c("time_snc_full_disturbance")   # both levels
  type = 'fixed'
)

plot(pred_hgt_level) + theme_classic()


# do i have any influential values???check raw data!
both_levels_re2 %>%
  dplyr::group_by(time_snc_full_disturbance) %>%
  dplyr::summarise(
    mean_hgt = mean(mean_hgt, na.rm = TRUE),
    n = dplyr::n()
  )

ggplot(both_levels_re2, aes(time_snc_full_disturbance, mean_hgt)) +
  geom_point(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 3, colour = "red") +
  theme_classic()


# test if level has any effect
gam_mean_hgt_both_level <- gam(
  mean_hgt ~ 
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, by = level, k = 4) +
    grndwrk_intensity+
    year_f +
    level +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)

gam_mean_hgt_sub <- gam(
  mean_hgt ~ 
    planting_intensity*anti_browsing_intensity +
     s(time_snc_full_disturbance, by = year_f, k = 4) +
    grndwrk_intensity+
    year_f +
    #level +
    s(plot_id, bs = "re"),
  data = df_sub_clean %>% dplyr::filter(mean_hgt < 6), #df_sub_clean,
  family = tw(link = "log"),
  method = "REML"
)

pred_hgt_level <- ggpredict(
  gam_mean_hgt_plot,
  terms = c("time_snc_full_disturbance")   # both levels
)

plot(pred_hgt_level) + theme_classic()


###### chceck for spatial autocorrelation  ------------------

# Extract residuals
df_sp <- df_plot_clean %>%
  #dplyr::filter(level == "plot") %>%
  dplyr::mutate(resid = residuals(gam_mean_hgt_plot))

# Create spatial object
coords <- as.matrix(df_sp[, c("x", "y")])

# k-nearest neighbours (choose k ~ 4–6 typically)
nb <- spdep::knn2nb(spdep::knearneigh(coords, k = 5))
lw <- spdep::nb2listw(nb, style = "W")

# Moran's I
spdep::moran.test(df_sp$resid, lw)






##### CV hgt bin  ----

# check raw data

both_levels_re2 %>%
  dplyr::group_by(time_snc_full_disturbance) %>%
  dplyr::summarise(
    cv_hgt = mean(cv_hgt, na.rm = TRUE),
    n = dplyr::n()
  )

ggplot(both_levels_re2, aes(time_snc_full_disturbance, cv_hgt)) +
  geom_point(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 3, colour = "red") +
  theme_classic()


# does CV increasese with increasing heights?
ggplot(both_levels_re2, aes(mean_hgt, cv_hgt_pos)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess") +
  theme_classic()


gam_cv_hgt_bin_plot <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    planting_intensity*anti_browsing_intensity +
    #ti(planting_intensity, anti_browsing_intensity) +
    grndwrk_intensity +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_plot_clean,
  method = "REML",
  family = binomial(link = "logit")
)

gam_cv_hgt_bin_sub <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    planting_intensity*anti_browsing_intensity +
    grndwrk_intensity +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_sub_clean,
  method = "REML",
  family = binomial(link = "logit")
)

gam_cv_hgt_bin_both <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    planting_intensity*anti_browsing_intensity +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = binomial(link = "logit")
)


pred <- ggpredict(
  gam_cv_hgt_bin_both,
  terms = c("time_snc_full_disturbance[0:8]", "level"),   # both levels
  #terms = c("time_snc_full_disturbance")   # both levels
  type = 'fixed'
)

plot(pred) + theme_classic()


gam_cv_hgt_bin_both_level <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, by = level, k = 7) +
    planting_intensity*anti_browsing_intensity +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = binomial(link = "logit")
)
AIC(gam_cv_hgt_bin_both_level, gam_cv_hgt_bin_both)

##### CV Height (Positive Part) -----------
gam_cv_hgt_pos_sub <- gam(
  cv_hgt_pos ~ planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity +
    #level +
    year_f +
   s(plot_id, bs = "re"),
  data = df_sub_clean,
  method = "REML",
  family = tw(link = "log")
)

gam_cv_hgt_pos_plot <- gam(
  cv_hgt_pos ~ planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_plot_clean,
  method = "REML",
  family = tw(link = "log")
)

gam_cv_hgt_pos_both <- gam(
  cv_hgt_pos ~ planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

# check teh effect of heigt - but not for plotting
gam_cv_hgt_pos_both_hgt <- gam(
  cv_hgt_pos ~ planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 5) +
    s(mean_hgt, k = 5) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)


gam_cv_hgt_pos_both_log <- mgcv::gam(
  log(cv_hgt_pos) ~ 
    planting_intensity * anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 5) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = gaussian()
)


pred_cv <- ggpredict(
  gam_cv_hgt_pos_both_log,
  terms = c("time_snc_full_disturbance[0:7]", "level"),
  type = "fixed"
)

pred_cv$predicted <- exp(pred_cv$predicted)
pred_cv$conf.low <- exp(pred_cv$conf.low)
pred_cv$conf.high <- exp(pred_cv$conf.high)

plot(pred_cv) + theme_classic()


pred_cv_log <- ggpredict(
  gam_cv_hgt_pos_both_log,
  terms = c("time_snc_full_disturbance [0:8]", "level"),
  type = "fixed"
)

plot(pred_cv_log) +
  theme_classic() +
  labs(
    y = "Predicted log(CV height)"
  )



pred_cv_time2 <- ggpredict(
  gam_cv_hgt_pos_both,
  terms = c("time_snc_full_disturbance[0:6]", "level"),
  type = "fixed"
)
plot(pred_cv_time2)




#fin.m.cv.pos <- gam_cv_hgt_pos_base_plot_lin #gam_cv_hgt_pos_base_plot #gam_cv_hgt_pos_re

##### Species Diversity (Effective Number) -----------


gam_eff_plot <- gam(
  effective_numbers ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    #planting_intensity, anti_browsing_intensity) +
    grndwrk_intensity +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_plot_clean, #both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_eff_sub <- gam(
  effective_numbers ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_sub_clean, #both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_eff_both <- gam(
  effective_numbers ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 3) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)


gam_eff_both_level <- gam(
  effective_numbers ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 5) +
    s(time_snc_full_disturbance, by = level, k = 5) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)
AIC(gam_eff_both_level, gam_eff_both)

pred_eff_time <- ggpredict(
  #gam_eff_both_level,
  gam_eff_both,
  terms = c("time_snc_full_disturbance", "level")
)
plot(pred_eff_time)

summary(gam_eff_both)

#fin.m.eff <- gam_effective_intensity_re_int1_plot_lin #gam_effective_intensity_re_int1_plot #gam_effective_intensity_re_lev_int1


##### Species Richness -----------
gam_rich_both <- gam(
  sp_richness ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = nb(link = "log")
)


gam_rich_both_level <- gam(
  sp_richness ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    s(time_snc_full_disturbance, by = level, k = 7) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = nb(link = "log")
)
summary(gam_rich_both)
summary(gam_rich_both_level)

AIC(gam_rich_both_level, gam_rich_both)

gam_rich_both_pois <- gam(
  sp_richness ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = poisson(link = "log"),
)


gam_rich_plot <- gam(
  sp_richness ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity +
   # level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_plot_clean,
  method = "REML",
  family = nb(link = "log")
)

gam_rich_sub <- gam(
  sp_richness ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_intensity +
    # level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_sub_clean,
  method = "REML",
  family = nb(link = "log")
)


pred_rich_time <- ggpredict(
  gam_rich_both,
  terms = c("time_snc_full_disturbance", "level")
)
plot(pred_rich_time)


#fin.m.rich <- gam_rich_both #gam_richness_intensity_nb


##### Species turnover: alpha, gamma and beta (species turnover) diversity ---------------------------
# Jaccard dissimilairty - == beta diversity/turnover
# long -> wide community matrix at subplot level
comm_sub <- dat_overlap_mng_upd2 %>%
  mutate(n = coalesce(n, 0L)) %>%
  filter(!is.na(species), species != "") %>%
  group_by(plot, year, subplot, species) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = n,
    values_fill = 0
  )

beta_within_plotyear <- function(df_plotyear) {
  meta <- df_plotyear %>% dplyr::select(plot, year, subplot)
  mat  <- df_plotyear %>% dplyr::select(-plot, -year, -subplot) %>% as.data.frame()
  
  # drop empty subplots (all zeros)
  keep <- rowSums(mat) > 0
  mat  <- mat[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  
  n_sub <- nrow(mat)
  
  # if <2 subplots with any stems, turnover is undefined
  if (n_sub < 2) {
    return(tibble::tibble(
      plot = unique(meta$plot),
      year = unique(meta$year),
      n_subplots_nonempty = n_sub,
      beta_jaccard_mean = NA_real_,
      beta_bray_mean = NA_real_,
      alpha_rich_mean = if (n_sub == 1) sum(mat[1, ] > 0) else NA_real_,
      gamma_rich = if (n_sub == 1) sum(colSums(mat) > 0) else NA_real_,
      beta_whittaker = NA_real_
    ))
  }
  
  # Presence/absence for Jaccard
  mat_pa <- (mat > 0) * 1
  
  # Pairwise dissimilarities
  d_jac  <- vegan::vegdist(mat_pa, method = "jaccard", binary = TRUE)
  d_bray <- vegan::vegdist(mat,    method = "bray")
  
  # Richness alpha/gamma and Whittaker beta
  alpha <- rowSums(mat_pa)
  gamma <- sum(colSums(mat_pa) > 0)
  betaW <- gamma / mean(alpha)
  
  tibble::tibble(
    plot = unique(meta$plot),
    year = unique(meta$year),
    n_subplots_nonempty = n_sub,
    beta_jaccard_mean = mean(as.numeric(d_jac),  na.rm = TRUE),   # 0..1 (higher = more turnover)
    beta_bray_mean    = mean(as.numeric(d_bray), na.rm = TRUE),   # 0..1 (higher = more turnover)
    alpha_rich_mean   = mean(alpha),
    gamma_rich        = gamma,
    beta_whittaker    = betaW                                     # >=1 (higher = more turnover)
  )
}


beta_plotyear <- comm_sub %>%
  group_by(plot, year) %>%
  group_split() %>%
  purrr::map_dfr(beta_within_plotyear)

beta_plotyear


#beta_jaccard_mean (0–1):
# ~0 = subplots share the same species (homogenized)
# ~1 = subplots share few/no species (high turnover)
# 
# beta_bray_mean (0–1):
#   like above, but accounts for abundance (stem counts)

# beta_whittaker (≥1):
#   1 = plot richness equals mean subplot richness (no turnover)
# 2.4 (like you estimated) = plot has ~2.4× the richness of an average subplot

plot_df_AEF2 <- plot_df_AEF %>% 
  left_join(beta_plotyear)

# small adjustment to use beta family
eps <- 0.001
plot_df_AEF2$beta_jaccard_mean <- 
  pmin(pmax(plot_df_AEF2$beta_jaccard_mean, eps), 1 - eps)


# quick check over time
ggplot(plot_df_AEF2, 
       aes(x = time_snc_full_disturbance, 
           y = beta_jaccard_mean)) +
  geom_jitter(alpha = 0.4) +
  geom_smooth(method = "gam", 
              formula = y ~ s(x, k = 7)) +
  theme_classic()


cor.test(plot_df_AEF2$time_snc_full_disturbance,
         plot_df_AEF2$beta_jaccard_mean,
         method = "spearman")

cor.test(plot_df_AEF2$planting_intensity,
         plot_df_AEF2$spruce_share,
         method = "spearman")

plot_df_AEF2 %>%
  dplyr::group_by(time_snc_full_disturbance) %>%
  dplyr::summarise(
    beta_jaccard_mean = mean(beta_jaccard_mean, na.rm = TRUE),
    n = dplyr::n()
  )

ggplot(plot_df_AEF2, aes(time_snc_full_disturbance, beta_jaccard_mean)) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 3, colour = "red") +
  theme_classic()


m_beta_add <- mgcv::gam(
  beta_jaccard_mean ~ 
    time_snc_full_disturbance +
    planting_intensity+anti_browsing_intensity +
    grndwrk_intensity+
    year_f +
    s(plot_id, bs = "re"),
  data = plot_df_AEF2,
  method = "REML"
)


m_beta_int <- mgcv::gam(
  beta_jaccard_mean ~ 
    time_snc_full_disturbance +
    planting_intensity*anti_browsing_intensity +
    grndwrk_intensity+
    year_f +
    s(plot_id, bs = "re"),
  data = plot_df_AEF2,
  method = "REML"
)

AIC(m_beta_add, m_beta_int)
appraise(m_beta_add)
summary(m_beta_add)

m_beta_betar <- gam(
  beta_jaccard_mean ~ 
    time_snc_full_disturbance +
    planting_intensity*anti_browsing_intensity +
    grndwrk_intensity+
    year_f +
    s(plot_id, bs = "re"),
  family = betar(link = "logit"),
  data = plot_df_AEF2
)
summary(m_beta_betar)
appraise(m_beta_betar)

pred <- ggpredict(
  m_beta_add,
  terms = c("time_snc_full_disturbance[0:8]"),   # both levels
  #terms = c("time_snc_full_disturbance")   # both levels
  type = 'fixed',
  exclude = "s(plot_id)"
)

plot(pred) + theme_classic()


pred_plant <- ggpredict(
  m_beta_add,
  terms = "planting_intensity [0:1 by=0.05]",
  exclude = "s(plot_id)"
)

pred_time <- ggpredict(
  m_beta_add,
  terms = "time_snc_full_disturbance [all]",
  exclude = "s(plot_id)"
)

p_time_beta <- ggplot(pred_time, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(x = "Time since full disturbance (years)",
       y = "Predicted Jaccard beta") +
  theme_classic() 



ggarrange(p_plant_beta, p_spruce_beta,  p_time_beta,
          ncol =3)
summary(m_beta_add)
AIC(m_beta_full, m_beta_add)
appraise(m_beta_full)
plot(m_beta_add, page = 1)

concurvity(m_beta_add, full = TRUE)



# check correlation between planting and anti-browsing
cor(plot_df_AEF2$planting_intensity,
    plot_df_AEF2$anti_browsing_intensity,
    use = "complete.obs")
#[1] 0.6491358

cor.test(plot_df_AEF2$planting_intensity,
         plot_df_AEF2$anti_browsing_intensity,
         use = "complete.obs")










##### Summaries and Export -----------
all.models <- list(
  #gam_mean_hgt_plot, 
  gam_mean_hgt_both, 
  #gam_mean_hgt_sub,
  # gam_cv_hgt_bin_plot, 
  # gam_cv_hgt_bin_sub, 
  # gam_cv_hgt_bin_both,
  # gam_cv_hgt_pos_sub, 
  # gam_cv_hgt_pos_plot, 
  gam_cv_hgt_pos_both,
  # gam_eff_plot, 
  # gam_eff_sub, 
  gam_eff_both,
  gam_rich_both,
  m_beta_add   # jaccard beta dissimilarity
  # gam_rich_plot, 
  # gam_rich_sub
)

lapply(all.models, summary)
#lapply(all.models, appraise)


fin.models <- list(
  hgt   = gam_mean_hgt_both,
 # cvbin = fin.m.cv.bin,
  cvpos = gam_cv_hgt_pos_both,
  eff   = gam_eff_both,
  rich  = gam_rich_both,
 beta = m_beta_add
)

lapply(fin.models, summary)

# Plots ----------------------------------------------------------------------------
##  Time since disturbnace cross-scale ---------------------------------------



pp <- function(model, terms, xlab = NULL, ylab = NULL,
               annot = NULL, scale_y = 1,
               annot_x = 3.5, annot_y = NULL) {
  
  pr <- ggpredict(model, 
                  terms = terms,
                  exclude = "s(plot_id)")
  pr <- as.data.frame(pr)
  
  pr$predicted <- pr$predicted * scale_y
  pr$conf.low  <- pr$conf.low  * scale_y
  pr$conf.high <- pr$conf.high * scale_y
  
  p <- ggplot(pr, aes(x = x, y = predicted, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.8) +
    labs(x = xlab, y = ylab, colour = "Scale") +
    scale_color_manual(values = c("subplot"="grey80","plot"="grey30")) +
    scale_fill_manual(values = c("subplot"="grey80","plot"="grey30"), guide="none") +
    theme_classic() +
    theme(axis.title = element_text(size = 7))
  
  if(!is.null(annot)){
    p <- p + annotate(
      "text",
      x = annot_x,
      y = annot_y,
      label = annot,
      size = 2.5,
      hjust = 0.5
    )
  }
  
  p
}


pp_inset_model <- function(model, scale_y = 1, p_lab = NULL,
                           annot_x = 1.5, annot_y = NULL){
  
  pr <- ggpredict(model, terms="level")
  pr <- as.data.frame(pr)
  
  pr$predicted <- pr$predicted * scale_y
  pr$conf.low  <- pr$conf.low  * scale_y
  pr$conf.high <- pr$conf.high * scale_y
  
  ggplot(pr, aes(x=x,y=predicted, colour=x))+
    geom_errorbar(aes(ymin=conf.low,ymax=conf.high), width=.08) +
    geom_point(size=2.4) +
    annotate(
      "text",
      x = annot_x,
      y = annot_y,
      label = p_lab,
      size = 2.5
    )+
    scale_color_manual(values=c("subplot"="grey80","plot"="grey30"),guide="none")+
    theme_classic(base_size=8)+
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(size=7),
      legend.position="none",
      plot.margin = margin(5,5,5,0)
    )
}

ylim_hgt  <- c(0.2,2.0)
ylim_cv   <- c(10,85)
ylim_eff  <- c(1,12)
ylim_rich <- c(0.8,7)

y_time_hgt  <- ylim_hgt[2]*0.98
y_time_cv   <- ylim_cv[2]*0.98
y_time_eff  <- ylim_eff[2]*0.98
y_time_rich <- ylim_rich[2]*0.98

p_hgt <- pp(
  gam_mean_hgt_base,
  terms=c("time_snc_full_disturbance[0:7]","level"),
  ylab="Mean height [m]",
  annot="Time:\np = 0.004",
  annot_y=y_time_hgt
) + coord_cartesian(ylim=ylim_hgt)

p_cv_pos <- pp(
  gam_cv_hgt_pos_both,
  terms=c("time_snc_full_disturbance[0:7]","level"),
  ylab="CV [%]",
  annot="Time:\np = 0.668",
  scale_y=100,
  annot_y=y_time_cv
) + coord_cartesian(ylim=ylim_cv)

p_eff <- pp(
  gam_eff_both,
  terms=c("time_snc_full_disturbance[0:7]","level"),
  xlab="Time since disturbance\n(years)",
  ylab="Effective species [#]",
  annot="Time:\np = 0.452",
  annot_y=y_time_eff
) + coord_cartesian(ylim=ylim_eff)

p_rich <- pp(
  gam_rich_both,
  terms=c("time_snc_full_disturbance[0:7]","level"),
  xlab="Time since disturbance\n(years)",
  ylab="Species richness [#]",
  annot="Time:\np = 0.896",
  annot_y=y_time_rich
) + coord_cartesian(ylim=ylim_rich)

p_time_since_all <- ggarrange(
  p_hgt,
  p_cv_pos,
  p_eff,
  p_rich,
  ncol = 2, nrow = 2,
  labels = c("[a]", "[b]", "[c]", "[d]"),
  font.label = list(size = 10, face = "plain"),
  common.legend = TRUE,
  legend = "bottom"
)

p_time_since_all

# ggsave('outFigs/p_time_since.png',
#        plot =  p_time_since_all, 
#        width = 4.5, height = 5,
#        bg = "white")



# combine time plots with erro plot from models  --------

p_inset_hgt  <- pp_inset_model(gam_mean_hgt_base, scale_y=1,
                               p_lab="Scale:\np < 0.001",
                               annot_y=y_time_hgt)

p_inset_cv   <- pp_inset_model(gam_cv_hgt_pos_both, scale_y=100,
                               p_lab="Scale:\np < 0.001",
                               annot_y=y_time_cv)

p_inset_eff  <- pp_inset_model(gam_eff_both, scale_y=1,
                               p_lab="Scale:\np < 0.001",
                               annot_y=y_time_eff)

p_inset_rich <- pp_inset_model(gam_rich_both, scale_y=1,
                               p_lab="Scale:\np < 0.001",
                               annot_y=y_time_rich)


# Apply same y-limits to both plots in the pair
pair_hgt <- (p_hgt + coord_cartesian(ylim = ylim_hgt)) +
  (p_inset_hgt + coord_cartesian(ylim = ylim_hgt)) +
  plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[a]")

pair_cv <- (p_cv_pos + coord_cartesian(ylim = ylim_cv)) +
  (p_inset_cv + coord_cartesian(ylim = ylim_cv)) +
  plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[b]")

pair_eff <- (p_eff + coord_cartesian(ylim = ylim_eff)) +
  (p_inset_eff + coord_cartesian(ylim = ylim_eff)) +
  plot_layout(widths = c(2.6, 1)) +
  plot_annotation(title = "[c]")

pair_rich <- (p_rich + coord_cartesian(ylim = ylim_rich)) +
  (p_inset_rich + coord_cartesian(ylim = ylim_rich)) +
  plot_layout(widths = c(2.6, 1)) +
  plot_annotation(title = "[d]")

p_final <- (pair_hgt | pair_cv) /
  (pair_eff | pair_rich) +
  plot_layout(guides = "collect") & #+  plot_annotation(tag_levels = list(c("[a]", "[b]", "[c]", "[d]")))  
#&
  theme(legend.position = "bottom")

p_final

ggsave('outFigs/p_time_scale.png',
       plot =  p_final, 
       width = 6, height = 5,
       bg = "white")



#### Effect of management -------------------------------------------------

# Extract parametric terms
model_intensity_all_df <- map_dfr(fin.models, tidy, parametric = TRUE, .id = "response") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    
    # Relabel terms nicely
    term = dplyr::recode(term,
                         "planting_intensity" = "Planting",
                         "anti_browsing_intensity" = "Browsing\nprotection",
                         "grndwrk_intensity" = "Soil\npreparation",
                         "planting_intensity:anti_browsing_intensity" = "Planting×Browsing\nprotection",
                         # "levelplot" = "Level: Plot",
                         #"year_f2025" = "Year: 2025",
                         "time_snc_full_disturbance" = "Time since disturbance",
                         .default = term
    ),
    
    response = dplyr::recode(response,
                             "eff" = "Effective species [#]",
                             "rich" = "Species richness [#]",
                             "hgt" = "Mean height [m]",
                             "cvpos" = "CV [%]",
                             "beta" = "Turnover [dim.]",
                             .default = response
    ),
    
    estimate_adj = ifelse(abs(estimate) < 0.01, sign(estimate + 1e-6) * 0.01, estimate),
    lower_adj = ifelse(abs(estimate) < 0.01, 0, lower),
    upper_adj = ifelse(abs(estimate) < 0.01, 0, upper),
    
    response = factor(response, 
                      levels = c(
                        "Mean height [m]",
     # "Presence Height\nvariability [CV, %]",
     "CV [%]",
     "Effective species [#]",
     "Species richness [#]",
     "Turnover [dim.]"
     
    )
    )
  )

# Convert to percentage change (on log-link scale)
model_intensity_all_df_pct <- model_intensity_all_df %>%
  mutate(
      estimate_pct = (exp(estimate) - 1) * 100,
    lower_pct = (exp(lower) - 1) * 100,
    upper_pct = (exp(upper) - 1) * 100,
    p_label = dplyr::case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ formatC(p.value, format = "f", digits = 3)
    ),
    sig_col = ifelse(p.value < 0.05, "sig", "n.s.")
  )

# Filter to include only management terms (not year or level)
management_terms <- c("Planting", 
                      "Browsing\nprotection", 
                      "Soil\npreparation"#, 
                      #"Planting×Browsing\nprotection"#,
                      #"Time since disturbance"
)

model_intensity_all_df_pct_mng <- model_intensity_all_df_pct %>%
  filter(term %in% management_terms)

names(model_intensity_all_df_pct_mng)

# subset only important columns to get specific results
model_intensity_all_df_pct_mng_clean <- model_intensity_all_df_pct_mng %>%
  dplyr::select(
    response,
    term,
    estimate_pct,
    lower_pct,
    upper_pct,
    p_label,
    sig_col
  )
#View(model_intensity_all_df_pct_mng_clean)


p_model_response <-ggplot(model_intensity_all_df_pct_mng, 
                          aes(x = response, y = estimate_pct, 
                              fill = response)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  geom_col(position = "dodge", width = 0.7, color = NA) +
  #geom_col(position = "dodge", width = 0.7) +
  # geom_errorbar(aes(ymin = lower_pct, ymax = upper_pct), 
  #               width = 0.15, linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower_pct, ymax = upper_pct), 
                width = 0.15, linewidth = 0.5, color = "grey30") +
  facet_wrap(~ term, ncol = 4) +
  theme_classic(base_size = 8) +
  scale_fill_brewer(palette = "Dark2") +
  geom_text(
    aes(
      label = p_label,
      y = ifelse(estimate_pct >= 0, 
                 upper_pct + 5, 
                 lower_pct - 14),
       fontface = ifelse(p.value < 0.05, "bold", "plain")
    ),
    vjust = ifelse(model_intensity_all_df_pct_mng$estimate_pct >= 1, 0, 1),
    size = 2.5
  )+
  # scale_color_manual(values = c("sig" = "black", "n.s." = "grey60")) +
  labs(
    x = "",
    #y = "Effect size (ratio)", # "Multiplicative effect on response"
    y = "Effect on response [%]",
    title = ""
  ) +
  #coord_cartesian(ylim = c(-80, 120)) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

p_model_response


ggsave('outFigs/p_model_response.png',
       plot =  p_model_response, 
       width = 7, height = 5)



##### Export summary table all models 
#### Export model tables  -----------------------------
library(tibble)


# Function to clean term names
clean_term <- function(x) {
  x %>%
    str_replace_all("[:()`,]", "_") %>%
    str_replace_all("__+", "_") %>%
    str_remove_all("^_|_$")
}

# Function to convert p-values to significance codes
p_to_signif <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    p <= 0.1   ~ ".",
    TRUE       ~ "n.s."
  )
}

# Format p-values to max 3 digits
format_pval <- function(p) {
  ifelse(is.na(p), NA, formatC(p, format = "f", digits = 3))
}

# Extractor function
# Extractor function
extract_gam_summary <- function(model, model_name) {
  s <- summary(model)
  
  # Parametric terms
  param <- as.data.frame(s$p.table)
  pval_col <- grep("Pr\\(>.*\\)", colnames(param), value = TRUE)
  param <- rownames_to_column(param, "term")
  
  # Extract intercept value separately
  intercept_val <- param %>%
    filter(term == "(Intercept)") %>%
    pull(Estimate) %>%
    exp()
  
  # Keep only non-intercept terms for p-values
  param_clean <- param %>%
    filter(term != "(Intercept)") %>%
    select(term, p.value = all_of(pval_col))
  
  # Smooth terms
  smooth <- if (!is.null(s$s.table)) {
    as.data.frame(s$s.table) %>%
      rownames_to_column("term") %>%
      select(term, p.value = `p-value`)
  } else {
    tibble(term = character(), p.value = numeric())
  }
  
  # Combine parametric + smooth terms
  all_terms <- bind_rows(param_clean, smooth) %>%
    mutate(
      term_clean = clean_term(term),
      signif = p_to_signif(p.value),
      pval_fmt = format_pval(p.value),
      model = model_name
    )
  
  # Model-level metrics
  model_metrics <- tibble(
    model = model_name,
    response_intercept = intercept_val,
    r_squared = s$r.sq,
    deviance_explained = s$dev.expl,
    n_samples = s$n
  )
  
  list(pvalues = all_terms, metrics = model_metrics)
}

# Apply to models
results <- map2(fin.models, names(fin.models), extract_gam_summary)

# Wide format for significance codes
# pvals_signif <- map_dfr(results, "pvalues") %>%
#   select(model, term_clean, signif) %>%
#   pivot_wider(names_from = term_clean, values_from = signif)


# Also wide format for p-values (optional)
pvals_fmt <- map_dfr(results, "pvalues") %>%
  select(model, term_clean, pval_fmt) %>%
  pivot_wider(names_from = term_clean, values_from = pval_fmt, names_prefix = "p_")


# Model metrics
metrics_df <- bind_rows(map(results, "metrics"))

# Join together
# final_results <- left_join(metrics_df, 
#                            pvals_signif, by = "model")
final_results <- left_join(metrics_df, 
                           pvals_fmt, 
                           by = "model")


# Export to Word
library(sjPlot)
tab_df(final_results,
       title = "Model Summary Table",
       file = "outTable/models_summary_final_results.doc")









#### share spruce by planting & anti-browsing ----------------------------------------------

# Plot for planting_intensity
p1 <- both_levels_re2 %>%
  filter(level == 'plot') %>%
  ggplot(aes(x = factor(planting_intensity),
             y = spruce_share*100,
             fill = factor(planting_intensity))) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_jitter(alpha = 0.5, size = 0.5) +
  stat_compare_means(method = "kruskal.test", 
                     label.x = 2,
                     label.y = 1.05 * max(both_levels_re2$spruce_share*100, na.rm = TRUE)) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "red", color = "white")+
  theme_classic2() +
  labs(x = "Planting Intensity\n", y = "Spruce share [%]") +
  theme(legend.position = "none")

# Plot for anti_browsing_intensity
p2 <- both_levels_re2 %>%
  filter(level == 'plot') %>%
  ggplot(aes(x = factor(anti_browsing_intensity),
             y = spruce_share*100,
             fill = factor(anti_browsing_intensity))) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_jitter(alpha = 0.5, size = 0.5) +
  stat_compare_means(method = "kruskal.test", 
                     label.y = 1.05 * max(both_levels_re2$spruce_share*100, na.rm = TRUE), 
                     label.x = 2) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "red", color = "white") +

  theme_classic2() +
  labs(x = "Browsing protection\nintensity", y = "Spruce share [%]") +
  theme(legend.position = "none")

# Combine with ggarrange
windows(width = 7, height = 3.2)
# Create the arranged plot object
p_spruce_shares_boxplot <- ggarrange(p1, p2, 
                           labels = c("[a]", "[b]"), 
                           ncol = 2, align = "v",
                           font.label = list(face = "plain"))

# Save to PNG
ggsave("outFigs/spruce_share_boxplots.png", 
       plot = p_spruce_shares_boxplot, 
       width = 7, height = 3.2, units = "in", dpi = 300)









### Convert to long format: plot and subplot levels ---------------
# List of response variables to pivot
response_vars <- c("mean_hgt",
  #"dens_m2", 
                   "cv_hgt", 
                     "effective_numbers" , 
                   "sp_richness"#, 
                   #"shannon_sp", 
                   #"evenness_sp",
                  
  #                 "dens_ha"
  )

# Pivot to long format
both_levels_long <- both_levels_re2 %>%
  pivot_longer(cols = all_of(response_vars),
               names_to = "variable",
               values_to = "value")

both_levels_long_plot <- both_levels_long %>% 
  filter(level == 'plot')


# !!!! 
# both_levels_long_plot <- both_levels_long %>% 
#   filter(level == "plot") %>% 

# capp all variables to 95%
# Cap values at the 95th percentile within each variable
both_levels_long_capped <- both_levels_long %>%
  group_by(variable) %>%
  mutate(
    value_capped = ifelse(
      value > quantile(value, 0.95, na.rm = TRUE),
      quantile(value, 0.95, na.rm = TRUE),
      value
    )
  ) %>%
  ungroup(.) %>% 
  mutate(value_capped = ifelse(variable == "cv_hgt", value_capped * 100, value_capped)) %>% 
  mutate(variable = factor(variable, levels = c('mean_hgt',
                                                "cv_hgt",
                                                'effective_numbers',
                                                'sp_richness')), 
         variable = recode(variable,
                    "effective_numbers" = "Species diversity\n[Effective #]",
                    "sp_richness"       = "Species richness\n[#]",
                    "mean_hgt"          = "Mean height\n[m]",
                    "cv_hgt"            = "Height variability\n[CV, %]"
         )) 



####  Cross-scales:  -------------------------------
##### Half-violin plot ------------------------------------------------

library(gghalves)
# Define which variables are discrete (integers)
discrete_vars <- c("Species richness\n[#]", "Species diversity\n[Effective #]")

# Create a flag in the dataframe - adjust violin plot for discrete vars
both_levels_long_capped <- both_levels_long_capped %>%
  mutate(is_discrete = variable %in% discrete_vars)

# Get min y per facet (variable), to mask below 0
facet_ymins <- both_levels_long_capped %>%
  distinct(variable) %>%
  mutate(
    ymin = c(-0.5, -20, -1.8, -0.9),  # your manual values here
    ymax = 0
  )

# add formal stats
library(dplyr)
library(purrr)
library(broom)

# Wilcoxon test per variable
wilcox_results <- both_levels_long_capped %>%
  group_by(variable) %>%
  summarise(
    test = list(wilcox.test(value_capped ~ level, data = cur_data())),
    .groups = "drop"
  ) %>%
  mutate(
    tidied = map(test, broom::tidy)
  ) %>%
  unnest(tidied) %>%
  mutate(
    p.value = round(p.value, 4),
    p.adj = round(p.adjust(p.value, method = "BH"), 4),
    p.display = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3)),
    signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      p.adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  select(variable, p.value, p.adj, p.display, signif)

# Get label y-position for each facet (e.g. just above Q3 max)
label_positions <- both_levels_long_capped %>%
  group_by(variable) %>%
  summarise(y = max(value_capped, na.rm = TRUE) * 1.20)

# Join with test results
wilcox_plot_labels <- wilcox_results %>%
  left_join(label_positions, by = "variable")



# Plot with two geom_half_violin calls (with filtering inside the geom)
p_violin <- ggplot() +
  # Continuous vars (default smoothing)
  gghalves::geom_half_violin(
    data = both_levels_long_capped %>% filter(!is_discrete),
    aes(x = level, y = value_capped, fill = level),
    side = "l",
    color = NA,
    trim = FALSE,
    scale = "width"
  ) +
  # Discrete vars (smoother KDE)
  gghalves::geom_half_violin(
    data = both_levels_long_capped %>% filter(is_discrete),
    aes(x = level, y = value_capped, fill = level),
    side = "l",
    color = NA,
    trim = FALSE,
    scale = "width",
    adjust = 2
  ) +
  # Shared boxplot layer
  geom_boxplot(
    data = both_levels_long_capped,
    aes(x = level, y = value_capped),
    width = 0.1,
    outlier.shape = NA,
    color = "black"
  ) +
  # White rectangle to hide KDE below zero
  geom_rect(
    data = facet_ymins,
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    fill = "white",
    color = NA
  ) +
  # Red mean dots
  stat_summary(
    data = both_levels_long_capped,
    aes(x = level, y = value_capped),
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2,
    fill = "red",
    color = "black"
  ) +
  scale_fill_manual(
    values = c("subplot" = "grey50", "plot" = "grey80")
  ) +
  geom_text(
    data = wilcox_plot_labels,
    inherit.aes = FALSE,
    aes(x = 1.5, y = y, label = p.display),  # x = 1.5 centers between subplot and plot
    size = 3
  )+
  facet_wrap(~ variable, scales = "free_y", ncol = 4, nrow = 1) +
  theme_classic(base_size = 10) +
  labs(x = NULL, y = NULL, fill = "Level") +
 
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 9, face = "plain"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

# Display the plot
p_violin


ggsave(
  filename = "outFigs/p_violin.png",
  plot = p_violin,
  width = 7,       # Adjust width as needed
  height = 2.5,      # Adjust height as needed
  units = "in",
  dpi = 300        # High resolution for publication
)

ggsave(
  filename = "outFigs/p_violin.svg",
  plot = p_violin,
  width = 7,       # Adjust width as needed
  height = 2.5,      # Adjust height as needed
  units = "in",
  dpi = 300        # High resolution for publication
)



###### Summary table ----------------------
# Step 1: Summarise values
summary_stats_violin <- both_levels_long_capped %>%
  group_by(variable, level) %>%
  summarise(
    n = n(),
    mean = mean(value_capped, na.rm = TRUE),
    median = median(value_capped, na.rm = TRUE),
    Q1 = quantile(value_capped, 0.25, na.rm = TRUE),
    Q3 = quantile(value_capped, 0.75, na.rm = TRUE),
    IQR = IQR(value_capped, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(variable, level)

# Step 2: Calculate difference and percent increase (rounded)
mean_diffs <- summary_stats_violin %>%
  select(variable, level, mean) %>%
  pivot_wider(names_from = level, values_from = mean) %>%
  mutate(
    mean_diff = round(plot - subplot, 2),
    percent_increase = round(100 * (plot - subplot) / subplot, 2) # change in % from subplot to plot level
  )

# Step 3: Join back to summary
summary_stats_violin_with_change <- summary_stats_violin %>%
  left_join(
    mean_diffs %>% select(variable, mean_diff, percent_increase),
    by = "variable"
  )

# View final table
print(summary_stats_violin_with_change)



# how many species per plot?

# Replace df and column names if needed
both_levels_long_capped %>%
  filter(variable == "Species richness\n[#]", level == "plot") %>%
  count(value_capped) %>%
  mutate(share = round(100 * n / sum(n), 1)) %>%
  arrange(value_capped)



### Plots:----------------
##### Height by species and management intensity ---------------
# what is the heigh distribution per species and per management level?

# Convert to tibble (safer than data.table auto behavior)
dat <- dat_overlap_mng_upd2 %>%
  as_tibble()

dat_clean <- dat %>%
  dplyr::filter(
    !is.na(hgt_est),
    !is.na(n),
    n > 0
  )

# Optional sanity check
summary(dat_clean$hgt_est)
summary(dat_clean$n)


# Calculate weighted mean height per species × management 

height_species_mng <- dat_clean %>%
  dplyr::group_by(
    species,
    planting_intensity,
    anti_browsing_intensity
  ) %>%
  dplyr::summarise(
    total_stems = sum(n, na.rm = TRUE),
    mean_height_weighted = weighted.mean(hgt_est, w = n, na.rm = TRUE),
    mean_height_unweighted = mean(hgt_est, na.rm = TRUE),
    .groups = "drop"
  )


# get average per year
height_species_mng_yr <- dat_clean %>%
  dplyr::group_by(
    year,
    species,
    planting_intensity,
    anti_browsing_intensity
  ) %>%
  dplyr::summarise(
    total_stems = sum(n, na.rm = TRUE),
    mean_height_weighted = weighted.mean(hgt_est, w = n, na.rm = TRUE),
    mean_height_unweighted = mean(hgt_est, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(
    seral_group = dplyr::if_else(species %in% earlyspecs_laura, "early", "late")
  )

# aggregate to: year × seral × planting × browsing
height_seral_mng_year <- height_species_mng_yr %>%
  dplyr::group_by(year, seral_group, planting_intensity, anti_browsing_intensity) %>%
  dplyr::summarise(
    mean_h = weighted.mean(mean_height_weighted, w = total_stems, na.rm = TRUE),
    stems  = sum(total_stems, na.rm = TRUE),
    .groups = "drop"
  )

# limit cells with low stem density
min_stems = 3
height_seral_mng_year_f <- height_seral_mng_year %>%
  dplyr::filter(stems >= min_stems)


ggplot(height_seral_mng_year_f,
       aes(x = planting_intensity,
           y = anti_browsing_intensity,
           fill = mean_h)) +
  geom_tile(color = "grey80") +
  
  facet_grid(seral_group ~ year) +
  
  scale_fill_viridis_c(limits = c(0, 2.5),
    option = "E",
    direction = -1,
    end = 0.95,
    name = "Mean height (m)"
  ) +
  
  geom_text(aes(label = round(mean_h, 2)),
            color = "white",
            size = 2.5,
            fontface = "bold") +
  
  coord_equal() +
  theme_classic() +
  labs(x = "Planting intensity",
       y = "Anti-browsing intensity") +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


# View result
height_species_mng


height_species_mng2 <- height_species_mng %>%
  dplyr::mutate(
    seral_group = dplyr::if_else(
      species %in% earlyspecs_laura,
      "early",
      "late"
    )
  )

# calculate height averageges based on seral groups
height_seral_mng <- height_species_mng2 %>%
  dplyr::group_by(
    seral_group,
    planting_intensity,
    anti_browsing_intensity
  ) %>%
  dplyr::summarise(
    mean_height_weighted =
      weighted.mean(mean_height_weighted,
                    w = total_stems,
                    na.rm = TRUE),
    total_stems = sum(total_stems),
    .groups = "drop"
  )



height_seral_mng %>%
 # dplyr::filter(total_stems >= 10) %>%
  ggplot(aes(x = planting_intensity,
             y = anti_browsing_intensity,
             fill = mean_height_weighted#,
             #alpha = total_stems
             )) +
  geom_tile(color = "grey80") +
  facet_wrap(~ seral_group) +
  scale_fill_viridis_c(limits = c(0, 2.5), 
                       oob = scales::squish, 
                       option = "E",
                       direction = -1,   # 🔹 reverse scale
                       end = 0.95,
                       name = "Mean height (m)") +
   geom_text(aes(label = round(mean_height_weighted, 1)),
            color = "white", size = 2.2, fontface = "bold") +
  
 # scale_alpha(range = c(0.4, 1), guide = "none") +
  coord_equal() +
  theme_classic() +
  labs(x = "Planting intensity",
       y = "Browsing protection intensity") + 
  theme(legend.position = "bottom")








height_species_mng %>%
  dplyr::filter(species != "") %>%   # keep only well-supported cells
  ggplot(aes(x = planting_intensity,
             y = anti_browsing_intensity,
             fill = mean_height_weighted,
             alpha = total_stems)) +
  geom_tile(color = "grey80") +
  facet_wrap(~ species) +
  scale_fill_gradientn(
    colours = c("#e5f5e0",  # very light green
                "#74c476",  # medium green
                "#00441b"), # dark forest green
    #limits = c(0, max(height_species_mng$mean_height_weighted, na.rm = TRUE)),
    limits = c(0, 2),
    oob = scales::squish,   # values >5 will be clipped to dark green
    name = "Mean height (m)"
  ) +
  scale_alpha(range = c(0.4, 1), guide = "none") +
  #scale_fill_viridis_c(name = "Mean height (m)") +
  coord_equal() +
  theme_classic() +
  labs(x = "Planting intensity",
       y = "Anti-browsing intensity")

height_species_mng %>%
  dplyr::filter(!is.na(species), species != "") %>%
  ggplot(aes(
    x = factor(planting_intensity),
    y = factor(anti_browsing_intensity),
    fill = mean_height_weighted
  )) +
  geom_tile(color = "grey80") +
  facet_wrap(~ species) +
  scale_fill_viridis_c(
    limits = c(0, 2.5),
    oob = scales::squish,
    option = "E",
    direction = -1,
    end = 0.95,
    name = "Mean height (m)"
  ) +
  # halo label (underlay)
  geom_text(aes(label = round(mean_height_weighted, 1)),
            color = "white", size = 2.2, fontface = "bold") +
  coord_equal() +
  theme_classic() +
  labs(x = "Planting intensity",
       y = "Anti-browsing intensity")



##### Disturbance characteristics (plot) --------------
plot_context_chars <- dat_overlap_mng_upd2%>% 
  dplyr::select(plot, year,
                status,
                disturbance_year, 
                forest_year, 
                disturbance_length, 
                time_snc_full_disturbance, 
                time_snc_part_disturbance,
                planting_intensity,
                clear_intensity,          
                grndwrk_intensity,         
                logging_trail_intensity,
                planting_intensity,
                anti_browsing_intensity
                #clear, grndwrk, logging_trail, planting, anti_browsing
  )  %>%  
  distinct()


plot_context_chars <- plot_context_chars %>%
  mutate(disturbance_year = case_when(
    disturbance_year < 2018 ~ 2018,
    disturbance_year > 2022 ~ 2022,
    TRUE ~ disturbance_year
  )) %>% 
  mutate(time_since_disturbance = year - disturbance_year)





p_hist_dist_year <- plot_context_chars %>%
  group_by(disturbance_year) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = disturbance_year, y = n)) +
  geom_col(fill = 'grey80', color = 'black',  width = 0.8) +
  #geom_histogram(binwidth = 0.8, ) +
  scale_x_continuous(breaks = seq(2018, 2022, 2)) +
  labs(x = "Disturbance Year\n", y = "Number of Plots [#]") +
  theme_classic2() +
  theme(
    axis.text = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 8)
  )

# Histogram: Time Since Disturbance
p_hist_time_since_dist <- plot_context_chars %>%
  group_by(time_since_disturbance) %>% 
  summarise(n = n()) %>% 
  
  ggplot(aes(x = time_since_disturbance, y = n)) +
  geom_col(fill = 'grey80', color = 'black', width = 0.8) +
  labs(x = "Time since disturbance\n[years]", 
       y = "Number of Plots [#]") +
  theme_classic2() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

# Combine
p_combined_disturb_fig <- ggarrange(
  p_hist_dist_year, p_hist_time_since_dist,
  # p_management_bin_plot,
  labels = c("[a]", "[b]"),
  font.label = list(size = 10, face = 'plain'),
  widths = c(0.9,1.1),
  ncol = 2, nrow = 1
)
p_combined_disturb_fig


windows(width = 6,2.5)
p_combined_disturb_fig
ggsave("outFigs/p_combined_disturb_fig.png", 
       plot = p_combined_disturb_fig, 
       width = 5, 
       height = 2.5, 
       units = "in", dpi = 300)

dev.off() 




# Summary: Disturbance Year
summary_disturbance_year <- plot_context_chars %>%
  count(disturbance_year, name = "n") %>%
  mutate(
    share = round(100 * n / sum(n), 1)
  )

# Time Since Disturbance
summary_time_since <- plot_context_chars %>%
  count(time_since_disturbance, name = "n") %>%
  mutate(
    share = round(100 * n / sum(n), 1)
  )



#  Count combinations
combo_counts <- both_levels_re2 %>%
  filter(level == "plot") %>%
  count(planting_intensity, anti_browsing_intensity)  # gives n per combo

# density 2D plot - number of plots
p_density_mng <- both_levels_re2 %>%
  filter(level == "plot") %>%
  ggplot(aes(x = planting_intensity,
             y = anti_browsing_intensity)) +
  geom_density_2d_filled(contour_var = "ndensity")+
  #facet_wrap(~year) +
  #scale_fill_viridis_d(option = "inferno") +
  #scale_fill_gradient(low = "grey90", high = "black")+
  scale_fill_viridis_d(option = "magma", 
                       begin = 0.2,
                       direction = -1) +
  geom_point(data = combo_counts, 
             aes(x = planting_intensity, y = anti_browsing_intensity, 
                 size = n), alpha = 0.9) +
  scale_size_continuous(range = c(1, 8)) +
  # geom_jitter(alpha = 0.5,
  #             width = 0.02, height = 0.02) +
  labs(
    x = "Planting intensity",
    y = "Browsing protection intensity",
    fill = "Density",
    size = "Plots [#]"
    #title = "Co-occurrence density of planting and anti-browsing"
  ) +
  theme_classic(base_size = 10) +
  guides(
    fill = guide_legend(ncol = 1),
    size = guide_legend(ncol = 1)
  ) +
  theme(
    legend.position = "right",
    legend.box = "horizontal"
  )
windows(width = 6,4)
p_density_mng
ggsave("outFigs/density_plot.png", 
       plot = p_density_mng, 
       width = 6, 
       height = 4, 
       units = "in", dpi = 300)

dev.off() 


# Combine
p_combined_management_bin <- ggarrange(
  #p_hist_dist_year, p_hist_time_since_dist,
  p_combined_disturb_fig,
  p_management_bin_plot,
  labels = c(" ", "[c]"),
  font.label = list(size = 10, face = 'plain'),
  ncol = 1, nrow = 2,
  # align = 'hv',
  widths = c(1,1.5),  
  heights = c(1.2,1)  
)

p_combined_management_fig
# Save as PNG
ggsave("outFigs/combined_management_bin.png", plot = combined_management_bin,
       width = 4, height = 4, units = "in", dpi = 300)

# # Save as SVG
# ggsave("outFigs/combined_management_fig.svg", plot = p_combined_management_fig,
#        width = 7, height = 3.5, units = "in", dpi = 300)



# Combine
p_combined_management_intens <- ggarrange(
  #p_hist_dist_year, p_hist_time_since_dist,
  p_combined_disturb_fig,
  p_management_intensity_plot,#p_density_mng,
  labels = c(" ", "[c]"),
  font.label = list(size = 10, face = 'plain'),
  ncol = 1, nrow = 2,
  # align = 'hv',
  widths = c(1,1.6),  
  heights = c(1.2,1.2)  
)

p_combined_management_intens

windows(6,5)
p_combined_management_intens
# Save as PNG
ggsave("outFigs/p_combined_management_intens.png", plot = p_combined_management_intens,
       width = 6, height = 5, units = "in", dpi = 300)







