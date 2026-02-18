
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





theme_set(theme_classic2(base_size = 10) +
            theme(axis.title = element_text(size = 10),
                  axis.text  = element_text(size = 10)))



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
# find species with the highest share of stems overall
# 0) Safe counts (treat NA counts as 0)
df <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>%
  filter(!is.na(species) & species != "")

# 1) Totals
overall_totals <- df %>%
  summarise(total_trees = sum(n),
            n_plots = n_distinct(plot),
            n_subplots = n_distinct(subplot))

n_trees_total <- overall_totals %>% 
  pull(total_trees)

# Species counts and presence - now without year! all pooled as composition have not changed much
species_stem_share <- 
  df %>%
  group_by(species) %>%
  summarise(
    stems = sum(n),                                   # number of trees
    #plots_present = n_distinct(plot[n > 0]),          # plots where species occurs
    .groups = "drop"
  ) %>%
  dplyr::arrange(species) %>% 
  mutate(#n_trees = n_trees,
        # trees25 = n_trees25,
         share = round(stems/n_trees_total*100,2)) 

# merge counst across 2 years:     
top_overall_stem_share <- 
  species_stem_share %>%
  #dplyr::select(species, total_stems, total_share) %>%
  dplyr::slice_max(order_by = share, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup() 

v_top_species_overall <- top_overall_stem_share %>% pull(species)

# Reverse the color palette and map to the species in the desired order
n_colors  <- length(v_top_species_overall)  # Number of species
my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n_colors)  # Generate colors
# make / use a palette of exactly the needed length
pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(length(v_top_species_overall))
pal <- rev(pal)  # start with dark green

# map colors to species automatically
species_colors <- setNames(pal, v_top_species_overall)
species_colors
#

# species_colors
# piab      besp      pisy      qusp      fasy      lade      saca      soau      acps      potr      absp      sasp 
#"#006837" "#17934D" "#58B65F" "#94D168" "#C6E77F" "#EDF7A7" "#FEF0A7" "#FDCD7B" "#FA9C58" "#EE613D" "#D22B26" "#A50026" 


# Print the color assignments for confirmation
print(species_colors)


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
  potr = "Populus tremula" #,
  #absp = "Abies sp.",
  #sasp = "Salix sp."
)

species_levels <- rev(names(species_labels))  # Custom order, matching color palette and labels

### How many plots per single species and what is the species? ----------------------

# Step 2: Summarize number of species per plot, only for n > 1 trees
plots_with_species_counts <- df %>%
  filter(n > 1) %>%
  group_by(plot) %>%
  summarise(
    n_species = n_distinct(species),
    .groups = "drop"
  )

# Step 3: Filter to plots with exactly 1 species
single_species_plots <- plots_with_species_counts %>%
  dplyr::filter(n_species == 1)

# Step 4: Join back to original to get the species name
species_per_single_species_plot <- df %>%
  filter(n > 1) %>%
  semi_join(single_species_plots, by = "plot") %>%
  group_by(plot, species) %>%
  summarise(n_trees = sum(n), .groups = "drop")  # just to show tree count if needed

# Step 5 (optional): Count how many plots per species
summary_species <- species_per_single_species_plot %>%
  count(species, name = "single_plots") %>%
  arrange(desc(single_plots)) %>% 
  mutate(share = single_plots/333)

# Print result
print(summary_species)





### Stem density per species per top 10 species --------------------------------
df_stem_dens_species <- df %>% 
  group_by(plot, year, species, n_subplots ) %>%
  summarize(sum_n = sum(n, na.rm =T)) %>% 
  mutate(scaling_factor = 10000/(n_subplots * 4),
         stem_dens = sum_n*scaling_factor) %>% 
  mutate(log_sum_stem_density = log10(stem_dens + 1)) #%>%  # Adding 1 to avoid log(0)
#ungroup()

# get total sum and calculate as average value over all sites 
df_stem_dens_species_sum <- 
  df_stem_dens_species %>% 
  group_by(species, year) %>% #, year
  summarise(stem_dens = sum(stem_dens, na.rm = T),
            log_sum_stem_density = sum(log_sum_stem_density, na.rm = T)) %>%
  mutate(stem_dens_avg = stem_dens/n_plots_total,
         log_sum_stem_density_avg = log_sum_stem_density/n_plots_total)


df_stem_dens_species <- df_stem_dens_species %>% 
  ungroup(.) %>% 
  filter(sum_n >0) %>% 
  filter(species %in% v_top_species_overall) %>% 
  dplyr::group_by(species, year) %>%
  dplyr::mutate(median_stem_density = median(stem_dens, na.rm = TRUE)) %>% 
  dplyr::ungroup(.) %>%
  mutate(species = factor(species, levels = rev(v_top_species_overall))) # Set custom order


# df_stem_dens_species2 <- df_stem_dens_species_year %>%
#   dplyr::filter(!is.na(log_sum_stem_density) & sum_n > 0) %>%
#   dplyr::mutate(year = factor(year, levels = c("2023","2025")),
#                 # order by mean log density (ascending → highest ends up at the TOP after coord_flip)
#                 species = forcats::fct_reorder(species, log_sum_stem_density, .fun = mean, na.rm = TRUE))

# 
p_density<-df_stem_dens_species_sum %>% #df_stem_dens_species_year2 %>% 
  mutate(species = factor(species, levels = species_levels)) %>% 
  filter(!is.na(species)) %>% 
  ggplot(aes(x = log_sum_stem_density, y = species,
             fill = species)) +
  geom_boxplot(
    #aes(group = interaction(species, year), 
    #   alpha = factor(year)),
    position = position_dodge(width = 0.6),
    outlier.shape = NA,
    width = 0.45#,
    # color = "black"
  ) +
  
  # coord_flip() +
  labs(
    x = "log(sum stem density)",
    y = "",  
    fill = "Year"
  ) +
  scale_fill_manual(values= species_colors) +
  theme_classic(base_size = 10) +
  scale_y_discrete(labels = species_labels) +
  theme(axis.text.y = element_text(face = "italic", size = 8),
        legend.position = 'none')


p_density




p_bar <- top_overall_stem_share %>%  
  mutate(species = factor(species, levels = species_levels)) %>% 
  ggplot(aes(x = share, y = species, fill = species)) +
 
  geom_col(aes(#group = interaction(species, year)#,
    #alpha = factor(year)
  ), 
  position = position_dodge(width = 0.7), width = 0.6) +
  # if 'share' is 0–100 already, just add a % suffix:
  scale_x_continuous(labels = label_number(accuracy = 1, 
                                           suffix = "")) +
  scale_y_discrete(
    limits = rev(names(species_labels)),
    labels = species_labels,
    drop = FALSE
  ) +
    labs(
    x = "Stems share [%]",
    y = ""
  ) +
  scale_fill_manual(values = species_colors) +
  theme_classic2(base_size = 10) +
  theme(axis.text.y = element_text(size = 8, face = "italic"),
        legend.position = "none")

p_bar

####  Get species occurence from total number of plots -------------
# Total number of unique plots
total_plots <- dat_overlap %>%
  pull(plot, year) %>%
  n_distinct()

# Share of plots per species (where species has non-zero stems)
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
species_order <- species_occurence %>%
  group_by(species) %>%
  summarise(max_share = max(share_of_plots)) %>%
  arrange(desc(max_share)) %>%
  pull(species)

species_plot_share <- species_occurence %>%
  group_by(species) %>% 
  summarize(share_of_plots_avg = mean(share_of_plots)) %>% 
  mutate(species = factor(species, levels = rev(species_order)))

# Plot
p_occurence <- species_plot_share %>% 
  filter(species %in% v_top_species_overall ) %>% 
  mutate(species = factor(species, levels = species_levels)) %>% 
  ggplot(aes(x = share_of_plots_avg, y = species, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Species occurence over plots [%]",
    y = "Species",
    fill = "Year"
  ) +
  scale_fill_manual(values = species_colors) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "none"
  )

p_occurence

# p_occurence with no y labels
p_occurence <- p_occurence +
  labs(y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# p_density with no y labels
p_density <- p_density +
  labs(y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


p_bar <- p_bar + theme(plot.margin = margin(t = 12, r = 5, b = 5, l = 5))
p_occurence <- p_occurence + theme(plot.margin = margin(t = 12, r = 5, b = 5, l = 5))


p_species_composition <- ggarrange(p_bar, p_occurence,# p_density,  
          ncol = 2, common.legend = F, 
          widths = c(1.5, 1),
          labels = c("[a]", "[b]"),
          font.label = list(size = 10, face = "plain"),
          label.x = 0.02,    # near left edge
          label.y = 1.01     # just above plot
)

p_species_composition
# Export to PNG
ggsave("outFigs/p_species_composition.png", 
       plot = p_species_composition,
       width = 7, height = 3, units = "in", dpi = 300)






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
  geom_col(width = 0.4, color = "black") +
  scale_fill_manual(values = fill_colors, name = "Intensity class",
                    breaks = intensity_levels) +
  geom_vline(xintercept = 0, color = "grey", 
             linewidth = 0.8, lty = 'dashed') +
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

sf::st_write(
  plot_sf_AEF,
  dsn = 'outDataShare/Karim_AEF/regeneration_chars_3035.gpkg',
  layer = "plot",
  delete_layer = TRUE  # overwrite layer if it already exists
)







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
  mutate(cv_hgt_present = as.integer(cv_hgt > 0))
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

both_levels_re2 %>% 
  filter(level == "plot") %>% 
  ggplot(aes(x = anti_brow_f,
             y = mean_hgt)) +
  geom_boxplot()



#### GAM management interaction + time since disturbance smooths --------------



##### Mean height -----------------------------------------------------
gam_mean_hgt_intensity_base <- gam(
  mean_hgt ~ 
    #s(planting_intensity, k = 3) +
    #s(anti_browsing_intensity, k = 3) +
    #ti(planting_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    s(time_snc_full_disturbance, k = 5) +
   # s(grndwrk_intensity, k = 3) +
    year_f +
    level +
    s(plot_id, bs = "re")
  ,
  data = both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)




gam_mean_hgt_intensity_base <- gam(
  mean_hgt ~ 
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    #ti(planting_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    s(time_snc_full_disturbance, k = 7) +
    s(grndwrk_intensity, k = 3) +
    year_f +
    level ,#+
  #  s(plot_id, bs = "re")
  ,
  data = both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)




# chewck the effect of k
gam_mean_hgt_intensity0 <- gam(
  mean_hgt ~ 
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    #ti(planting_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    s(time_snc_full_disturbance, k = 7) +
    s(grndwrk_intensity, k = 3) +
    year_f +
    level +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)

AIC(gam_mean_hgt_intensity0, gam_mean_hgt_intensity_base)

# test also only for plot vs subplot!
gam_mean_hgt_intensity02_int1_plot <- gam(
  mean_hgt ~ 
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    ti(planting_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    #ti(grndwrk_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    
    s(time_snc_full_disturbance, by = year_f, k = 4) +
    s(grndwrk_intensity, k = 3) +
    year_f +
   # level +
    s(plot_id, bs = "re"),
  data = filter(both_levels_re2, level == 'plot'),
  family = tw(link = "log"),
  method = "REML"
)


# test linear model
gam_mean_hgt_intensity02_int1_plot_lin <- gam(
  mean_hgt ~ 
    planting_intensity*anti_browsing_intensity +
    #ti(planting_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    #ti(grndwrk_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    
    s(time_snc_full_disturbance, by = year_f, k = 4) +
    grndwrk_intensity+
    year_f +
    # level +
    s(plot_id, bs = "re"),
  data = filter(both_levels_re2, level == 'plot'),
  family = tw(link = "log"),
  method = "REML"
)
AIC(gam_mean_hgt_intensity02_int1_plot_lin, gam_mean_hgt_intensity02_int1_plot)

gam_mean_hgt_intensity02_int1_plot2 <- gam(
  mean_hgt ~ 
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    ti(planting_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    #ti(grndwrk_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    
    s(time_snc_full_disturbance, by = year_f, k = 4) +
    s(grndwrk_intensity, k = 3) +
    year_f +
    # level +
    s(plot_id, bs = "re"),
  data = filter(both_levels_re2, level == 'plot'),
  family = tw(link = "log"),
  method = "REML"
)

pp <- ggpredict(gam_mean_hgt_intensity02_int1_plot, 
          'time_snc_full_disturbance')

plot(pp)

gam_mean_hgt_intensity02_int1_sub <- gam(
  mean_hgt ~ 
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    ti(planting_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    #ti(grndwrk_intensity, anti_browsing_intensity, k = 3) +      # nonlinear interaction
    
    s(time_snc_full_disturbance, k = 7) +
    s(grndwrk_intensity, k = 3) +
    year_f +
    # level +
    s(plot_id, bs = "re"),
  data = filter(both_levels_re2, level == 'subplot'),
  family = tw(link = "log"),
  method = "REML"
)


gam_mean_hgt_intensity03_int2 <- gam(
  mean_hgt ~ 
    s(planting_intensity, k = 6) +
    s(anti_browsing_intensity, k = 3) +
    #ti(planting_intensity, anti_browsing_intensity) +      # nonlinear interaction
    ti(planting_intensity, grndwrk_intensity) +      # nonlinear interaction
    s(time_snc_full_disturbance, k = 7) +
    s(grndwrk_intensity, k = 3) +
    year_f +
    #level +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)

gam_mean_hgt_intensity03_int3 <- gam(
  mean_hgt ~ 
    s(planting_intensity, k = 6) +
    s(anti_browsing_intensity, k = 3) +
    #ti(planting_intensity, anti_browsing_intensity) +      # nonlinear interaction
   # ti(planting_intensity, grndwrk_intensity) +      # nonlinear interaction
    ti(anti_browsing_intensity, grndwrk_intensity) +      # nonlinear interaction
    
     s(time_snc_full_disturbance, k = 7) +
    s(grndwrk_intensity, k = 3) +
    year_f +
    #level +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  family = tw(link = "log"),
  method = "REML"
)

AIC(gam_mean_hgt_intensity02_int1,
  gam_mean_hgt_intensity03_int2,
    gam_mean_hgt_intensity03_int3,
  gam_mean_hgt_intensity0,
  gam_mean_hgt_intensity_base)

fin.m.hgt <- gam_mean_hgt_intensity02_int1_plot_lin #gam_mean_hgt_intensity02_int1_plot#gam_mean_hgt_intensity02_int1

vis.gam(
  fin.m.hgt,
  view = c("planting_intensity", "anti_browsing_intensity"), 
  plot.type = "contour",
  theta = 30, phi = 30, # viewing angle
  color = "terrain",
  main = "Interaction: Planting x Anti-browsing",
  xlab = "Planting intensity",
  ylab = "Anti-browsing intensity",
  zlab = "Partial effect"
)
out_pred <- ggpredict(fin.m.hgt ,
          terms = 'planting_intensity',
          condition = list(level = "plot")) 

plot(out_pred)

# height 

plt <- ggpredict(
  fin.m.hgt,
  terms = c("time_snc_full_disturbance [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)
plt <- ggpredict(
  fin.m.hgt,
  terms = c("planting_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)

plt <- ggpredict(
  fin.m.hgt,
  terms = c("anti_browsing_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)

plt <- ggpredict(
  fin.m.hgt,
  terms = c("grndwrk_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)








plt <- ggpredict(
  fin.m.rich,
  terms = c("planting_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)

plt <- ggpredict(
  fin.m.rich,
  terms = c("anti_browsing_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)

plt <- ggpredict(
  fin.m.rich,
  terms = c("grndwrk_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)

# effiecinet 

plt <- ggpredict(
  fin.m.eff,
  terms = c("planting_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)

plt <- ggpredict(
  fin.m.eff,
  terms = c("anti_browsing_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)

plt <- ggpredict(
  fin.m.eff,
  terms = c("grndwrk_intensity [all]"),
  condition = list(level = "plot")
) #%>%
plot(plt)



##### CV hgt bin  ----
gam_cv_hgt_bin_intensity_base <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    #ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f,# +
   # s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = binomial(link = "logit")
)

gam_cv_hgt_bin_intensity_re <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    #ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = binomial(link = "logit")
)
# 
k.check(gam_cv_hgt_bin_intensity_re)

gam_cv_hgt_bin_intensity_re_int1 <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = binomial(link = "logit")
)

gam_cv_hgt_bin_intensity_re_int1_plot <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = filter(both_levels_re2, level == 'plot'),
  method = "REML",
  family = binomial(link = "logit")
)


gam_cv_hgt_bin_intensity_re_int1_plot_lin <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    planting_intensity*anti_browsing_intensity +
    #ti(planting_intensity, anti_browsing_intensity) +
    grndwrk_intensity +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = filter(both_levels_re2, level == 'plot'),
  method = "REML",
  family = binomial(link = "logit")
)
AIC(gam_cv_hgt_bin_intensity_re_int1_plot_lin, gam_cv_hgt_bin_intensity_re_int1_sub, gam_cv_hgt_bin_intensity_re_int1_plot)
gam_cv_hgt_bin_intensity_re_int1_sub <- gam(
  cv_hgt_present ~
    s(time_snc_full_disturbance, k = 7) +
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = filter(both_levels_re2, level == 'subplot'),
  method = "REML",
  family = binomial(link = "logit")
)


fin.m.cv.bin <- gam_cv_hgt_bin_intensity_re_int1_plot_lin #gam_cv_hgt_bin_intensity_re_int1_plot#gam_cv_hgt_bin_intensity_re_int1

AIC(gam_cv_hgt_bin_intensity_re, gam_cv_hgt_bin_intensity_base, gam_cv_hgt_bin_intensity_re_int1)
draw(gam_cv_hgt_bin_intensity_re_int1)


##### CV Height (Positive Part) -----------


gam_cv_hgt_pos_base <- gam(
  cv_hgt_pos ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
   # ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f , #+
   # s(plot_id, bs = "re"),
  data = subset(both_levels_re2, cv_hgt > 0),
  method = "REML",
  family = tw(link = "log")
)

gam_cv_hgt_pos_base_plot <- gam(
  cv_hgt_pos ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    # ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    #level +
    year_f , #+
  # s(plot_id, bs = "re"),
  data = df_plot_clean,
  method = "REML",
  family = tw(link = "log")
)

gam_cv_hgt_pos_base_plot_lin <- gam(
  cv_hgt_pos ~
    planting_intensity*anti_browsing_intensity +
    time_snc_full_disturbance +
    # ti(planting_intensity, anti_browsing_intensity) +
    grndwrk_intensity +
    #level +
    year_f ,# +
   #s(plot_id, bs = "re"),
  data = df_plot_clean,
  method = "REML",
  family = tw(link = "log")
)
AIC(gam_cv_hgt_pos_base_plot_lin, gam_cv_hgt_pos_base_plot)
gam_cv_hgt_pos_base_sub <- gam(
  cv_hgt_pos ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    # ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    #level +
    year_f , #+
  # s(plot_id, bs = "re"),
  data = df_sub_clean,
  method = "REML",
  family = tw(link = "log")
)

gam_cv_hgt_pos_re <- gam(
  cv_hgt_pos ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    # ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f +
   s(plot_id, bs = "re"),
  data = subset(both_levels_re2, cv_hgt > 0, sp_richness > 1),
  method = "REML",
  family = tw(link = "log")
)
AIC(gam_cv_hgt_pos_base, gam_cv_hgt_pos_re)


gam_cv_hgt_pos_re_plot <- gam(
  cv_hgt_pos ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    # ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
   # level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_plot_clean,
  method = "REML",
  family = tw(link = "log")
)


gam_cv_hgt_pos_re_sub <- gam(
  cv_hgt_pos ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    # ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    # level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_sub_clean,
  method = "REML",
  family = tw(link = "log")
)

AIC(gam_cv_hgt_pos_re, gam_cv_hgt_pos_re_lev)

gam_cv_hgt_pos_re_lev <- gam(
  cv_hgt_pos ~ 
    s(time_snc_full_disturbance, by = level, k = 7) +
    s(planting_intensity,  k = 3) +  
    s(anti_browsing_intensity, k = 3) +  
    s(grndwrk_intensity, k = 3) + 
    level +  # include the parametric effect
    year_f +
    s(plot_id, bs = "re"),
  data = subset(both_levels_re2, cv_hgt > 0, sp_richness > 1),
  family = tw(link = "log"),
  method = "REML"
)


gam_cv_hgt_pos_re_lev_int1 <- gam(
  cv_hgt_pos ~ 
    s(time_snc_full_disturbance, by = level, k = 7) +
    s(planting_intensity,  k = 3) +  
    s(anti_browsing_intensity, k = 3) +  
    s(grndwrk_intensity, k = 3) + 
    ti(planting_intensity,anti_browsing_intensity) + 
    level +  # include the parametric effect
    year_f +
    s(plot_id, bs = "re"),
  data = subset(both_levels_re2, cv_hgt > 0, sp_richness > 1),
  family = tw(link = "log"),
  method = "REML"
)
AIC(gam_cv_hgt_pos_re_lev, gam_cv_hgt_pos_re, gam_cv_hgt_pos_re_lev_int1)
draw(gam_cv_hgt_pos_re)

fin.m.cv.pos <- gam_cv_hgt_pos_base_plot_lin #gam_cv_hgt_pos_base_plot #gam_cv_hgt_pos_re

##### Species Diversity (Effective Number) -----------

gam_effective_intensity_base <- gam(
  effective_numbers ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    #ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f ,#+
    #s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw() #Gamma(link = "log")
)


gam_effective_intensity_re <- gam(
  effective_numbers ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    #ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f +
  s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_effective_intensity_re_lev <- gam(
  effective_numbers ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, by = level, k = 7) +
    #ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)


gam_effective_intensity_re_int1 <- gam(
  effective_numbers ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_effective_intensity_re_int1_plot <- gam(
  effective_numbers ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_plot_clean, #both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_effective_intensity_re_int1_plot_lin <- gam(
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
AIC(gam_effective_intensity_re_int1_plot_lin, gam_effective_intensity_re_int1_plot
    )
gam_effective_intensity_re_int1_sub <- gam(
  effective_numbers ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    #level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_sub_clean, #both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_effective_intensity_re_int2 <- gam(
  effective_numbers ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    ti(planting_intensity, grndwrk_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)
AIC(gam_effective_intensity_re, gam_effective_intensity_base, 
    gam_effective_intensity_re_lev_int1,
    gam_effective_intensity_re_lev_int2)

summary(gam_effective_intensity_re)

fin.m.eff <- gam_effective_intensity_re_int1_plot_lin #gam_effective_intensity_re_int1_plot #gam_effective_intensity_re_lev_int1


##### Species Richness -----------
gam_richness_intensity <- gam(
  sp_richness ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 3) +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)


gam_richness_intensity_nb <- gam(
  sp_richness ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    #s(time_snc_full_disturbance, k = 7) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 4) +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = nb(link = "log")
)

gam_richness_intensity_nb_lin <- gam(
  sp_richness ~
    planting_intensity*anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 7) +
    #ti(planting_intensity, anti_browsing_intensity) +
    grndwrk_intensity +
    level +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = nb(link = "log")
)

AIC(gam_richness_intensity_nb_lin, gam_richness_intensity_nb)
vis.gam(gam_richness_intensity_nb)
gam_richness_intensity_nb_plot <- gam(
  sp_richness ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 4) +
   # level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_plot_clean,
  method = "REML",
  family = nb(link = "log")
)

gam_richness_intensity_nb_sub <- gam(
  sp_richness ~
    s(planting_intensity, k = 3) +
    s(anti_browsing_intensity, k = 3) +
    s(time_snc_full_disturbance, k = 7) +
    ti(planting_intensity, anti_browsing_intensity) +
    s(grndwrk_intensity, k = 4) +
    # level +
    year_f +
    s(plot_id, bs = "re"),
  data = df_sub_clean,
  method = "REML",
  family = nb(link = "log")
)


fin.m.rich <- gam_richness_intensity_nb_lin #gam_richness_intensity_nb


##### CV hurdle ------------------------------------------------------

# Load the package
library(glmmTMB)

# Full hurdle model on CV data
m_cv_hurdle <- glmmTMB(
  cv_hgt_pos ~ planting_intensity * anti_browsing_intensity +
    time_snc_full_disturbance +
    grndwrk_intensity +
    year_f +
    (1 | plot_id),
  
  ziformula = ~ planting_intensity * anti_browsing_intensity +
    time_snc_full_disturbance +
    grndwrk_intensity +
    year_f +
    (1 | plot_id),
  
  family = Gamma(link = "log"),
  data = df_plot_clean
)

m_cv_zigamma_simple <- glmmTMB(
  cv_hgt_pos ~ planting_intensity * anti_browsing_intensity +
    time_snc_full_disturbance +
    grndwrk_intensity +
    year_f ,#+
   # (1 | plot_id),
  
  ziformula = ~1,  # constant zero-inflation probability
  
  family = Gamma(link = "log"),
  data = df_plot_clean
)

library(performance)

r2(m_cv_zigamma_simple)

##### Summaries and Export -----------
summary(fin.m.hgt)
appraise(fin.m.hgt)

summary(fin.m.cv.bin)
appraise(fin.m.cv.bin)

summary(fin.m.cv.pos)
appraise(fin.m.cv.pos)

summary(fin.m.eff)
appraise(fin.m.eff)

summary(fin.m.rich)
appraise(fin.m.rich)



fin.models <- list(
  hgt   = fin.m.hgt,
 # cvbin = fin.m.cv.bin,
  cvpos = fin.m.cv.pos,
  eff   = fin.m.eff,
  rich  = fin.m.rich
)


##### Get bar plot of effect sizes -------------------------------------------------

models_intensity_all <- fin.models

# Extract parametric terms
model_intensity_all_df <- map_dfr(models_intensity_all, tidy, parametric = TRUE, .id = "response") %>%
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
                             "eff" = "Species diversity\n[Effective #]",
                             "rich" = "Species richness\n[#]",
                             "hgt" = "Mean height\n[m]",
                             "cvpos" = "Height variability\n[CV, %]",
                             .default = response
    ),
    
    estimate_adj = ifelse(abs(estimate) < 0.01, sign(estimate + 1e-6) * 0.01, estimate),
    lower_adj = ifelse(abs(estimate) < 0.01, 0, lower),
    upper_adj = ifelse(abs(estimate) < 0.01, 0, upper),
    
    response = factor(response, levels = c(
      "Mean height\n[m]",
     # "Presence Height\nvariability [CV, %]",
      "Height variability\n[CV, %]",
      "Species diversity\n[Effective #]",
      "Species richness\n[#]"
    )
    )
  )

# Convert to percentage change (on log-link scale)
model_intensity_all_df_pct <- model_intensity_all_df %>%
  mutate(
    estimate_pct = (exp(estimate) - 1) * 100,
    lower_pct = (exp(lower) - 1) * 100,
    upper_pct = (exp(upper) - 1) * 100,
    
    # p_signif = case_when(
    #   p.value <= 0.001 ~ "***",
    #   p.value <= 0.01  ~ "**",
    #   p.value <= 0.05  ~ "*",
    #   TRUE             ~ "n.s."
    # ),
    # format p-values: max 3 decimals, always shown
    p_label = paste0(
      "",
      formatC(p.value, format = "f", digits = 3)
    ),
    sig_col = ifelse(p.value < 0.05, "sig", "n.s.")
  )

# Filter to include only management terms (not year or level)
management_terms <- c("Planting", 
                      "Browsing\nprotection", 
                      "Soil\npreparation", 
                      "Planting×Browsing\nprotection"#,
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
View(model_intensity_all_df_pct_mng_clean)

model_intensity_all_df_pct_mng_clean %>%
  mutate(
    ci_width = upper_pct - lower_pct,
    ci_to_estimate_ratio = ci_width / abs(estimate_pct)
  ) %>%
  summarise(
    mean_ratio = mean(ci_to_estimate_ratio),
    median_ratio = median(ci_to_estimate_ratio)
  )


p_model_response <-ggplot(model_intensity_all_df_pct_mng, # model_intensity_all_df_pct_mng
                          aes(x = term, y = estimate_pct, fill = term)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = lower_pct, ymax = upper_pct), 
                width = 0.15, linewidth = 0.5) +
  facet_wrap(~ response, ncol = 4) +
  theme_classic(base_size = 9) +
  scale_fill_brewer(palette = "Dark2") +
  # NEW: Add stars above bars
  geom_text(
    aes(label = p_label, #p_signif, 
        y = ifelse(estimate_pct >= 0, 
                   upper_pct + 5, 
                   lower_pct - 14),
        #color = sig_col,
        fontface = ifelse(p.value < 0.05, "bold", "plain")),
    vjust = 0,
    size = 2.5
  ) +
  # scale_color_manual(values = c("sig" = "black", "n.s." = "grey60")) +
  labs(
    x = "Management Intensity",
    y = "Effect on response [%]",
    title = ""
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1, face = 'italic'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

p_model_response

ggsave('outFigs/p_model_response.png',
       plot =  p_model_response, 
       width = 7, height = 5)

appraise(models_intensity_all[[4]])
draw(models_intensity[[4]])
plot.gam(models_intensity[[4]], page = 1)







### Make plots : 
##### Time since disturbance  ----------------------------
#  Mean height
preds_hgt <- ggpredict(
  fin.m.hgt,
  terms = c("time_snc_full_disturbance"),
  condition = list(level = "plot")
)

p_hgt <- plot(preds_hgt) +
  labs(y = "", x = "", title = "[a] Mean\nheight [m]") +
  theme_classic(base_size = 9) +
  coord_cartesian(ylim = c(0, 5)) +  # preserve data but fix visible y-range
  scale_y_continuous(breaks = seq(0, 5, by = 2.5)) +
  
  annotate(
    "text",
    x = 4,
    y = Inf,  # small padding above plo
    label = "0.000",       # your p-value here
    #hjust = 1.1,               # center from right
    vjust = 1.5,               # slightly below top
    size = 3
  )
p_hgt

# 2️ CV binary
pred_bin <- ggpredict(
  fin.m.cv.bin,
  terms = c("time_snc_full_disturbance"),
  condition = list(level = "plot")
) %>%
  rename(prob = predicted, prob_low = conf.low, prob_high = conf.high)

#  CV positive
pred_pos <- ggpredict(
  fin.m.cv.pos,
  terms = c("time_snc_full_disturbance"),
  condition = list(level = "plot")
) %>%
  rename(mean_pos = predicted, pos_low = conf.low, pos_high = conf.high)

# 4️⃣ Combine CV
pred_cv_combined <- left_join(pred_bin, pred_pos, by = c("x")) %>%
  mutate(
    expected_cv = prob * mean_pos,
    lower = prob_low * pos_low,
    upper = prob_high * pos_high
  )

p_cv <- ggplot(pred_cv_combined,
               aes(x = x, y = expected_cv * 100)) +
  geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.7) +
  labs(y = "", x = "", title = "[b] Height\nvariability [CV, %]") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 100)) +  # preserve data but fix visible y-range
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  
  annotate(
    "text",
    x = 4,
    y = Inf,  # small padding above plo
    label = "0.876",       # your p-value here
    #hjust = 1.1,               # center from right
    vjust = 1.5,               # slightly below top
    size = 3
  )
p_cv

#  Effective species number
preds_eff <- ggpredict(
  fin.m.eff,
  terms = c("time_snc_full_disturbance"),
  condition = list(level = "plot")
)

p_eff <- plot(preds_eff) +
  labs(y = "", x = "", title = "[c] Effective\nspecies [#]") +
  scale_y_continuous(breaks = seq(0, 15, by = 2)) +
  theme_classic(base_size = 9) +
  annotate(
    "text",
    x = 4,
    y = Inf,  # small padding above plo
    label = "0.489",       # your p-value here
    #hjust = 1.1,               # center from right
    vjust = 1.5,               # slightly below top
    size = 3
  )


# Species richness
preds_rich <- ggpredict(
  fin.m.rich,
  terms = c("time_snc_full_disturbance"),
  condition = list(level = "plot")
)

p_rich <- plot(preds_rich) +
  labs(y = "", x = "", title = "[d] Species\nrichness [#]") +
  coord_cartesian(ylim = c(0, 10)) +  # preserve data but fix visible y-range
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  
  theme_classic(base_size = 9) +
  annotate(
    "text",
    x = 4,
    y = Inf,  # small padding above plo
    label = "0.927",       # your p-value here
    #hjust = 1.1,               # center from right
    vjust = 1.5,               # slightly below top
    size = 3
  )
p_rich

combined_plot <- ggarrange(
  p_hgt, p_cv,
  p_eff, p_rich,
#  labels = c("[a]", "[b]", "[c]", "[d]"),
  ncol = 4, nrow = 1,
  align = "hv",
  common.legend = TRUE,
  legend = "top",
  font.label = list(size = 9, face = "plain", color = "black"),
  label.x = 0.1,
  label.y = 0.95
)

# Optional: annotate bottom
p_time <- annotate_figure(
  combined_plot,
  bottom = text_grob("Years since disturbance", size = 10, face = "plain")
)

p_time

ggsave('outFigs/time_since.png',
       p_time, 
       width = 7, height = 2.5)



##### Export summary table all models 
###### Collect all values from models -----------------------------
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
# View
View(final_results)



# Export to Word
tab_df(final_results,
       title = "Model Summary Table",
       file = "outTable/models_summary_final_results.doc")




# Export model summaries to Word
library(sjPlot)
tab_model(
  unlist(fin.models),
  show.icc = FALSE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  file = "outTable/forest_model_summaries_intensity.doc"
)





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





# check richness by planting?
both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance),
             y = sp_richness,
             fill = plant_f))  +
  geom_boxplot() #+
  geom_smooth(method = "lm") +
  #geom_boxplot() + 
  facet_grid(.~plant_f)
  #geom_jitter() +
  

# check species richness if planted or not --------
  both_levels_re2 %>% 
    filter(level == 'plot') %>% 
    ggplot(aes(x = time_snc_full_disturbance,
               y = sp_richness,
               fill = plant_f,
               col = plant_f))  +
   # geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm")
  #geom_jitter() +
    
    
both_levels_re2 %>% 
  filter(level == 'plot') %>% 
      ggplot(aes(x = disturbance_length,
                 y = sp_richness,
                 fill = plant_f,
                 col = plant_f))  +
     # geom_jitter(alpha = 0.5) +
      geom_smooth(method = "lm")
  


# effective species if planted or not -----------
both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = time_snc_full_disturbance,
             y = effective_numbers,
             fill = plant_f,
             col = plant_f))  +
  # geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm")
#geom_jitter() +


both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = disturbance_length,
             y = effective_numbers,
             fill = plant_f,
             col = plant_f))  +
  geom_jitter() +
  # geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm")# +

# CV vs time since -----------------------------------------
both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  filter(time_snc_full_disturbance < 8) %>% 
  ggplot(aes(x = time_snc_full_disturbance,
             y = cv_hgt_pos,
             fill = plant_f,
             col = plant_f))  +
  #geom_jitter() +
  # geom_jitter(alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) 

# mean height vs time since -----------------------
both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  filter(time_snc_full_disturbance < 8) %>% 
  ggplot(aes(x = time_snc_full_disturbance,
             y = mean_hgt,
             fill = plant_f,
             col = plant_f))  +
  #geom_jitter() +
  # geom_jitter(alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) 


# check richness by planting?
both_levels_re2 %>% 
  ggplot(aes(x = planting_intensity,
             y = sp_richness)) + 
  geom_jitter() +
  geom_smooth()

both_levels_re2 %>% 
  ggplot(aes(x = time_snc_full_disturbance,
             y = sp_richness)) + 
  geom_jitter() +
  geom_smooth()


both_levels_re2 %>% 
  filter(level == 'subplot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance),
             y = cv_hgt_pos,
             fill = early_class)) + 
  geom_boxplot() + 
  facet_grid(.~early_class)

both_levels_re2 %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance), y = cv_hgt)) +
  geom_boxplot() +
  geom_jitter() + facet_grid(.~level)

hist(both_levels_re2$area_m2, breaks = 50)








# Convert to long format: plot and subplot levels ---------------
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



## All indicators:  Create half-violin plot ------------------------------------------------

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



### Graphics: disturbance chars, species composition and management

## Get context and disturbance characteristics (plot) --------------
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


# 1. Summary: Disturbance Year
summary_disturbance_year <- plot_context_chars %>%
  count(disturbance_year, name = "n") %>%
  mutate(
    share = round(100 * n / sum(n), 1)
  )

# 2. Summary: Time Since Disturbance
summary_time_since <- plot_context_chars %>%
  count(time_since_disturbance, name = "n") %>%
  mutate(
    share = round(100 * n / sum(n), 1)
  )








### Get context, disturbnace characteristics --------------------------------------

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







