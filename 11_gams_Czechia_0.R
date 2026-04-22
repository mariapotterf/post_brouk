
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
library(broom);      library(emmeans); 
library(ggalluvial)

library(tibble)
library(spdep)
library(sjPlot)
library(patchwork)
library(flextable)

## Project variables (palettes, labels, species lists, species_class, earlyspecs_cz)
source('my_variables.R')



# Filter data from both level: plot and subplot from EU  
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

# data with individual species per subplot and plots each species,
# contains ALL plots: overlapping and no overlapping ones
dat_full <- fread("outData/dat_full_species_all_countries.csv") %>%
  filter(country_name == "Czechia") %>%
  mutate(
    plot      = factor(plot)
  )


dat_full <- dat_full %>% 
  mutate(
    time_since_cap = case_when(
      time_snc_full_disturbance == 0  ~ 1L,
      time_snc_full_disturbance == 7  ~ 6L,
      time_snc_full_disturbance == 11 ~ 6L,
      time_snc_full_disturbance == 13 ~ 6L,
      TRUE ~ time_snc_full_disturbance
    )
  ) 


dat_overlap <- dat_overlap %>% 
  mutate(
    time_since_cap = case_when(
      time_snc_full_disturbance == 0  ~ 1L,
      time_snc_full_disturbance == 7  ~ 6L,
      time_snc_full_disturbance == 11 ~ 6L,
      time_snc_full_disturbance == 13 ~ 6L,
      TRUE ~ time_snc_full_disturbance
    )
  ) 


# get management values by year - planting is teh most important
mng_year <-  fread("outData/mng_intensity_year_CZ.csv") 

# Clean up input data -----------------------

# add a unique plot and subplot id per year - so i can easily calculate how many unique plots i have
dat_full <- dat_full %>% 
  mutate(plot_year = paste0(plot, year, '_'),
         subplot_year = paste0(subplot, year, '_'))

# select only overlapping plots
dat_overlap <- dat_full %>%
  filter(status == "both")

plots_per_year <- dat_overlap %>%
  distinct(plot, year) %>%
  count(plot) %>%
  count(n, name = "n_plots") %>%
  rename(n_years = n)

plots_per_year


# get plots with overlapping years  
plots_with_both <- dat_overlap %>%
  filter(status == "both") %>%
  distinct(plot, year) %>%
  count(plot) %>%
  filter(n == 2) %>%
  pull(plot)


dat_overlap <- dat_overlap %>% 
  filter(plot %in% plots_with_both)


n_overlap_plots    <- length(plots_with_both)  #  125
n_overlap_subplots <- length(unique(dat_overlap$subplot))  #  1250

n_total_plots   <- length(unique(dat_full$plot_year))  #  333
n_subplots_total <-length(unique(dat_full$subplot_year))  # 1665


# convert management by year into long format
mng_year_long_all <- mng_year %>%
  select(plot, ends_with("_2023"), ends_with("_2025"),
         -starts_with("delta_")) %>%
  pivot_longer(
    cols         = -plot,
    names_to     = c("variable", "year"),
    names_pattern = "(.+)_(\\d{4})",
    values_to    = "intensity"
  ) %>%
  mutate(year = as.integer(year))

mng_year_wide_clean <- mng_year_long_all %>%
  pivot_wider(
    names_from  = variable,
    values_from = intensity
  )


# Data overview: landscape ------------------
## ALL plots ----------------------
# get master table, having all unique plots and subplots - even teh empty ones
dat_master_subplot <- dat_full %>% 
  dplyr::select(plot, subplot, year, status) %>% 
  distinct()

table(dat_master_subplot$year)  
table(dat_master_subplot$status)  


## ── Total stems and plots: overall and per year ───────────────────────────
n_trees_total <- sum(dat_full$n, na.rm = T)
n_trees_overlap <- sum(dat_overlap$n, na.rm = T)

n_trees_total
n_trees_overlap

# per year
trees_per_year <- dat_full %>%
  group_by(year) %>%
  summarise(
    total_stems    = sum(n, na.rm = T),
    .groups = "drop"
  )





## Overlapping plots  -------------------------
### Species composition --------------------------------------------------

### Get summary across all stems and study sites  
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

# Species x year summary + year-specific share
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
v_top_species <- species_levels # c('other', species_levels)
v_top_species_other <- c('other', species_levels)

# add 'other' to rename spacies and keep them there
species_stem_share_year_other <- 
  species_stem_share_year %>%
  mutate(species_other = ifelse(species %in% v_top_species, species, "other")) %>%
  group_by(year, species_other) %>%
  #  View()
  summarise(
    stems         = sum(stems, na.rm = TRUE),
    plots_present = sum(plots_present),  
    total_trees   = first(total_trees),
    .groups = "drop"
  ) %>%
  mutate(share = round(100 * stems / total_trees, 2))



#### Single species plots chars  -----------

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
  mutate(share = single_plots/n_overlap_plots )

# Print result
print(summary_species)


#Make sure every species appears in both years (missing gets 0)
top10_by_year_clean <- species_stem_share_year_other %>%
  dplyr::select(year, species_other, share) %>%
  tidyr::complete(year, species_other = v_top_species_other, fill = list(share = 0)) %>%
  mutate(
    species = factor(species_other, levels = v_top_species_other),
    year = factor(year)  # for alpha legend control
  )

# 3) Plot: same hue per species, lighter/darker by year
p_bar <- 
  top10_by_year_clean %>%
  ggplot(aes(x = share, y = species, fill = species_other, alpha = year)) +
  geom_col(position = position_dodge(width = 0.7), 
           width = 0.6, 
           aes(colour = factor(year))) + # "grey50"
  scale_x_continuous(labels = label_number(accuracy = 1, suffix = ""),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_discrete(
    limits = v_top_species_other,
    labels = species_labels_other,   # keep your mapping
    drop = FALSE
  ) +
  scale_colour_manual(
    values = c("2023" = "grey50",
               "2025" = "black"),
    guide = "none"   # prevents second legend
  )+
  scale_fill_manual(values = species_colors_other, guide = "none") +
  scale_alpha_manual(
    name = "Year of inventory",
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

species_occurence <- 
  df %>%
  mutate(species_other = ifelse(species %in% v_top_species, species, "other")) %>%
  ungroup(.) %>% 
  dplyr::filter(n > 0) %>%                 # Only where species occurred
  distinct(species_other, year, plot) %>%    #year,       # Unique species × plot combos
  count(species_other, year, name = "n_plots")  %>% #year, %>%  # Count number of plots per species
  mutate(share_of_plots = n_plots / n_overlap_plots  *100,
         species_other  = factor(species_other, levels = v_top_species_other)) %>% 
  arrange(species_other)


species_occurence

# Plot
p_occurence <- 
  species_occurence %>% 
  filter(species_other %in% v_top_species_other ) %>% 
  mutate(species = factor(species_other  , levels = v_top_species_other)) %>% 
  ggplot(aes(x = share_of_plots, 
             y = species_other, 
             fill = species_other,
             alpha = factor(year))) +
  geom_col(position = position_dodge(width = 0.7), 
           width = 0.6,
           aes(colour = factor(year))) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), 
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_manual(values = species_colors_other, 
                    guide = "none") +
  scale_alpha_manual(
    name = "Year of inventory",
    values = c("2025" = 1,
               "2023" = 0.45),
    breaks = c("2025", "2023")
  ) +
  scale_colour_manual(
    name = "Year of inventory",
    values = c("2025" = "black",
               "2023" = "grey50"),
    breaks = c("2025", "2023")
  ) +
  labs(
    x = "Plot share [%]",
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




### Summary table  -----------------------
# ── bind top species + other ───────────────────────────────────────────────
species_table <- 
  species_stem_share_year_other %>%
  filter(species_other %in% v_top_species_other) %>%
  select(!plots_present) %>% # can be removed and teh number of plots is calculated in plot occurence
  left_join(
    species_occurence,
    by = c("species_other", "year")
  )%>%
  #head()
  rename(share_stems = share ) %>%
  mutate(species = factor(species_other, levels = rev(v_top_species_other)))

# ── wide format, order by stems_2025, other last ──────────────────────────
species_table_wide <- 
  species_table %>%
    select(-total_trees) %>% 
  pivot_wider(
    names_from  = year,
    values_from = c(stems, share_stems, n_plots, share_of_plots),
    names_glue  = "{.value}_{year}"
  ) %>%
  mutate(is_other = species_other == "other") %>%
  arrange(is_other, desc(stems_2025)) %>%
  select(-is_other)


# ── format and export ─────────────────────────────────────────────────────

species_table_wide %>%
  mutate(species_other = recode(as.character(species_other),
                                "piab"  = "Picea abies",
                                "besp"  = "Betula sp.",
                                "pisy"  = "Pinus sylvestris",
                                "qusp"  = "Quercus sp.",
                                "fasy"  = "Fagus sylvatica",
                                "potr"  = "Populus tremula",
                                "absp"  = "Abies sp.",
                                "acps"  = "Acer pseudoplatanus",
                                "lade"  = "Larix decidua",
                                "saca"  = "Salix caprea",
                                "soau"  = "Sorbus aucuparia",
                                "other" = "Other species"
  )) %>%
  mutate(
    stems_2023_lab = paste0(stems_2023, " (", round(share_stems_2023, 1), "%)"),
    stems_2025_lab = paste0(stems_2025, " (", round(share_stems_2025, 1), "%)"),
    plots_2023_lab = paste0(n_plots_2023, " (", round(share_of_plots_2023, 1), "%)"),
    plots_2025_lab = paste0(n_plots_2025, " (", round(share_of_plots_2025, 1), "%)"),
    delta_stems = paste0(ifelse(stems_2025 - stems_2023 > 0, "+", ""),
                         stems_2025 - stems_2023),
    delta_share = paste0(ifelse(share_stems_2025 - share_stems_2023 > 0, "+", ""),
                         round(share_stems_2025 - share_stems_2023, 1), "%"),
    delta_plots = paste0(ifelse(n_plots_2025 - n_plots_2023 > 0, "+", ""),
                         n_plots_2025 - n_plots_2023)
  ) %>%
  select(species_other,
         stems_2023_lab, plots_2023_lab,
         stems_2025_lab, plots_2025_lab,
         delta_stems, delta_share, delta_plots) %>%
  flextable() %>%
  set_header_labels(
    species_other  = "Species",
    stems_2023_lab = "Stems n (%)",
    plots_2023_lab = "Plots n (%)",
    stems_2025_lab = "Stems n (%)",
    plots_2025_lab = "Plots n (%)",
    delta_stems    = "Stems",
    delta_share    = "Share",
    delta_plots    = "Plots"
  ) %>%
  add_header_row(
    values    = c("", "2023", "2025", "\u0394 2023\u21922025"),
    colwidths = c(1, 2, 2, 3)
  ) %>%
  italic(j = 1, i = 1:11) %>%
  align(align = "right", j = 2:8, part = "all") %>%
  align(align = "left",  j = 1,   part = "all") %>%
  hline(i = 11) %>%
  autofit() %>%
  save_as_docx(path = "outTable/species_table.docx")







### Drought tolerance classes ------------------------------


# Recode TSD outliers — applied once, used everywhere
dat_overlap_recoded <- dat_overlap %>%
  mutate(
    time_since_cap = case_when(
      time_snc_full_disturbance == 0  ~ 1L,
      time_snc_full_disturbance == 7  ~ 6L,
      time_snc_full_disturbance == 11 ~ 6L,
      time_snc_full_disturbance == 13 ~ 6L,
      TRUE ~ time_snc_full_disturbance
    )
  ) %>%
  left_join(species_functional %>% select(species, func_group), by = "species") %>%
  mutate(
    n          = if_else(is.na(n), 0L, n),
    func_group = replace_na(func_group, "other")
  )


# Get species traits -------------------------
species_traits <- dat_overlap_recoded %>%
  filter(!is.na(species), species != "",
         !is.na(Drought_tolerance)) %>%
  group_by(species) %>%
  summarise(
    drought_tol = mean(Drought_tolerance, na.rm = TRUE),
    shade_tol   = mean(Shade_tolerance,   na.rm = TRUE),
    n_obs       = n(),
    .groups = "drop"
  ) %>%
  arrange(drought_tol)

# ── Look at the full distribution first ───────────────────────────────────────
print(species_traits, n = Inf)
summary(species_traits$drought_tol)
hist(species_traits$drought_tol, breaks = 15,
     main = "Species drought tolerance (Niinemets)",
     xlab = "Drought tolerance (higher = more tolerant)")


# ── Data-driven functional classification using Niinemets drought tolerance ───
species_functional_v2 <- species_traits %>%
  mutate(
    func_group_drought = case_when(
      # Norway spruce gets its own group — pre-disturbance dominant,
      # ecologically distinct regardless of trait value
      species == "piab"        ~ "Norway_spruce",
      drought_tol < 2.2        ~ "drought_sensitive",   # clear gap before frex
      drought_tol > 3.3        ~ "drought_tolerant",    # clear gap after osca
      TRUE                     ~ "intermediate"
    ),
    func_group_drought = factor(func_group_drought,
                                levels = c("Norway_spruce", "drought_sensitive",
                                           "intermediate", "drought_tolerant"))
  )

# verify — check the assignment
species_functional_v2 %>%
  select(species, drought_tol, shade_tol, func_group_drought) %>%
  arrange(func_group_drought, drought_tol) %>%
  print(n = Inf)

# ── Update the trait space plot with group colors ─────────────────────────────
drought_labels <- c(
  "Norway_spruce"     = "Norway spruce",
  "drought_sensitive" = "Drought-sensitive",
  "intermediate"      = "Intermediate",
  "drought_tolerant"  = "Drought-tolerant",
  "no_regeneration"   = "No regeneration"
)


drought_class_title = "Drought tolerance\nclasses"

drought_cl_levels <- c("Norway_spruce", "drought_sensitive",
                       "intermediate", "drought_tolerant",
                       "no_regeneration")

drought_colors <- c(
  "Norway_spruce"     = "#1a5c1a",   # dark forest green — matches piab in species fig
  "drought_sensitive" = "#a8d9a8",   # light green — low drought, shade-tolerant
  "intermediate"      = "#f4a736",   # warm orange — matches mid species colors
  "drought_tolerant"  = "#d73027",   # deep red — high drought, matches Abies/rare end
  "no_regeneration"   = "#d9d9d9"    # grey
)


shared_fill <- scale_fill_manual(
  values = drought_colors,
  name   =  drought_class_title,# "Drought tolerance\nclass",
  breaks = drought_cl_levels,
  labels = drought_labels
)



top11_species <- unique(top10_by_year$species)  # gives your 11 species across both years

species_functional_v2 <- species_functional_v2 %>%
  mutate(is_top11 = species %in% top11_species)

p_function_drought <- ggplot(species_functional_v2,
                             aes(x = drought_tol, y = shade_tol,
                                 color = func_group_drought, label = species)) +
  geom_point(aes(size = is_top11, alpha = is_top11)) +
  ggrepel::geom_text_repel(
    data = . %>% filter(is_top11),
    aes(label = species),
    size = 3.5, fontface = "bold"
  ) +
  ggrepel::geom_text_repel(
    data = . %>% filter(!is_top11),
    aes(label = species),
    size = 2.5, color = "grey60", alpha = 0.7
  ) +
  scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3.5), guide = "none") +
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none") +
  geom_vline(xintercept = c(2.2, 3.3), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = drought_colors,
                     name = drought_class_title,
                     labels = drought_labels) +
  guides(color = guide_legend(override.aes = list(label = ""))) +
  labs(x = "Drought tolerance (Niinemets)",
       y = "Shade tolerance (Niinemets)") +
  theme_classic()


p_function_drought



# ── Join new classification ───────────────────────────────────────────────────
dat_overlap_recoded2 <- dat_overlap_recoded %>%
  left_join(
    species_functional_v2 %>% select(species, func_group_drought),
    by = "species"
  ) %>%
  mutate(
    func_group_drought = replace_na(as.character(func_group_drought), "other")
  )


hist(dat_overlap_recoded2$Drought_tolerance)

# Functional groups based on drought tolerance

# ── Update compute_func_shares for new group names ────────────────────────────
compute_func_shares_v2 <- function(df) {
  df %>%
    mutate(
      total_stems          = Norway_spruce + drought_sensitive + 
        intermediate + drought_tolerant + other,
      share_spruce         = if_else(total_stems == 0, 0, Norway_spruce    / total_stems),
      share_drought_sens   = if_else(total_stems == 0, 0, drought_sensitive / total_stems),
      share_intermediate   = if_else(total_stems == 0, 0, intermediate      / total_stems),
      share_drought_tol    = if_else(total_stems == 0, 0, drought_tolerant  / total_stems),
      share_other          = if_else(total_stems == 0, 0, other             / total_stems)
    )
}

# ── classify_dom_func for new groups ───────────────────────────────────
classify_dom_func_v2 <- function(share_spruce, share_drought_sens,
                                 share_intermediate, share_drought_tol,
                                 total_stems) {
  case_when(
    total_stems == 0 ~ "no_regeneration",
    share_spruce       == pmax(share_spruce, share_drought_sens,
                               share_intermediate, share_drought_tol) ~ "Norway_spruce",
    share_drought_tol  == pmax(share_spruce, share_drought_sens,
                               share_intermediate, share_drought_tol) ~ "drought_tolerant",
    share_intermediate == pmax(share_spruce, share_drought_sens,
                               share_intermediate, share_drought_tol) ~ "intermediate",
    TRUE                                                               ~ "drought_sensitive"
  )
}


# ── Rebuild base stem table with new classification ───────────────────────────

dat_full_func <- dat_full %>%
  left_join(
    species_functional_v2 %>% select(species, func_group_drought),
    by = "species"
  ) %>%
  mutate(
    func_group_drought = replace_na(as.character(func_group_drought), "other")
  )



# get only for overlapping sites -> for alluvial
func_stems_base_v2 <- dat_overlap_recoded2 %>%
 # filter(status == "both") %>%
  mutate(func_group_drought = replace_na(as.character(func_group_drought), 
                                         "other")) %>%
  group_by(plot, year, time_since_cap, func_group_drought) %>%
  summarise(stems = sum(n, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = func_group_drought,
              values_from = stems,
              values_fill = 0) %>%
  compute_func_shares_v2()

# from all sites -> for time plots
func_stems_base_full <- dat_full_func %>%
  # filter(status == "both") %>%
  mutate(func_group_drought = replace_na(as.character(func_group_drought), 
                                         "other")) %>%
  group_by(plot, year, time_since_cap, func_group_drought) %>%
  summarise(stems = sum(n, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = func_group_drought,
              values_from = stems,
              values_fill = 0) %>%
  compute_func_shares_v2()



func_tsd_bar_v2 <- func_stems_base_full %>%
  group_by(plot, time_since_cap) %>%
  summarise(
    across(starts_with("share_"), ~ mean(.x, na.rm = TRUE)),
    total_stems = mean(total_stems, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dom_func = factor(
      classify_dom_func_v2(share_spruce, share_drought_sens,
                           share_intermediate, share_drought_tol,
                           total_stems),
      levels = drought_cl_levels
    ),
    tsd     = factor(time_since_cap),
    tsd_num = time_since_cap
  ) %>%
  group_by(tsd) %>%
  mutate(n_plots_tsd = n_distinct(plot)) %>%
  ungroup() %>%
  count(tsd, tsd_num, dom_func, n_plots_tsd) %>%
  mutate(pct = n / n_plots_tsd * 100)

# get text for classes
func_tsd_n_v2 <- func_tsd_bar_v2 %>% distinct(tsd, tsd_num, n_plots_tsd)

func_tsd_bar_v2 <- func_tsd_bar_v2 %>%
  mutate(dom_func = factor(dom_func, levels = drought_cl_levels))


# Option B: stacked bar 
p_bar_TSD <- ggplot(func_tsd_bar_v2,
                    aes(x = tsd_num, y = pct, fill = dom_func)) +
  geom_col(color = "black", linewidth = 0.3) +
  shared_fill +
  scale_x_continuous(breaks = 1:6,
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(limits = c(0, 110),          # bars capped at 100%
                     expand = expansion(mult = c(0, 0))) +
  geom_text(data = func_tsd_n_v2,
            aes(x = tsd_num, y = 104, label = n_plots_tsd),
            inherit.aes = FALSE, size = 2.8, color = "grey50") +
  coord_cartesian(ylim = c(0, 100), clip = "off") +   # clip="off" here +                        # allow text outside plot area
  labs(x = "Time since disturbance\n(years)", y = "Share of plots [%]") +
  theme(legend.position  = "right",
        plot.margin      = margin(t = 40, r = 5, b = 5, l = 5))  # room for n labels
p_bar_TSD


# ── 2. Aggregate stems per plot × TSD × drought group ────────────────────────
func_tsd_stems <- dat_full_func %>%
  filter(!is.na(n), n > 0) %>%
  group_by(plot, time_since_cap, func_group_drought) %>%
  summarise(stems = sum(n, na.rm = TRUE), .groups = "drop") %>%
  group_by(time_since_cap, func_group_drought) %>%
  summarise(stems   = sum(stems),
            n_plots = n_distinct(plot),
            .groups = "drop") %>%
  group_by(time_since_cap) %>%
  mutate(share = stems / sum(stems) * 100) %>%
  ungroup() %>% 
  mutate(func_group_drought = factor(func_group_drought, 
                                     levels = drought_cl_levels)) 


# ── 3. Plot ───────────────────────────────────────────────────────────────────
p_bar_TSD_stems_share <- ggplot(func_tsd_stems,
                          aes(x = time_since_cap,
                              y = share,
                              fill = func_group_drought)) +
  geom_col(color = "black", linewidth = 0.3) +
  shared_fill +
  scale_x_continuous(breaks = 1:6,
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  # geom_text(data = func_tsd_stems_n,
  #           aes(x = time_since_cap, y = 104, label = n_plots),
  #           inherit.aes = FALSE, size = 2.8, color = "grey50") +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(x = "Time since disturbance\n(years)", 
       y = "Stem share [%]") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right",
        plot.margin     = margin(t = 40, r = 5, b = 5, l = 5))

p_bar_TSD_stems_share

p_tsd_merged <- ggarrange( p_bar_TSD_stems_share, 
                           p_bar_TSD,
                           common.legend = T,
                           labels = c("[a]", "[b]"),
                           legend = 'right',
                           font.label = list(size = 10, face = "plain"))
p_tsd_merged


pdf("outFigsCZ/p_tsd_merged_corel.pdf", width = 7, height = 4)
print(p_tsd_merged)
dev.off()





##  Functional group alluvial: drought classification 2023 → 2025 --------------
func_alluvial_v2 <- func_stems_base_v2 %>%
  mutate(
    dom_func = classify_dom_func_v2(share_spruce, share_drought_sens,
                                    share_intermediate, share_drought_tol,
                                    total_stems)
  ) %>%
  select(plot, year, dom_func) %>%
  pivot_wider(names_from = year, values_from = dom_func, names_prefix = "func_") %>%
  filter(!is.na(func_2023), !is.na(func_2025)) %>%
  mutate(
    func_2023 = factor(func_2023, levels = drought_cl_levels),  # codes only
    func_2025 = factor(func_2025, levels = drought_cl_levels)   # codes only
  ) %>%
  count(func_2023, func_2025, name = "n")


p_func_alluvial_v2 <- ggplot(func_alluvial_v2,
                             aes(axis1 = func_2023,
                                 axis2 = func_2025,
                                 y = n)) +
  geom_alluvium(aes(fill = func_2023),
                width = 0.4, alpha = 0.7, knot.pos = 0.4) +
  geom_stratum(aes(fill = after_stat(stratum)),
               width = 0.4, color = "black", linewidth = 0.5) +
  geom_text(stat = "stratum",
            aes(label = ifelse(after_stat(prop) >= 0.02,
                               paste0(round(after_stat(prop) * 100, 1)), # , "%"
                               "")),
            size = 3, color = "grey20") +
  scale_x_discrete(limits = c("2023", "2025"),
                   name = 'Year of inventory'
                   #expand = expansion(mult = c(0.25, 0.25)
                   #                  )
  ) +
  scale_fill_manual(values = drought_colors, guide = "none") +  # no legend
  labs(x = NULL, y = "Number of plots") +
  theme(
    axis.text.x  = element_text(),
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank()
  )
p_func_alluvial_v2


# verify levels are codes not display labels
levels(func_alluvial_v2$func_2023)
# should show: "Norway_spruce" "drought_sensitive" "intermediate" "drought_tolerant" "no_regeneration"


p_combined_func <- (p_func_alluvial_v2 | p_bar_TSD) +
  plot_layout(widths = c(1.2, 1), heights = c(1.2,1)) +
  plot_annotation(tag_levels = list(c("[a]", "[b]")))# &
p_combined_func


# export again


pdf("outFigsCZ/p_combined_func_laura.pdf", width = 7, height = 3.5)
print(p_combined_func)
dev.off()




# # ── Data ──────────────────────────────────────────────────────────────────────
func_bar_data <- func_stems_base_v2 %>%
  filter(!is.na(year)) %>%
  group_by(year) %>%
  summarise(
    Norway_spruce     = sum(Norway_spruce,     na.rm = TRUE),
    drought_sensitive = sum(drought_sensitive, na.rm = TRUE),
    intermediate      = sum(intermediate,      na.rm = TRUE),
    drought_tolerant  = sum(drought_tolerant,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(total = rowSums(across(c(Norway_spruce, drought_sensitive,
                                  intermediate, drought_tolerant))),
         across(c(Norway_spruce, drought_sensitive,
                  intermediate, drought_tolerant),
                ~ .x / total * 100)) %>%
  select(-total) %>%
  pivot_longer(cols      = -year,
               names_to  = "group",
               values_to = "share") %>%
  mutate(
    group = factor(group, levels = rev(drought_cl_levels)),  # rev for y-axis order
    year  = factor(year)
  )

# ── Plot ──────────────────────────────────────────────────────────────────────

# get order from 2025 stem shares
drought_group_order <- func_bar_data %>%
  filter(year == "2025") %>%
  arrange(share) %>%          # ascending so top of y-axis = highest
  pull(group) %>%
  as.character()

func_bar_data <- func_bar_data %>%
  mutate(group = factor(group, levels = drought_group_order))


p_bar_drought <- ggplot(func_bar_data,
                        aes(x = share, y = group,
                            fill = group, alpha = year)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           aes(colour = factor(year))) +
  scale_fill_manual(values = drought_colors, guide = "none") +
  scale_alpha_manual(
    name   = "Year of inventory",
    values = c("2023" = 0.45, "2025" = 1.00),
    breaks = c("2025", "2023")
  ) +
  scale_colour_manual(
    name   = "Year of inventory",
    values = c("2023" = "grey50", "2025" = "black"),
    breaks = c("2025", "2023")
  ) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_discrete(
    labels = drought_labels
  ) +
  labs(x = "Stem share [%]", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y      = element_text(size = 9),
    legend.position = "none"
  )

p_bar_drought

# p occurence drought
func_occurence_data <- func_stems_base_v2 %>%
  filter(!is.na(year)) %>%
  mutate(
    Norway_spruce_occ     = Norway_spruce     > 0,
    drought_sensitive_occ = drought_sensitive > 0,
    intermediate_occ      = intermediate      > 0,
    drought_tolerant_occ  = drought_tolerant  > 0
  ) %>%
  group_by(year) %>%
  summarise(
    Norway_spruce     = sum(Norway_spruce_occ,     na.rm = TRUE),
    drought_sensitive = sum(drought_sensitive_occ, na.rm = TRUE),
    intermediate      = sum(intermediate_occ,      na.rm = TRUE),
    drought_tolerant  = sum(drought_tolerant_occ,  na.rm = TRUE),
    n_plots           = n(),
    .groups = "drop"
  ) %>%
  mutate(across(c(Norway_spruce, drought_sensitive,
                  intermediate, drought_tolerant),
                ~ .x / n_plots * 100)) %>%
  select(-n_plots) %>%
  pivot_longer(cols     = -year,
               names_to  = "group",
               values_to = "share_plots") %>%
  mutate(
    group = factor(group, levels = rev(drought_cl_levels)),
    year  = factor(year)
  )

# ── Plot ──────────────────────────────────────────────────────────────────────


# apply to both plots

func_occurence_data <- func_occurence_data %>%
  mutate(group = factor(group, levels = drought_group_order))


p_occurence_drought <- ggplot(func_occurence_data,
                              aes(x = share_plots, y = group,
                                  fill = group, alpha = year)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           aes(colour = factor(year))) +
  scale_fill_manual(values = drought_colors, guide = "none") +
  scale_alpha_manual(
    name   = "Year of inventory",
    values = c("2023" = 0.45, "2025" = 1.00),
    breaks = c("2025", "2023")
  ) +
  scale_colour_manual(
    name   = "Year of inventory",
    values = c("2023" = "grey50", "2025" = "black"),
    breaks = c("2025", "2023")
  ) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_discrete(labels = drought_labels) +
  labs(x = "Plot share [%]", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = 'none',
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank()
  #   axis.text.y          = element_blank(),
  #   legend.position      = c(0.98, 0.02),
  #   legend.justification = c(1, 0),
  #   legend.key           = element_rect(fill = NA),
  #   legend.title         = element_text(size = 9),
  #   legend.text          = element_text(size = 8),
  #   legend.key.width     = unit(1,   "lines"),
  #   legend.key.height    = unit(0.4, "lines")
   )

p_occurence_drought





p_combined_drought <- ggarrange(
  #p_bar, p_occurence,
  p_bar_drought, p_occurence_drought,
  #p_func_alluvial_v2,
  ncol = 2, common.legend = FALSE,
  align = "h",
  widths = c(1.5, 1),
  labels = c("[c]", "[d]"),
  font.label = list(size = 10, face = "plain"),
  label.x = 0.02,
  label.y = 1.01
)

p_combined_drought

p_species_composition



p_fig1_combined <- ggarrange(
  p_species_composition,
         
          p_func_alluvial_v2,
  p_combined_drought,
  align = "hv",
          widths = c(3, 1)
  )

p_fig1_combined

pdf("outFigsCZ/p_fig1_combined_corel.pdf", width = 7, height = 6)
print(p_fig1_combined)
dev.off()


pdf("outFigsCZ/p_fig1_alluvial_corel.pdf", width = 4, height = 3)
print(p_func_alluvial_v2)
dev.off()

### ── Climate adaptation score per plot -----------------
# higher = more climate adapted
# use 2025 composition as the outcome (current state)

climate_adapt <- func_stems_base_v2 %>%
  filter(year == 2025) %>%
  mutate(
    # climate adapted = drought tolerant + intermediate
    # climate maladapted = Norway spruce + drought sensitive
    share_adapted    = share_drought_tol + share_intermediate,
    share_maladapted = share_spruce + share_drought_sens
  ) %>%
  left_join(
    mng_year_wide_clean %>%
      select(plot, year, planting_intensity, anti_browsing_intensity),
    by = c("plot", "year")
  ) %>%
  mutate(
    planted = planting_intensity > 0
  )

# ── Test: is climate adaptation higher on planted plots? ─────────────────────
wilcox.test(share_adapted ~ planted, data = climate_adapt)

# ── Summary ───────────────────────────────────────────────────────────────────
climate_adapt %>%
  group_by(planted) %>%
  summarise(
    n                 = n(),
    mean_adapted      = round(mean(share_adapted,    na.rm = TRUE), 2),
    mean_maladapted   = round(mean(share_maladapted, na.rm = TRUE), 2),
    mean_spruce       = round(mean(share_spruce,     na.rm = TRUE), 2),
    .groups = "drop"
  )

# ── Continuous: planting intensity ~ share adapted ────────────────────────────
cor.test(~ planting_intensity + share_adapted,
         data = climate_adapt, method = "spearman")

# ── Visualize ─────────────────────────────────────────────────────────────────
ggplot(climate_adapt,
       aes(x = planting_intensity, y = share_adapted)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#1a5c1a") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  labs(x = "Planting intensity ",
       y = "Share of climate-adapted stems\n(drought-tolerant + intermediate)") +
  theme_classic()




## 6b. Management intensity plot (simpler stacked bar)
fill_colors <- brewer.pal(length(intensity_levels), "YlOrRd")



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




# make management intensity plot simpler
p_management_intensity_plot_simpler <- 
  mng_shifted %>% 
  filter(activity != "logging_trail_intensity") %>% 
  droplevels(.) %>%
  mutate(intensity_class_plot = factor(intensity_class, 
                                       levels = rev(intensity_levels))) %>%  # ← fix here
  ggplot(aes(x = proportion , 
             y = activity,
             fill = intensity_class_plot)) +
  geom_vline(xintercept = c(25,50,75), color = "grey80", linewidth = 0.5, lty = "dashed")+ # ← add
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






## Disturbance history ---------------------------------------
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



#  Models (GAMs) -----------------------------------------------
# All GAMs: REML, random effect s(plot_id, bs="re")
# Families: tw(log) continuous, nb(log) counts, binomial(logit) binary

# make trully cross scale: get management intensity on pot level, 
# management binary on subplot level

## ── Prepare cross-scale data  --------------- ─
# subplot level: use subplot binary management
# plot level: use plot intensity management
# both in one dataset with level as factor

# ── Step 1: get subplot-level response variables + binary management ──────────
# both_levels_cz subplot rows already have the response variables (mean_hgt etc.)
# but we need to add binary management from dat_overlap

# get binary management per subplot
mng_subplot_binary <- dat_full %>% # from full data, not only overlaps
  #filter(status == "both") %>%
  distinct(plot, subplot, year, planting, anti_browsing, grndwrk) %>%
  rename(
    plot_id           = plot,
    planting_bin      = planting,
    anti_browsing_bin = anti_browsing,
    grndwrk_bin       = grndwrk
  ) %>%
  mutate(plot_id = factor(plot_id))

# ── Step 2: join binary management to subplot rows in both_levels_cz ──────────
both_levels_crossscale <- both_levels_cz %>%
  rename(subplot = ID) %>%
  left_join(mng_subplot_binary,
            by = c("plot_id", 
                   "subplot",
                   "year")) %>%
  mutate(
    # subplot level: use binary (0/1) from dat_overlap
    # plot level: use intensity (0-1) already in both_levels_cz
    # for plot rows: planting_bin will be NA — use intensity instead
    planting_pred = case_when(
      level == "subplot" & !is.na(planting_bin) ~ as.numeric(planting_bin),
      level == "plot"                            ~ planting_intensity,
      TRUE                                       ~ planting_intensity  # fallback
    ),
    browsing_pred = case_when(
      level == "subplot" & !is.na(anti_browsing_bin) ~ as.numeric(anti_browsing_bin),
      level == "plot"                                 ~ anti_browsing_intensity,
      TRUE                                            ~ anti_browsing_intensity
    ),
    grndwrk_pred = case_when(
      level == "subplot" & !is.na(grndwrk_bin) ~ as.numeric(grndwrk_bin),
      level == "plot"                           ~ grndwrk_intensity,
      TRUE                                      ~ grndwrk_intensity
    )
  )





##  Share adapted at plot level -----------

# Test share adapted
func_stems_base_all <- dat_full_func %>%
  mutate(func_group_drought = replace_na(as.character(func_group_drought), 
                                         "other")) %>%
  group_by(plot, year, time_snc_full_disturbance, func_group_drought) %>%
  summarise(stems = sum(n, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = func_group_drought,
              values_from = stems,
              values_fill = 0) %>%
  compute_func_shares_v2()


# calculate climate adaptation on plot level (as having 1-2 species per subplot level will be not representative for community level)
share_adapted_plot <- func_stems_base_all %>%
  mutate(
    share_adapted    = share_drought_tol + share_intermediate,
    share_maladapted = share_spruce + share_drought_sens
  ) %>%
  select(plot, year, time_snc_full_disturbance, 
         share_adapted, share_maladapted, share_spruce) %>%
  mutate(
    plot_id = factor(plot),
    year_f  = factor(year)
  ) %>%
  left_join(
    both_levels_cz %>%
      filter(level == "plot") %>%
      select(plot_id, year, planting_intensity, anti_browsing_intensity,
             grndwrk_intensity) %>%
      distinct(),
    by = c("plot_id", "year")
  ) %>% 
  mutate(  # bound between 0-1
    share_adapted_adj = case_when(
      share_adapted == 0 ~ 0.001,
      share_adapted == 1 ~ 0.999,
      TRUE               ~ share_adapted
    )
  )



## Spruce shares at plot level -------------------------
spruce_share_plot <- both_levels_cz %>%
  filter(level == "plot") %>%
  select(plot_id, year, year_f, time_snc_full_disturbance,
         spruce_share, planting_intensity, anti_browsing_intensity,
         grndwrk_intensity) %>%
  distinct() %>%
  mutate(
    plot = as.character(plot_id),
    # squeeze for beta family
    spruce_share_adj = case_when(
      spruce_share == 0  ~ 0.001,
      spruce_share == 1  ~ 0.999,
      is.na(spruce_share) ~ NA_real_,
      TRUE               ~ spruce_share
    )
  ) %>%
  filter(!is.na(spruce_share_adj))



## collinearity check -----------------------
# ── Plot level ────────────────────────────────────────────────────────────────
cor.test(~ planting_intensity + anti_browsing_intensity,
         data = both_levels_cz %>% filter(level == "plot"),
         method = "spearman")

# ── Subplot level ─────────────────────────────────────────────────────────────
cor.test(~ planting_pred + browsing_pred,
         data = both_levels_crossscale %>% filter(level == "subplot"),
         method = "spearman")

# ── Full dataset (both levels) ────────────────────────────────────────────────
cor.test(~ planting_pred + browsing_pred,
         data = both_levels_crossscale,
         method = "spearman")





### Management combination composite?  -------------------

# ── Option 1: mean of planting and browsing 
# interpretation: average management intensity across the two key interventions
both_levels_crossscale <- both_levels_crossscale %>%
  mutate(
    mng_composite = (planting_pred + browsing_pred) / 2
  )

share_adapted_plot <- share_adapted_plot %>%
  mutate(
    mng_composite = (planting_intensity + anti_browsing_intensity) / 2
  )

spruce_share_plot <- spruce_share_plot %>%
  mutate(
    mng_composite = (planting_intensity + anti_browsing_intensity) / 2
  )

plot_df_cz <- plot_df_cz %>% 
  mutate(
    mng_composite = (planting_intensity + anti_browsing_intensity) / 2
  )



# Model comparison: browsing include or not? --------------------------------
# check wheather including browisng protection along with planing is actualy meaningful
# create models - only planting, only browsing, both effects and and tehir interaction 
# Yes, I should include browing on cross scale analysis!!


# ── Refit with ML for fair AIC comparison ────────────────────────────────────
compare_mng_aic <- function(response, data, family, k_tsd = 4) {
  
  base_formula <- as.formula(paste0(
    response, " ~ grndwrk_pred + year_f + level +",
    "s(time_snc_full_disturbance, k =", k_tsd, ") +",
    "s(plot_id, bs = 're')"
  ))
  
  m_both <- gam(update(base_formula, . ~ . + planting_pred + browsing_pred),
                data = data, family = family, method = "ML")
  m_inter <- gam(update(base_formula, . ~ . + planting_pred * browsing_pred),
                 data = data, family = family, method = "ML")
  m_plant <- gam(update(base_formula, . ~ . + planting_pred),
                 data = data, family = family, method = "ML")
  m_browse <- gam(update(base_formula, . ~ . + browsing_pred),
                  data = data, family = family, method = "ML")
  m_none  <- gam(base_formula,
                 data = data, family = family, method = "ML")
  
  tibble(
    response = response,
    model    = c("both", "interaction", "plant_only", 
                 "browse_only", "none"),
    AIC      = c(AIC(m_both), AIC(m_inter), AIC(m_plant),
                 AIC(m_browse), AIC(m_none)),
  ) %>%
    mutate(
      delta_AIC = AIC - min(AIC),
      best      = delta_AIC < 2
    ) %>%
    arrange(AIC)
}

### ── Run for each response Plot and subplot level ─────────────────────────────────────────────────────
aic_results_cross <- bind_rows(
  compare_mng_aic("mean_hgt",        
                  both_levels_crossscale %>% filter(mean_hgt < 6),
                  tw(link = "log"), k_tsd = 4),
  compare_mng_aic("cv_hgt_pos",      
                  both_levels_crossscale,
                  tw(link = "log"), k_tsd = 7),
  compare_mng_aic("effective_numbers",
                  both_levels_crossscale,
                  tw(link = "log"), k_tsd = 3),
  compare_mng_aic("sp_richness",     
                  both_levels_crossscale,
                  nb(link = "log"),  k_tsd = 7)
)

aic_results_cross %>%
  select(response, model, AIC, delta_AIC, best) %>%
  print(n = Inf)

  distinct()

#####  AIC comparison - plot level ---------------------
compare_mng_aic_plot <- function(response, data, family, k_tsd = 4) {
    base_formula <- as.formula(paste0(
      response, " ~ year_f +",
      "s(time_snc_full_disturbance, k =", k_tsd, ") +",
      "s(plot_id, bs = 're')"
    ))
    
    m_inter  <- gam(update(base_formula, . ~ . + planting_intensity * anti_browsing_intensity + grndwrk_intensity),
                    data = data, family = family, method = "ML")
    m_both   <- gam(update(base_formula, . ~ . + planting_intensity + anti_browsing_intensity + grndwrk_intensity),
                    data = data, family = family, method = "ML")
    m_plant  <- gam(update(base_formula, . ~ . + planting_intensity + grndwrk_intensity),
                    data = data, family = family, method = "ML")
    m_browse <- gam(update(base_formula, . ~ . + anti_browsing_intensity + grndwrk_intensity),
                    data = data, family = family, method = "ML")
    m_grndwrk <- gam(update(base_formula, . ~ . + grndwrk_intensity),
                     data = data, family = family, method = "ML")
    m_none   <- gam(base_formula, data = data, family = family, method = "ML")
    
    tibble(
      model = c("interaction", "both", "plant_only", 
                "browse_only", "grndwrk_only", "none"),
      AIC   = c(AIC(m_inter), AIC(m_both), AIC(m_plant),
                AIC(m_browse), AIC(m_grndwrk), AIC(m_none))
    ) %>%
      mutate(
        delta_AIC = AIC - min(AIC),
        best      = delta_AIC < 2
      ) %>%
      arrange(AIC)
  }
# run
# AIC for all three plot-level responses
aic_adapted <- compare_mng_aic_plot("share_adapted_adj",  share_adapted_plot, betar(),   k_tsd = 4)
aic_spruce  <- compare_mng_aic_plot("spruce_share_adj",   spruce_share_plot,  betar(),   k_tsd = 4)
aic_beta    <- compare_mng_aic_plot("beta_jaccard_mean",  plot_df_cz,         gaussian(), k_tsd = 3)
  
aic_adapted
aic_spruce
aic_beta


###  spruce
gam_spruce <- gam(
  spruce_share_adj ~
    planting_intensity+anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 4) +
    grndwrk_intensity + year_f +
    s(plot_id, bs = "re"),
  data   = spruce_share_plot,
  family = betar(),
  method = "REML"
)

# 
# update this based on AIC result — likely plant_only
gam_adapted_final <- gam(
  share_adapted_adj ~
    planting_intensity + anti_browsing_intensity +
    s(time_snc_full_disturbance, k = 6) +
    grndwrk_intensity + year_f +
    s(plot_id, bs = "re"),
  data   = share_adapted_plot,
  family = betar(),
  method = "REML"
)

summary(gam_adapted_final)
gratia::draw(gam_adapted_final, select = 1)





summary(gam_spruce)
gratia::draw(gam_spruce, select = 1)

# ── Effect size: planting 0 vs 1 ─────────────────────────────────────────────
emmeans(gam_spruce, ~ planting_intensity,
        at   = list(planting_intensity = c(0, 1)),
        type = "response") %>%
  summary(infer = TRUE)



# ── Collinearity check ────────────────────────────────────────────────────────
cor.test(~ planting_intensity + anti_browsing_intensity,
         data = share_adapted_plot, method = "spearman")

# ── Effect size at planting 0 vs 1 ───────────────────────────────────────────
emmeans(gam_adapted_final, ~ planting_intensity,
        at   = list(planting_intensity = c(0, 1)),
        type = "response") %>%
  summary(infer = TRUE)




### Refit GAMs with cross-scale management predictors  ---------------------------
gam_mean_hgt_cross <- gam(
  mean_hgt ~
    planting_pred * browsing_pred +
    s(time_snc_full_disturbance, k = 4) +
    grndwrk_pred + year_f + level +
    s(plot_id, bs = "re"),
  data   = both_levels_crossscale %>% filter(mean_hgt < 6),
  family = tw(link = "log"),
  method = "REML"
)

plot(gam_mean_hgt_cross, page = 1)

gam_cv_hgt_pos_cross <- gam(
  cv_hgt_pos ~
    planting_pred * browsing_pred +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_pred + year_f + level +
    s(plot_id, bs = "re"),
  data   = both_levels_crossscale,
  family = tw(link = "log"),
  method = "REML"
)

gam_eff_cross <- gam(
  effective_numbers ~
    planting_pred * browsing_pred +
    s(time_snc_full_disturbance, k = 3) +
    grndwrk_pred + year_f + level +
    s(plot_id, bs = "re"),
  data   = both_levels_crossscale,
  family = tw(link = "log"),
  method = "REML"
)

gam_rich_cross <- gam(
  sp_richness ~
    planting_pred * browsing_pred +
    s(time_snc_full_disturbance, k = 7) +
    grndwrk_pred + year_f + level +
    s(plot_id, bs = "re"),
  data   = both_levels_crossscale,
  family = nb(link = "log"),
  method = "REML"
)


## Plot level: Beta diversity (Jaccard)
m_beta_add <- mgcv::gam(
  beta_jaccard_mean ~
    time_snc_full_disturbance +
    planting_intensity + anti_browsing_intensity +
    grndwrk_intensity + year_f +
    s(plot_id, bs = "re"),
  data   = plot_df_cz,
  method = "REML"
)

# ## Final model list
fin.models.all <- list(
  hgt   = gam_mean_hgt_cross,
  cvpos = gam_cv_hgt_pos_cross,
  eff   = gam_eff_cross,
  rich  = gam_rich_cross,
  beta  = m_beta_add,
  adapt = gam_adapted_final,
  spruce = gam_spruce
)

lapply(fin.models.all, summary)



# ── Response labels ───────────────────────────────────────────────────────────
response_labels <- c(
  hgt    = "Mean height [m]",
  cvpos  = "CV [%]",
  eff    = "Effective species [#]",
  rich   = "Species richness [#]",
  beta   = "Dissimilarity [dim.]",
  adapt  = "Climate-adapted share [%] ",
  spruce = "Norway spruce share [%]"
)

# ── Term labels ───────────────────────────────────────────────────────────────
term_labels <- c(
  planting_intensity      = "Planting",
  anti_browsing_intensity = "Browsing\nprotection",
  grndwrk_intensity       = "Soil\npreparation"
)


# ── Rerun pvals_mng ───────────────────────────────────────────────────────────
pvals_mng <- map_dfr(fin.models.all,
                     ~ broom::tidy(.x, parametric = TRUE),
                     .id = "model") %>%
  filter(term %in% c("planting_intensity", "anti_browsing_intensity",
                     "grndwrk_intensity",
                     "planting_pred", "browsing_pred", "grndwrk_pred")) %>%
  mutate(
    focal_term = recode(term,
                        "planting_pred"          = "planting_intensity",
                        "browsing_pred"          = "anti_browsing_intensity",
                        "grndwrk_pred"           = "grndwrk_intensity",
                        "planting_intensity"     = "planting_intensity",
                        "anti_browsing_intensity"= "anti_browsing_intensity",
                        "grndwrk_intensity"      = "grndwrk_intensity")
  ) %>%
  transmute(model, focal_term, p.value,
            p_label = case_when(p.value < 0.001 ~ "<0.001",
                                TRUE ~ formatC(p.value, format = "f", digits = 3)))





# management present vs no present
get_mng_emm <- function(model, model_name, focal_term) {
  
  # ── detect internal term names ──────────────────────────────────────────────
  model_terms <- names(model$coefficients)
  uses_pred   <- any(str_detect(model_terms, "_pred"))
  
  # ── map focal_term to actual term in model ───────────────────────────────── 
  actual_term <- if (uses_pred) {
    recode(focal_term,
           "planting_intensity"      = "planting_pred",
           "anti_browsing_intensity" = "browsing_pred",
           "grndwrk_intensity"       = "grndwrk_pred")
  } else {
    focal_term
  }
  
  # ── build at_list with actual term name ──────────────────────────────────── 
  at_list <- list(c(0, 1))
  names(at_list) <- actual_term
  
  # ── models without level term ────────────────────────────────────────────── 
  no_level <- model_name %in% c("beta", "adapt", "spruce")
  
  emm_out <- if (no_level) {
    emmeans::emmeans(model,
                     specs = stats::as.formula(paste0("~ ", actual_term)),
                     at    = at_list,
                     type  = "response")
  } else {
    emmeans::emmeans(model,
                     specs = stats::as.formula(paste0("~ ", actual_term)),
                     at    = at_list,
                     type  = "response")
  }
  
  emm_out %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    # ── standardize column names back to _intensity ──────────────────────────
    rename_with(~ recode(.x,
                         "planting_pred"  = "planting_intensity",
                         "browsing_pred"  = "anti_browsing_intensity",
                         "grndwrk_pred"   = "grndwrk_intensity")) %>%
    mutate(model = model_name, focal_term = focal_term)
}

# ── Run pipeline ──────────────────────────────────────────────────────────────
emm_mng <- map_dfr(
  names(fin.models.all),
  function(nm) {
    bind_rows(
      get_mng_emm(fin.models.all[[nm]], nm, "planting_intensity"),
      get_mng_emm(fin.models.all[[nm]], nm, "anti_browsing_intensity"),
      get_mng_emm(fin.models.all[[nm]], nm, "grndwrk_intensity")
    )
  }
)

# verify — should have planting_intensity column, not planting_pred
names(emm_mng)
emm_mng %>% count(model, focal_term)

# ── Step 2: build emm_mng2 ────────────────────────────────────────────────────
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


emm_mng2 %>%
  select(model, focal_term, intensity, response_plot, lower.CL, upper.CL) %>%
  filter(focal_term %in% c("planting_intensity", "anti_browsing_intensity")) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  arrange(model, focal_term, intensity) #%>%
#  print(n = Inf)
  


# Get estimates from the teh model: campre by % change from  ---------------

# ── Build the summary table ───────────────────────────────────────────────────

emm_table <- emm_mng2 %>%
  select(model, focal_term, term, response_lab,
         intensity, response_plot, lower.CL, upper.CL) %>%
  rename(response = response_plot) %>%
  # scale proportions to % for readable output
  mutate(
    response  = case_when(model %in% c("adapt", "spruce") ~ response * 100, TRUE ~ response),
    lower.CL  = case_when(model %in% c("adapt", "spruce") ~ lower.CL * 100, TRUE ~ lower.CL),
    upper.CL  = case_when(model %in% c("adapt", "spruce") ~ upper.CL * 100, TRUE ~ upper.CL),
    response  = round(response, 2),
    lower.CL  = round(lower.CL, 2),
    upper.CL  = round(upper.CL, 2)
  ) %>%
  pivot_wider(
    names_from  = intensity,
    values_from = c(response, lower.CL, upper.CL),
    names_sep   = "_"
  ) %>%
  mutate(
    estimate_pct = round(100 * (response_1 - response_0) / response_0, 1),
    lower_pct    = round(100 * (lower.CL_1 - upper.CL_0) / response_0, 1),
    upper_pct    = round(100 * (upper.CL_1 - lower.CL_0) / response_0, 1),
    CI_95        = paste0("[", lower_pct, "%, ", upper_pct, "%]"),
    change_pct   = paste0(ifelse(estimate_pct > 0, "+", ""), estimate_pct, "%")
  ) %>%
  left_join(pvals_mng, by = c("model", "focal_term")) %>%
  select(
    Management    = term,
    Response      = response_lab,
    `Intensity 0` = response_0,
    `Intensity 1` = response_1,
    `Change (%)`  = change_pct,
    `95% CI`      = CI_95,
    `p-value`     = p_label
  ) %>%
  arrange(Management, Response)

# ── Export to Word ────────────────────────────────────────────────────────────
sjPlot::tab_df(
  emm_table,
  title  = "Predicted marginal means at management intensity 0 vs 1",
  footnote = "Predicted marginal means from GAMs, all other variables held at their means. Proportions (climate-adapted share, Norway spruce share) scaled to %. 95% CI on the difference computed as lower.CL_1 − upper.CL_0 and upper.CL_1 − lower.CL_0.",
  file   = "outTable/management_effects_emmeans.doc"
)



# ── Step 3: build pct change table ───────────────────────────────────────────
model_intensity_all_df_pct_mng <- emm_mng2 %>%
  select(model, focal_term, term, response_lab,
         intensity, response_plot, lower.CL, upper.CL) %>%
  rename(response = response_plot) %>%
  pivot_wider(names_from  = intensity,
              values_from = c(response, lower.CL, upper.CL),
              names_sep   = "_") %>%
  left_join(pvals_mng, by = c("model", "focal_term")) %>%
  mutate(
    estimate_pct = 100 * (response_1 - response_0) / response_0,
    lower_pct    = 100 * (lower.CL_1  - upper.CL_0) / response_0,  # <-- lower of diff
    upper_pct    = 100 * (upper.CL_1  - lower.CL_0) / response_0,  # <-- upper of diff
    sig_col      = ifelse(p.value < 0.05, "sig", "n.s."),
    response     = factor(response_lab,
                          levels = c("Mean height [m]", "CV [%]",
                                     "Effective species [#]",
                                     "Species richness [#]",
                                     "Dissimilarity [dim.]",
                                     "Climate-adapted share [%] ",
                                     "Norway spruce share [%]"))
  )

# ── Verify ────────────────────────────────────────────────────────────────────
model_intensity_all_df_pct_mng %>%
  select(model, focal_term, estimate_pct, p_label , response) %>%
  filter(!is.na(estimate_pct)) %>%
  print(n = Inf)

model_intensity_all_df_pct_mng <- model_intensity_all_df_pct_mng %>%
  mutate(
    plot_only    = model %in% c("beta", "adapt", "spruce"),
    bar_color    = if_else(!plot_only, "grey20", NA_character_),  # border on cross-scale
    bar_linewidth = if_else(!plot_only, 0.8, 0)                   # no border on plot-only
  )

# ── Check no NAs ──────────────────────────────────────────────────────────────
model_intensity_all_df_pct_mng %>%
  filter(is.na(estimate_pct) | is.na(p.value) | is.na(response)) %>%
  select(model, focal_term, response_lab)



## 6e. Time since disturbance x scale (panel)
ylim_hgt  <- c(0.2, 2.0)
ylim_cv   <- c(10, 85)
ylim_eff  <- c(1, 12)
ylim_rich <- c(0.8, 7)

y_time_hgt  <- ylim_hgt[2]  * 0.98
y_time_cv   <- ylim_cv[2]   * 0.98
y_time_eff  <- ylim_eff[2]  * 0.98
y_time_rich <- ylim_rich[2] * 0.98

p_hgt <- pp(fin.models.all$hgt,
            terms   = c("time_snc_full_disturbance[0:7]", "level"),
            ylab    = "Mean height [m]",
            annot   = "Time:\np = 0.004",
            ylim    = ylim_hgt)

p_cv_pos <- pp(fin.models.all$cvpos,
               terms   = c("time_snc_full_disturbance[0:7]", "level"),
               ylab    = "CV [%]",
               scale_y = 100,
               annot   = "Time:\np = 0.763",
               ylim    = ylim_cv)

p_eff <- pp(fin.models.all$eff,
            terms   = c("time_snc_full_disturbance[0:7]", "level"),
            xlab    = "Time since disturbance\n(years)",
            ylab    = "Effective species [#]",
            annot   = "Time:\np = 0.571",
            ylim    = ylim_eff)

p_rich <- pp(fin.models.all$rich,
             terms   = c("time_snc_full_disturbance[0:7]", "level"),
             xlab    = "Time since disturbance\n(years)",
             ylab    = "Species richness [#]",
             annot   = "Time:\np = 0.840",
             ylim    = ylim_rich)

p_inset_hgt  <- pp_inset_model(fin.models.all$hgt,   scale_y = 1,
                               p_lab = "Scale:\np < 0.001", ylim = ylim_hgt)
p_inset_cv   <- pp_inset_model(fin.models.all$cvpos, scale_y = 100,
                               p_lab = "Scale:\np < 0.001", ylim = ylim_cv)
p_inset_eff  <- pp_inset_model(fin.models.all$eff,   scale_y = 1,
                               p_lab = "Scale:\np < 0.001", ylim = ylim_eff)
p_inset_rich <- pp_inset_model(fin.models.all$rich,  scale_y = 1,
                               p_lab = "Scale:\np < 0.001", ylim = ylim_rich)

pair_hgt  <- (p_hgt  + p_inset_hgt)  + plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[a]")
pair_cv   <- (p_cv_pos + p_inset_cv) + plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[b]")
pair_eff  <- (p_eff  + p_inset_eff)  + plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[c]")
pair_rich <- (p_rich + p_inset_rich) + plot_layout(widths = c(2.6, 1)) + plot_annotation(title = "[d]")

p_final <- (pair_hgt | pair_cv) / (pair_eff | pair_rich) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_final





## Management effect (% change 0 -> 1) ----------------------------------
p_model_response <- ggplot(
  model_intensity_all_df_pct_mng,
  aes(x = response, y = estimate_pct, fill = response)
) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  #geom_col(width = 0.7, color = NA) +
  geom_col(aes(color = plot_only, linewidth = plot_only),
           width = 0.7) +
  scale_color_manual(
    values = c("TRUE" = "grey20", "FALSE" = NA),
    name = "Measurement scale",
    labels = c("TRUE" = "Plot level", "FALSE" = "Cross-scale"),
    guide = guide_legend(override.aes = list(fill = "grey70"))
  ) +
  scale_linewidth_manual(
    values = c("TRUE" = 0.6, "FALSE" = 0),
    guide = "none"
  ) +
  geom_errorbar(aes(ymin = lower_pct, ymax = upper_pct),
                width = 0.15, linewidth = 0.5, color = "grey30") +
  # geom_text(aes(label    = p_label,
  #               y        = upper_pct + 12,
  #               fontface = ifelse(p.value < 0.05, "bold", "plain")),
  #           size = 2.5) +
  
  geom_text(aes(label    = p_label,
                y        = ifelse(estimate_pct >= 0, 
                                  upper_pct + 15, 
                                  lower_pct - 22),
                fontface = ifelse(p.value < 0.05, "bold", "plain"),
                vjust    = ifelse(estimate_pct >= 0, 0, 1)),
            size = 2.5)+
  
  facet_wrap(~ term, ncol = 3) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "", y = "Effect on response\n
       [% change relative to unmanaged baseline]") +
  theme_classic(base_size = 8) +
  theme(legend.position    = "none",
        strip.text         = element_text(size = 8),
        strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.6),
        axis.text.x        = element_text(angle = 30, hjust = 1),
        panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3))


p_model_response

# plot effect sizes as betas 

# START
library(tidyverse)
library(broom)

# ── 1. Extract coefficients with CI ──────────────────────────────────────────
focal_terms_map <- c(
  "planting_pred"          = "planting_intensity",
  "browsing_pred"          = "anti_browsing_intensity",
  "grndwrk_pred"           = "grndwrk_intensity",
  "planting_intensity"     = "planting_intensity",
  "anti_browsing_intensity"= "anti_browsing_intensity",
  "grndwrk_intensity"      = "grndwrk_intensity"
)

# Model families for transformation decision
log_link_models  <- c("hgt", "cvpos", "eff", "rich")   # exp(β) → IRR / multiplier
identity_models  <- c("beta")                            # β on response scale
logit_models     <- c("adapt", "spruce")                 # exp(β) → OR, or marginal Δp

coef_df <- map_dfr(names(fin.models.all), function(nm) {
  broom::tidy(fin.models.all[[nm]], parametric = TRUE, conf.int = TRUE) %>%
    filter(term %in% names(focal_terms_map)) %>%
    mutate(
      model      = nm,
      focal_term = recode(term, !!!focal_terms_map)
    )
})

# ── 2. Transform to interpretable effect size ─────────────────────────────────
coef_df2 <- coef_df %>%
  mutate(
    link = case_when(
      model %in% log_link_models  ~ "log",
      model %in% identity_models  ~ "identity",
      model %in% logit_models     ~ "logit"
    ),
    # For log-link: exponentiate → "response multiplier" (IRR/rate ratio)
    # For identity: keep as-is (Jaccard units)
    # For logit: exponentiate → odds ratio (or convert to Δp below)
    effect     = case_when(link == "identity" ~ estimate,
                           TRUE               ~ exp(estimate)),
    effect_lo  = case_when(link == "identity" ~ conf.low,
                           TRUE               ~ exp(conf.low)),
    effect_hi  = case_when(link == "identity" ~ conf.high,
                           TRUE               ~ exp(conf.high)),
    # Human-readable: for log-link, subtract 1 → % change from baseline
    # (now baseline-independent, because exp(β)-1 is pure multiplicative effect)
    effect_pct      = case_when(link == "log"      ~ (effect - 1) * 100,
                                link == "identity"  ~ effect * 100,   # Jaccard is 0-1, ×100 for readability
                                link == "logit"     ~ (effect - 1) * 100),  # OR expressed as % change in odds
    effect_pct_lo   = case_when(link == "log"      ~ (effect_lo - 1) * 100,
                                link == "identity"  ~ effect_lo * 100,
                                link == "logit"     ~ (effect_lo - 1) * 100),
    effect_pct_hi   = case_when(link == "log"      ~ (effect_hi - 1) * 100,
                                link == "identity"  ~ effect_hi * 100,
                                link == "logit"     ~ (effect_hi - 1) * 100),
    sig_col    = ifelse(p.value < 0.05, "sig", "n.s."),
    response   = factor(recode(model, !!!response_labels),
                        levels = c("Mean height [m]", "CV [%]",
                                   "Effective species [#]", "Species richness [#]",
                                   "Dissimilarity [dim.]",
                                   "Climate-adapted share [%] ",
                                   "Norway spruce share [%]")),
    term_lab   = recode(focal_term, !!!term_labels),
    p_label    = case_when(p.value < 0.001 ~ "<0.001",
                           TRUE ~ formatC(p.value, format = "f", digits = 3)),
    plot_only  = model %in% c("beta", "adapt", "spruce")
  )

# ── 3. Plot ───────────────────────────────────────────────────────────────────
p_coef <- ggplot(coef_df2,
                 aes(x = response, y = effect_pct, fill = response)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  geom_col(aes(color = plot_only, linewidth = plot_only), width = 0.7) +
  scale_color_manual(values = c("TRUE" = "grey20", "FALSE" = NA),
                     name   = "Measurement scale",
                     labels = c("TRUE" = "Plot level", "FALSE" = "Cross-scale"),
                     guide  = guide_legend(override.aes = list(fill = "grey70"))) +
  scale_linewidth_manual(values = c("TRUE" = 0.6, "FALSE" = 0), guide = "none") +
  geom_errorbar(aes(ymin = effect_pct_lo, ymax = effect_pct_hi),
                width = 0.15, linewidth = 0.5, color = "grey30") +
  geom_text(aes(label    = p_label,
                y        = ifelse(effect_pct >= 0, effect_pct_hi + 5, effect_pct_lo - 5),
                fontface = ifelse(p.value < 0.05, "bold", "plain"),
                vjust    = ifelse(effect_pct >= 0, 0, 1)),
            size = 2.5) +
  facet_wrap(~ term_lab, ncol = 3) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "", y = "Effect size [% change in response or % change in odds]") +
  theme_classic(base_size = 8) +
  theme(legend.position    = "none",
        strip.text         = element_text(size = 8),
        strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.6),
        axis.text.x        = element_text(angle = 30, hjust = 1),
        panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3))

p_coef

# END








# does browsing protection correlate with lower spruce density
# or just lower spruce SHARE (because other species increase)?
both_levels_cz %>%
  filter(level == "plot") %>%
  group_by(anti_browsing_intensity > 0) %>%
  summarise(
    mean_spruce_share  = mean(spruce_share, na.rm = TRUE),
    mean_spruce_dens   = mean(spruce_share * dens_ha, na.rm = TRUE),
    n = n()
  )




# Tables & export ---------------------------------------------


ggsave("outFigsCZ/p_function_drought.png",
       p_function_drought, width = 6, height = 4, dpi = 300)


## Save figures
ggsave("outFigsCZ/p_species_composition.png",
       p_species_composition, width = 7, height = 3, dpi = 300)

ggsave("outFigsCZ/p_time_scale.png",
       p_final, width = 6, height = 5, dpi = 300, bg = "white")

ggsave("outFigsCZ/p_time_scale.svg",
       p_final, width = 6, height = 5, dpi = 300, bg = "white")


ggsave("outFigsCZ/p_model_response.png",
       p_model_response, width = 7, height = 4, dpi = 300)

ggsave("outFigsCZ/p_management_intensity_plot_simpler.png",# - this one cuts
       p_management_intensity_plot_simpler, width = 5, height = 2.8, dpi = 300)

ggsave("outFigsCZ/p_combined_disturb_fig.png",
       p_combined_disturb_fig, width = 5, height = 2.5, dpi = 300)

ggsave("outFigsCZ/density_plot.png",
       p_height_seral_mng, width = 6, height = 4, dpi = 300)


ggsave("outFigsCZ/p_func_tsd_col_v2.png",
       p_bar_TSD, width = 6, height = 3, dpi = 300, bg = "white")

# alternative: svg is also fully editable in Corel
ggsave("outFigsCZ/p_combined_function.svg",
       p_combined_func,
       width  = 7,
       height = 3.5,
       device = "svg")

ggsave("outFigsCZ/p_combined_function.pdf",
       p_combined_func, 
       width  = 7, 
       height = 3.5, 
       device = cairo_pdf)

ggsave("outFigsCZ/p_combined_function.png",
       p_combined_func, width = 7, 
       height = 3.5, dpi = 300, bg = "white")








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

results <- map2(fin.models.all, 
                names(fin.models.all), extract_gam_summary)

pvals_fmt <- map_dfr(results, "pvalues") %>%
  select(model, term_clean, pval_fmt) %>%
  pivot_wider(names_from = term_clean, values_from = pval_fmt, names_prefix = "p_")

metrics_df    <- bind_rows(map(results, "metrics"))
final_results <- left_join(metrics_df, pvals_fmt, by = "model")

tab_df(final_results,
       title = "Model Summary Table",
       file  = "outTable/models_summary_final_results.doc")


# ## ── AIC table: cross-scale models -------------------------------
aic_cross_table <- aic_results %>%
  mutate(
    response = recode(response,
                      "mean_hgt"         = "Mean height",
                      "cv_hgt_pos"       = "CV height",
                      "effective_numbers" = "Effective species",
                      "sp_richness"      = "Species richness"),
    model = recode(model,
                   "interaction" = "Planting × Browsing",
                   "both"        = "Planting + Browsing",
                   "plant_only"  = "Planting only",
                   "browse_only" = "Browsing only",
                   "none"        = "Null"),
    AIC       = round(AIC, 1),
    delta_AIC = round(delta_AIC, 2),
    best      = ifelse(best, "✓", "")
  ) %>%
  select(Response = response,
         `Management structure` = model,
         AIC,
         `ΔAIC` = delta_AIC,
         `Best (ΔAIC<2)` = best)

tab_df(aic_cross_table,
       title   = "AIC comparison of management structures — cross-scale models (subplot + plot level)",
       file    = "outTable/aic_crossscale_models.doc")

### ── AIC table: plot-level adapted model ------─────
aic_adapted_table <- compare_mng_aic_plot(
  "share_adapted_adj", share_adapted_plot, betar(), k_tsd = 4
) %>%
  mutate(
    model = recode(model,
                   "interaction" = "Planting × Browsing",
                   "both"        = "Planting + Browsing",
                   "plant_only"  = "Planting only",
                   "browse_only" = "Browsing only",
                   "none"        = "Null"),
    AIC       = round(AIC, 1),
    delta_AIC = round(delta_AIC, 2),
    best      = ifelse(best, "✓", "")
  ) %>%
  mutate(Response = "Drought-tolerant species share") %>%
  select(Response,
         `Management structure` = model,
         AIC,
         `ΔAIC` = delta_AIC,
         `Best (ΔAIC<2)` = best)

tab_df(aic_adapted_table,
       title   = "AIC comparison of management structures — plot-level climate adaptation model",
       file    = "outTable/aic_adapted_model.doc")

# combine aic into one table
aic_combined <- bind_rows(aic_cross_table, aic_adapted_table)

tab_df(aic_combined,
       title   = "AIC comparison of management structures across all response variables",
       file    = "outTable/aic_all_models.doc")





