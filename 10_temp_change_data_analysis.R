
# Analyze data temporal change
#   

# read overlapping data: from 223 and 2025
# investigate how they change with time since disturbnace
# structure (height, cv, stem density), composition (shannon, richness, eveness)
# run models:


gc()

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




theme_set(theme_classic2(base_size = 10) +
            theme(axis.title = element_text(size = 10),
                  axis.text  = element_text(size = 10)))


# Read data -----------------------------
dat_overlap  <- fread('outData/full_table_overlap_23_25.csv')


# get master table, having all unique plots and subplots
dat_master_subplot <- dat_overlap %>% 
  dplyr::select(plot, subplot, year) %>% 
  distinct()

table(dat_overlap$year)  

n_plots_total <- length(unique(dat_master_subplot$plot))     # 126
n_plots_total

n_subplots_total <-length(unique(dat_master_subplot$subplot))  # 1250
n_subplots_total

# histogram of stem denisty per vertcal class:
dat_overlap %>% 
  dplyr::filter(n > 0) %>% 
  ggplot(aes(n , fill = year)) +
  geom_histogram() + 
  facet_grid(year~vegtype, scales = 'free')


# Clean up dat on subplot & plot level ---------------------------------
## Get disturbance characteristics (plot) --------------
plot_disturb_chars <- dat_overlap %>% 
  dplyr::select(plot, year, 
                disturbance_year, 
                forest_year, 
                disturbance_length, 
                time_snc_full_disturbance, 
                time_snc_part_disturbance #,
                #clear, grndwrk, logging_trail, planting, anti_browsing
  )  %>%  
  distinct()

plot_disturb_chars %>% 
  ggplot(aes(time_snc_full_disturbance)) + 
  geom_histogram() + 
  facet_grid(.~year)


p_hist_dist_length <- plot_disturb_chars %>% 
 # filter(year == 2023) %>%  # to keep half of values
  ggplot(aes(disturbance_length)) + 
  geom_histogram(fill = 'grey90', color = 'black') 

p_hist_dist_year <- plot_disturb_chars %>%
  #filter(year == 2023) %>%  # to keep half of values
  filter(disturbance_year > 2012) %>% 
  ggplot(aes(disturbance_year)) + 
  geom_histogram(fill = 'grey90', color = 'black') 


p_hist_time_since_dist <- plot_disturb_chars %>%
 # filter(year == 2023) %>%  # to keep half of values
  filter(disturbance_year > 2012) %>% 
  ggplot(aes(time_snc_full_disturbance)) + 
  geom_histogram(fill = 'grey90', color = 'black') 

ggarrange(p_hist_dist_year, p_hist_dist_length, p_hist_time_since_dist,
          ncol = 3)


# Clean up management - if site has been mamaged in 2023, it needs to be managemt (cleared)
# in 2025 as well!!
df_master_mng <- dat_overlap %>% 
  dplyr::filter(year == "2023") %>% # keep management oionly frm 2023 for consistency
  distinct(plot, subplot,year,
           clear,
           grndwrk,
           logging_trail,
           planting,
           anti_browsing)

prop.table(table(df_master_mng$clear ))
prop.table(table(df_master_mng$grndwrk))
prop.table(table(df_master_mng$logging_trail))
prop.table(table(df_master_mng$planting))
prop.table(table(df_master_mng$anti_browsing))

# > table(df_master_mng$clear )
# 0   1 
# 7 618 
# > table(df_master_mng$grndwrk)
# 0   1 
# 79 546 

# > table(df_master_mng$logging_trail)
# 0   1 
# 522 103 

# > table(df_master_mng$planting)
# 0   1 
# 297 328 

# > table(df_master_mng$anti_browsing)
# 0   1 
# 377 248

# get management characteristics only from 2023
# clear,
# grndwrk,
# logging_trail,
# planting,
# anti_browsing


## Early vs late (landscape, plot, subplot) ---------------
# Summarize total number of trees per year and recovery type
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

print(tree_summary)

# Stacked barplot
tree_summary %>% 
  dplyr::filter(seral_stage != 'other') %>% 
  ggplot( aes(x = year, y = share, fill = seral_stage)) +
  geom_bar(stat = "identity", color = "black") +
  labs(
    title = "Tree Counts by Recovery Type per Year",
    x = "Year",
    y = "Share [%]",
    fill = "Recovery Type"
  ) +
  scale_fill_manual(values = c("early" = "#66c2a5", "late" = "#fc8d62")) +
  theme_classic2(base_size = 8)




dat_sum_recovery <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>% 
  # mutate(year = as.factor(year)) %>% 
  group_by(seral_stage, year) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE))

dat_sum_recovery


### Early vs. late species over time --------------------

dat_overlap %>% 
  group_by(plot, seral_stage) %>%
  summarise(n_stems = n(), .groups = "drop") %>%
  pivot_wider(names_from = seral_stage,
              values_from = n_stems,
              values_fill = 0)  # fill missing types with 0


# check up development of early vs ;late shares givet time since disturbance

# Calculate stem counts by recovery type at the plot level
share_early_vs_late <- 
  dat_overlap %>%
  filter(!is.na(n)) %>%  # Optional: remove NAs if present
  group_by(plot, seral_stage, time_snc_full_disturbance) %>%
  summarise(n_stems = n(), .groups = "drop") %>%
  pivot_wider(names_from = seral_stage,
              values_from = n_stems,
              values_fill = 0) %>% 
  mutate(total = early + late,
         share_early = early/total*100,
         share_late = late/total*100) %>% 
  select(plot, time_snc_full_disturbance, share_early, share_late) %>%
  pivot_longer(cols = starts_with("share_"),
               names_to = "seral_stage",
               values_to = "share")# %>%

### Early vs late: plot level 
df_plot_share_early <-  
  dat_overlap %>%
  filter(!is.na(n)) %>%  # Optional: remove NAs if present
  group_by(plot,year, seral_stage, time_snc_full_disturbance) %>%
  summarise(n_stems = n(), .groups = "drop") %>%
  pivot_wider(names_from = seral_stage,
              values_from = n_stems,
              values_fill = 0) %>% 
  mutate(total = early + late,
         share_early = early/total*100,
         share_late = late/total*100) %>% 
  select(plot, year, share_early, share_late,time_snc_full_disturbance) 

### early vs late : subplot level 
df_sub_share_early <-  
  dat_overlap %>%
  filter(!is.na(n)) %>%  # Optional: remove NAs if present
  group_by(subplot, plot,year, seral_stage, time_snc_full_disturbance) %>%
  summarise(n_stems = n(), .groups = "drop") %>%
  pivot_wider(names_from = seral_stage,
              values_from = n_stems,
              values_fill = 0) %>% 
  mutate(total = early + late,
         share_early = early/total*100,
         share_late = late/total*100) %>% 
  select(subplot, plot, year, share_early, share_late,time_snc_full_disturbance) 


# Plot
ggplot(share_early_vs_late, 
       aes(x = factor(time_snc_full_disturbance),
           y = share,
           fill = seral_stage)) +
  geom_boxplot() +
  #geom_jitter() +
  #geom_bar(stat = "identity", position = "stack") +
  labs(x = "Time since stand replacing\ndisturbance (years)",
       y = "Share of stems (%)",
       fill = "Seral stage") 




## Traits: Weighted community mean (subplot & plot) ------------------

# choose  weighting : by stems (for early communities), by structure?
# Option A (default): weight by stem counts
wvar <- "n"

# Option B: weight by structure (uncomment ONE)
# wvar <- "basal_area_cm2"
# wvar <- "hgt_est"

# helper to pull a numeric weight safely
wfun <- function(x) ifelse(is.na(x) | x < 0, 0, x)

# Subplot × year CWMs 
cwm_subplot <- dat_overlap %>%
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
  )

# 2) Plot × year CWMs
cwm_plot <- dat_overlap %>%
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
  )

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
cwm_all_long %>% 
  ggplot(aes(x = time_snc_full_disturbance, 
             y = CWM, 
             fill = factor(time_snc_full_disturbance))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
 # geom_jitter(width = 0.1, alpha = 0.25, size = 0.7) +
  facet_grid(level ~ trait, 
             labeller = labeller(trait = trait_labs),
             scales = 'free_y') +
  labs(x = "Year", y = "Community-weighted mean (CWM)",
       title = "Trait CWMs",
       subtitle = "by level and time since disturbance") +
  theme_bw(base_size = 10) + 
  theme(legend.position = "none")


## Field data summary: subplot metrics ----------------------------------------------------------
field_sub_summ <- dat_overlap %>%
  mutate(n = coalesce(n, 0L)) %>% # change NAs to 0 (0L = literrary 0 = integer (not double))
  # mutate(year = as.factor(year)) %>% 
  group_by(plot, subplot, year, 
           time_snc_full_disturbance, 
           time_snc_part_disturbance,
           disturbance_year, 
           forest_year, disturbance_length,
           clear,
           grndwrk,
           logging_trail,
           planting,
           anti_browsing,
           management_intensity) %>%
  summarise(
    stems_total = sum(n, na.rm = TRUE),
    sp_richness = if (stems_total > 0) n_distinct(species[n > 0]) else 0,
    shannon_sp  = if (stems_total > 0) vegan::diversity(n, index = "shannon") else 0,
    evenness_sp = if (sp_richness > 1) shannon_sp / log(sp_richness) else 0,
    effective_numbers = ifelse(stems_total > 0, exp(shannon_sp), NA_real_),
    # weighted mean height using only present stems
    mean_hgt    = if (stems_total > 0) weighted.mean(hgt_est[n > 0], n[n > 0]) else NA_real_,
    
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

# mutate(cv_hgt = ifelse(is.na(cv_hgt), 0L, cv_hgt),
#       mean_hgt = ifelse(is.na(mean_hgt), 0L, mean_hgt)) # replace NA by 0 if stems are missing


# is teh change in community shading/drought tolerance driven by planting????
# Plot: Shade ~ Time since disturbance, by planting
x_lab_time_snc_full_dist = "Time since stand\nreplacing disturbance (years)"

#### Effect of planting? ---------------------------------------------
p_shade_planting <- field_sub_summ %>% 
  ggplot(aes(x = as.factor(time_snc_full_disturbance),
             y = CWM_shade,
             fill = factor(planting))) +
 #
  geom_boxplot() +
  labs(x = x_lab_time_snc_full_dist,
       y = "CWM shade",
       fill = "Planting") +
  theme_classic2() +
  theme(text  = element_text(size = 10))

# Plot: Drought ~ Time since disturbance, by planting
p_shade_drought <- field_sub_summ %>% 
  ggplot(aes(x = as.factor(time_snc_full_disturbance),
             y = CWM_drought,
             fill = factor(planting))) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = x_lab_time_snc_full_dist,
       y = "CWM drought",
       fill = "Planting") +
  theme_classic2()+
  theme(text  = element_text(size = 10))

# Plot: Shade ~ Planting
p_shade_total <- field_sub_summ %>% 
  ggplot(aes(x = factor(planting),
             y = CWM_shade,
             fill = factor(planting))) +
  #
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.y = 4, size = 3) +
  labs(x = "Planting",
       y = "CWM shade",
       fill = "Planting") +
  theme_classic2() +
  theme(text  = element_text(size = 10))

# Plot: Drought ~ Planting
p_drought_total <- field_sub_summ %>% 
  ggplot(aes(x = factor(planting),
             y = CWM_drought,
             fill = factor(planting))) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.y = 4, size = 3) +
  labs(x = "Planting",
       y = "CWM drought",
       fill = "Planting") +
  theme_classic2() +
  theme(text  = element_text(size = 10))

# Arrange all plots
annotate_figure(
  ggarrange(p_shade_planting, p_shade_total,
            p_shade_drought, p_drought_total,
            ncol = 2, nrow = 2,
            common.legend = TRUE, legend = "bottom"),
  top = text_grob("Subplot level", 
                  face = "bold", size = 12)
)



#### Subplot quick plotting: all vars ---------------------
df_sub_long <- field_sub_summ %>%
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
         management_intensity,
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
                  anti_browsing,
                  management_intensity),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value))# %>%



### Subplot: Time since disturbnace ----------------------------------------

# see CV with time since disturbnace : poartial disturbance
p_partial_disturbance <- df_sub_long %>% 
  ggplot(aes(x = time_snc_part_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_part_disturbance ), outlier.shape = NA) +
  
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Subplot: Partial disturbance') +
  theme_classic2()


# see CV with time since disturbnace : full disturbance
p_full_disturbance <- df_sub_long %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Subplot: Full disturbance') +
  theme_classic2()

ggarrange(p_partial_disturbance, p_full_disturbance, ncol = 2)



subplot_summary_tbl <- field_sub_summ %>%
  ungroup() %>%
  select(year, mean_hgt, cv_hgt, shannon_sp, sp_richness) %>%
  pivot_longer(-year, names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  group_by(metric, year) %>%
  summarise(
    n_total  = n(),
    n_non_na = sum(!is.na(value)),
    n_na     = sum(is.na(value)),
    mean     = mean(value, na.rm = TRUE),
    sd       = sd(value, na.rm = TRUE),
    se       = sd/sqrt(n_non_na),
    median   = median(value, na.rm = TRUE),
    p25      = quantile(value, 0.25, na.rm = TRUE),
    p75      = quantile(value, 0.75, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(metric, year)

subplot_summary_tbl


# summarize subplot information 
# how many trees?
sum(field_sub_summ$stems_total)            # 3658 stems
length(unique(field_sub_summ$subplot))     # 1250
sum(field_sub_summ$stems_total == 0)       # 271 
sum(field_sub_summ$cv_hgt > 1, na.rm = T)  # 20


field_sub_summ_filt <- field_sub_summ %>%
  filter(stems_total > 0 & cv_hgt > 0)

p1 <- ggplot(field_sub_summ_filt, aes(x = stems_total, y = cv_hgt, color = year,
                                      fill = year)) +
  geom_point(alpha = 0.4, size = 2 ) +
  # geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3),
              se = TRUE,  linetype = "dashed") +
  labs(
    x = "Stem count per subplot",
    y = "CV of tree height",
    title = "Stem density vs. Vertical Structural Variation",
    subtitle = "Empty subplots excluded"
  ) +
  theme_classic2(base_size = 8)


p2 <- ggplot(field_sub_summ_filt, aes(x = shannon_sp, y = cv_hgt, color = year,
                                      fill = year)) +
  geom_point(alpha = 0.4, size = 2) +
  # geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3),
              se = TRUE, linetype = "dashed") +
  labs(
    x = "Shannon per subplot",
    y = "CV of tree height",
    title = "Shannon vs. Vertical Structural Variation",
    subtitle = "Empty subplots excluded"
  ) +
  theme_classic2(base_size = 8)

ggarrange(p1, p2)
# ssame analysis on both scales? 

# Species composition --------------------------------------------------


## Get summary across all trees and study sites --------------------------
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

by_year_totals <- df %>%
  group_by(year) %>%
  summarise(total_trees = sum(n),
            #n_plots = n_distinct(plot),
            # n_subplots = n_distinct(subplot),
            .groups = "drop")

n_trees23 <- by_year_totals %>% 
  filter(year == 2023) %>% 
  pull()

n_trees25 <- by_year_totals %>% 
  filter(year == 2025) %>% 
  pull()

# 2) Species × year counts and presence
species_stem_share_year <- 
  df %>%
  group_by(year, species) %>%
  summarise(
    stems = sum(n),                                   # number of trees
    #plots_present = n_distinct(plot[n > 0]),          # plots where species occurs
    .groups = "drop"
  ) %>% 
  tidyr::pivot_wider(
    names_from  = year,
    values_from = c(stems),
    names_glue  = "{.value}_{year}",
    values_fill = list(stems = 0L)
  ) %>%
  dplyr::arrange(species) %>% 
  mutate(trees23 = n_trees23,
         trees25 = n_trees25,
         share23 = round(stems_2023/trees23*100,2),
         share25 = round(stems_2025/trees25*100,2),
         total_stems = stems_2023 + stems_2025,
         total_trees = n_trees23 + n_trees25,
         total_share = round(total_stems/total_trees * 100,2)) 

# merge counst across 2 years:     
top_overall_stem_share <- 
  species_stem_share_year %>%
  dplyr::select(species, total_stems, total_share) %>%
  dplyr::slice_max(order_by = total_share, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup() 

# get top 10 per each year (the order changes a bit)
top_stems_by_year <- species_stem_share_year %>%
  dplyr::select(species, share23, share25) %>%
  tidyr::pivot_longer(
    dplyr::starts_with("share"),
    names_to = "year",
    names_prefix = "share",
    values_to = "share"
  ) %>%
  dplyr::mutate(year = factor(paste0("20", year), levels = c("2023", "2025"))) %>%
  dplyr::group_by(year) %>%
  dplyr::slice_max(order_by = share, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

top_2023 <- top_stems_by_year %>% dplyr::filter(year == "2023") %>% dplyr::pull(species)
top_2025 <- top_stems_by_year %>% dplyr::filter(year == "2025") %>% dplyr::pull(species)

length(union(top_2023, top_2025))      # 11 (your result)
length(intersect(top_2023, top_2025))  # 9  (overlap size)

setdiff(top_2025, top_2023)  # species only in 2025's top10
setdiff(top_2023, top_2025)  # species only in 2023's top10

unique(top_stems_by_year$species)

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
# "#006837" "#17934D" "#58B65F" "#94D168" "#C6E77F" "#EDF7A7" "#FEF0A7" "#FDCD7B" "#FA9C58" "#EE613D" "#D22B26" "#A50026" 


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
  potr = "Populus tremula",
  absp = "Abies sp.",
  sasp = "Salix sp."
)

### Overall stats --------------------------------------------------------
#### Identify plots without any stems present ---------------------------------------------------------------------------------

# Step 1: Summarise total stems per plot and year
plot_year_summ <- df %>%
  group_by(plot, year) %>%
  summarise(total_stems = sum(n, na.rm = TRUE), .groups = "drop")

# Step 2: Reshape to wide format
plot_year_wide <- plot_year_summ %>%
  tidyr::pivot_wider(names_from = year, values_from = total_stems, values_fill = 0)

# Step 3: Filter based on presence in each year
plots_empty23    <- plot_year_wide %>% filter(`2023` == 0) %>% pull(plot) # 6 ~ 4.7%
plots_empty25    <- plot_year_wide %>% filter(`2025` == 0) %>% pull(plot) # 2 ~ 1.5%
plots_empty_both <- plot_year_wide %>% filter(`2023` == 0 & `2025` == 0)  # zero

plots_empty23
plots_empty25



### get average stem density per species per top 10 species --------------------------------
df_stem_dens_species <- df %>% 
  group_by(plot, species, year, n_subplots ) %>%
  summarize(sum_n = sum(n, na.rm =T)) %>% 
  mutate(scaling_factor = 10000/(n_subplots * 4),
         stem_dens = sum_n*scaling_factor) %>% 
  mutate(log_sum_stem_density = log10(stem_dens + 1)) #%>%  # Adding 1 to avoid log(0)
#ungroup()

# get total sum and calculate as average value over all sites 
df_stem_dens_species_sum <- 
  df_stem_dens_species %>% 
  group_by(species, year) %>% 
  summarise(stem_dens = sum(stem_dens, na.rm = T),
            log_sum_stem_density = sum(log_sum_stem_density, na.rm = T)) %>%
  mutate(stem_dens_avg = stem_dens/n_plots_total,
         log_sum_stem_density_avg = log_sum_stem_density/n_plots_total)


df_stem_dens_species_year <- df_stem_dens_species %>% 
  ungroup(.) %>% 
  filter(sum_n >0) %>% 
  filter(species %in% v_top_species_overall) %>% 
  dplyr::group_by(species, year) %>%
  dplyr::mutate(median_stem_density = median(stem_dens, na.rm = TRUE)) %>% 
  dplyr::ungroup(.) %>%
  mutate(species = factor(species, levels = rev(v_top_species_overall))) # Set custom order



### Boxplot for stem density -------------
df_stem_dens_species_year2 <- df_stem_dens_species_year %>%
  dplyr::filter(!is.na(log_sum_stem_density) & sum_n > 0) %>%
  dplyr::mutate(year = factor(year, levels = c("2023","2025")),
                # order by mean log density (ascending → highest ends up at the TOP after coord_flip)
                species = forcats::fct_reorder(species, log_sum_stem_density, .fun = mean, na.rm = TRUE))

# 
p_density<-df_stem_dens_species_year2 %>% 
  filter(!is.na(species)) %>% 
  ggplot(aes(x = log_sum_stem_density, y = species,
             fill = year)) +
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
  # scale_fill_manual(values= species_colors) +
  # scale_alpha_manual(
  #   values = c("2023" = 0.5, "2025" = 1)#,  # 2023 = lighter
  #   #guide = "none"  # hides alpha legend
  # ) +
  theme_classic(base_size = 10) +
  scale_y_discrete(labels = species_labels) +
  theme(axis.text.y = element_text(face = "italic", size = 8))


p_density


# barplot of change over years

# Order species by their maximum share across years (nice stable ordering)
order_levels <- top_stems_by_year %>%
  group_by(species) %>%
  summarise(max_share = max(share, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_share)) %>%
  pull(species)

top_stems_by_year <- top_stems_by_year %>%
  mutate(
    species = factor(species, levels = order_levels),
    year    = factor(year, levels = c("2023","2025"))
  )

# master mapping: code -> Latin name
species_labels_all <- c(
  piab = "Picea abies",
  besp = "Betula sp.",
  pisy = "Pinus sylvestris",
  qusp = "Quercus sp.",
  fasy = "Fagus sylvatica",
  saca = "Salix caprea",
  lade = "Larix decidua",
  soau = "Sorbus aucuparia",
  acps = "Acer pseudoplatanus",
  potr = "Populus tremula",
  absp = "Abies sp.",
  sasp = "Salix sp.",
  cabe = "Carpinus betulus"
)

p_bar <- ggplot(top_stems_by_year, aes(x = share, y = species, fill = year)) +
  geom_col(aes(#group = interaction(species, year)#,
    #alpha = factor(year)
  ), 
  position = position_dodge(width = 0.7), width = 0.6) +
  # if 'share' is 0–100 already, just add a % suffix:
  scale_x_continuous(labels = label_number(accuracy = 0.1, suffix = "")) +
  scale_y_discrete(
    limits = rev(names(species_labels)),
    labels = species_labels,
    drop = FALSE
  ) +
  # scale_alpha_manual(
  #   values = c("2023" = 0.5, "2025" = 1.0)#,  # 2023 = lighter
  #   #guide = "none"  # hides alpha legend
  # ) +
  # if you prefer proportions (0–1), use:
  # scale_x_continuous(labels = label_percent()) 
  labs(
    x = "Share of stems",
    y = "Species",
    fill = "Year"#,
    #title = "Top 12 species by share, by year"
  ) +
  # scale_fill_manual(values = species_colors) +
  theme_classic2(base_size = 10) +
  theme(axis.text.y = element_text(size = 8, face = "italic"))

p_bar

# Get species occurence from total number of plots 
# Total number of unique plots
total_plots <- df_stem_dens_species %>%
  pull(plot) %>%
  n_distinct()

# Share of plots per species (where species has non-zero stems)
species_occurence <- 
  df_stem_dens_species %>%
  ungroup(.) %>% 
  dplyr::filter(sum_n > 0) %>%                 # Only where species occurred
  distinct(year, species, plot) %>%           # Unique species × plot combos
  count(year, species, name = "n_plots") %>%  # Count number of plots per species
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
  mutate(species = factor(species, levels = rev(species_order)))

# Plot
p_occurence <- species_plot_share %>% 
  filter(species %in% v_top_species_overall ) %>% 
  ggplot(aes(x = share_of_plots, y = species, fill = year)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Species occurence over plots (%)",
    y = "Species",
    fill = "Year"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "right"
  )

p_occurence



ggarrange(p_bar, p_occurence, p_density,  
          ncol = 3, common.legend = T)

### heights by species -------------------------------------------------------
# get heights by species
library(forcats)
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


# Drivers -----------------------------------------------------------------

## Summary stats on plot level stem density per m² --------------------
area_subplot_m2 <- 4      # 4 m²
area_plot_m2    <- 5*4    # 20 m²



## --- Plot-level metrics (aggregate over subplots) ---
plot_metrics_mean <- field_sub_summ %>%
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
  left_join(cwm_plot)


# --- pooled CV directly from dat23_subplot_recode ---
plot_metrics_pooled  <- dat_overlap %>%
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
  left_join(cwm_plot)
#mutate(cv_hgt = ifelse(is.na(cv_hgt), 0L, cv_hgt)) # replace NA by 0 if stems are missing

# get only context and disturbance information on plot level
df_plot_context <- dat_overlap %>%
  dplyr::select(plot,year,  
                pre_dist_trees_n,  area_m2, 
                pre_dist_dens_ha, 
                time_snc_full_disturbance,
                time_snc_part_disturbance,
                disturbance_year, forest_year, disturbance_length,
                protection_intensity, management_intensity) %>% 
  distinct() 


## create final table for both levels -----------------------------------------------------
# 1) Subplot table (has subplot mean_hgt and stems_total as weights)
sub_df <- field_sub_summ %>%
  filter(stems_total > 0) %>% #, cv_hgt > 0
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
                                "time_snc_full_disturbance" = "time_snc_full_disturbance") )

nrow(sub_df)
hist(sub_df$cv_hgt, breaks = 80)

# i have some many NA in cv_hgt?
# check
dat_overlap %>% 
  filter(subplot == "514_T4_TP_20250827") %>% 
  arrange(n)



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
                             "time_snc_full_disturbance" = "time_snc_full_disturbance") )


#  3) Bind & clean final table with both levels
both_levels_re2 <- bind_rows(sub_df, plot_df) %>%
  left_join(df_plot_share_early, 
            by = c('plot_id' = 'plot',
                   "year" = "year",
                   "time_snc_full_disturbance" = "time_snc_full_disturbance")) %>% 
  mutate(early_class = case_when(
    share_early >= 60 ~ "early_dom",
    share_early <= 20 ~ "late_dom",
    TRUE              ~ "mixed"
  )) |> 
  mutate(trait_class = factor(early_class, 
                              levels = c("late_dom",
                                         "mixed",
                                         "early_dom"))) %>% 
  mutate(
    level   = factor(level, levels = c("subplot","plot")),
    plot_id = factor(plot_id),
    w       = pmin(pmax(w, 1), 50)   # cap weights so a few dense plots don't dominate
  ) %>% 
  mutate(time_since_f = ifelse(time_snc_full_disturbance <= 2, "early", "later") %>% 
           factor(levels = c("early", "later")))

table(both_levels_re2$trait_class)

hist(both_levels_re2$stems_total, breaks = 150)
range(both_levels_re2$stems_total, na.rm = T)



# add small value to CV = 0
both_levels_re2 <- both_levels_re2 %>%
  mutate(
    cv_hgt_pos = ifelse(cv_hgt <= 0, 1e-4, cv_hgt)#is.na(cv_hgt) | 
  ) %>% 
  filter(!is.na(mean_hgt), !is.na(w), !is.na(cv_hgt), !is.na(dens_m2))



# check the distribution oof weights
ggplot(both_levels_re2, aes(x = w)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  labs(x = "Weight (capped stem count)", y = "Frequency") +
  theme_minimal()




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

# create species data only for Cv_heights testing - becaue i have many zeros 
# if i have only one tree species present!
dat_cv_hgt <- both_levels_re2 %>% 
  filter(cv_hgt>0) %>% 
  filter(time_snc_full_disturbance>0)

##### CV_hgt vs disturbance length ----------------
dat_cv_hgt %>%
  ggplot(aes(x = disturbance_length,
             y = cv_hgt_pos)) +
  geom_point() +
  geom_smooth(method = 'loess')




dat_cv_hgt %>%
  filter(time_snc_full_disturbance>0) %>% 
   ggplot(aes(x = time_snc_full_disturbance,
             y = cv_hgt_pos)) +
  geom_point() +
  geom_smooth(method = 'loess')


##### CV x Disturbnace length  

# cv_hgt vs disturbance length 
# Fit GAMM: smooth term + random effect of plot_id
m0 <- gam(cv_hgt_pos ~ disturbance_length +
            #s(disturbance_length, by = level, k = 3) + 
            #level + 
            s(plot_id, bs = "re"),
          #family = tw(),
          #family = Gamma(link = "log"),
          data = both_levels_re2 %>% filter(stems_total > 0),
          method = "REML")

summary(m0)

# add time since disturbance as a fector
m1 <- gam(cv_hgt_pos ~ s(disturbance_length, by = level, k = 3) + 
            level + time_since_f +
            s(plot_id, bs = "re"),
          #family = tw(),
          #family = Gamma(link = "log"),
          data = both_levels_re2 %>% filter(stems_total > 0),
          method = "REML")

summary(m1)

# add interaction between level and time since disturbance
m2 <- gam(cv_hgt_pos ~ s(disturbance_length, by = level, k = 3) + 
            level*time_since_f +
            s(plot_id, bs = "re"),
          #family = tw(),
          #family = Gamma(link = "log"),
          data = both_levels_re2 %>% filter(stems_total > 0),
          method = "REML")

summary(m2)

plot(m2, page = 1)



# remove 

# nicely plot the predicted values 
p <- predict_response(m_cv_cont, terms = c("time_snc_full_disturbance [all]"))

plot(p, one_plot = TRUE)

##### CV x time since disturbnace ------------------------------------

# add interaction between level and time since disturbance
m_cv_hgt1 <- gam(cv_hgt_pos ~ 
                   s(time_snc_full_disturbance , k = 3) +
                   s(plot_id, bs = "re"),
                 
                 #family = tw(),
                 family = Gamma(link = "log"),
                 data = dat_cv_hgt, #%>% filter(stems_total > 0),
                 method = "REML")

dat_cv_hgt %>% 
  filter(stems_total > 0) %>% 
  filter(cv_hgt > 0) %>% 
  ggplot(aes(x = time_snc_part_disturbance,
             y = cv_hgt_pos,
             color = level)) +
  geom_point() + 
  geom_smooth(method = 'lm') + facet_grid(.~level)


summary(m_cv_hgt1)
appraise(m_cv_hgt1)

m_cv_hgt1 <- gam(cv_hgt_pos ~ 
                     s(time_snc_full_disturbance , k = 3) +
                     s(plot_id, bs = "re"),
                 family = tw(link = "log"),
                   #family = Gamma(link = "log"),
                 weights = w,
                   data = dat_cv_hgt, #%>% filter(stems_total > 0 &cv_hgt >0),
                   method = "REML")

summary(m_cv_hgt1)
appraise(m_cv_hgt1)
hist(both_levels_re2$cv_hgt_pos)
plot(m_cv_hgt1, page = 1)

m_cv_hgt_gamma <- gam(cv_hgt_pos ~ 
                     s(time_snc_full_disturbance , by = level, k = 3) + 
                     level +
                     s(plot_id, bs = "re"),
                   #family = tw(),
                   #family = tw(link = "log"),
                   #family = Gamma(link = "log"),
                   data = dat_cv_hgt,
                   weights = w,
                   method = "REML")

summary(m_cv_hgt_gamma)
appraise(m_cv_hgt_gamma)

AIC(m_cv_hgt_gamma, m_cv_hgt_tw)

m_cv_hgt_tw <- gam(cv_hgt_pos ~ 
                     s(time_snc_full_disturbance , by = level, k = 3) + 
            level +
            s(plot_id, bs = "re"),
          #family = tw(),
          family = tw(link = "log"),
          #family = Gamma(link = "log"),
          data = dat_cv_hgt,
          weights = w,
          method = "REML")

summary(m_cv_hgt_tw)
appraise(m_cv_hgt_tw)

plot(m_cv_hgt_tw, page = 1)

p <- predict_response(
  m_cv_hgt_tw,
  terms = c("time_snc_full_disturbance [all]",    # x-sequence
            "level [all]"
            #"early_class")                    # fix level (or drop to average)
  ))
p_cv_hgt <- plot(p, one_plot = TRUE)
#p_mean_hgt <- plot(p, one_plot = TRUE)


# remove 

# nicely plot the predicted values 
p <- predict_response(m_cv_cont, terms = c("time_snc_full_disturbance [all]"))

plot(p, one_plot = TRUE)




# get test scater plots ---------------------------------
p_early_cls_plot <- both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance),
             y = cv_hgt,
             color = early_class,
             fill = early_class,
             group = early_class)) +
  geom_point(alpha = 0.1) +
  ggtitle('plot') +
  geom_smooth(method = 'lm') +
  coord_cartesian(y = c(0,1))
#facet_grid(.~early_class)

p_early_cls_subplot <- both_levels_re2 %>% 
  filter(level == 'subplot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance),
             y = cv_hgt,
             color = early_class,
             fill = early_class,
             group = early_class)) +
  geom_point(alpha = 0.1) +
  ggtitle('subplot') +
  geom_smooth(method = 'lm') +
  coord_cartesian(y = c(0,1))
#facet_grid(.~early_class)
ggarrange(p_early_cls_plot, p_early_cls_subplot, common.legend = T)

# test effect of early seral dominance ---------------------
m_cv_cont <- gam(
  cv_hgt ~ s(share_early, k=5) + s(time_snc_full_disturbance, k=5) +
    te(share_early, time_snc_full_disturbance, k=c(4,4)) +
    level + s(plot_id, bs="re"),
  data = both_levels_re2, 
  family = gaussian()
)

summary(m_cv_cont)
appraise(m_cv_cont)
plot.gam(m_cv_cont, page = 1)


# nicely plot the predicted values 
p <- predict_response(m_cv_cont, terms = c("time_snc_full_disturbance [all]"))

plot(p, one_plot = TRUE)


summary(m_cv_cont)
draw(m_cv_cont)


# share classes --------------------
m_cv <- gam(
  cv_hgt ~ early_class + level +
    s(time_snc_full_disturbance, by = early_class, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  data = both_levels_re2, family = gaussian()
)
summary(m_cv)


# nicely plot the predicted values 
p <- predict_response(
  m_cv,
  terms = c("time_snc_full_disturbance [all]",    # x-sequence
            "early_class",                        # 3 groups
            "level [subplot]")                    # fix level (or drop to average)
)
plot(p, one_plot = TRUE)


boxplot(cv_hgt ~ early_class, 
        data = both_levels_re2, 
        col = "grey90")
boxplot(mean_hgt ~ early_class, 
        data = both_levels_re2, 
        col = "grey90")

###### mean height: time since disturbnace  ---------------------
library(mgcv)
library(gratia)
m_hgt_tw <- gam(
  mean_hgt ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2
)

p <- predict_response(
  m_hgt_tw,
  terms = c("time_snc_full_disturbance [all]",    # x-sequence
            "level [all]"
            #"early_class")                    # fix level (or drop to average)
))
p_mean_hgt <- plot(p, one_plot = TRUE)
ggarrange(p_mean_hgt, p_cv_hgt, common.legend = TRUE, 
          labels = "auto")


summary(m_hgt_tw)
appraise(m_hgt_tw)
draw(m_hgt_tw, select = c(1,2))
plot(m_hgt_tw, page = 1)

m_stems_tw <- gam(
  stems_total ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2,
)
summary(m_stems_tw)
appraise(m_stems_tw)
draw(m_stems_tw, select = c(1,2))

p_stems_tw <- predict_response(
  m_stems_tw,
  terms = c("time_snc_full_disturbance [all]",    # x-sequence
            "level [all]"
            #"early_class")                    # fix level (or drop to average)
  ))
p_stems_tw <- plot(p_stems_tw, one_plot = TRUE)
p_stems_tw



m_share_early_tw <- gam(
  share_early ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2,
)
summary(m_share_early_tw)
appraise(m_share_early_tw)
draw(m_share_early_tw, select = c(1,2))



m_shade_tw <- gam(
  CWM_shade ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2,
)
summary(m_shade_tw)
appraise(m_shade_tw)
draw(m_shade_tw, select = c(1,2))


hist(both_levels_re2$share_early)

m_cv_tw <- gam(
  cv_hgt ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  #family = tw(link = "log"),  # Tweedie with log-link
  data = filter(both_levels_re2, cv_hgt >0)
)
summary(m_cv_tw)
appraise(m_cv_tw)
draw(m_cv_tw, select = c(1,2))

m_shan_tw <- gam(
  shannon_sp ~ level + 
    s(time_snc_full_disturbance, by = level, k = 5, bs = "cs") +
    s(plot_id, bs = "re"),
  #family = tw(link = "log"),  # Tweedie with log-link
  data = both_levels_re2
)
summary(m_shan_tw)
appraise(m_shan_tw)
draw(m_shan_tw, select = c(1,2))


# nicely plot the predicted values 
p <- predict_response(m_hgt, terms = c("time_snc_full_disturbance [all]"))

plot(p, one_plot = TRUE)









p_plot <- plot_df %>%
  #ungroup() %>%
  select(year, time_snc_full_disturbance, time_snc_part_disturbance, 
         mean_hgt, cv_hgt, shannon_sp, sp_richness,  CWM_shade ,
         CWM_drought) %>%
  pivot_longer(c(-year, 
                 - time_snc_full_disturbance, 
                 - time_snc_part_disturbance),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Plot level') +
  theme_classic2()


p_subplot <- sub_df %>%
  #ungroup() %>%
  select(year, time_snc_full_disturbance, time_snc_part_disturbance, 
         mean_hgt, cv_hgt, shannon_sp, sp_richness,  CWM_shade ,
         CWM_drought) %>%
  pivot_longer(c(-year, 
                 - time_snc_full_disturbance, 
                 - time_snc_part_disturbance),
               names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, leave others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time_snc_full_disturbance, y = value)) +
  geom_boxplot(aes(group = time_snc_full_disturbance ), outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, linewidth = 0.2, col = 'red') +
  stat_summary(fun = \(y) mean(y, na.rm = TRUE), geom = "point", size = 2, col = 'red') +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL, title = 'Subplot level') +
  theme_classic2()

ggarrange(p_plot, p_subplot)
# add comparison between years at plot and subplot levels

make_violin_per_year <- function(df, y,
                                 ylim = NULL,
                                 drop_zeros = FALSE,
                                 # p_y = NULL,           # e.g., p_y = 5  (data units)
                                 p_y_npc = 0.9,       # or p_y_npc = 0.9 (90% up)
                                 p_size = 3,
                                 p_method = "wilcox.test") {
  
  pd <- position_dodge(0.9)
  
  d <- df %>% filter(!is.na({{y}}))
  if (drop_zeros) d <- d %>% filter({{y}} > 0)
  
  p <- ggplot(d, aes(x = year, y = {{y}}, fill = year, color = year)) +
    geom_violin(alpha = 0.5, trim = TRUE, width = 0.8, position = pd) +
    geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.5,
                 width = 0.2, position = pd) +
    theme_grey(base_size = 8)
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  
  p + ggpubr::stat_compare_means(method = p_method, label = "p.format",
                                 size = p_size,
                                 #label.y = p_y,       # use either…
                                 label.y.npc = p_y_npc)  # …or this (0–1)
}


# porovnanie cez subplot
p.dens     <- make_violin_per_year(plot_df, dens_ha, ylim = c(0, 45000), p_y_npc = 0.80) #,   
p.dens
p.height   <- make_violin_per_year(plot_df, mean_hgt, drop_zeros = TRUE, ylim = c(0, 6), 
                                   p_y_npc = 0.15) #,  
p.cv       <- make_violin_per_year(plot_df, cv_hgt)
p.shannon  <- make_violin_per_year(plot_df, shannon_sp)
p.richness <- make_violin_per_year(plot_df, sp_richness)
p.eveness  <- make_violin_per_year(plot_df, evenness_sp)
p.eff      <- make_violin_per_year(plot_df, effective_numbers) #, ylim = c(0, 10)

# Arrange with a shared legend
out_plot <- ggarrange(p.dens, p.height, p.cv, p.shannon, p.richness,p.eveness,p.eff,
                      common.legend = TRUE, legend = "bottom")

annotate_figure(out_plot, top = text_grob("Plot level", 
                                          color = "black", face = "bold", size = 14))


# get summary statistics
out_summary_full <- plot_df %>%
  ungroup() %>%
  select(year, level, mean_hgt, cv_hgt, shannon_sp, sp_richness, evenness_sp,effective_numbers) %>%
  pivot_longer(-c(year, level), names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  group_by(metric, year, level) %>%
  summarise(
    n_total  = n(),
    n_non_na = sum(!is.na(value)),
    n_na     = sum(is.na(value)),
    mean     = mean(value, na.rm = TRUE),
    sd       = sd(value, na.rm = TRUE),
    se       = sd/sqrt(n_non_na),
    median   = median(value, na.rm = TRUE),
    p25      = quantile(value, 0.25, na.rm = TRUE),
    p75      = quantile(value, 0.75, na.rm = TRUE),
    .groups  = "drop"
  ) 

out_summary_full

out_summary_years <- both_levels_re2 %>%
  ungroup() %>%
  select(year,  mean_hgt, cv_hgt, shannon_sp, sp_richness, evenness_sp,effective_numbers) %>%
  pivot_longer(-c(year), names_to = "metric", values_to = "value") %>%
  # keep mean_hgt > 0, others as-is
  filter(!(metric == "mean_hgt" & (is.na(value) | value <= 0))) %>%
  group_by(metric, year) %>%
  summarise(
    n_total  = n(),
    n_non_na = sum(!is.na(value)),
    n_na     = sum(is.na(value)),
    mean     = mean(value, na.rm = TRUE),
    sd       = sd(value, na.rm = TRUE),
    se       = sd/sqrt(n_non_na),
    median   = median(value, na.rm = TRUE),
    p25      = quantile(value, 0.25, na.rm = TRUE),
    p75      = quantile(value, 0.75, na.rm = TRUE),
    .groups  = "drop"
  ) 

out_summary_years


# how does the diversity and composition changes over years?
make_box <- function(df, y, x = "time_snc_full_disturbance",
                     ylim = NULL, fill = "grey95") {
  x_sym <- sym(x); y_sym <- sym(y)
  
  p <- ggplot(df, aes(x = factor(!!x_sym), y = !!y_sym)) +
    geom_boxplot(outlier.shape = NA, fill = fill) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 0.6, color = "grey30") +
    labs(x = "Years since dist.", y = gsub("_", " ", y)) +
    theme()
  theme_classic2()
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}
# !!!

# Example usage ---------------------------------------------------------------
ys <- c("mean_hgt", "cv_hgt", "dens_ha", "sp_richness", "shannon_sp", "CWM_shade")
ylims <- list(mean_hgt = c(0, 6), cv_hgt = c(0, 1.5))  # others left NULL

plots <- lapply(ys, function(v) make_box(filter(plot_df, !time_snc_full_disturbance %in% c(0, 11, 13)),
                                         y = v, ylim = ylims[[v]]))

ggarrange(plotlist = plots, ncol = 3, nrow = 2, labels = 'auto')






both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(y = share_early, x = factor(time_snc_full_disturbance))) +
  geom_boxplot(outlier.shape = NA) +geom_jitter(alpha = 0.5, size = 0.5)

# change over time




# Step 1: Reduce to one row per plot × year
plot_wide <- plot_df %>%
  select(plot_id, year, mean_hgt, shannon_sp, cv_hgt) %>%
  distinct() %>%  # make sure it's one row per plot-year
  pivot_wider(
    names_from = year,
    values_from = c(mean_hgt, shannon_sp,cv_hgt),
    names_sep = "_"
  ) %>%
  mutate(
    delta_hgt = mean_hgt_2025 - mean_hgt_2023,
    delta_div = shannon_sp_2025 - shannon_sp_2023,
    delta_cv = cv_hgt_2025 - cv_hgt_2023
  )

# View result
head(plot_wide)

ggplot(plot_wide, aes(x = delta_hgt, y = delta_div)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = "Δ Mean Height", y = "Δ Shannon Diversity") +
  theme_classic()

ggplot(plot_wide, aes(x = delta_hgt, y = delta_cv)) +
  geom_point() +
  geom_smooth(method = "gam", se = TRUE) +
  labs(x = "Δ Mean CV", y = "Δ Shannon Diversity") +
  theme_classic()


# Make a threshold for legacy effects: > 4 m



#