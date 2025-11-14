
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

library(forcats) # order factors





theme_set(theme_classic2(base_size = 10) +
            theme(axis.title = element_text(size = 10),
                  axis.text  = element_text(size = 10)))


# Read data -----------------------------
#dat_overlap  <- fread('outData/full_table_overlap_23_25.csv')

# test 2025/11/03 -> then needs to rename the layer to 'dat'!
dat_overlap  <- fread('outData/full_table_23_25.csv')  # accound for all data points, not just the ovelapping ones
# this can help me to benefit from all disturbnace history sites

# Summary --------------------------------------
## Analyze data: first check up ----------------------
# get master table, having all unique plots and subplots - ven teh empty ones
dat_master_subplot <- dat_overlap %>% 
  dplyr::select(plot, subplot, year) %>% 
  distinct()

table(dat_overlap$year)  

n_plots_total <- length(unique(dat_master_subplot$plot))     # 126
n_plots_total # 208

n_subplots_total <-length(unique(dat_master_subplot$subplot))  # 1250
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
  mutate(n_trees = n_trees,
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


### Stem density per species per top 10 species --------------------------------
df_stem_dens_species <- df %>% 
  group_by(plot, species, n_subplots ) %>%
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


df_stem_dens_species <- df_stem_dens_species %>% 
  ungroup(.) %>% 
  filter(sum_n >0) %>% 
  filter(species %in% v_top_species_overall) %>% 
  dplyr::group_by(species) %>%
  dplyr::mutate(median_stem_density = median(stem_dens, na.rm = TRUE)) %>% 
  dplyr::ungroup(.) %>%
  mutate(species = factor(species, levels = rev(v_top_species_overall))) # Set custom order


df_stem_dens_species2 <- df_stem_dens_species_year %>%
  dplyr::filter(!is.na(log_sum_stem_density) & sum_n > 0) %>%
  dplyr::mutate(year = factor(year, levels = c("2023","2025")),
                # order by mean log density (ascending → highest ends up at the TOP after coord_flip)
                species = forcats::fct_reorder(species, log_sum_stem_density, .fun = mean, na.rm = TRUE))

# 
p_density<-df_stem_dens_species_year2 %>% 
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
  distinct(species, year, plot) %>%           # Unique species × plot combos
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

### Richness: Subplot: Species richness per planting -----------------------
richness_per_subplot <- dat_overlap_n0 %>%
  group_by(plot, subplot, year, planting) %>% #
  summarize(
    species_richness = n_distinct(species[n > 0 & !is.na(n)])#,
    #planting = first(planting)   # first(planting) assuming consistent per subplot/year
  ) %>%
  ungroup()

# Convert planting to factor for clarity
richness_per_subplot <- richness_per_subplot %>%
  mutate(planting = factor(planting, levels = c(0,1), 
                           labels = c("no planting", "planted")))

# Quick plot
ggplot(richness_per_subplot, 
       aes(x = planting, y = species_richness)) +
  geom_boxplot() +
  labs(title = "Species Richness by Planting",
       x = "Planting",
       y = "Species Richness") +
  theme_classic2()


### Richnesss: Plot - per planting intensity ---------------------------------------

### Calculate species richness per subplot -----------------------
richness_per_plot <- dat_overlap_n0 %>%
  group_by(plot, year, planting_intensity) %>% #
  summarize(
    species_richness = n_distinct(species[n > 0 & !is.na(n)])#,
   # planting = first(planting)   # first(planting) assuming consistent per subplot/year
  ) %>%
  ungroup()


# Quick plot
ggplot(richness_per_plot, 
       aes(x = planting_intensity, 
           y = species_richness,
           group = planting_intensity)) +
  geom_boxplot() +
  labs(title = "Plot: Species Richness by Planting",
       x = "Planting",
       y = "Species Richness") +
  theme_classic2()


ggplot(richness_per_plot, 
       aes(x = planting_intensity, 
           y = species_richness)) +
  geom_jitter(size = 0.5) +
  geom_smooth()


#### Get single tree species per subplot/plot ----------------------------------------

# check 613_T2_AH_20250827 - seems to have very many trees!! ()

# identify subplots where i have  single tree species but but several stems
# what happend at the subplot/plot scale? are they heterogenous/homogenous from photos?
# how many subplots have only 1 tree species?
# Identify subplots with only one species present (with at least 1 individual)
single_species_subplots <- dat_overlap %>%
  filter(!is.na(n), n > 0) %>%
  group_by(subplot, year) %>%
  summarize(n_species = n_distinct(species), .groups = "drop") %>%
  filter(n_species == 1)

# View detailed records from those subplots
single_species_sub_details <- dat_overlap %>%
  filter(subplot %in% single_species_subplots$subplot,
         !is.na(n), n > 0) %>%
  select(plot, subplot, species, n) %>%
  arrange(subplot)

# Print or View
#print(single_species_sub_details)
length(unique(single_species_sub_details$subplot))  # 709!!


#### single species per plot 
single_species_plots <- dat_overlap %>%
  filter(!is.na(n), n > 0) %>%
  group_by(plot, year, planting_intensity) %>%
  summarize(n_species = n_distinct(species), .groups = "drop") %>%
  filter(n_species == 1)


# Print or View
#print(single_species_plot_details)
length(unique(single_species_plots$plot))  # 30
hist(single_species_plots$planting_intensity)


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



## Get context and disturbance characteristics (plot) --------------
plot_context_chars <- dat_overlap %>% 
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
                anti_browsing_intensity,
                salvage_intensity,
                protection_intensity,
                management_intensity  
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


# p_hist_dist_length <- plot_context_chars %>% 
#  # filter(year == 2023) %>%  # to keep half of values
#   ggplot(aes(disturbance_length)) + 
#   geom_histogram(fill = 'grey90', color = 'black') 

p_hist_dist_year <- plot_context_chars %>%
  #filter(year == 2023) %>%  # to keep half of values
  #filter(disturbance_year > 2012) %>% 
  ggplot(aes(disturbance_year)) + 
  geom_histogram(fill = 'grey90', color = 'black') 


p_hist_time_since_dist <- plot_context_chars %>%
 # filter(year == 2023) %>%  # to keep half of values
  #filter(disturbance_year > 2012) %>% 
  ggplot(aes(time_since_disturbance)) + 
  geom_histogram(fill = 'grey90', color = 'black') 

ggarrange(p_hist_dist_year, p_hist_time_since_dist,
          ncol = 2)


p_hist_dist_year <- plot_context_chars %>%
  group_by(disturbance_year) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = disturbance_year, y = n)) +
  geom_col(fill = 'grey80', color = 'black') +
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
  geom_col(fill = 'grey80', color = 'black') +
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
ggarrange(
  p_hist_dist_year, p_hist_time_since_dist,
  labels = c("[a]", "[b]"),
  font.label = list(size = 10, face = 'plain'),
  ncol = 2,
  align = 'hv'
)



# Graphics: share of management per subplot and plot level --------------
## Subplot level -----------------------------------------------------------
# Clean up management - if site has been mamaged in 2023, it needs to be managemt (cleared)
# in 2025 as well!!
df_mng_sub <- dat_overlap %>% 
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
           anti_browsing_intensity,   
           salvage_intensity,         
           protection_intensity,      
           management_intensity
           
           )


df_mng_plot <- df_mng_sub %>% 
  select(plot, year, clear_intensity,           
         grndwrk_intensity,         
         logging_trail_intensity,  
         planting_intensity ,       
         anti_browsing_intensity) %>% 
  distinct()


mng_sub_props <- df_mng_sub %>%
  # group by year???
  pivot_longer(cols = c(clear, grndwrk, logging_trail, planting, anti_browsing),
               names_to = "activity",
               values_to = "applied") %>%
  group_by(activity, applied, year) %>%
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
  "grndwrk" = "Soil preparation",
  "planting" = "Planting",
  "anti_browsing" = "Browsing protection",
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
           anti_browsing_intensity,   
           salvage_intensity,         
           protection_intensity,      
           management_intensity
           
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

# igroup into classes
mng_sub_conv <- mng_intensity_props %>%
  mutate(
    intensity_binary = factor(intensity_binary, levels = c('no', 'yes')),
    intensity_class = factor(intensity_class, levels = intensity_levels),
    activity = factor(activity, levels = applied_mng_intens_order)
  ) %>% 
  group_by(activity, intensity_binary) %>% 
  summarise(proportion = sum(proportion, na.rm = T))

# define colors 
fill_colors <- brewer.pal(length(intensity_levels), "YlOrRd")

# Flip "no" values to negative
mng_sub_conv_plot <- mng_sub_conv %>%
  mutate(proportion_plot = if_else(intensity_binary == "no", -proportion, proportion))


activity_intens_labels <- c(
  "clear_intensity" = "Salvage\nlogging",
  "grndwrk_intensity" = "Soil\npreparation",
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


### Graphics: disturbance chars, species composition and management
# Combine
p_combined_disturb_fig <- ggarrange(
  p_hist_dist_year, p_hist_time_since_dist,
 # p_management_bin_plot,
  labels = c("[a]", "[b]"),
  font.label = list(size = 10, face = 'plain'),
  ncol = 2, nrow = 1
)

# Combine
p_combined_management_fig <- ggarrange(
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
ggsave("outFigs/combined_management_fig.png", plot = p_combined_management_fig,
       width = 4, height = 4, units = "in", dpi = 300)

# # Save as SVG
# ggsave("outFigs/combined_management_fig.svg", plot = p_combined_management_fig,
#        width = 7, height = 3.5, units = "in", dpi = 300)




### Plot level: intensities - summary  ------------------------------------

mng_intensity_props <- mng_intensity_props %>%
  mutate(
    intensity_class = cut(
      applied,
      breaks = c(-Inf, 0.19, 0.39, 0.59, 0.79, Inf),
      labels = c("no", "yes", "yes", "yes", "yes"),
      right = TRUE
    )
  )





# Define the factor levels in correct order
intensity_levels <-  c("0–19", "20–39", "40–59", "60–79", "80–100")
low_classes <- c("0–19", "20–39")

# igroup into classes
mng_sub_conv <- mng_intensity_props %>%
  mutate(
    intensity_class = factor(intensity_class, levels = intensity_levels),
    activity = factor(activity, levels = applied_mng_intens_order)
  ) %>% 
  group_by(activity, intensity_class) %>% 
  summarise(proportion = sum(proportion, na.rm = T))

# define colors 
fill_colors <- brewer.pal(length(intensity_levels), "YlOrRd")

mng_shifted <- mng_sub_conv %>%
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
  "clear_intensity" = "Salvage logging",
  "grndwrk_intensity" = "Soil preparation",
  "planting_intensity" = "Planting",
  "anti_browsing_intensity" = "Browsing protection",
  "logging_trail_intensity" = "Logging trail"
)

ggplot(mng_shifted, aes(x = proportion_shifted, y = activity,
                        fill = intensity_class_plot)) +
  geom_col(width = 0.4, color = "black") +
  scale_fill_manual(values = fill_colors, name = "Intensity class",
                    breaks = intensity_levels) +
  geom_vline(xintercept = 0, color = "grey", 
             linewidth = 0.8, lty = 'dashed') +
  ylab('') +
  scale_y_discrete(labels = activity_intens_labels) +   # 👈 this does the relabeling
  scale_x_continuous(labels = abs, name = "Plots share [%]") +
  annotate("text", x = -80, y = 5.5, label = "Low intensity", hjust = 0, size = 2.5, fontface = "bold") +
  annotate("text", x =  80, y = 5.5, label = "High intensity",     hjust = 1, size = 2.5, fontface = "bold") +
  theme_classic2() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )




# Generate new data -----------------------------------------------------------
## Early vs late ( plot, subplot) ---------------

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
df_sub_share_early <-  dat_overlap %>%
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


## Spruce share ( plot, subplot) ----------------------------------------------------------------------


# Calculate spruce share at the plot level
spruce_share_plot <- dat_overlap %>%
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
spruce_share_sub <- dat_overlap %>%
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
  ) %>% 
  select(-stems_total, - stems_with_traits, -trait_coverage)

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


# Create table on subplot and plot level ----------------------------------------
## Field data summary: subplot metrics ----------------------------------------------------------
field_sub_summ <- dat_overlap %>%
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


# is teh change in community shading/drought tolerance driven by planting????
# Plot: Shade ~ Time since disturbance, by planting
x_lab_time_snc_full_dist = "Time since stand\nreplacing disturbance (years)"

#### Effect of planting? 
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



#### Subplot: Time since disturbnace ----------------------------------------

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
  left_join(cwm_plot, by = join_by(plot, year)) %>% 
  left_join(df_mng_sub,by = join_by(plot, year)) %>% 
  left_join(spruce_share_plot, by = join_by(plot, year)) #%>% 


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
  left_join(cwm_plot, by = join_by(plot, year)) %>% 
  left_join(df_mng_sub,by = join_by(plot, year)) %>% 
  left_join(spruce_share_plot, by = join_by(plot, year)) #%>% 

#mutate(cv_hgt = ifelse(is.na(cv_hgt), 0L, cv_hgt)) # replace NA by 0 if stems are missing

# get only context and disturbance information on plot level
df_plot_context <- dat_overlap %>%
  dplyr::select(plot,year,  
                pre_dist_trees_n,  
                area_m2, 
                pre_dist_dens_ha, 
                time_snc_full_disturbance,
                time_snc_part_disturbance,
                disturbance_year, 
                forest_year, 
                disturbance_length,
                protection_intensity, 
                management_intensity) %>% 
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
                                "time_snc_full_disturbance" = "time_snc_full_disturbance") ) %>% 
  left_join(df_mng_plot, by = c("plot_id" = "plot",
                                       "year" = "year")) %>%  # need to separate here by year?
  left_join(spruce_share_sub, by = c("plot_id" = "plot",
                                     "year" = "year",
                                     "ID" = 'subplot'))# %>%

nrow(sub_df)
hist(sub_df$cv_hgt, breaks = 80)

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
                                       "year" = "year"))#


names(plot_df)


#### Bind & clean final table with both levels ----------
both_levels_re2 <- bind_rows(sub_df, plot_df) %>%
    mutate(
      # cap time since disturbance
      time_snc_full_disturbance = pmin(time_snc_full_disturbance, 8),
      level   = factor(level, levels = c("subplot","plot")),
      plot_id = factor(plot_id),
      w       = pmin(pmax(w, 1), 50)   # cap weights so a few dense plots don't dominate
  ) %>% 
  mutate(time_since_f = ifelse(time_snc_full_disturbance <= 2, "early", "later") %>% 
           factor(levels = c("early", "later")),
         dist_length_f = ifelse(disturbance_length <= 2, "abrupt", "continuous") %>% 
           factor(levels = c( "abrupt", "continuous")),
         plant_f = ifelse(planting_intensity < 0.2, "no", "yes") %>% 
           factor(levels = c("no", "yes")), 
         anti_brow_f = ifelse(anti_browsing_intensity < 0.2, "no", "yes") %>% 
           factor(levels = c("no", "yes")),
         grndwrk_f = ifelse(grndwrk_intensity < 0.2, "no", "yes") %>% 
           factor(levels = c("no", "yes")),
         year_f = factor(year)
           )


# add small value to CV = 0
both_levels_re2 <- both_levels_re2 %>%
  mutate(
    cv_hgt_pos = ifelse(cv_hgt <= 0, 1e-4, cv_hgt)#is.na(cv_hgt) | 
  ) %>% 
  filter(!is.na(mean_hgt), !is.na(w), !is.na(cv_hgt), !is.na(dens_m2))


both_levels_re2 <- both_levels_re2 %>%
  mutate(plant_browse_f = interaction(plant_f, anti_brow_f, sep = "_")) %>% 
  mutate(treatment_grp = case_when(
    plant_f == "yes" & anti_brow_f == "yes" ~ "yes",
    TRUE ~ "no"
  ) %>% factor(levels = c("yes", "no")))


both_levels_re2$mng_f <- dplyr::case_when(
  both_levels_re2$plant_browse_f == "no_no"   ~ "none",
  both_levels_re2$plant_browse_f == "yes_yes" ~ "both",
  TRUE                                        ~ "some"
)


both_levels_re2$mng_f <- factor(
  both_levels_re2$mng_f,
  levels = c("none", "some", "both")
)

table(both_levels_re2$mng_f)

##### Boxplot: which management is most important?  --------------------
# Create long table for plotting across factors
both_levels_long <- both_levels_re2 %>%
  filter(level == 'plot') %>% 
  select(plot_id, level,
         dens_m2, cv_hgt, mean_hgt, sp_richness, shannon_sp,
         time_since_f, dist_length_f, plant_f, 
         anti_brow_f,
         grndwrk_f) %>%
  pivot_longer(
    cols = c(dens_m2, cv_hgt, mean_hgt, sp_richness, shannon_sp),
    names_to = "variable",
    values_to = "value"
  )

# check which one of management activities is more important for species richness 
# and vertical structure?
p1<-ggplot(both_levels_long %>% filter(level == "plot"),
       aes(x = plant_f, y = value, fill = plant_f)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),  # nicely formatted p-values
    method = "wilcox.test",                     # or "t.test" depending on data
    label.y.npc = 0.95,                         # position near top of panel
    size = 3
  ) +
  facet_wrap(~ variable, scales = "free_y", ncol = 5) +
  theme_bw() + 
  theme(legend.position = 'none')

p2<-ggplot(both_levels_long %>% filter(level == "plot"),
       aes(x = anti_brow_f, y = value,
           fill = anti_brow_f)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),  # nicely formatted p-values
    method = "wilcox.test",                     # or "t.test" depending on data
    label.y.npc = 0.95,                         # position near top of panel
    size = 3
  ) +
  facet_wrap(~ variable, scales = "free_y", ncol = 5) +  theme_bw() +
  theme(legend.position = 'none')

p3<-ggplot(both_levels_long %>% filter(level == "plot"),
       aes(x = grndwrk_f, y = value,
           fill = grndwrk_f)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),  # nicely formatted p-values
    method = "wilcox.test",                     # or "t.test" depending on data
    label.y.npc = 0.95,                         # position near top of panel
    size = 3
  ) +
  facet_wrap(~ variable, scales = "free_y", ncol = 5) +
  theme_bw() +
  theme(legend.position = 'none')


ggarrange(p1,p2,p3, ncol = 1, nrow = 3)

both_levels_re2 %>% 
  filter(level == "plot") %>% 
  ggplot(aes(x = anti_brow_f,
             y = mean_hgt)) +
  geom_boxplot()

#### GAM: test which management activity is most important ----------------------------
library(broom)
library(dplyr)
library(purrr)

models <- both_levels_re2 %>%
  filter(level == "plot") %>%
  {
    list(
      # shannon_sp   = lm(shannon_sp ~ plant_f + anti_brow_f + grndwrk_f, data = .),
      # sp_richness  = lm(sp_richness ~ plant_f + anti_brow_f + grndwrk_f, data = .),
      # mean_hgt     = lm(mean_hgt ~ plant_f + anti_brow_f + grndwrk_f, data = .),
      # cv_hgt       = lm(cv_hgt ~ plant_f + anti_brow_f + grndwrk_f, data = .)
      effective_numbers  = gam(effective_numbers ~ plant_f + anti_brow_f + grndwrk_f +
                          s(plot_id, bs = "re"),
                        data = ., family = tw(link = "log"), method = "REML"),
      
      sp_richness = gam(sp_richness ~ plant_f + anti_brow_f + grndwrk_f +
                          s(plot_id, bs = "re"),
                        data = ., family = nb(link = "log"), method = "REML"),
      
      mean_hgt    = gam(mean_hgt ~ plant_f + anti_brow_f + grndwrk_f +
                          s(plot_id, bs = "re"),
                        data = ., family = tw(link = "log"), method = "REML"),
      
      cv_hgt      = gam(cv_hgt_pos ~ plant_f + anti_brow_f + grndwrk_f +
                          s(plot_id, bs = "re"),
                        data = ., family = tw(link = "log"), method = "REML")
    )
  }

##### check model assumptions  - from gam model

# Extract model summaries
model_stats <- map_dfr(names(models), function(nm) {
  m <- models[[nm]]
  sm <- summary(m)
  
  tibble(
    response = nm,
    r2_adj = sm$r.sq,                     # adjusted R² for GAM
    dev_expl = sm$dev.expl * 100          # explained deviance in %
  )
})

model_stats


model_results_df <- map_dfr(models, tidy, parametric = TRUE, .id = "response") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    term = recode(term,
                  "plant_fyes" = "Planting",
                  "anti_brow_fyes" = "Anti-browsing",
                  "grndwrk_fyes" = "Groundwork"),
    # Add a tiny offset so very small bars are still visible
    estimate_adj = ifelse(abs(estimate) < 0.01, sign(estimate + 1e-6) * 0.01, estimate),
    lower_adj = ifelse(abs(estimate) < 0.01, 0, lower),
    upper_adj = ifelse(abs(estimate) < 0.01, 0, upper),
    response = recode(response,
                      "effective_numbers" = "Species diversity\n[Effective #]",
                      "sp_richness"       = "Species richness\n[#]",
                      "mean_hgt"          = "Mean height\n[m]",
                      "cv_hgt"            = "Height variability\n[CV, %]"
                      ),
    # order facets manually
    response = factor(response,
                      levels = c("Mean height\n[m]",
                                 "Height variability\n[CV, %]",
                                 "Species diversity\n[Effective #]",
                                 "Species richness\n[#]"))
  )

ggplot(model_results_df, aes(x = term, y = estimate, fill = term)) +
  geom_hline(yintercept = 0, col = 'grey', lwd = 0.5, lty = 'solid') +
  geom_col(position = "dodge", width = 0.7) +
   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 0.5) +
  facet_wrap(~ response, ncol = 4) + #scales = "free_y"
  theme_classic2(base_size = 8) +
  labs(
    x = "",
    y = "Estimated effect (parametric term)",
    title = ""
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text( size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )


# how much is the change in %? 

effects_df <- model_results_df %>%
  select(response, term, estimate) %>%
  mutate(
    perc_change = (exp(estimate) - 1) * 100
  ) %>%
  group_by(response, term) %>%
  summarise(mean_effect = mean(perc_change), .groups = "drop")
effects_df


table(both_levels_re2$treatment_grp)
table(both_levels_re2$plant_browse_f )







### GAM Planting -------------------------------------
##### GAM: height x time_since   ------------------
# test by presence_absence of planting
gam_mean_hgt <- gam(
  mean_hgt  ~ s(planting_intensity,  k = 3) +
    s(grndwrk_intensity        ,  k = 3) +
    s(anti_browsing_intensity  ,  k = 3) +
    s(time_snc_full_disturbance, by = plant_f,  k = 3)+
    level + 
    plant_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

# test for present/absent groundwork
gam_mean_hgt_comb <- gam(
  mean_hgt  ~ #s(planting_intensity,  k = 3) +
   # s(grndwrk_intensity        ,  k = 3) +
    #s(anti_browsing_intensity  ,  k = 3) +
    s(time_snc_full_disturbance, by = treatment_grp ,  k = 3)+
    level + 
   # plant_f +
    treatment_grp +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_mean_hgt_plt <- gam(
  mean_hgt  ~ 
    s(time_snc_full_disturbance, by = plant_f ,  k = 3)+
    level + 
    plant_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_mean_hgt_mng <- gam(
  mean_hgt  ~ 
    s(time_snc_full_disturbance, by = mng_f ,  k = 3)+
    level + 
    mng_f +
    year_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

# Generate predictions
df <- ggpredict(gam_mean_hgt_mng, terms = c("time_snc_full_disturbance [0:8]", "mng_f"))

# Quick plot
plot(df)



AIC(gam_mean_hgt_comb, gam_mean_hgt_plt, gam_mean_hgt_mng)
summary(gam_mean_hgt)
summary(gam_mean_hgt_mng)
plot(gam_mean_hgt_mng, shade = TRUE, pages = 1)
plot(gam_mean_hgt_comb, shade = TRUE, pages = 1)
appraise(gam_mean_hgt_plt)
appraise(gam_mean_hgt_comb)

pred_smooths <- ggpredict(
  gam_mean_hgt,
  terms = c("time_snc_full_disturbance [1:8]",
            "plant_f [all]"),  # variable and grouping factor
  condition = list(level = "plot")
  )

p_gam_height <- plot(pred_smooths)
p_gam_height


# pred_smooths <- ggpredict(
#   gam_mean_hgt2,
#   terms = c("time_snc_full_disturbance [0:8]",
#             "grndwrk_f [all]")  # variable and grouping factor
# )

#p_gam_height <- plot(pred_smooths)
#p_gam_height


# Convert to data frame for ggplot
pred_df <- as.data.frame(pred_smooths)

# Custom ggplot (same structure as your CV plot)
p_gam_height <- ggplot(pred_df,
                       aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.1) +
  labs(
    x = "Years since disturbance",
    y = "Mean tree height [m]",
    colour = "Groundwork",
    fill = "Groundwork"
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.title = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

p_gam_height


# is teh pattern real? eg planted sites are slightly taller then no planted ones?

both_levels_re2 %>%
  filter(level == "plot") %>%
  group_by(plant_f, time_snc_full_disturbance) %>%
  summarise(n = n(), mean_hgt = mean(mean_hgt, na.rm = TRUE))



##### GAM: CV_height x time_since   ------------------

# as binary/continuous part

#  1. Binary part: presence/absence 
both_levels_re2 <- both_levels_re2 %>%
  mutate(cv_hgt_present = as.integer(cv_hgt > 0))

table(both_levels_re2$cv_hgt_present)

gam_cv_hgt_bin <- gam(
  cv_hgt_present ~   s(planting_intensity,  k = 3) +
    s(grndwrk_intensity        ,  k = 3) +
    s(anti_browsing_intensity  ,  k = 3) +
    s(time_snc_full_disturbance, by = plant_f,  k = 3)+
    level + 
    plant_f +
    s(plot_id, bs = "re"),
  data   = both_levels_re2,
  method = "REML",
  family = binomial(link = "logit")
)

# 2. Positive part: only where cv_hgt_pos > 0
gam_cv_hgt_pos <- gam(
  cv_hgt_pos ~ s(planting_intensity,  k = 3) +
    s(grndwrk_intensity        ,  k = 3) +
    s(anti_browsing_intensity  ,  k = 3) +
    s(time_snc_full_disturbance, by = plant_f,  k = 3)+
    level + 
    plant_f +
    s(plot_id, bs = "re"),
  data   = subset(both_levels_re2, cv_hgt > 0),
  method = "REML",
  family = tw(link = "log")   # or tw(link="log") if still skewed
)

summary(gam_cv_hgt_pos)
plot(gam_cv_hgt_pos, page = 1)
summary(gam_cv_hgt_bin)
plot(gam_cv_hgt_bin, page = 1)
appraise(gam_cv_hgt_pos)

pred_bin <- ggpredict(
  gam_cv_hgt_bin,
  terms = c("time_snc_full_disturbance [0:8]", "plant_f"),
  condition = list(level = "plot")
)

plot(pred_bin) + 
  labs(y = "Probability of CV height > 0", x = "Years since disturbance")


pred_pos <- ggpredict(
  gam_cv_hgt_pos,
  terms = c("time_snc_full_disturbance [0:8]", "plant_f"),
  condition = list(level = "plot")
)

plot(pred_pos) +
  labs(y = "Mean CV height (given >0)", x = "Years since disturbance")


# 1️⃣ Predict probability (binary model)
pred_bin <- ggpredict(
  gam_cv_hgt_bin,
  terms = c("time_snc_full_disturbance [0:8]", "plant_f")
) %>%
  rename(prob = predicted,
         prob_low = conf.low,
         prob_high = conf.high)

# 2️⃣ Predict magnitude (positive model)
pred_pos <- ggpredict(
  gam_cv_hgt_pos,
  terms = c("time_snc_full_disturbance [0:8]", "plant_f")
) %>%
  rename(mean_pos = predicted,
         pos_low = conf.low,
         pos_high = conf.high)

# 3️⃣ Merge and compute expected CV = probability × conditional mean
pred_combined <- left_join(pred_bin, pred_pos, by = c("x", "group")) %>%
  mutate(
    expected_cv = prob * mean_pos,
    lower = prob_low * pos_low,
    upper = prob_high * pos_high
  )

# 4️⃣ Plot combined expected CV
p_gam_cv <- ggplot(pred_combined,
       aes(x = x, y = expected_cv*100, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lower*100, ymax = upper*100), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.1) +
  labs(
    x = "Years since disturbance",
    y = "Height variability [CV, %]",
    colour = "Planting",
    fill = "Planting"
  ) +
   theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )
p_gam_cv





##### GAM: Effective number of species x time_since   ------------------
gam_effective_both <- gam(
  effective_numbers ~ s(planting_intensity,  k = 3) +
    s(grndwrk_intensity        ,  k = 3) +
    s(anti_browsing_intensity  ,  k = 3) +
    s(time_snc_full_disturbance, by = plant_f,  k = 3)+
    level + 
    plant_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = Gamma(link = "log")
)
summary(gam_effective_both)
plot(gam_effective_both, shade = TRUE, pages = 1)
mgcv::gam.check(gam_effective_both)
appraise(gam_effective_both)


pred_smooths <- ggpredict(
  gam_effective_both,
  terms = c("time_snc_full_disturbance [0:8]", "plant_f"),  # variable and grouping factor
  condition = list(level = "plot")
  )

#p_gam_eff <- plot(pred_smooths)
# Convert predictions to data frame
pred_df <- as.data.frame(pred_smooths)

# Build the ggplot manually
p_gam_eff <- ggplot(pred_df,
                    aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.1) +
  labs(
    x = "Years since disturbance",
    y = "Effective species number",
    colour = "Planting",
    fill = "Planting"
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

p_gam_eff



##### GAM: Richness  x time_since   ------------------
gam_richness <- gam(
  sp_richness  ~ s(planting_intensity,  k = 3) +
    s(grndwrk_intensity        ,  k = 3) +
    s(anti_browsing_intensity  ,  k = 3) +
    s(time_snc_full_disturbance, by = plant_f,  k = 3)+
    level + 
    plant_f +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)
summary(gam_richness)
plot(gam_richness, shade = TRUE, pages = 1)
mgcv::gam.check(gam_richness)
appraise(gam_richness)


pred_smooths <- ggpredict(
  gam_richness,
  terms = c("time_snc_full_disturbance [0:8]", "plant_f"),
  condition = list(level = "plot")# variable and grouping factor
)

#p_gam_richness <- plot(pred_smooths)

# Convert to data frame for ggplot
pred_df <- as.data.frame(pred_smooths)

# Build the plot manually with consistent style
p_gam_richness <- ggplot(pred_df,
                         aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.1) +
  labs(
    x = "Years since disturbance",
    y = "Species richness [#]",
    colour = "Planting",
    fill = "Planting"
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

p_gam_richness


####### Plot all GAMs planting  ------------------------------------------------------------
p_gam_combined <- ggarrange(p_gam_height, p_gam_cv, p_gam_eff,p_gam_richness, common.legend = TRUE,
          ncol = 2, nrow = 2,
          align = "hv", legend = 'bottom',
          labels = c("[a]", "[b]", "[c]", "[d]"),
          font.label = list(size = 11, face = "plain", color ="black"),
          label.y = 0.95,
          label.x = 0.15)

annotate_figure(
  p_gam_combined,
  bottom = text_grob("Planting", face = "plain", size = 11)
)


#### GAM interaction plantoing&browsing  ------------------------------------
##### Height & time since by planting&anti-browsing -----------
gam_mean_hgt_int <- gam(
  mean_hgt ~ 0 + treatment_grp +  # removes global intercept
    s(time_snc_full_disturbance, by = treatment_grp, k = 3) +
    level +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)

gam_mean_hgt0 <- gam(
  mean_hgt ~ 0 + treatment_grp + 
    s(time_snc_full_disturbance, by = treatment_grp, k = 3) +
    level +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)
table(both_levels_re2$treatment_grp)
table(both_levels_re2$plant_browse_f )

AIC(gam_mean_hgt_int, gam_mean_hgt0)

summary(gam_mean_hgt_int)
summary(gam_mean_hgt0)

plot(gam_mean_hgt_int)

df_hgt <- ggpredict(gam_mean_hgt_int, 
                    terms = c("time_snc_full_disturbance [1:7]", "treatment_grp"))
plot(df_hgt)

##### CV_hgt vs time  ----------------------------------------
gam_cv_hgt_bin0 <- gam(
  cv_hgt_present ~   0 + treatment_grp + 
    s(time_snc_full_disturbance, by = treatment_grp, k = 3) +
    level +
    s(plot_id, bs = "re"),,
  data   = both_levels_re2,
  method = "REML",
  family = binomial(link = "logit")
)

# 2. Positive part: only where cv_hgt_pos > 0
gam_cv_hgt_pos0 <- gam(
  cv_hgt_pos ~ 0 + treatment_grp + 
    s(time_snc_full_disturbance, by = treatment_grp, k = 3) +
    level +
    s(plot_id, bs = "re"),
  data   = subset(both_levels_re2, cv_hgt > 0),
  method = "REML",
  family = tw(link = "log")   # or tw(link="log") if still skewed
)


##### EFficient number vs time ----------------------------------
gam_effective_both0 <- gam(
  effective_numbers ~ 0 + treatment_grp + 
    s(time_snc_full_disturbance, by = treatment_grp, k = 3) +
    level +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = Gamma(link = "log")
)

##### GAM: Richness  x time_since   ------------------
gam_richness0 <- gam(
  sp_richness  ~  0 + treatment_grp + 
    s(time_snc_full_disturbance, by = treatment_grp, k = 3) +
    level +
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = tw(link = "log")
)



###### Gam plotting  function ---------------------

plot_gam_smooth <- function(model, y_lab, group_var, level_value = "plot") {
  # Get predictions
  pred <- ggpredict(
    model,
    terms = c("time_snc_full_disturbance [0:8]", group_var),
    condition = list(level = level_value)
  )
  
  # Convert to data frame
  df <- as.data.frame(pred)
  
  # Make the plot
  ggplot(df, aes(x = x, y = predicted, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1.1) +
    labs(
      x = "Years since disturbance",
      y = y_lab,
      colour = group_var,
      fill = group_var
    ) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    )
}


# Mean height
p_hgt <- plot_gam_smooth(gam_mean_hgt0, 
                         y_lab = "Mean height [m]", group_var = "treatment_grp")

# CV height — binary part
df_bin <- ggpredict(gam_cv_hgt_bin0, terms = c("time_snc_full_disturbance [1:7]", "treatment_grp"),
                    condition = list(level = "plot"))

# CV height — positive part
df_pos <- ggpredict(gam_cv_hgt_pos0, terms = c("time_snc_full_disturbance [1:7]", "treatment_grp"),
                    condition = list(level = "plot"))

# Combine binary and positive parts
df_cv <- merge(df_bin, df_pos, by = c("x", "group"))
df_cv$predicted <- df_cv$predicted.x * df_cv$predicted.y
df_cv$conf.low  <- df_cv$conf.low.x * df_cv$conf.low.y
df_cv$conf.high <- df_cv$conf.high.x * df_cv$conf.high.y

# CV plot
p_cv <- ggplot(df_cv, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.1) +
  labs(
    x = "Years since disturbance",
    y = "Height variability [CV, %]",
    colour = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )

# Effective number
p_eff <- plot_gam_smooth(gam_effective_both0, y_lab = "Effective species number", group_var = "treatment_grp")

# Richness
p_rich <- plot_gam_smooth(gam_richness0, 
                          y_lab = "Species richness [#]", 
                          group_var = "treatment_grp")









## make categories for planting or browsing and both/none -----------

# Step 1: Create the treatment group at plot level
plot_treatment <- both_levels_re2 %>%
  filter(level == "plot") %>%
  mutate(
    treatment = case_when(
      plant_f == "no" & anti_brow_f == "no" ~ "none",
      plant_f == "yes" & anti_brow_f == "no" ~ "planting",
      plant_f == "no" & anti_brow_f == "yes" ~ "anti_brw",
      plant_f == "yes" & anti_brow_f == "yes" ~ "both"
    ),
    treatment = factor(treatment, levels = c("none", 
                                             "planting", 
                                             "anti_brw", 
                                             "both"))
  )

# Step 2: Check summary
summary(plot_treatment$treatment)

#
p_height <- ggplot(plot_treatment, aes(x = treatment, y = mean_hgt, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
 # geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.y = max(plot_treatment$mean_hgt, na.rm = TRUE) * 1.05) +
  labs(y = "Mean Height [m]", x = NULL) +
  theme_bw() +
  theme(legend.position = "none")

# 
p_cv <- ggplot(plot_treatment, aes(x = treatment, y = cv_hgt*100, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
 # geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.y = max(plot_treatment$cv_hgt*100, na.rm = TRUE) * 1.05) +
  labs(y = "Height CV [%]", x = NULL) +
  theme_bw() +
  theme(legend.position = "none")




# Step 3: Plot richness and diversity (Shannon or effective_numbers)
p_richness <- ggplot(plot_treatment, aes(x = treatment, y = sp_richness, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
 # geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  stat_compare_means(method = "kruskal.test", 
                     label = "p.format",
                     label.y = max(plot_treatment$sp_richness, na.rm = TRUE) * 1.05) +
  labs(y = "Species Richness [#]", x = NULL) +
  theme_bw() +
  theme(legend.position = "none")

p_diversity <- ggplot(plot_treatment, aes(x = treatment, y = effective_numbers, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
 # geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  stat_compare_means(method = "kruskal.test", 
                     label = "p.format",
                     label.y = max(plot_treatment$effective_numbers, na.rm = TRUE) * 1.05) +
  labs(y = "Effective Species [#]", x = NULL) +
  theme_bw() +
  theme(legend.position = "none")

# Step 4: Display both plots
ggarrange(
   p_height, p_cv,p_richness, p_diversity,
  labels = c("[a]", "[b]", "[c]", "[d]"),
  font.label = list(size = 10, face = "plain"),
  ncol = 2, nrow = 2,
  align = "hv",
  common.legend = FALSE,
  legend = "none"
)


# Dunn test to for pairwise group comparison
library(rstatix)
# Dunn test for mean height
# Full pairwise Dunn test (all groups)
dunn_height <- plot_treatment %>%
  dunn_test(mean_hgt ~ treatment, p.adjust.method = "BH")

dunn_cv <- plot_treatment %>%
  dunn_test(cv_hgt ~ treatment, p.adjust.method = "BH")

# View results
dunn_height
dunn_cv


#### share spruce by planting ----------------------------------------------
both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = plant_f,
             y = spruce_share,
             fill = plant_f)) + 
  geom_boxplot(notch = T) +
  stat_compare_means(method = "wilcox.test", 
                     label.y = max(both_levels_re2$spruce_share, na.rm = TRUE) * 1.05) + 
  geom_jitter(alpha = 0.2) +
  facet_grid(.~level) +
  theme(legend.position = 'none')


both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(planting_intensity),
             y = spruce_share,
             fill = factor(planting_intensity))) + 
  geom_boxplot(notch = T)  + 
  theme(legend.position = 'none')



both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = planting_intensity,
             y = spruce_share)) + 
  geom_point() +
  geom_smooth(method = 'lm')
  




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








# how did disturbnace lenght affected post-dsiturbance recovery? ---------------
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


library(gghalves)

# All indicators:  Create half-violin plot ------------------------------------------------
p_violin <- 
  ggplot(both_levels_long_capped, 
         aes(x = level,
             y = value_capped,
             fill = level)) +
    gghalves::geom_half_violin(
      side = "l",  # use "r" for right
      #fill = 'grey', 
      color = NA, 
      trim = FALSE,
      scale = 'width'
    ) +
  geom_boxplot(width = 0.1, 
               #fill = 'grey',
               outlier.shape = NA, 
               color = "black") +
    scale_fill_manual(
      values = c("subplot" = "grey50", 
                 "plot" = "grey80")
    ) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4, nrow = 1) +
  #scale_fill_manual(values = c("2023" = "skyblue", "2025" = "tomato")) +
  theme_classic(base_size = 10) +
  labs(x = NULL, y = NULL, fill = "Level") +
  theme(legend.position = 'none',
    strip.text = element_text(size = 9, 
                              face = "plain"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

# Print plot
p_violin

ggsave(
  filename = "outFigs/p_violin.png",
  plot = p_violin,
  width = 7,       # Adjust width as needed
  height = 2.5,      # Adjust height as needed
  units = "in",
  dpi = 300        # High resolution for publication
)
  


# make lollipop plot: two levels
plot_level_means <- both_levels_re2 %>%
  group_by( level) %>%
  summarize(
    shannon_sp = mean(shannon_sp, na.rm = TRUE),
    effect_n_sp = mean(effective_numbers, na.rm = TRUE),
    richness = mean(sp_richness, na.rm = TRUE), # assuming w = effective species number
    dens_m2 = mean(dens_m2, na.rm = TRUE),
    cv_hgt = mean(cv_hgt, na.rm = TRUE),
    mean_hgt = mean(mean_hgt, na.rm = TRUE),
    .groups = "drop"
  )

plot_level_long <- plot_level_means %>%
  pivot_longer(cols = c(shannon_sp, 
                        effect_n_sp,
                        richness, 
                        dens_m2, 
                        cv_hgt, 
                        mean_hgt),
               names_to = "metric",
               values_to = "value")



# Convert subplot values to negative
level_means_diverged <- plot_level_long %>%
  mutate(value = ifelse(level == "subplot", -value, value),
         level = fct_relevel(level, "subplot", "plot"))


level_means_diverged <- level_means_diverged %>%
  mutate(metric = fct_reorder(metric, abs(value), .desc = F))


### Plot vs subplot: diverging lollipop  --------------------------------
ggplot(level_means_diverged, aes(x = value, y = fct_rev(metric), fill = level,
                                 color = level)) +
  geom_segment(aes(x = 0, xend = value, y = metric, yend = metric),
               linewidth = 0.8, 
               color = "grey",
               lty = 'dashed'
               ) +
  scale_y_discrete(labels = c(
    "cv_hgt"      = "CV height [dim.]",
    "dens_m2"     = "Stem density [m²]",
    "effect_n_sp" = "Effective species [#]",
    "mean_hgt"    = "Mean height [m]",
    "richness"    = "Species richness [#]",
    "shannon_sp"  = "Shannon diversity [dim.]"
  )) +
  geom_vline(xintercept = 0) +
  geom_point(size = 4) +
  scale_x_continuous(labels = abs) +  # show positive axis labels
  labs(x = "Mean Value", y = "",
       title = "") +
  annotate("text", x = -2, y = 6.5, label = "Subplot", hjust = 0, size = 2.5, fontface = "bold") +
  annotate("text", x =  3, y = 6.5, label = "Plot",     hjust = 1, size = 2.5, fontface = "bold") +
  theme_classic2() +
  theme(legend.position = "none")

# Boxplot - all vars vs time since, disturb. length ----------------------------
# how does the disturbnac elength affects structure and composition?
# how does time since stand replacing disturbnace affects str & composition?
both_levels_long %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance), y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.3) +
  facet_wrap(variable ~ . , scales = "free_y") +
  ggtitle("Time since stand replacing disturb.")

# time since disturbance start 
both_levels_long %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(time_snc_part_disturbance), y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.3) +
  facet_wrap(variable ~ . , scales = "free_y") +
  ggtitle("Time since disturbnace start")




# do time over axis, length as point size
both_levels_long %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(time_snc_full_disturbance), 
             y = factor(time_snc_part_disturbance),
             size = disturbance_length)) +
 geom_point()


