
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



### Stem density per species -------------
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


p_hist_dist_length <- plot_context_chars %>% 
 # filter(year == 2023) %>%  # to keep half of values
  ggplot(aes(disturbance_length)) + 
  geom_histogram(fill = 'grey90', color = 'black') 

p_hist_dist_year <- plot_context_chars %>%
  #filter(year == 2023) %>%  # to keep half of values
  filter(disturbance_year > 2012) %>% 
  ggplot(aes(disturbance_year)) + 
  geom_histogram(fill = 'grey90', color = 'black') 


p_hist_time_since_dist <- plot_context_chars %>%
 # filter(year == 2023) %>%  # to keep half of values
  filter(disturbance_year > 2012) %>% 
  ggplot(aes(time_snc_full_disturbance)) + 
  geom_histogram(fill = 'grey90', color = 'black') 

ggarrange(p_hist_dist_year, p_hist_dist_length, p_hist_time_since_dist,
          ncol = 3)



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
### Plot level: intensities  ------------------------------------
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
    )
  )



# First extract proportion for applied == 1 per activity
applied_mng_intens_order <- c(
  "clear_intensity",           
  "grndwrk_intensity",         
  "logging_trail_intensity",  
  "planting_intensity",       
  "anti_browsing_intensity"#,   
  #"salvage_intensity",         
  #"protection_intensity",      
  #"management_intensity"
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


#  3) Bind & clean final table with both levels
both_levels_re2 <- bind_rows(sub_df, plot_df) %>%
    mutate(
    level   = factor(level, levels = c("subplot","plot")),
    plot_id = factor(plot_id),
    w       = pmin(pmax(w, 1), 50)   # cap weights so a few dense plots don't dominate
  ) %>% 
  mutate(time_since_f = ifelse(time_snc_full_disturbance <= 2, "early", "later") %>% 
           factor(levels = c("early", "later")),
         dist_length_f = ifelse(disturbance_length <= 2, "abrupt", "continuous") %>% 
           factor(levels = c( "abrupt", "continuous")),
         plant_f = ifelse(planting_intensity < 0.2, "no", "yes") %>% 
           factor(levels = c("no", "yes")))


# add small value to CV = 0
both_levels_re2 <- both_levels_re2 %>%
  mutate(
    cv_hgt_pos = ifelse(cv_hgt <= 0, 1e-4, cv_hgt)#is.na(cv_hgt) | 
  ) %>% 
  filter(!is.na(mean_hgt), !is.na(w), !is.na(cv_hgt), !is.na(dens_m2))


### GAM: Richness by planting?  ------------------
both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = plant_f,
             y = sp_richness,
             fill = plant_f)) + 
  geom_boxplot(notch = T) +
  stat_compare_means(method = "wilcox.test", label.y = max(both_levels_re2$sp_richness, na.rm = TRUE) * 1.05) + 
  geom_jitter(alpha = 0.2) +
  facet_grid(.~level) +
  theme(legend.position = 'none')



p1<-both_levels_re2 %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(planting_intensity),
             y = shannon_sp)) + 
  geom_boxplot(notch = T) 


# effective number sof species
p2<-  both_levels_re2 %>% 
    filter(level == 'plot') %>% 
    ggplot(aes(x = factor(planting_intensity),
               y = effective_numbers)) + 
    geom_boxplot(notch = T) 

p3<-  both_levels_re2 %>% 
    filter(level == 'plot') %>% 
    ggplot(aes(x = factor(planting_intensity),
               y = sp_richness)) + 
    geom_boxplot(notch = T) 
  
  
ggarrange(p1, p2, p3, ncol = 1)
  
hist(both_levels_re2$planting_intensity)


gam_shannon_both <- gam(
  shannon_sp ~ s(planting_intensity, by = level, k = 5) +
    level + year + 
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = gaussian()
)

gam_shannon_both_k6 <- gam(
  shannon_sp ~ s(planting_intensity, by = level, k = 6) +
    level + year + s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML"
)

mgcv::gam.check(gam_shannon_both_k6)

summary(gam_shannon_both_k6)
appraise(gam_shannon_both)
plot(gam_shannon_both, page = 1)
mgcv::gam.check(gam_shannon_both)

draw(gam_shannon_both, select = 1:2)  # shows subplot vs plot smooths


# Predict partial effects for both levels -------------------------
pred_smooths <- ggpredict(
  gam_shannon_both_k6,
  terms = c("planting_intensity [all]", "level")  # variable and grouping factor
)

# Check structure
head(pred_smooths)

# Plot both smooths in one figure ---------------------------------
ggplot(pred_smooths, aes(x = x, y = predicted,
                         color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("#1b9e77", "#7570b3"),
                     labels = c("Subplot level", "Plot level")) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"),
                    labels = c("Subplot level", "Plot level")) +
  labs(
    x = "Planting intensity (0–1)",
    y = "Predicted Shannon diversity",
    color = "Sampling level",
    fill = "Sampling level",
    title = "Effect of planting intensity on species diversity",
    subtitle = "Predicted smooths from GAM with random plot effects"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )


# test GAm: effective # by planting intensity  ------------------------

gam_effective_both <- gam(
  effective_numbers ~ s(planting_intensity, by = level, k = 5) +
    level + year + 
    s(plot_id, bs = "re"),
  data = both_levels_re2,
  method = "REML",
  family = Gamma(link = "log")
)
summary(gam_effective_both)
plot(gam_effective_both, shade = TRUE, pages = 1)
mgcv::gam.check(gam_effective_both)
appraise(gam_effective_both)


# share spruce by planting ----------------------------------------------
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
response_vars <- c(#"dens_m2", 
                   "cv_hgt", "mean_hgt", 
                   "sp_richness", 
                   "shannon_sp", 
                   "evenness_sp",
                   "effective_numbers", 
                   "dens_ha")

# Pivot to long format
both_levels_long <- both_levels_re2 %>%
  pivot_longer(cols = all_of(response_vars),
               names_to = "variable",
               values_to = "value")


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
#Plot diverging lollipops
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
both_levels_long %>% 
  filter(level == 'plot') %>% 
  ggplot(aes(x = factor(disturbance_length), y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.3) +
  facet_wrap(variable ~ . , scales = "free_y") +
  ggtitle("Disturbance length")

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






#