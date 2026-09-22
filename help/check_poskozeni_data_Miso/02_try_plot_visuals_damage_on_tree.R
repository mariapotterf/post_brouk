# test plots visualization example

gc()

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(stringr)
library(tibble)
library(tidyr)

library(forcats)
library(RColorBrewer)


# get list of the icdentified patogenes from Michal
df <- read_excel("inData_Michal/Patogeny_hmyz_postbrouk_1_clean.xlsx",
                 sheet = 1,        # or a sheet name, e.g. "Sheet1"
                 na = c("", "NA"))

# need to read orifinal data to properly link the pathogenes with tree species via sample number 
df_trees  <- fread("inData_Michal/poskoz_jedinci_clean.csv") %>%
  mutate(individual_id = row_number())   # one row = one surveyed individual (confirm 'n' isn't pre-aggregating multiple stems per row before relying on this)


# sum up the get damage counts and % of damage -----------------
summary_damaged <- df_trees %>%
  group_by(species) %>%
  summarise(
    n_examined         = sum(n, na.rm = TRUE),
    n_damaged_Terminal = sum(tot_dmg_terminal, na.rm = TRUE),
    n_damaged_Foliage  = sum(tot_dmg_foliage, na.rm = TRUE),
    n_damaged_Stem     = sum(dmg_stem_numeric * n, na.rm = TRUE),  # CHECK: does a 1 here mean all `n` individuals in the row had stem damage, or just ">=1 of them"?
    .groups = "drop"
  ) %>%
  pivot_longer(starts_with("n_damaged_"), names_to = "damage_type", values_to = "n_damaged",
               names_prefix = "n_damaged_") %>%
  mutate(pct_damaged = 100 * n_damaged / n_examined,
         species = factor(species, levels = rev(sort(unique(species)))))



# ---- proper species order: by sample size, largest first ----
# species_rank <- summary_damaged %>%
#   group_by(species) %>%
#   summarise(total_damaged = sum(n_damaged), .groups = "drop") %>%
#   arrange(desc(total_damaged)) %>%
#   pull(species)

# teh foliage damage is prevalent, so order accordingly; check Terminal as well
species_rank <- summary_damaged %>%
  #filter(damage_type == "Terminal") %>%
  filter(damage_type == "Foliage") %>%
  arrange(desc(n_damaged)) %>%
  pull(species)


# ---- keep only the top N species, for a quick look ----
top_n_species <- 10   # try 5, then 7

summary_damage_df <- summary_damaged %>%
  filter(species %in% head(species_rank, top_n_species)) %>%
  mutate(species = factor(species, levels = rev(head(species_rank, top_n_species))))
# rev() again, same reason as before: first species in the ranking ends up plotted at the top

# get colors per species
# ---- fixed species -> color mapping, same order as species_rank, reused across all 6 panels ----
species_levels <- levels(summary_damage_df$species)
library(viridisLite)

species_colors <- setNames(
  viridis(length(species_levels), option = "D", begin = 0.05, end = 0.90, alpha = 0.85),
  species_levels
)
species_colors

# > species_colors
# VRX          JD          KL          JR 
# "#481467D9" "#46337FD9" "#3C4F8AD9" "#31688ED9" 
# MD          BO          BR          BK 
# "#277F8ED9" "#1F958BD9" "#25AC82D9" "#49C16ED9" 
# DB          SM 
# "#7DD34FD9" "#BBDF27D9" 


# ---- 2. one panel-builder, reused for every column ----
# ---- 2. one panel-builder, reused for every column ----
make_panel <- function(data, damage_type_i, value_col, x_lab, x_limits,
                       show_y_labels = FALSE, title = NULL, ref_line = NULL) {
  d <- data %>% filter(damage_type == damage_type_i)
  p <- ggplot(d, aes(x = .data[[value_col]], y = species, fill = species))
  
  if (!is.null(ref_line)) {
    p <- p + geom_vline(xintercept = ref_line, linetype = "dashed",
                        color = "grey60", linewidth = 0.4)
  }
  
  p +
    geom_col(width = 0.8, col = 'black') +
    scale_x_continuous(
      limits = x_limits,
      breaks = scales::breaks_extended(n = 3)(x_limits),   # <- nice round min/mid/max instead of the exact mean
      labels = scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.08))
    ) +
    scale_fill_manual(values = species_colors, guide = "none") +
    labs(x = x_lab, y = NULL, title = title) +
    theme_classic2(base_size = 10) +
    theme(
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y  = if (show_y_labels) element_text(size = 9) else element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title   = element_text(size = 10, face = "bold", hjust = 0),
      plot.margin  = margin(t = 5, r = 2, b = 5, l = 2)
    )
}

make_panel_lollipop <- function(data, damage_type_i, value_col, x_lab, x_limits,
                                show_y_labels = FALSE, title = NULL,
                                ref_line = NULL) {
  d <- data %>% filter(damage_type == damage_type_i)
  ggplot(d, aes(x = .data[[value_col]], y = species, color = species)) +
    geom_segment(aes(x = 0, xend = .data[[value_col]], yend = species), linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_x_continuous(
      limits = x_limits,
      breaks = scales::breaks_extended(n = 3)(x_limits),
      labels = scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.08))
    ) +
    scale_color_manual(values = species_colors, guide = "none") +   # note: color, not fill, for points/lines
    labs(x = x_lab, y = NULL, title = title) +
    theme_classic2(base_size = 10) +
    theme(
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y  = if (show_y_labels) element_text(size = 9) else element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title   = element_text(size = 10, face = "bold", hjust = 0),
      plot.margin  = margin(t = 5, r = 2, b = 5, l = 2)
    )
}
count_max <- max(summary_damage_df$n_damaged) * 1.05
pct_max   <- max(summary_damage_df$pct_damaged) * 1.05

p_term_n   <- make_panel(summary_damage_df, "Terminal", "n_damaged",   "#\nCount",    c(0, count_max), show_y_labels = TRUE, title = "Terminal",
                         ref_line = 50)
#p_term_pct <- make_panel(summary_damage_df, "Terminal", "pct_damaged", "%\ndamaged", c(0, 40), title = "")
p_fol_n    <- make_panel(summary_damage_df, "Foliage",  "n_damaged",   "#\nCount",    c(0, count_max), title = "Foliage",
                         ref_line = 50)
#p_fol_pct  <- make_panel(summary_damage_df, "Foliage",  "pct_damaged", "%\ndamaged", c(0, 40), title = "")
p_stem_n   <- make_panel(summary_damage_df, "Stem",     "n_damaged",   "#\nCount",    c(0, count_max), title = "Stem",
                         ref_line = 50)
#p_stem_pct <- make_panel(summary_damage_df, "Stem",     "pct_damaged", "%\ndamaged", c(0, 40), title = "")


p_term_pct <- make_panel_lollipop(summary_damage_df, "Terminal", "pct_damaged","%\ndamaged", c(0, pct_max),title = "")
p_fol_pct  <- make_panel_lollipop(summary_damage_df, "Foliage",  "pct_damaged", "%\ndamaged", c(0, pct_max), title = "")
p_stem_pct <- make_panel_lollipop(summary_damage_df, "Stem",     "pct_damaged", "%\ndamaged", c(0, pct_max), title = "")


ggarrange(
  p_term_n, p_term_pct,
  p_fol_n,  p_fol_pct,
  p_stem_n, p_stem_pct,
  ncol = 6, nrow = 1,
  widths = c(2.5, 1, 2, 1, 2, 1)
)








# Clean up pathogens data --------------------------------------------


# teh sample name is not corretly allocated to subplot!!! need to get back the subplot number
df <- df %>%
  separate_wider_delim(
    Subplot, delim = "_",
    names = c("tablet_subplot", "who_collected", "date_raw", "subplot_number"),
    cols_remove = FALSE
  ) %>%
  mutate(
    tablet_subplot = str_remove(tablet_subplot, "^T") %>% as.integer(),
    subplot_ch = as.character(subplot_number),
   # date_subplot   = ymd(date_raw)     # turns "20250827" into a proper Date, useful later
  ) %>%
  select(-date_raw, -tablet_subplot, -who_collected)


# reshape the three sample-code columns into long format, keep only real vzorek codes
df_trees_samples <- df_trees %>%
  select(species, plot, subplot, vegtype, n,
         dmg_term_sample, dmg_stem_sample, dmg_foliage_sample) %>%
  pivot_longer(cols = starts_with("dmg_") & ends_with("_sample"),
               names_to = "sample_location", values_to = "vzorek") %>%
  filter(str_detect(vzorek, "^T\\d")) |>     # drops "Ožnut", "Okus zvěří" etc - those are browsing notes, not pathogen codes
 select(-sample_location)


#View(df)

# check if patogenes names are correct
unique(df$patogen)
unique(df$vzorek)

# one sample has two records: "Lophodermium pinastri, Dothideomycetes indet"
# need to split it in two
# also - decode the sample name back: T1-XXX-O-SP-1 = T1-ploska-vyska_obnovy-poskozeni_kde 

df_pathogen <- df |> 
  filter(vzorek != "Blank")

df_blank <- df |> 
  filter(vzorek == "Blank") |> 
  select(-vzorek, -Species) |> 
  select(subplot_number, Subplot, patogen) |> 
  distinct()


# Vrstva
# O – obnova (<2 m)
# P – pokrocila obnova
# (< 2 m - < 10 DBH)
# D – dospele stromy
# (< 10 cm DBH)
# 
# Drevina
# BK
# JD
# SM
# 
# Poskozeni
# 1 – terminal (do 2 m)
# 2 – kmen
# 3 – báze kmene
# 4 – olistění

# reconstruct the species from sample name
# --- species abbreviation lookup, from "Seznam zkratek dřevin" ---
species_lookup <- tribble(
  ~code, ~czech_name,                ~latin_name,
  "SM",  "smrk ztepilý",             "Picea abies",
  "JD",  "jedle bělokorá",           "Abies alba",
  "DG",  "douglaska tisolistá",      "Pseudotsuga menziesii",
  "BO",  "borovice lesní",           "Pinus sylvestris",
  "MD",  "modřín opadavý",           "Larix decidua",
  "TS",  "tis červený",              "Taxus baccata",
  "DB",  "dub letní",                "Quercus sp.",
  "DBZ", "dub zimní",                "Quercus sp.",
  "BK",  "buk lesní",                "Fagus sylvatica",
  "HB",  "habr obecný",              "Carpinus betulus",
  "JV",  "javor mléč",               "Acer platanoides",
  "KL",  "javor klen",               "Acer pseudoplatanus",
  "BB",  "javor babyka",             "Acer campestre",
  "JS",  "jasan ztepilý",            "Fraxinus excelsior",
  "JL",  "jilm habrolistý",          "Ulmus minor",
  "AK",  "trnovník akát",            "Robinia pseudoacacia",
  "BR",  "bříza bělokorá",           "Betula sp.",
  "JR",  "jeřáb ptačí",              "Sorbus aucuparia",
  "OR",  "ořešák královský",         "Juglans regia",
  "TR",  "třešeň ptačí",             "Prunus avium",
  "HR",  "hrušeň polnička",          "Pyrus pyraster",
  "JB",  "jabloň lesní",             "Malus sylvestris",
  "LP",  "lípa srdčitá",             "Tilia cordata",
  "LO",  "líska obecná",             "Corylus avellana",
  "OL",  "olše lepkavá",             "Alnus glutinosa",
  "OS",  "topol osika",              "Populus tremula",
  "TP",  "topol bílý",               "Populus alba"
  # extend with more rows from the PDF if other codes turn up
)

layer_lookup <- tribble(
  ~layer, ~layer_label,
  "O", "small",
  "P", "advanced",
  "D", "mature"     # NB: the source slide says "<10 cm DBH" for D,
  # which looks like it should read ">10 cm DBH" -
  # worth double-checking with whoever made the key
)

location_lookup <- tribble(
  ~location, ~location_label,
  "1", "terminal",
  "2", "trunk",
  "3", "trunk base",
  "4", "leaves"
)

# --- parse vzorek into components ---
vzorek_pattern <- "^T(\\d)(\\d+)([OPD])([A-Za-z]+)(\\d)$"

df_parsed <- df_pathogen %>%
  mutate(
    vzorek_match       = str_match(vzorek, vzorek_pattern),
    tablet             = vzorek_match[, 2],
    subplot_vzorek     = vzorek_match[, 3],
    layer              = vzorek_match[, 4],
    tree_species_code  = vzorek_match[, 5],
    location           = vzorek_match[, 6]#,
   # parse_ok           = !is.na(vzorek_match[, 1]) | vzorek == "Blank"
  ) %>%
  select(-vzorek_match) %>%
  left_join(species_lookup, by = c("tree_species_code" = "code")) %>%
  left_join(layer_lookup, by = "layer") %>%
  left_join(location_lookup, by = "location") |> 
  select(-Subplot, - Species,  - czech_name, -location,
         -subplot_number) 


nrow(df_parsed)
head(df_parsed)
head(df_blank)

df_parsed_merged <- df_parsed |> 
  full_join(df_blank, by = c("subplot_vzorek" = "subplot_number"))
  
# i have missing data - Michal need to check it up, and we will meet on wednesday to plan the paper ahead.



# get summary tables -------------------------------------------------------------

# --- pathogen occurrence by species (long) ---
table(df_clean$patogen, df_clean$tree_species_code) 

# --- pathogen occurrence by regeneration layer (long) ---
table(df_clean$patogen, df_clean$location_label)


# --- overall total per pathogen, for context / ordering ---
table(df_clean$patogen) 



# attach plot / subplot / species context to every pathogen record
patho_linked <- df_clean %>%
  inner_join(df_trees_samples, by = "vzorek")

# the ones that still don't match - same handful needing a manual fix at the source
anti_join(df_clean, df_trees_samples, by = "vzorek") %>% distinct(vzorek, patogen)

head(patho_linked)



plot_df <- df_clean %>%
  count(patogen, tree_species_code, name = "n") %>%
  filter(n > 0) %>%                                  # drop empty cells - the whole point of a dot matrix
  group_by(patogen) %>%
  mutate(patogen_total = sum(n)) %>%
  ungroup() %>%
  mutate(patogen = fct_reorder(patogen, patogen_total))  # rarest pathogens at bottom, most common at top

ggplot(plot_df, aes(x = tree_species_code, y = patogen)) +
  geom_point(aes(size = n, color = n)) +
  scale_size_continuous(range = c(3, 12), name = "occurrences") +
  scale_color_viridis_c(name = "occurrences", option = "D") +
  labs(x = "Tree species", y = NULL,
       title = "Pathogen / insect occurrence by tree species") +
  theme_classic2(base_size = 12) +
  theme(panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank())

