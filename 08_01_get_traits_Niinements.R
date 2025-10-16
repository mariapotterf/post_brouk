# Get functional traits database ---------------------------------------

# tolerance scales range from 0 (no tolerance) to 5 (maximal tolerance)

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)

#  Input
niin_path   <- "raw/traits_database/Niinemets_2006.xls"  # adjust if needed
dat_overlap <- fread('outData/full_table_overlap_23_25.csv')

# Your species codes from dat_subplot_mng$species
sp_codes <- unique(dat_overlap$species)
sp_codes <- sp_codes[!sp_codes %in% c("", "ots1")]

# Read Niinemets + keep essential columns 
eco_traits_raw <- read_excel(
  path  = niin_path,
  sheet = "Niinemets_2006_appendix",
  skip  = 3,
  .name_repair = function(x) gsub("\\s+", "_", x)
)

# filter out only teh important columns
eco_traits <- eco_traits_raw %>%
  select(Species, Shade_tolerance, Drought_tolerance) %>%
  mutate(Species = str_squish(Species))


# Map my acronyms to latin ones, indicate reason fror merge
# 'species_latin' must either be an exact species in Niinemets OR a genus label
# for which we will compute a genus mean (e.g., "Quercus", "Salix", "Betula").
sp_map <- tribble(
  ~species_code, ~species_latin,            ~notes,
  "piab",        "Picea abies",             "",
  "pisy",        "Pinus sylvestris",        "",
  "lade",        "Larix decidua",           "",
  "absp",        "Abies alba",              "Abies sp. -> proxy A. alba",
  "psme",        "Pseudotsuga menziesii",   "Fix spelling if needed (not 'mensiesii')",
  "taba",        "Taxus baccata",           "",
  "fasy",        "Fagus sylvatica",         "",
  "qusp",        "Quercus",                 "Quercus spp. -> genus mean (Q. robur & Q. petraea)",
  "acca",        "Acer campestre",          "",
  "acpl",        "Acer platanoides",        "",
  "acps",        "Acer pseudoplatanus",     "",
  "algl",        "Alnus glutinosa",         "",
  "alin",        "Alnus incana",            "",
  "alvi",        "Alnus viridis",           "Close to A. alnobetula complex",
  "potr",        "Populus tremula",         "",
  "posp",        "Populus tremula",         "Populus spp. -> use P. tremula proxy (keep consistent)",
  "besp",        "Betula",                  "Betula spp. -> genus mean (B. pendula & B. pubescens)",
  "frex",        "Fraxinus excelsior",      "",
  "tisp",        "Tilia cordata",           "Tilia spp. -> proxy T. cordata",
  "prav",        "Prunus avium",            "",
  "soau",        "Sorbus aucuparia",        "",
  "soto",        "Sorbus torminalis",       "",
  "soar",        "Sorbus aria",             "",
  "casa",        "Castanea sativa",         "",
  "aehi",        "Aesculus hippocastanum",  "",
  "cabe",        "Carpinus betulus",        "",
  "ulsp",        "Ulmus glabra",            "Ulmus spp. -> proxy U. glabra (switch to U. minor if that’s your system)",
  "rops",        "Robinia pseudoacacia",    "",
  "saca",        "Salix caprea",            "",
  "juni",        "Juniperus communis",      "",
  "jure",        "Juglans regia",           "",
  "sasp",        "Salix",                   "Salix spp. -> genus mean (S. caprea & S. alba)",
  "aial",        "Ailanthus altissima",        "Map to A. viridis in Niinemets if needed",
  "osca",        "Ostrya carpinifolia",     "If missing in Niinemets, proxy with Carpinus betulus",
  "fror",        "Fraxinus ornus",          ""
  
)

# get compleyte binomial latin names for my acronyms
df_filter_binomial <- sp_map %>% 
  filter(!species_latin %in% c("Quercus", "Salix", "Betula")) %>%  # filter out only Genus names
  select(-notes)

if (!"Alnus alnobetula" %in% eco_traits$Species && "Alnus viridis" %in% eco_traits$Species) {
  eco_traits <- eco_traits %>%
    add_row(
      Species = "Alnus alnobetula",
      Shade_tolerance   = eco_traits$Shade_tolerance[eco_traits$Species == "Alnus viridis"],
      Drought_tolerance = eco_traits$Drought_tolerance[eco_traits$Species == "Alnus viridis"]
    )
}

# (2) fix Pseudotsuga spelling if your inputs sometimes carry a typo
if (!"Pseudotsuga mensiesii" %in% eco_traits$Species && "Pseudotsuga menziesii" %in% eco_traits$Species) {
  eco_traits <- eco_traits %>%
    add_row(
      Species = "Pseudotsuga mensiesii",
      Shade_tolerance   = eco_traits$Shade_tolerance[eco_traits$Species == "Pseudotsuga menziesii"],
      Drought_tolerance = eco_traits$Drought_tolerance[eco_traits$Species == "Pseudotsuga menziesii"]
    )
}



# filter complete records 
eco_traits_binomial <- eco_traits %>% 
  filter(Species %in% df_filter_binomial$species_latin) %>% 
  rename(species = Species)



#  3) Create genus means for Quercus / Salix / Betula 
# Helper function: mean traits for a set of species -> single row with genus name
genus_mean <- function(species_vec, genus_label) {
  eco_traits %>%
    filter(Species %in% species_vec) %>%
    summarise(
      species = genus_label,
      Shade_tolerance   = mean(Shade_tolerance, na.rm = TRUE),
      Drought_tolerance = mean(Drought_tolerance, na.rm = TRUE)#,
      #n_species_used    = dplyr::n()
    )
}

# Define which species to use per genus (adjust if your region differs)
q_species <- c("Quercus robur", "Quercus petraea")
s_species <- c("Salix caprea", "Salix alba")  #  Salix alba is not present 
b_species <- c("Betula pendula", "Betula pubescens") # if only one is present, mean() still works

traits_quercus <- genus_mean(q_species, "Quercus")
traits_salix   <- genus_mean(s_species, "Salix")
traits_betula  <- genus_mean(b_species, "Betula")

# Bind genus means
traits_full <- bind_rows(
  eco_traits_binomial,
  traits_quercus,
  traits_salix,
  traits_betula
) 


# add back acronyms
traits_full <- traits_full %>% 
  left_join(sp_map, by = c("species" = "species_latin")) %>% 
  select(-species, -notes) %>% 
  rename(species = species_code)

fwrite(traits_full, 'outData/my_species_traits.csv')

library(ggrepel)


traits_full %>%
  ggplot(aes(x = Shade_tolerance,
             y = Drought_tolerance,
             label = species,
             color = species)) +         # <--- map color here
   geom_hline(yintercept = 3, linetype = "dashed", color = "grey90") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey90") +
geom_point() +
  geom_text_repel(size = 3, show.legend = T) +   # hide duplicated legend
   theme(legend.position = "none")


print(sort(traits_full$Shade_tolerance))

traits_full %>% 
  arrange(Shade_tolerance) %>% 
  print(n = 40)

