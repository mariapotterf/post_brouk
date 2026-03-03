

# hold on all variables regarding species classifications



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

# species_colors
# piab      besp      pisy      qusp      fasy      lade      saca      soau      acps      potr      absp      sasp 
#"#006837" "#17934D" "#58B65F" "#94D168" "#C6E77F" "#EDF7A7" "#FEF0A7" "#FDCD7B" "#FA9C58" "#EE613D" "#D22B26" "#A50026" 





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

