
# Get ecological traits per species from Ellenberg/Niinements
# have not been run yet! 

# should work for identification of teh early vs late seral species 
# based on Shade tolerance


# packages
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# --- your species list as a tibble (paste your CSV if stored elsewhere) ----
species_csv <- 'species_code,species_latin,seral_class,primary_basis,notes
piab,Picea abies,late,Niinemets/EUFORGEN,shade-tolerant conifer
pisy,Pinus sylvestris,early,Niinemets/EUFORGEN,light-demanding pioneer on poor/dry sites
lade,Larix decidua,early,Niinemets/EUFORGEN,very light-demanding deciduous conifer
absp,Abies alba,late,Niinemets/EUFORGEN,very shade-tolerant
psme,Pseudotsuga menziesii,intermediate,EUFORGEN/plantation trials,fast juvenile but needs light; not truly tolerant
taba,Taxus baccata,late,Niinemets/EUFORGEN,very shade-tolerant understory conifer
fasy,Fagus sylvatica,late,Niinemets/Ellenberg,very shade-tolerant broadleaf
acca,Acer campestre,intermediate,Niinemets/Ellenberg,moderately shade-tolerant
acpl,Acer platanoides,intermediate,Niinemets/Ellenberg,intermediate–tolerant
acps,Acer pseudoplatanus,intermediate,Niinemets/EUFORGEN,gap-recruiting; moderate tolerance
algl,Alnus glutinosa,early,Niinemets/EUFORGEN,N-fixing pioneer on wet soils
alin,Alnus incana,early,Niinemets/EUFORGEN,pioneer on riparian/gravels
alvi,Alnus viridis,early,EUFORGEN,montane/subalpine pioneer shrub
potr,Populus tremula,early,Niinemets/EUFORGEN,fast pioneer; strong suckering
posp,Populus tremula,early,Niinemets/EUFORGEN,alias of Populus tremula (same class)
frex,Fraxinus excelsior,intermediate,Niinemets/EUFORGEN,light-demanding seedling but tolerant in shade as sapling
tisp,Tilia cordata,late,Niinemets/EUFORGEN,shade-tolerant
prav,Prunus avium,early,EUFORGEN,light-demanding gap/edge species
soau,Sorbus aucuparia,intermediate,EUFORGEN,gap species; survives some shade
soto,Sorbus torminalis,early,EUFORGEN,light-demanding thermophilous
soar,Sorbus aria,early,EUFORGEN,light-demanding on dry/rocky sites
casa,Castanea sativa,intermediate,EUFORGEN,needs light/heat; not very tolerant
aehi,Aesculus hippocastanum,intermediate,EUFORGEN,moderate tolerance when young
cabe,Carpinus betulus,late,Niinemets/EUFORGEN,shade-tolerant
ulsp,Ulmus glabra,intermediate,EUFORGEN,intermediate–tolerant mesic broadleaf
rops,Robinia pseudoacacia,early,EUFORGEN/BIOLFLOR,very fast pioneer; N-fixer
saca,Salix caprea,early,EUFORGEN,classic pioneer willow
juni,Juniperus communis,early,EUFORGEN,light-demanding (slow but not tolerant)
jure,Juglans regia,intermediate,EUFORGEN,light-demanding tree; intermediate overall
aial,Alnus alnobetula,early,EUFORGEN,subalpine pioneer shrub
osca,Ostrya carpinifolia,late,EUFORGEN,shade-tolerant thermophilous
fror,Fraxinus ornus,intermediate,EUFORGEN,thermophilous; intermediate/light-demanding
,NA,NA,NA,missing species name/code
ots1,NA,NA,NA,unmatched placeholder
'
sp <- read_csv(I(species_csv), show_col_types = FALSE) %>%
  mutate(species_binom = str_squish(str_replace_all(species_latin, "\\.", "")) %>%
           str_replace("^([A-Z][a-z]+)\\s+(\\S+).*", "\\1 \\2"))

# --- paths to your files ----------------------------------------------------
path_niin <- "Niinemets_2006.xls"  # Niinemets & Valladares 2006 (XLS)
path_eive <- "Indicator_values_Tichy_et_al 2022-11-29 (1).xlsx"  # EIVE (XLSX)

# --- read Niinemets: find sheet & columns with shadow/shade tolerance -------
# readxl handles .xls too
niin_sheetnames <- excel_sheets(path_niin)

read_niin <- function(sheet) {
  df <- read_excel(path_niin, sheet = sheet)
  nm <- tolower(names(df))
  df %>%
    setNames(nm) %>%
    mutate(across(everything(), ~.x)) %>%
    # guess species & shade columns
    { 
      spec_col <- names(.)[str_detect(names(.), "species|latin|name|taxon")][1]
      shade_col <- names(.)[str_detect(names(.), "shade")][1]
      if (is.na(spec_col) || is.na(shade_col)) return(tibble())
      select(., species_raw = all_of(spec_col), shade_raw = all_of(shade_col))
    } %>%
    mutate(species_binom = species_raw %>%
             as.character() %>%
             str_squish() %>%
             str_replace("^([A-Z][a-z]+)\\s+(\\S+).*", "\\1 \\2"),
           niin_shade = readr::parse_number(as.character(shade_raw)))
}

niin_all <- map_dfr(niin_sheetnames, read_niin) %>%
  distinct(species_binom, .keep_all = TRUE) %>%
  filter(!is.na(species_binom))

# --- read EIVE/Tichý: get Ellenberg Light L ---------------------------------
eive_sheetnames <- excel_sheets(path_eive)

read_eive <- function(sheet) {
  df <- read_excel(path_eive, sheet = sheet)
  nm <- tolower(names(df))
  df <- setNames(df, nm)
  spec_col <- names(df)[str_detect(names(df), "species|scientific|latin|taxon|name")][1]
  light_col <- names(df)[names(df) %in% c("light","ellenberg light","l","eiv l","eiv_l","eiv_light")][1]
  if (is.na(spec_col) || is.na(light_col)) return(tibble())
  df %>%
    transmute(species_raw = .data[[spec_col]],
              ellenberg_l = suppressWarnings(as.numeric(.data[[light_col]]))) %>%
    mutate(species_binom = species_raw %>%
             as.character() %>%
             str_squish() %>%
             str_replace("^([A-Z][a-z]+)\\s+(\\S+).*", "\\1 \\2"))
}

eive_all <- map_dfr(eive_sheetnames, read_eive) %>%
  filter(!is.na(species_binom)) %>%
  distinct(species_binom, .keep_all = TRUE)

# --- join everything for your species ---------------------------------------
out <- sp %>%
  left_join(select(niin_all, species_binom, niin_shade), by = "species_binom") %>%
  left_join(select(eive_all, species_binom, ellenberg_l), by = "species_binom") %>%
  mutate(
    niin_class = case_when(
      !is.na(niin_shade) & niin_shade <= 2 ~ "early/light-demanding",
      !is.na(niin_shade) & niin_shade >= 4 ~ "late/shade-tolerant",
      !is.na(niin_shade)                  ~ "intermediate",
      TRUE ~ NA_character_
    ),
    # Ellenberg: high L = light-demanding (early)
    eive_class = case_when(
      !is.na(ellenberg_l) & ellenberg_l >= 7 ~ "early/light-demanding",
      !is.na(ellenberg_l) & ellenberg_l <= 4 ~ "late/shade-tolerant",
      !is.na(ellenberg_l)                    ~ "intermediate",
      TRUE ~ NA_character_
    )
  ) %>%
  select(species_code, species_latin, species_binom,
         niin_shade, niin_class, ellenberg_l, eive_class,
         seral_class, primary_basis, notes)

# write out (so you can inspect/share)
write_csv(out, "shade_tolerance_lookup_joined.csv")

# quick peek
print(out, n = nrow(out))
