# ==============================================================================
# 01_build_subplot_species_table.R
# Damage data (QField survey 2025) -> one row per SPECIES x SUBPLOT (x vegtype)
#
# Decisions agreed so far
#   * unit of analysis  : species within subplot (plot kept as grouping variable)
#   * damage            : kept BY TYPE (terminal / foliage / stem), counts of individuals
#                         "any damage" only as bounds (lower / upper), no overlap assumed
#   * management        : planting & anti_browsing  NA -> 0   (subplot level)
#   * protection type   : fence / repellent / coating / other, parsed from `comments`
#   * dmg_jedinec       : 1 = damage on tallest individual, 2 = on a lower individual
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# --- paths (edit) -------------------------------------------------------------
path_csv <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else
  "inData_Michal/poskoz_jedinci_clean.csv"
path_out <- "outDataMichal/subplot_species_damage.csv"

raw <- fread(path_csv, na = c("", "NA"), encoding = "UTF-8") %>% as_tibble()

# --- 1. basic cleaning ---------------------------------------------------------
d <- raw %>%
  filter(species != "ots1") %>%                         # unidentified species
  mutate(
    planting      = replace_na(planting, 0),            # NA ("not sure"/none) -> 0
    anti_browsing = replace_na(anti_browsing, 0),
    comments_lc   = str_to_lower(replace_na(comments, "")),
    tablet        = str_extract(subplot, "T\\d"),
    # the date in the subplot ID (20250717 / 20250827) is the export batch, NOT the survey day.
    # `time` is the tablet clock: 42 subplots of tablet T1 show 2025-02-01 (clock error)
    # -> set to NA instead of using a wrong date (needed for the frost-damage argument)
    date          = as.Date(time),
    date          = if_else(format(date, "%m") == "02", as.Date(NA), date)
  )

# --- 2. protection type from comments (subplot level) ---------------------------
# fence    : oplocenka / oplocení / oploceka / opolocenka / "plot" ...
# repellent: repelent ; coating: nátěr ; wrap: "v obalu" (individual wrap)
# failed   : fallen / removed / former fence  (anti_browsing is already 0 there)
# partial  : only part of the subplot or only some species protected
subplot_info <- d %>%
  distinct(plot, subplot, planting, anti_browsing, comments, comments_lc,
           stump, clear, grndwrk, logging_trail, windthrow, standing_deadwood,
           date, tablet, x, y) %>%
  mutate(
    prot_fence      = anti_browsing == 1 & str_detect(comments_lc, "oplo|opolo|^plot$"),
    prot_repellent  = anti_browsing == 1 & str_detect(comments_lc, "repelent"),
    prot_coating    = anti_browsing == 1 & str_detect(comments_lc, "n[aá]t[eě]r"),
    prot_wrap       = anti_browsing == 1 & str_detect(comments_lc, "obal"),
    prot_fence_failed = anti_browsing == 0 &                   # fallen/removed fence
      str_detect(comments_lc, "oplo|opolo"),
    prot_partial    = anti_browsing == 1 &
      str_detect(comments_lc, "část|cast|1/2|50 ?%|polovin|pouze"),
    prot_n_types    = prot_fence + prot_repellent + prot_coating + prot_wrap,
    prot_type = case_when(
      anti_browsing == 0                 ~ "none",
      prot_n_types > 1                   ~ "mixed",
      prot_fence                         ~ "fence",
      prot_repellent                     ~ "repellent",
      prot_coating                       ~ "coating",
      prot_wrap                          ~ "wrap",
      TRUE                               ~ "unknown"      # protected, type not stated
    ),
    # fence = collective protection (all species); the rest = individual protection
    prot_class = case_when(
      prot_type == "none"  ~ "none",
      prot_type == "fence" ~ "fence",
      prot_type == "unknown" ~ "unknown",
      TRUE                 ~ "individual"
    )
  ) %>%
  select(-comments_lc)

stopifnot(!anyDuplicated(subplot_info$subplot))   # management must be unique per subplot

# --- 3. aggregate rows to species x subplot x vegtype --------------------------
# (some species have several rows within a subplot, e.g. different height classes)
sp <- d %>%
  group_by(plot, subplot, species, vegtype) %>%
  summarise(
    n_rows        = n(),
    n_censored    = any(n == 17),                        # ">16" code, true n >= 17
    n             = sum(n),                              # (must come after n_censored)
    T_raw         = sum(tot_dmg_terminal),               # individuals, terminal damage
    F_raw         = sum(tot_dmg_foliage),                # individuals, foliage damage
    S             = max(dmg_stem_numeric),               # stem: presence only (0/1)
    dmg_on_lower  = any(dmg_jedinec == 2, na.rm = TRUE), # damage on a lower individual
    T_missing_rec = any(dmg_terminal == "chybajuci", na.rm = TRUE),
    T_dead_rec    = any(dmg_terminal == "odumrety",  na.rm = TRUE),
    stem_cause    = first(na.omit(dmg_stem_cause)),
    hgt_tallest   = suppressWarnings(max(hgt_est, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    hgt_tallest   = ifelse(is.finite(hgt_tallest), hgt_tallest, NA_real_),
    stem_cause    = replace_na(stem_cause, "none"),
    stem_game     = stem_cause == "zver",
    stem_mech     = stem_cause == "mechanizace",
    cnt_gt_n      = T_raw > n | F_raw > n,               # counts exceeding n -> check
    T             = pmin(T_raw, n),                      # capped at n
    F             = pmin(F_raw, n),
    # "any damage": bounds only (overlap between types is not identifiable)
    any_lo        = pmin(n, pmax(T, F, S)),
    any_hi        = pmin(n, T + F + S),
    any_ambiguous = any_lo != any_hi
  ) %>%
  left_join(subplot_info, by = c("plot", "subplot"))

# --- 4. checks (compare with the manuscript) -----------------------------------
cat("\n--- CHECKS ---------------------------------------------------------\n")
cat("records (species x subplot x vegtype):", nrow(sp), "\n")
cat("plots:", n_distinct(sp$plot), "| subplots:", n_distinct(sp$subplot), "\n")
cat("individuals (sum n)       :", sum(sp$n), "   [manuscript: 2986]\n")
cat("terminal damaged (T)      :", sum(sp$T), "   [172]\n")
cat("foliage damaged (F, capped):", sum(sp$F), "   [501 uncapped:", sum(sp$F_raw), "]\n")
cat("stem records (S)          :", sum(sp$S), "   [129]\n")
cat("any damage, upper bound   :", sum(sp$any_hi), "   [784]\n")
cat("any damage, lower bound   :", sum(sp$any_lo), "\n")
cat("records with ambiguous 'any':", sum(sp$any_ambiguous), "\n")
cat("records with counts > n   :", sum(sp$cnt_gt_n), "\n")
cat("records with n = 17 (>16) :", sum(sp$n_censored), "\n")
cat("records with damage on a lower individual:", sum(sp$dmg_on_lower), "\n")

cat("\nSubplots with regeneration per plot (missing subplots = no regeneration or not surveyed):\n")
print(table(n_subplots = sp %>% distinct(plot, subplot) %>% count(plot) %>% pull(n)))

cat("\nplanting x anti_browsing (subplots):\n")
print(table(planting = subplot_info$planting, anti_browsing = subplot_info$anti_browsing))
cat("\nprotection type (subplots):\n")
print(table(subplot_info$prot_type))
cat("\nfailed fences coded 0:", sum(subplot_info$prot_fence_failed),
    "| partial protection:", sum(subplot_info$prot_partial), "\n")

cat("\nProtected subplots with 'unknown' type (check comments):\n")
print(subplot_info %>% filter(prot_type == "unknown") %>% count(comments, sort = TRUE) %>% head(10))
cat("\nAll non-standard comments classified as 'mixed' or 'wrap':\n")
print(subplot_info %>% filter(prot_type %in% c("mixed", "wrap")) %>% distinct(comments, prot_type))

# --- 5. save ------------------------------------------------------------------
dir.create(dirname(path_out), showWarnings = FALSE, recursive = TRUE)
fwrite(sp, path_out)
cat("\nSaved:", path_out, "\n")
