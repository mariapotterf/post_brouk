# ==============================================================================
# 02_descriptives.R
# Descriptive overview of damage by species and management (no models yet)
#
# Input : outData/subplot_species_damage.csv   (from 01_build_subplot_species_table.R)
# Unit  : species x subplot (records); regeneration = small + advanced
#         terminal damage -> small class only (protocol records it only there)
# Output: outData/descriptives/*.csv  and  *.png
#
# NOTE: CIs are pooled Wilson intervals - they IGNORE clustering (species within
#       subplot within plot), so they are too narrow. They are for orientation only;
#       the mixed models in step 3 handle the clustering.
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

path_tab <- "outData/subplot_species_damage.csv"
out_dir  <- "outData/descriptives"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sp <- fread(path_tab) %>% as_tibble() %>%
  filter(vegtype %in% c("small", "advanced")) %>%            # regeneration only
  mutate(
    prot_class = factor(prot_class, levels = c("none", "fence", "individual", "unknown")),
    protected  = factor(anti_browsing, levels = c(0, 1), labels = c("not protected", "protected")),
    planted    = factor(planting,      levels = c(0, 1), labels = c("not planted", "planted")),
    hgt_class  = cut(hgt_tallest, breaks = c(0, 0.4, 0.8, 1.3, 2.1),
                     labels = c("0.2-0.4 m", "0.4-0.8 m", "0.8-1.3 m", "1.3-2.0 m"))
  )

# species with >= 100 individuals shown separately, the rest pooled
sp_order <- sp %>% group_by(species) %>% summarise(N = sum(n)) %>% arrange(desc(N))
major    <- sp_order$species[sp_order$N >= 100]
sp <- sp %>% mutate(species_g = factor(if_else(species %in% major, species, "other"),
                                       levels = c(major, "other")))

# --- long format: one row per record x outcome ----------------------------------
mk <- function(d, outcome, x, den) {
  tibble(plot = d$plot, subplot = d$subplot, species_g = d$species_g,
         vegtype = d$vegtype, protected = d$protected, planted = d$planted,
         prot_class = d$prot_class, hgt_class = d$hgt_class,
         outcome = outcome, num = x, den = den)
}
small <- filter(sp, vegtype == "small")
long <- bind_rows(
  mk(small, "terminal",     small$T,                         small$n),   # individuals
  mk(small, "term_missing", as.integer(small$T_missing_rec), 1L),        # records
  mk(sp,    "foliage",      sp$F,                            sp$n),      # individuals
  mk(sp,    "stem",         sp$S,                            1L),        # records
  mk(sp,    "stem_game",    as.integer(sp$stem_game),        1L)         # records
) %>%
  mutate(outcome = factor(outcome,
    levels = c("terminal", "term_missing", "foliage", "stem", "stem_game"),
    labels = c("Terminal damage\n(% individuals, small)",
               "Missing terminal\n(% records, small)",
               "Foliage damage\n(% individuals)  [negative control]",
               "Stem damage\n(% records)",
               "Stem damage by game\n(% records)")))

# --- Wilson interval ------------------------------------------------------------
wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  tibble(lwr = 100 * (ctr - hw), upr = 100 * (ctr + hw))
}
rates <- function(d, ...) {
  d %>% group_by(outcome, ...) %>%
    summarise(records = n(), num = sum(num), den = sum(den), .groups = "drop") %>%
    mutate(pct = 100 * num / den, wilson(num, den))
}

# ==============================================================================
# A. Management balance (subplot level and species level)
# ==============================================================================
sub <- sp %>% distinct(plot, subplot, planted, protected, prot_class)
cat("\n=== A1. Subplots: planting x protection class ===\n")
print(table(sub$planted, sub$prot_class))

bal <- sp %>% group_by(species_g, planted, prot_class) %>%
  summarise(records = n(), individuals = sum(n), .groups = "drop")
write.csv(bal, file.path(out_dir, "A_balance_species_planting_protection.csv"), row.names = FALSE)
cat("\n=== A2. Records per species x protection class ===\n")
print(sp %>% count(species_g, prot_class) %>% pivot_wider(names_from = prot_class,
      values_from = n, values_fill = 0))
cat("\n=== A3. Records per species x planting x protected ===\n")
print(sp %>% count(species_g, planted, protected) %>%
      unite("grp", planted, protected, sep = " / ") %>%
      pivot_wider(names_from = grp, values_from = n, values_fill = 0))

# ==============================================================================
# B. Rates
# ==============================================================================
r_overall <- rates(long)
r_prot    <- rates(long, protected)
r_pp      <- rates(long, planted, protected)
r_class   <- rates(long, prot_class)
r_sp_prot <- rates(long, species_g, protected)
r_sp_cls  <- rates(long, species_g, prot_class)
r_hgt     <- rates(filter(long, outcome == levels(long$outcome)[1]), hgt_class, protected)
r_sp      <- rates(long, species_g)

for (nm in c("overall", "prot", "pp", "class", "sp_prot", "sp_cls", "hgt", "sp")) {
  write.csv(get(paste0("r_", nm)), file.path(out_dir, paste0("B_rates_", nm, ".csv")),
            row.names = FALSE)
}

cat("\n=== B1. Pooled rates, all species ===\n")
print(r_overall %>% mutate(outcome = gsub("\n", " ", outcome)) %>%
      select(outcome, records, num, den, pct, lwr, upr), width = 140)
cat("\n=== B2. By protection (pooled over species and planting) ===\n")
print(r_prot %>% mutate(outcome = gsub("\n", " ", outcome)) %>%
      select(outcome, protected, records, den, pct, lwr, upr), width = 140)
cat("\n=== B3. By planting x protection ===\n")
print(r_pp %>% mutate(outcome = gsub("\n", " ", outcome)) %>%
      select(outcome, planted, protected, records, den, pct, lwr, upr), n = 40, width = 140)
cat("\n=== B4. By protection class ===\n")
print(r_class %>% mutate(outcome = gsub("\n", " ", outcome)) %>%
      select(outcome, prot_class, records, den, pct, lwr, upr), n = 40, width = 140)

# ==============================================================================
# C. Figures
# ==============================================================================
col_prot <- c("not protected" = "#8f8e88", "protected" = "#2a78d6")
col_cls  <- c(none = "#8f8e88", fence = "#2a78d6", individual = "#eb6834", unknown = "#1baf7a")
th <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom", panel.grid.major.x = element_line(colour = "grey92"),
        axis.line = element_line(colour = "grey60"), axis.ticks = element_line(colour = "grey60"))
dodge <- position_dodge(width = 0.6)
cap <- "Pooled Wilson 95% CI (ignores plot/subplot clustering); hollow = fewer than 20 denominators."

# C1: planting x protection, all outcomes
p1 <- r_pp %>%
  mutate(fillcol = ifelse(den < 20, "white", unname(col_prot[as.character(protected)]))) %>%
  mutate(grp = factor(paste(planted, "/", protected),
         levels = rev(c("not planted / not protected", "not planted / protected",
                        "planted / not protected", "planted / protected")))) %>%
  ggplot(aes(pct, grp, colour = protected, shape = planted)) +
  geom_linerange(aes(xmin = lwr, xmax = upr), linewidth = 0.6) +
  geom_point(aes(fill = fillcol), size = 2.8, stroke = 0.9) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = col_prot, name = NULL) + scale_fill_identity() +
  scale_shape_manual(values = c("not planted" = 21, "planted" = 24), name = NULL) +
  labs(x = "% damaged", y = NULL, title = "Damage by planting and browsing protection",
       subtitle = "all species pooled; foliage = negative control (protection should not matter)",
       caption = cap) + th
ggsave(file.path(out_dir, "C1_planting_x_protection.png"), p1, width = 7, height = 10, dpi = 200)

# C2: protection class (fence vs individual vs unknown)
p2 <- r_class %>%
  mutate(fillcol = ifelse(den < 20, "white", unname(col_cls[as.character(prot_class)]))) %>%
  mutate(prot_class = factor(prot_class, levels = rev(levels(prot_class)))) %>%
  ggplot(aes(pct, prot_class, colour = prot_class)) +
  geom_linerange(aes(xmin = lwr, xmax = upr), linewidth = 0.6) +
  geom_point(aes(fill = fillcol), shape = 21, size = 2.8, stroke = 0.9) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = col_cls, guide = "none") + scale_fill_identity() +
  labs(x = "% damaged", y = NULL, title = "Damage by type of protection",
       subtitle = "fence = collective; individual = repellent / coating / wrap; unknown = no comment",
       caption = cap) + th
ggsave(file.path(out_dir, "C2_protection_class.png"), p2, width = 7, height = 9, dpi = 200)

# C3: species x protection (terminal, missing terminal, foliage, stem)
p3 <- r_sp_prot %>%
  filter(outcome %in% levels(long$outcome)[1:4]) %>%
  mutate(fillcol = ifelse(den < 20, "white", unname(col_prot[as.character(protected)]))) %>%
  mutate(species_g = factor(species_g, levels = rev(levels(species_g)))) %>%
  ggplot(aes(pct, species_g, colour = protected)) +
  geom_linerange(aes(xmin = lwr, xmax = upr), position = dodge, linewidth = 0.5) +
  geom_point(aes(fill = fillcol), shape = 21,
             position = dodge, size = 2.2, stroke = 0.8) +
  facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = col_prot, name = NULL) + scale_fill_identity() +
  labs(x = "% damaged", y = NULL, title = "Damage by species and browsing protection",
       subtitle = "planting not separated here (see C1); species ordered by abundance", caption = cap) + th
ggsave(file.path(out_dir, "C3_species_x_protection.png"), p3, width = 13, height = 5.5, dpi = 200)

# C4: terminal damage by height class
p4 <- r_hgt %>% filter(!is.na(hgt_class)) %>%
  mutate(fillcol = ifelse(den < 20, "white", unname(col_prot[as.character(protected)]))) %>%
  ggplot(aes(hgt_class, pct, colour = protected, group = protected)) +
  geom_linerange(aes(ymin = lwr, ymax = upr), position = position_dodge(width = 0.4), linewidth = 0.6) +
  geom_line(position = position_dodge(width = 0.4), linewidth = 0.5) +
  geom_point(aes(fill = fillcol), shape = 21,
             position = position_dodge(width = 0.4), size = 2.8, stroke = 0.9) +
  scale_colour_manual(values = col_prot, name = NULL) + scale_fill_identity() +
  labs(x = "height class of the tallest individual in the record (small class)",
       y = "% individuals with terminal damage",
       title = "Terminal damage by height class and protection", caption = cap) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "C4_terminal_by_height.png"), p4, width = 7, height = 4.5, dpi = 200)

cat("\nSaved tables and figures to:", out_dir, "\n")
