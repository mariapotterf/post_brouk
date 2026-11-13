

# Read data from Miso - poskozeni
# identify the tallest stem per species
# check wheather the taller/shorter individual has damaged terminal, or more of them
# decide whatther we work on plot or siubplot level
# make summary boxplots

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)


df <- fread("inData_Michal/poskoz_jedinci_clean.csv",
               na = c("", "NA"))

View(df)

head(df)

table(df$dmg_stem_cause)

table(df$dmg_foliage )



# test dummy approaches START !!!

library(dplyr)

# dummy raw records. n = individuals, T/F = individuals with terminal/foliage damage,
# S = stem damage (0/1), jedinec = dmg_jedinec (1 = tallest, 2 = only a lower one, NA = none)
toy <- tribble(
  ~plot, ~subplot, ~species, ~vegtype,   ~n, ~T, ~F, ~S, ~jed,
  1,       11, "piab",   "small",     5,  0,  5,  0,    1,
  1,       11, "fasy",   "small",     3,  1,  2,  0,    1,
  1,       11, "betu",   "small",    10,  0,  0,  0,   NA,
  1,       12, "piab",   "small",     1,  1,  0,  0,    1,
  1,       12, "fasy",   "small",     4,  1,  0,  0,    1,
  2,       21, "piab",   "small",     6,  1,  0,  0,    2,   # only a lower one damaged
  2,       21, "betu",   "small",     8,  0,  3,  0,    1,
  2,       22, "fasy",   "small",     2,  0,  0,  1,    1,
  2,       22, "piab",   "small",     7,  0,  0,  0,   NA,
  2,       23, "piab",   "advanced",  3,  0,  1,  0,    1,   # terminal not measured in advanced
  2,       23, "betu",   "advanced",  2,  0,  0,  0,   NA
)

toy
# VERSION 1: by damage type (share of individuals)
v1 <- toy %>%
  mutate(term_dmg = ifelse(vegtype == "small", T, NA),          # terminal: small only (NA, not 0)
         term_n   = ifelse(vegtype == "small", n, NA),
         term_perc = round(100 * term_dmg / term_n),
         fol_dmg  = F,                                          # foliage: small + advanced
         fol_perc = round(100 * F / n),
         stem_pres = S) %>%                                     # stem: presence only
  select(plot, subplot, species, vegtype, n, term_dmg, term_perc, fol_dmg, fol_perc, stem_pres)

v1
# VERSION 2: damage on the tallest individual (cumulative, any type)
v2 <- toy %>%
  mutate(tallest_any = as.integer(!is.na(jed) & jed == 1),
         lower_only  = as.integer(!is.na(jed) & jed == 2),
         n_types     = (T > 0) + (F > 0) + (S > 0),
         type_known  = n_types <= 1) %>%                        # FALSE = type per tree ambiguous
  select(plot, subplot, species, vegtype, n, tallest_any, lower_only, n_types, type_known)
v2
# !!! END










# # get total sum of indivudual per plots
# get total sum of damaged terminals per plot
# quantify the anti_browing intensity
df_small <- df %>% 
  filter(vegtype == "small", species != "ots1")   # for foliage: vegtype %in% c("small","advanced")

df_small %>%
  mutate(dmg     = tot_dmg_terminal > 0,
         tallest = dmg & dmg_jedinec == 1,
         lower   = dmg & dmg_jedinec == 2) %>%
  summarise(records         = n(),
            pct_individuals = 100 * sum(tot_dmg_terminal) / sum(n),
            pct_tallest     = 100 * mean(tallest),
            pct_lower_only  = 100 * mean(lower),
            pct_any         = 100 * mean(dmg))                       # 3. any individual in the record

# 1. total individuals per plot
tot_ind <- df_small %>%
  group_by(plot) %>%
  summarise(n_ind = sum(n), .groups = "drop")

# 2. damaged terminals per plot (tallest individual only, terminal damage only)
tot_dmg <- df_small %>%
  filter(dmg_type == "terminal") %>%
  group_by(plot) %>%
  summarise(n_dmg_term = sum(tot_dmg_terminal), .groups = "drop")

# 3. anti-browsing intensity = share of subplots with protection
df_protected <- df %>%
  distinct(plot, subplot, anti_browsing) %>%
  group_by(plot) %>%
  summarise(n_sub_prot = sum(anti_browsing), .groups = "drop") %>%
  mutate(n_sub = 5,
         protection_intensity = n_sub_prot / n_sub) %>%
  select(plot, n_sub, n_sub_prot, protection_intensity)


# 4. put together + %
plot_lvl <- tot_ind %>%
  left_join(tot_dmg, by = "plot") %>%
  left_join(df_protected,    by = "plot") %>%
  mutate(n_dmg_term = tidyr::replace_na(n_dmg_term, 0),
         pct_dmg    = 100 * n_dmg_term / n_ind)

plot_lvl

ggplot(plot_lvl, 
       aes(x = factor(protection_intensity), 
           y = pct_dmg)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", width = 0.6) +
  geom_jitter(aes(size = n_ind), width = 0.15, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
  scale_size_continuous(range = c(0.8, 4), name = "individuals\nper plot") +
  labs(x = "anti-browsing intensity (share of subplots protected)",
       y = "% individuals with terminal damage",
       title = "Terminal damage vs. browsing protection intensity",
       subtitle = "one point = one plot; red diamond = mean") +
  theme_classic(base_size = 12)



# test on species levels -----------------------------------------

target_sp <- c("SM", "BK", "BR", "DB", "BO")

# 1. individuals per plot x species
tot_ind_sp <- df_small %>%
  filter(species %in% target_sp) %>%
  group_by(plot, species) %>%
  summarise(n_ind = sum(n), .groups = "drop")

# 2. damaged terminals per plot x species
tot_dmg_sp <- df_small %>%
  filter(species %in% target_sp, dmg_type == "terminal") %>%
  group_by(plot, species) %>%
  summarise(n_dmg_term = sum(tot_dmg_terminal), .groups = "drop")

# 3. protection intensity is a PLOT property - unchanged, joined on plot only
# (df_protected from before)

# 4. combine
plot_sp <- tot_ind_sp %>%
  left_join(tot_dmg_sp,   by = c("plot", "species")) %>%
  left_join(df_protected, by = "plot") %>%
  mutate(n_dmg_term = tidyr::replace_na(n_dmg_term, 0),
         pct_dmg    = 100 * n_dmg_term / n_ind,
         species    = factor(species, levels = target_sp))

# how many plots per facet x protection level - check before trusting the boxes
plot_sp %>% count(species, protection_intensity) %>% tidyr::pivot_wider(
  names_from = protection_intensity, values_from = n, values_fill = 0)


# filter only if more individuals are present
plot_sp %>%
  filter(n_ind >= 3) %>%
ggplot(aes(factor(protection_intensity), pct_dmg)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", width = 0.6) +
  geom_jitter(aes(size = n_ind), width = 0.15, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "red") +
  facet_wrap(~ species, nrow = 2) +
  scale_size_continuous(range = c(0.8, 4), name = "individuals\nper plot") +
  ggtitle("Only sites > 3 individuals/species") +
  labs(x = "anti-browsing intensity (share of subplots protected)",
       y = "% individuals with terminal damage",
       title = "Terminal damage vs. protection intensity, by species, > 3",
       subtitle = "one point = one plot x species; red diamond = mean") +
  theme_classic(base_size = 12)



# counts number of individuals per protection level
plot_sp %>%
  group_by(species, protection_intensity) %>%
  summarise(plots = n(), n_ind = sum(n_ind), n_dmg = sum(n_dmg_term),
            pct = 100 * n_dmg / n_ind,
            lwr = 100 * binom.test(n_dmg, n_ind)$conf.int[1],
            upr = 100 * binom.test(n_dmg, n_ind)$conf.int[2],
            .groups = "drop") %>%
  ggplot(aes(x = factor(protection_intensity), pct)) +
 # geom_boxplot() + 
  geom_pointrange(aes(ymin = lwr, ymax = upr, size = n_ind)) +
 # scale_size_continuous(range = c(0.3, 1), guide = "none") +
 
  facet_wrap(~ species) +
  labs(x = "anti-browsing intensity", y = "% terminal damage (95% CI)") +
  theme_classic(base_size = 12)



# binomial glmm 

library(lme4)

m <- glmer(cbind(n_dmg_term, n_ind - n_dmg_term) ~ protection_intensity * species + (1 | plot),
           data = plot_sp, family = binomial)
summary(m)

# 1. is it a false positive? recompute gradient more carefully
relgrad <- with(m@optinfo$derivs, solve(Hessian, gradient))
max(abs(relgrad))   # if < 0.001, the warning is spurious

# 2. try a different optimizer
m2 <- update(m, control = glmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 2e5)))

# 3. does the interaction earn its keep?
m_add <- glmer(cbind(n_dmg_term, n_ind - n_dmg_term) ~ protection_intensity + species + (1 | plot),
               data = plot_sp, family = binomial,
               control = glmerControl(optimizer = "bobyqa"))
anova(m_add, m)


summary(m_add)

# odds ratios with CIs
exp(cbind(OR = fixef(m_add), confint(m_add, method = "Wald", parm = "beta_")))

# marginal means on the response scale
library(emmeans)
emmeans(m_add, ~ species, type = "response")
emmeans(m_add, ~ protection_intensity,
        at = list(protection_intensity = c(0, 0.5, 1)), type = "response")


library(DHARMa)
simulateResiduals(m_add, plot = TRUE)



emm_sp <- emmeans(m_add, ~ species, type = "response") %>% as.data.frame()

ggplot(emm_sp, aes(reorder(species, prob), prob * 100)) +
  geom_pointrange(aes(ymin = asymp.LCL * 100, ymax = asymp.UCL * 100)) +
  coord_flip() +
  labs(x = NULL, y = "predicted probability of terminal damage (%)",
       title = "Probability a terminal shoot is damaged, by species",
       subtitle = "GLMM estimate for an average plot, 95% CI") +
  theme_classic(base_size = 12)


# show with raw data

raw_sp <- plot_sp %>%
  group_by(species) %>%
  summarise(raw = 100 * sum(n_dmg_term) / sum(n_ind), n_ind = sum(n_ind))

ggplot(emm_sp, aes(reorder(species, prob), prob * 100)) +
  geom_pointrange(aes(ymin = asymp.LCL * 100, ymax = asymp.UCL * 100)) +
  geom_point(data = raw_sp, aes(species, raw),
             shape = 17, size = 3, colour = "red", 
             inherit.aes = FALSE) +
  geom_text(data = raw_sp, aes(species, y = 0, label = paste0("n=", n_ind)),
            hjust = 1.2, size = 3, inherit.aes = FALSE) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "% terminal damage",
       subtitle = "black = model estimate (average plot); red x = raw pooled rate") +
  theme_classic(base_size = 12)



# dummy test on real data -----------------------------------------------------------
library(data.table); library(dplyr); library(tidyr); library(lme4)

df <- fread("inData_Michal/poskoz_jedinci_clean.csv",
            na = c("", "NA"))



# ---- 1. one input table for ALL models ------------------------------------------------
d <- df %>%
  filter(species != "ots1", plot != 117, vegtype %in% c("small", "advanced")) %>%
  mutate(planting      = replace_na(planting, 0),
         anti_browsing = replace_na(anti_browsing, 0),
         # V1: counts of damaged individuals (capped at n)
         term_dmg = pmin(tot_dmg_terminal, n),                              # terminal (small only)
         fol_dmg  = pmin(tot_dmg_foliage, n),                               # foliage
         any_dmg  = pmin(n, tot_dmg_terminal + tot_dmg_foliage + dmg_stem_numeric),  # all types (upper bound)
         # V2: damage on the tallest individual, any type (frost leftovers with no counts -> 0)
         tallest_dmg = as.integer(!is.na(dmg_jedinec) & dmg_jedinec == 1 & any_dmg > 0))

# 5 main species as in your script, the rest pooled
d <- d %>% mutate(species = relevel(factor(ifelse(species %in% c("SM","BK","BR","DB","BO"), species, "other")), ref = "SM"))

small <- filter(d, vegtype == "small")      # terminal was measured only here
reg   <- d                                  # small + advanced

# ---- 2. models: identical fixed and random effects ----------------------------------------
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
fit <- function(resp, data) glmer(as.formula(paste(resp, "~ species + planting + anti_browsing + (1 | plot)")),
                                  data, binomial, control = ctl)

mods <- list(
  "V1 terminal  (small)"      = fit("cbind(term_dmg, n - term_dmg)", small),
  "V1 foliage   (small)"      = fit("cbind(fol_dmg, n - fol_dmg)",   small),
  "V1 any type  (small)"      = fit("cbind(any_dmg, n - any_dmg)",   small),
  "V2 tallest   (small)"      = fit("tallest_dmg",                    small),
  "V1 foliage   (small+adv)"  = fit("cbind(fol_dmg, n - fol_dmg)",   reg),
  "V1 any type  (small+adv)"  = fit("cbind(any_dmg, n - any_dmg)",   reg),
  "V2 tallest   (small+adv)"  = fit("tallest_dmg",                    reg)
)

# ---- 3. odds ratios (planting, anti_browsing) ----------------------------------------------
or_tab <- function(m, label) {
  s <- summary(m)$coefficients
  data.frame(model = label, term = rownames(s), OR = exp(s[, 1]),
             lwr = exp(s[, 1] - 1.96 * s[, 2]), upr = exp(s[, 1] + 1.96 * s[, 2]), p = s[, 4], row.names = NULL)
}
res <- bind_rows(Map(or_tab, mods, names(mods))) %>%
  filter(term %in% c("planting", "anti_browsing")) %>%
  mutate(across(c(OR, lwr, upr), ~ round(.x, 2)), p = round(p, 3))
write.csv(res, "outDataMichal/cmp_res.csv", row.names = FALSE)

cat("\nplot variance / records / damaged:\n")
for (nm in names(mods)) cat(sprintf("%-26s plot var = %.2f | records = %d\n", nm, as.data.frame(VarCorr(mods[[nm]]))$vcov[1], nobs(mods[[nm]])))
cat("\nraw: terminal", sum(small$term_dmg), "ind | foliage", sum(small$fol_dmg), "| any (V1)", sum(small$any_dmg),
    "ind | tallest (V2)", sum(small$tallest_dmg), "records of", nrow(small), "\n")
cat("\nconvergence messages per model:\n")
for (nm in names(mods)) cat(sprintf("%-26s %s\n", nm, if (length(mods[[nm]]@optinfo$conv$lme4$messages)) paste(mods[[nm]]@optinfo$conv$lme4$messages, collapse = " | ") else "none"))


# make plots ----------------------------------------------------------------
# needs the list `mods` of fitted models (created by cmp_v1_v2.R) - run that script first
library(dplyr); library(ggplot2); library(MASS, exclude = "select")   # MASS only for mvrnorm

use <- c("V1 terminal  (small)", "V1 foliage   (small+adv)", "V2 tallest   (small+adv)")
lab <- c("V1 terminal (small)",  "V1 foliage (small+adv)",  "V2 tallest (small+adv)")

# ---- predicted probabilities averaged over plots (marginal), with simulation CIs -------------
marg_draws <- function(m, newdata, nsim = 1000) {
  X    <- model.matrix(delete.response(terms(formula(m, fixed.only = TRUE))), newdata)
  B    <- MASS::mvrnorm(nsim, fixef(m), as.matrix(vcov(m)))          # uncertainty of fixed effects
  sd_u <- sqrt(as.data.frame(VarCorr(m))$vcov[1])                    # between-plot sd (logit scale)
  u    <- qnorm(ppoints(60)) * sd_u                                  # plot effects to average over
  est  <- rowMeans(plogis(outer(drop(X %*% fixef(m)), u, "+")))      # point estimate
  dr   <- sapply(1:nsim, function(i) rowMeans(plogis(outer(drop(X %*% B[i, ]), u, "+"))))
  list(est = est, draws = dr)
}

# ---- 1. predicted % damage by species x management variable (var = "anti_browsing" or "planting") ----
pred_by_species <- function(var, var_labels) {
  bind_rows(Map(function(nm, lb) {
    m  <- mods[[nm]]; mf <- model.frame(m)
    nd <- expand.grid(species = levels(mf$species), v = c(0, 1))
    nd$planting <- mean(mf$planting); nd$anti_browsing <- mean(mf$anti_browsing)   # other variable at its mean
    nd[[var]] <- nd$v
    r  <- marg_draws(m, nd)
    nd %>% mutate(model = lb, pct = 100 * r$est,
                  lwr = 100 * apply(r$draws, 1, quantile, 0.025), upr = 100 * apply(r$draws, 1, quantile, 0.975),
                  group = factor(v, labels = var_labels))
  }, use, lab)) %>% mutate(model = factor(model, levels = lab))
}

plot_pred <- function(d, title) {
  ord <- d %>% group_by(species) %>% summarise(m = mean(pct)) %>% arrange(m) %>% pull(species)
  d %>% mutate(species = factor(species, levels = ord)) %>%
    ggplot(aes(pct, species, colour = group)) +
    geom_linerange(aes(xmin = lwr, xmax = upr), position = position_dodge(0.6), linewidth = 0.6) +
    geom_point(position = position_dodge(0.6), size = 2.4) +
    facet_wrap(~ model, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = c("#8f8e88", "#2a78d6"), name = NULL) +
    labs(x = "predicted % damaged (average over plots, 95% CI)", y = NULL, title = title,
         caption = "other management variable held at its mean; CI = fixed-effect uncertainty only") +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
          legend.position = "bottom", panel.grid.major.x = element_line(colour = "grey92"))
}

d_prot <- pred_by_species("anti_browsing", c("not protected", "protected"))
d_plan <- pred_by_species("planting",      c("not planted", "planted"))
ggsave("outDataMichal/prob_species_protection.png", plot_pred(d_prot, "Predicted damage by species and browsing protection"), width = 7, height = 4.2, dpi = 200, bg = 'white')
ggsave("outDataMichal/prob_species_planting.png",   plot_pred(d_plan, "Predicted damage by species and planting"),            width = 7, height = 4.2, dpi = 200, bg = "white")

# ---- 2. management effect on the probability scale: percentage-point difference per species ----------
pp_diff <- function(var) {
  bind_rows(Map(function(nm, lb) {
    m  <- mods[[nm]]; mf <- model.frame(m)
    nd <- expand.grid(species = levels(mf$species), v = c(0, 1))
    nd$planting <- mean(mf$planting); nd$anti_browsing <- mean(mf$anti_browsing); nd[[var]] <- nd$v
    r  <- marg_draws(m, nd); i1 <- nd$v == 1; i0 <- nd$v == 0
    dd <- 100 * (r$draws[i1, ] - r$draws[i0, ])
    data.frame(species = nd$species[i1], model = lb, diff = 100 * (r$est[i1] - r$est[i0]),
               lwr = apply(dd, 1, quantile, 0.025), upr = apply(dd, 1, quantile, 0.975))
  }, use, lab)) %>% mutate(model = factor(model, levels = lab))
}
p_diff <- pp_diff("anti_browsing") %>%
  ggplot(aes(diff, species)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_linerange(aes(xmin = lwr, xmax = upr), colour = "#2a78d6", linewidth = 0.7) +
  geom_point(colour = "#2a78d6", size = 2.6) +
  facet_wrap(~ model, nrow = 1, scales = "free_x") +
  labs(x = "change in % damaged, protected minus not protected (percentage points)", y = NULL,
       title = "Effect of browsing protection on the probability scale") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_line(colour = "grey92"))
ggsave("outDataMichal/prob_protection_difference.png", p_diff, width = 7, height = 4, dpi = 300)

# ---- calibration check: averaged prediction vs raw rate -----------------------------------------
for (k in seq_along(use)) {
  m <- mods[[use[k]]]; mf <- model.frame(m); r <- marg_draws(m, mf[, c("species", "planting", "anti_browsing")])
  y <- mf[[1]]; if (is.matrix(y)) { n <- rowSums(y); dm <- y[, 1] } else { n <- rep(1, nrow(mf)); dm <- y }
  cat(sprintf("%-26s raw %.1f%% | model average %.1f%%\n", lab[k], 100 * sum(dm) / sum(n), 100 * sum(n * r$est) / sum(n)))
}
print(d_prot %>% filter(species %in% c("SM", "BK", "DB"), model == lab[1]) %>% mutate(across(c(pct, lwr, upr), ~ round(.x, 1))) %>% select(species, group, pct, lwr, upr))
