

# Read adat from Miso - poskozeni
# identify the tallest stem per species
# check wheather the taller/shorter individual has damaged terminal, or more of them
# decide whatther we work on plot or siubplot level
# make summary boxplots

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)


df <- fread("from_Miso_poskozeni/poskoz_jedinci_clean.csv",
               na = c("", "NA"))

View(df)

# # get total sum of indivudual per plots
# get total sum of damaged terminals per plot
# quantify the anti_browing intensity
df_small <- df %>% 
  filter(vegtype == "small") 

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
