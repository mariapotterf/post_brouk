

# spatial correlation in damage types?

#separate per terminal/foliage/stem


library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(sf)
library(spdep)
library(tidyr)


df <- fread("inData_Michal/poskoz_jedinci_clean.csv",
            na = c("", "NA"))


plot_coords <- df %>%
  distinct(plot, subplot, x, y) %>%
  group_by(plot) %>%
  summarise(x = mean(x, na.rm = TRUE),
            y = mean(y, na.rm = TRUE),
            n_subplots = n(),
            .groups = "drop")

# convert to sf object
plot_terminal_sf <- st_as_sf(plot_terminal, coords = c("x", "y"), crs = 4326) %>%
  st_transform(3035)

plot(st_geometry(plot_terminal_sf))

#View(df)

head(df)

table(df$dmg_stem_cause)

table(df$dmg_foliage )

# sanity check first - should be 0
df %>% filter(tot_dmg_terminal > n) %>% nrow()

# decide, at which level to do the analysis? subplot = can be too noisy or plot? ------

## damage clustering for Terminal + terminal similar -----------------------------


subplot_terminal <- df %>%
  group_by(plot, subplot) %>%
  summarise(n_ind        = sum(n, na.rm = TRUE),
            n_terminal   = sum(tot_dmg_terminal, na.rm = TRUE),
            pct_terminal = 100 * n_terminal / n_ind,
            .groups = "drop")

plot_terminal <- df %>%
  group_by(plot) %>%
  summarise(n_ind        = sum(n, na.rm = TRUE),
            n_terminal   = sum(tot_dmg_terminal, na.rm = TRUE),
            pct_terminal = 100 * n_terminal / n_ind,
            .groups = "drop")

summary(subplot_terminal$n_ind)
summary(plot_terminal$n_ind)


# same thing by height class - where does terminal damage actually sit?
ggplot(df, aes(x = hgt, y = as.numeric(tot_dmg_terminal))) +
  stat_summary(fun = mean, geom = "col") +
  labs(x = "Height class", y = "Proportion with terminal damage",
       title = "Terminal damage by height class, before collapsing to tallest") +
  theme_minimal()


# how many plots are we trusting on very few individuals?
plot_terminal %>% filter(n_ind < 5) %>% nrow()

ggplot(plot_terminal, aes(x = n_ind, y = pct_terminal)) +
  geom_point(alpha = 0.5) +
  labs(x = "N individuals in plot", y = "% terminal damage",
       title = "Plot-level terminal damage rate vs. sample size backing it") +
  theme_minimal()


# add coordinates 

plot_terminal <- plot_terminal %>%
  left_join(plot_coords, by = "plot")

# sanity check - any plot missing a coordinate?
plot_terminal %>% filter(is.na(x) | is.na(y)) %>% nrow()

## make up neighbours list k = 6 ------------------

coords_m <- st_coordinates(plot_terminal_sf)

nb <- knn2nb(knearneigh(coords_m, k = 6))
summary(nb)          # check link count, look for anything odd

lw <- nb2listw(nb, style = "W")   # row-standardized weights

# waht are the groups with far away neighbors?
comp <- n.comp.nb(nb)
table(comp$comp.id)

plot_terminal_sf$component <- factor(comp$comp.id)

ggplot(plot_terminal_sf) +
  geom_sf(aes(color = component), size = 2) +
  theme_minimal() +
  labs(title = "Plot groups implied by the k=6 neighbor graph")


nbdists(nb, coords_m)

# if my plots have exactly teh same numbers of neighbours, those 'neighbors' 
# can be very far away from each other - 25 km - hence, maybe not representative for 
# local deer pressure

# maybe define it by distance range? chenck first how far apart my plots are
nn1 <- knn2nb(knearneigh(coords_m, k = 1))
d1  <- unlist(nbdists(nn1, coords_m))

summary(d1)
hist(d1, breaks = 40)

# define teh max distance between teh neighbours
d_max <- 1500  # metres - sits in the gap between the dense local cluster
# (median 420 m, Q3 633 m) and the next distance tier (~3500 m+),
# so the exact value here doesn't change which plots connect

nb_dist <- dnearneigh(coords_m, d1 = 0, d2 = d_max)
summary(nb_dist)

# show tehm n the map
plot_terminal_sf$isolated <- card(nb_dist) == 0

plot(st_geometry(plot_terminal_sf),
     col = ifelse(plot_terminal_sf$isolated, "red", "black"),
     pch = 16, main = "Neighbour links at d_max = 1500 m (red = isolated)")
plot(nb_dist, coords_m, add = TRUE, col = "steelblue", lwd = 0.5)


# compare each plot to its neighbours - empirical Bayes smooth
eb <- EBlocal(ri = plot_terminal$n_terminal,
              ni = plot_terminal$n_ind,
              nb = nb_dist,
              zero.policy = TRUE)

plot_terminal$pct_terminal_smoothed <- eb$est * 100

# isolated plots have no neighbours to borrow strength from - see what
# EBlocal did with them before trusting this column
sum(is.na(plot_terminal$pct_terminal_smoothed))

plot_terminal %>%
  filter(is.na(pct_terminal_smoothed)) %>%
  select(plot, n_ind, pct_terminal, pct_terminal_smoothed) |> 
  print(n = 50)

# in isolated cases (34 locations) cases, pct_terminal should be  0 (as tere is no neighbors to compare with)

plot_terminal <- plot_terminal %>%
  mutate(pct_terminal_smoothed = ifelse(is.na(pct_terminal_smoothed),
                                        pct_terminal, pct_terminal_smoothed))

sum(is.na(plot_terminal$pct_terminal_smoothed))   # should be 0 now

set.seed(1)

lw_dist <- nb2listw(nb_dist, style = "W", zero.policy = TRUE)
moran.mc(plot_terminal$pct_terminal_smoothed, listw = lw_dist, nsim = 999,
         zero.policy = TRUE)
# ok, Morasn's foudn thare are differences, but now i need t find hotspots

nb_star <- include.self(nb_dist)
lw_star <- nb2listw(nb_star, style = "W", zero.policy = TRUE)

gi <- localG(plot_terminal$pct_terminal_smoothed, lw_star, zero.policy = TRUE)
plot_terminal$gi_z <- as.numeric(gi)

summary(plot_terminal$gi_z)

# check for isolated cases (have no neighbours to compare with)
plot_terminal %>%
  filter(is.na(gi_z) | is.nan(gi_z)) %>%
  select(plot, n_ind, pct_terminal_smoothed, gi_z)

# classify teh hotspots
plot_terminal <- plot_terminal %>%
  mutate(hotspot = case_when(
    gi_z >=  1.96 ~ "Hotspot (95%)",
    gi_z <= -1.96 ~ "Coldspot (95%)",
    is.na(gi_z)   ~ "Isolated",
    TRUE          ~ "ns"
  ))

table(plot_terminal$hotspot, useNA = "always")

isolated_lookup <- tibble(plot = plot_terminal_sf$plot, isolated = card(nb_dist) == 0)

plot_terminal <- plot_terminal %>%
  left_join(isolated_lookup, by = "plot") %>%
  mutate(hotspot = ifelse(isolated, "No neighbours", hotspot))

table(plot_terminal$hotspot, useNA = "always")

plot_terminal_sf <- plot_terminal_sf %>%
  left_join(plot_terminal %>% select(plot, gi_z, hotspot), by = "plot")

ggplot(plot_terminal_sf) +
  geom_sf(aes(color = hotspot), size = 2.5) +
  scale_color_manual(values = c("Hotspot (95%)"  = "red",
                                "Coldspot (95%)" = "blue",
                                "ns" = "grey70",
                                "Isolated"   = "grey90")) +
  theme_minimal() +
  labs(title = "Terminal damage hot/cold spots (Gi*, plot level)")

# do hotspots have anything in common in term of management?
site_by_plot <- df %>%
  distinct(plot, subplot, stump, clear, grndwrk, logging_trail,
           windthrow, standing_deadwood, anti_browsing, planting) %>%
  group_by(plot) %>%
  summarise(across(c(stump, clear, grndwrk, logging_trail, windthrow,
                     standing_deadwood, anti_browsing, planting),
                   ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

plot_terminal %>%
  select(plot, hotspot) %>%
  left_join(site_by_plot, by = "plot") %>%
  group_by(hotspot) %>%
  summarise(across(c(stump, clear, grndwrk, logging_trail, windthrow,
                     standing_deadwood, anti_browsing, planting),
                   ~ mean(.x, na.rm = TRUE)),
            n = n())


# Foliage - also recorded on tallest and similar -----------------
plot_foliage <- df %>%
  group_by(plot) %>%
  summarise(n_ind       = sum(n, na.rm = TRUE),
            n_foliage   = sum(tot_dmg_foliage, na.rm = TRUE),
            pct_foliage = 100 * n_foliage / n_ind,
            .groups = "drop") %>%
  left_join(plot_coords, by = "plot")

summary(plot_foliage$n_foliage)

# Stem damage by cause ---------------------------------------------
# subplots that had any trees recorded at all - the real "opportunity"
# to have observed game-caused stem damage
subplots_with_trees <- df %>%
  filter(!is.na(n) & n > 0) %>%
  distinct(plot, subplot)

# subplots where at least one row shows zver-caused stem damage
subplots_with_zver <- df %>%
  filter(dmg_stem_cause == "zver") %>%
  distinct(plot, subplot)

plot_zver <- subplots_with_trees %>%
  group_by(plot) %>%
  summarise(n_subplots = n(), .groups = "drop") %>%
  left_join(
    subplots_with_zver %>% group_by(plot) %>% summarise(n_zver_subplots = n(), .groups = "drop"),
    by = "plot"
  ) %>%
  mutate(n_zver_subplots = replace_na(n_zver_subplots, 0),
         pct_zver = 100 * n_zver_subplots / n_subplots) %>%
  left_join(plot_coords, by = "plot")

summary(plot_zver$n_subplots)
table(plot_zver$n_zver_subplots)
sum(plot_zver$n_zver_subplots > 0)   # how many plots have ANY zver-positive subplot
