

# spatial correlation in damage types?

#separate per terminal/foliage/stem


library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(sf)
library(spdep)


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


# add cpoordinates 

plot_terminal <- plot_terminal %>%
  left_join(plot_coords, by = "plot")

# sanity check - any plot missing a coordinate?
plot_terminal %>% filter(is.na(x) | is.na(y)) %>% nrow()

## make up neighbours list k = 6 ------------------

coords_m <- st_coordinates(plot_terminal_sf)

nb <- knn2nb(knearneigh(coords_m, k = 6))
summary(nb)          # check link count, look for anything odd

lw <- nb2listw(nb, style = "W")   # row-standardized weights


