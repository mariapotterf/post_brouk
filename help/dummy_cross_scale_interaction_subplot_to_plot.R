# cross scale interaction: dummy example


# processes: cancellation, reinforcement, persistence

# dummy example: cross scale interaction
# Two subplots in one plot: perfectly uniform inside each subplot,
# but different means between subplots
h1 <- c(1,1,1,1)   # subplot 1, plot 1
h2 <- c(3,3,3,3)   # subplot 2

h3 <- c(1,2,3,4)   # subplot 3
h4 <- c(1,2,3,4)   # subplot 3


cv1 <- sd(h1)/mean(h1)         # 0
cv2 <- sd(h2)/mean(h2)         # 0
cv3 <- sd(h3)/mean(h3)         # 0.51
cv4 <- sd(h4)/mean(h4)         # 0.51

mean_cvA <- mean(c(cv1, cv2))  # 0
mean_cvB <- mean(c(cv3, cv4))  # 0.51

pooledA <- c(h1, h2)
pooled_cvA <- sd(pooledA)/mean(pooledA)  # > 0 (here = 0.5)

pooledB <- c(h3, h4)
pooled_cvB <- sd(pooledB)/mean(pooledB)  # > 0 (here = 0.5)

c(mean_cv = mean_cvA, pooled_cv = pooled_cvA, delta = pooled_cvA - mean_cvA)
c(mean_cv = mean_cvB, pooled_cv = pooled_cvB, delta = pooled_cvB - mean_cvB)


# dummay example - continue 

library(tidyverse)

# Input vectors

# A: different mean, SD = 0 -> Increase in avg CV per plot
h1 <- c(1,1,1,1)   
h2 <- c(3,3,3,3)   

# B: double SD, double mean -> same avg CV per plot
h3 <- c(2,4,6,8)   
h4 <- c(4,8,12,16)   

# C: mean is about the same, SD is very different -> negative
h5 <- c(1,2,3,100)   # High variance, high CV
h6 <- c(51,52,51,52) # Low variance, low CV

# D: mean and sd are the same -> variance remains the same
h7 <- c(1,1,1,1)
h8 <- c(1,1,1,1)

# Create raw data table with values expanded
cv_data_long <- tibble(
  plot = c("A", "A", "B", "B", "C", "C", "D", "D"),
  subplot = c("h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8"),
  values = list(h1, h2, h3, h4, h5, h6, h7, h8)
) %>%
  unnest(values) %>%
  group_by(plot, subplot) %>%
  mutate(
    mean = mean(values),
    sd = sd(values),
    cv = sd / mean
  ) %>%
  ungroup()

# Per-subplot summary (optional)
cv_data <- cv_data_long %>%
  group_by(plot, subplot) %>%
  summarise(
    mean = mean(values),
    sd = sd(values),
    cv = sd / mean,
    .groups = "drop"
  )

# Correct summary table with proper pooling per plot
summary_cv <- cv_data_long %>%
  group_by(plot) %>%
  summarise(
    pooled_mean = mean(values),
    pooled_sd = sd(values), 
    mean_cv = mean(cv),
    pooled_cv = pooled_sd / pooled_mean,
    delta = pooled_cv - mean_cv,
    .groups = "drop"
  )

# View outputs
#cv_data
summary_cv


