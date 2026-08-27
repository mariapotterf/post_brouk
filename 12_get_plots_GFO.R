# ============================================================
# Final model plots — Czechia
# Lightweight script: loads exported final models + summaries only.
# No data prep, no model fitting, no diagnostics code needed here.
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(mgcv)
library(ggplot2)
library(marginaleffects)   # for back-transforming coefficients to % change later

# ---- load exported model bundle ----

fin.models.all           <- readRDS("outModel/final_gam_models.rds")
models_with_cluster      <- readRDS("outModel/final_gam_models_with_cluster.rds")

final_coefs              <- read.csv("outModel/final_gam_coefficients.csv")
cluster_sensitivity_wide <- read.csv("outModel/cluster_sensitivity_wide.csv")

plot_coords              <- readRDS("outModel/plot_coords.rds")

# sanity check everything loaded as expected
names(fin.models.all)
names(models_with_cluster)
str(final_coefs)
str(cluster_sensitivity_wide)