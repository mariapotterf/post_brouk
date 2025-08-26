

# HELPS: check up climatic data summaries

# get climatic data on yearly bases, and as averages across several seasons;

# data stahed from Stepanek, Czech Globe

# Load required libraries
library(terra)
library(dplyr)
library(tidyr)
library(stringr)

# Set path to the folder
data_path <- "raw/clim_data_CZ_annual"

# Get all .tif files recursively
raster_files <- list.files(data_path, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# Extract stats
results <- lapply(raster_files, function(file) {
  r <- try(rast(file), silent = TRUE)
  
  if (inherits(r, "try-error")) {
    return(data.frame(
      file = basename(file),
      min = NA,
      median = NA,
      max = NA
    ))
  }
  
  r_min <- round(global(r, fun = "min", na.rm = TRUE)[1,1], 2)
  r_max <- round(global(r, fun = "max", na.rm = TRUE)[1,1], 2)
  
  # Extract values and compute median manually
  r_vals <- values(r, na.rm = TRUE)
  r_median <- round(median(r_vals, na.rm = TRUE), 2)
  
  data.frame(
    file = basename(file),
    min = r_min,
    median = r_median,
    max = r_max
  )
})

# Combine into one dataframe
summary_df <- bind_rows(results)

# Parse filename into variable, year, stat
summary_df <- summary_df %>%
  mutate(file = str_remove(file, "\\.tif$")) %>%
  separate(file, into = c("variable", "year", "period_stat"), sep = "__") %>%
  mutate(stat = str_extract(period_stat, "(mean|sum)")) %>%
  select(variable, year, stat, min, median, max)

# Display result
print(summary_df)

write.csv(summary_df, "outTable/annual_climate_summary_stats.csv", row.names = FALSE)
