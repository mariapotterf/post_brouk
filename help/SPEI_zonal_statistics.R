# --- PACKAGES ---------------------------------------------------------------
library(terra)
library(dplyr)
library(stringr)
library(tibble)
library(readr)
library(fs)
library(tidyr)
library(data.table)

library(writexl)  # export xlsx for Roman

# --- INPUTS -----------------------------------------------------------------
spei_folder <- "raw/SPEI12/spei12_harg1_split"  # 
okres_shp   <- terra::vect("raw/CR_administrativa/OKRESY_P.shp")

# --- READ FILES ---
files <- list.files(spei_folder, pattern = "\\.tif$", full.names = TRUE)
meta <- tibble(file = files) |>
  mutate(
    name = basename(file),
    year = as.integer(str_extract(name, "(?<=_)(\\d{4})(?=_)")),
    month = as.integer(str_extract(name, "(?<=_\\d{4}_)(\\d{2})(?=_)")),
    day = as.integer(str_extract(name, "(?<=_\\d{4}_\\d{2}_)(\\d{2})(?=\\.)"))
  ) |>
  dplyr::filter(!is.na(year), !is.na(month)) |>
  arrange(year, month, day)

# --- PREPARE VECTORS & RASTER TEMPLATE ---
# Convert SpatVector attribute table to tibble
okres_attrib <- as_tibble(okres_shp) |>
  mutate(okres_id = row_number())

# Assign new attribute table back to SpatVector
okres_shp$okres_id <- okres_attrib$okres_id

# Now it's safe to use `okres_shp` with the new field
v_ok <- okres_shp

r_template <- rast(meta$file[1])

if (!terra::same.crs(r_template, v_ok)) {
  v_ok <- terra::project(v_ok, r_template)
}


# --- RASTERIZE OKRES POLYGONS ---
r_okres <- terra::rasterize(v_ok, r_template, field = "okres_id")

# --- STACK RASTERS & EXTRACT VALUES ---
r_stack <- rast(meta$file)
names(r_stack) <- paste0(meta$year, "_", meta$month)

# combine raster stack with rasterized okres
r_all <- c(r_okres, r_stack)
v_all <- values(r_all)

# Convert matrix to data.table directly
dt <- as.data.table(v_all)

# Name columns - columns are properly named
#setnames(dt, 1, "okres_id")
#date_names <- paste0(meta$year, "_", meta$month)
#setnames(dt, 2:ncol(dt), date_names)

# Filter NA and melt to long format
dt_filt <- dt[!is.na(okres_id)]  # remove NA values
dt_long <- melt(dt_filt,
                id.vars = "okres_id",
                variable.name = "date",
                value.name = "value"
)

# Split date into year and month
#dt_long[, c("year", "month") := tstrsplit(date, "_", fixed = TRUE)]
#dt_long[, `:=`(
#  year = as.integer(year),
#  month = as.integer(month)
#)]

# Fast zonal summary
zonal_stats <- dt_long[
  !is.na(value),
  .(
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    min = min(value),
    max = max(value)
  ),
  by = .(okres_id, date)
]

# add new columns
zonal_stats[, c("year", "month") := tstrsplit(date, "_", fixed = TRUE)]

# --- JOIN ATTRIBUTES ---
okres_attrib <- as_tibble(v_ok) |> dplyr::select(okres_id, NAZEV)
df_final <- left_join(zonal_stats, okres_attrib, by = "okres_id") |>
  relocate(NAZEV, .after = okres_id) |>
  arrange(year, month, okres_id)

# --- SAVE ---
write_csv(df_final, "outTable/okres_SPEI12_monthly_stats.csv")
write_xlsx(df_final, "outTable/okres_SPEI12_monthly_stats.xlsx")
