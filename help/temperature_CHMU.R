# ---- Packages ----
library(googlesheets4)   # install.packages("googlesheets4")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)

# ---- Source: ČHMÚ territorial monthly means (via Fakta o klimatu sheet that cites ČHMÚ) ----
# Sheet: "Teplota ČR od roku 1961 (faktaoklimatu.cz, verze 2022-03-14)"
# Link from https://faktaoklimatu.cz/infografiky/teplota-cr-mesice  -> "Naše tabulka s daty"
sheet_url <- "https://docs.google.com/spreadsheets/d/1uxIsR2OcKT1rF7cvLiT_QrezAXN63ikOIfpef7TZvMc"

# The sheet has a tab with monthly temperatures by year (territorial average, °C)
# We’ll read all tabs and guess the one that has the monthly table
sheets <- gs4_get_sheets(sheet_url)$name
tab_name <- sheets[grepl("měsíce|mesice|měsíční|monthly|tabulka", tolower(sheets))][1]
if (is.na(tab_name)) tab_name <- sheets[1]  # fallback

raw <- read_sheet(sheet_url, sheet = tab_name)

# Expect either a wide shape (columns: Rok, Leden, Únor, ..., Prosinec) or a long table.
# Try to detect/normalize:
nm <- tolower(names(raw))
if (!("rok" %in% nm) && "year" %in% nm) names(raw)[match("year", nm)] <- "rok"
names(raw) <- tolower(names(raw))

# Map month names to numbers
cz_months <- c("leden"=1,"únor"=2,"unor"=2,"březen"=3,"brezen"=3,"duben"=4,"květen"=5,"kveten"=5,
               "červen"=6,"cerven"=6,"červenec"=7,"cervenec"=7,"srpen"=8,"září"=9,"zari"=9,
               "říjen"=10,"rijen"=10,"listopad"=11,"prosinec"=12)

# If wide: pivot to long
if ("rok" %in% names(raw) && any(names(raw) %in% names(cz_months))) {
  dat_long <- raw |>
    rename(year = rok) |>
    pivot_longer(cols = any_of(names(cz_months)), names_to = "month_name", values_to = "tmean_c") |>
    mutate(month = cz_months[month_name]) |>
    select(year, month, tmean_c)
} else {
  # If already long, try to standardize columns
  # Expect columns like: year, month (1–12 or CZ names), value
  mcol <- if ("month" %in% names(raw)) "month" else names(raw)[grepl("měs|mes|month", names(raw))][1]
  vcol <- names(raw)[grepl("tepl|temp|tmean|hodnota|value", names(raw))][1]
  dat_long <- raw |>
    rename(year = any_of(c("rok","year")),
           month_raw = all_of(mcol),
           tmean_c = all_of(vcol)) |>
    mutate(month = ifelse(suppressWarnings(!is.na(as.numeric(month_raw))),
                          as.integer(month_raw),
                          cz_months[tolower(month_raw)])) |>
    select(year, month, tmean_c)
}

# Keep 2000–2025 and months 4–9 (IV–IX)
season <- dat_long |>
  mutate(year = as.integer(year)) |>
  filter(year >= 2000, year <= 2025, month %in% 4:9) |>
  group_by(year) |>
  summarise(
    months_present = sum(!is.na(tmean_c)),
    mean_temp_iv_ix = ifelse(months_present == 6, mean(tmean_c, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  arrange(year) |>
  mutate(source = "ČHMÚ (territorial monthly means) via Fakta o klimatu sheet",
         note_2025 = ifelse(year == 2025 & is.na(mean_temp_iv_ix),
                            "Season incomplete as of 2025-08-29", NA_character_))

# Write CSV
out_csv <- "cz_iv-ix_mean_temp_2000-2025_CHMU.csv"
write_csv(season |> select(year, mean_temp_iv_ix), out_csv)

# Plot
p <- ggplot(season, aes(x = year, y = mean_temp_iv_ix)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Czechia — Mean Air Temperature (°C), Apr–Sep (IV–IX), 2000–2025",
    subtitle = "Territorial monthly means from ČHMÚ (via Fakta o klimatu dataset). 2025 shown only if complete.",
    x = "Year", y = "Mean temperature (°C)"
  ) +
  theme_bw()

ggsave("cz_iv-ix_mean_temp_2000-2025_CHMU.png", p, width = 8, height = 4.5, dpi = 150)

season
