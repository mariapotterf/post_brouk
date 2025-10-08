
# ---------------------------------
# make plots for each SPEI
# ---------------------------------

# SPEI multi-window summarizer (folders: spei1_new, spei3_new, spei6_new, spei12_new)
# Dependencies
library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(data.table)
library(ggplot2)
library(scales)

# ---- CONFIG -----------------------------------------------------------------
# Root folder that contains spei*_new subfolders (adjust to your project)
root_dir <- "raw"   # e.g., "." or "raw" depending on where those folders live

# discover spei*_new folders (case-insensitive)
spei_roots <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE) |>
  (\(x) x[grepl("(?i)/SPEI(1|3|6|12)_new$", x)])()
print(spei_roots)

# ---- HELPERS ----------------------------------------------------------------
# Get SPEI window (1/3/6/12) from a path; returns integer or NA
spei_from_path <- function(path) {
  m <- stringr::str_match(basename(path), "(?i)SPEI\\s*[-_ ]?(\\d{1,2})")[,2]
  suppressWarnings(as.integer(m))
}

# Parse YYYY, MM, DD from filename
get_ymd <- function(x) {
  m <- stringr::str_match(basename(x), "((19|20)\\d{2})[-_ ]?(\\d{1,2})[-_ ]?(\\d{1,2})")
  tibble::tibble(
    year  = suppressWarnings(as.integer(m[,2])),
    month = suppressWarnings(as.integer(m[,4])),
    day   = suppressWarnings(as.integer(m[,5]))
  )
}

# Summarise a single scene folder (e.g., SPEI1_harg1_split)
summarize_scene <- function(scene_path, spei_hint = NA_integer_) {
  files <- list.files(scene_path, pattern = "\\.(tif|tiff)$",
                      full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    warning("No GeoTIFFs in: ", scene_path, " — skipping.")
    return(list(monthly = tibble(), yearly = tibble(), inventory = tibble()))
  }
  
  # metadata (SPEI from filename; fall back to hint from parent spei*_new)
  meta_date <- get_ymd(files)
  spei_file <- vapply(files, spei_from_path, integer(1))
  spei_vec  <- ifelse(is.na(spei_file), spei_hint, spei_file)
  
  keep <- which(!is.na(meta_date$year) & !is.na(meta_date$month) & !is.na(spei_vec))
  if (!length(keep)) {
    warning("No parsable SPEI+date in: ", scene_path, " — skipping.")
    return(list(monthly = tibble(), yearly = tibble(), inventory = tibble()))
  }
  files     <- files[keep]
  meta_date <- meta_date[keep, , drop = FALSE]
  spei_vec  <- spei_vec[keep]
  
  # quick inventory of what's found
  inventory <- tibble::tibble(
    scene = basename(scene_path),
    spei  = factor(spei_vec, levels = c(1,3,6,12), ordered = TRUE),
    year  = meta_date$year,
    month = meta_date$month
  ) |>
    dplyr::count(scene, spei, year, month, name = "n_files") |>
    dplyr::arrange(spei, year, month)
  
  # load rasters
  r <- terra::rast(files)
  names(r) <- basename(files)
  
  # ---------- (1) MONTHLY: per-layer stats across cells ----------------------
  g <- terra::global(r, fun = c("mean", "sd"), na.rm = TRUE)
  monthly_df <- tibble::as_tibble(g) |>
    dplyr::bind_cols(inventory |> dplyr::select(scene, spei, year, month)) |>
    dplyr::mutate(
      date = as.Date(sprintf("%04d-%02d-01", year, month))
    ) |>
    dplyr::relocate(scene, spei, year, month, date, mean, sd) |>
    dplyr::arrange(spei, year, month)
  
  # ---------- (2) YEARLY: across months × cells for year×SPEI ----------------
  idx_tbl <- tibble::tibble(idx = seq_len(terra::nlyr(r)), year = meta_date$year, spei = spei_vec)
  
  yearly_df <- idx_tbl |>
    dplyr::group_by(year, spei) |>
    dplyr::summarise(
      mean = {
        rr <- r[[idx]]
        vals <- terra::values(rr, mat = TRUE)
        vals <- as.vector(vals)
        vals <- vals[is.finite(vals)]
        if (!length(vals)) NA_real_ else mean(vals, na.rm = TRUE)
      },
      sd = {
        rr <- r[[idx]]
        vals <- terra::values(rr, mat = TRUE)
        vals <- as.vector(vals)
        vals <- vals[is.finite(vals)]
        if (!length(vals)) NA_real_ else sd(vals, na.rm = TRUE)
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      scene = basename(scene_path),
      spei  = factor(spei, levels = c(1,3,6,12), ordered = TRUE)
    ) |>
    dplyr::select(scene, spei, year, mean, sd) |>
    dplyr::arrange(spei, year)
  
  list(monthly = monthly_df, yearly = yearly_df, inventory = inventory)
}

# ---- WALK: for each spei*_new, summarise its scene subfolders ---------------
res_all <- purrr::map(spei_roots, function(s_root) {
  # deduce SPEI window from folder name (used as hint)
  spei_hint <- spei_from_path(s_root)
  scenes <- list.dirs(s_root, full.names = TRUE, recursive = FALSE)
  purrr::map(scenes, ~ summarize_scene(.x, spei_hint = spei_hint))
})

# bind results
summary_monthly <- dplyr::bind_rows(purrr::flatten(res_all) |> purrr::map("monthly"))
summary_yearly  <- dplyr::bind_rows(purrr::flatten(res_all) |> purrr::map("yearly"))
inventory_all   <- dplyr::bind_rows(purrr::flatten(res_all) |> purrr::map("inventory"))

# ensure clean types
summary_monthly <- summary_monthly |> dplyr::mutate(year = as.integer(year))
summary_yearly  <- summary_yearly  |> dplyr::mutate(year = as.integer(year))

# write out
dir.create("outTable", showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(inventory_all,   "outTable/SPEI_inventory.csv")
data.table::fwrite(summary_monthly, "outTable/SPEI_monthly_summary.csv")
data.table::fwrite(summary_yearly,  "outTable/SPEI_yearly_summary.csv")

# ---- PLOTS ------------------------------------------------------------------
# Monthly: mean ± SD over time, faceted by SPEI window
year_bands <- summary_monthly |>
  dplyr::distinct(year) |>
  dplyr::transmute(
    year, xmin = as.Date(sprintf("%04d-01-01", year)),
    xmax = as.Date(sprintf("%04d-12-31", year)),
    is_even = (year %% 2) == 0
  ) |>
  dplyr::filter(is_even)

# ---- PLOTS ------------------------------------------------------------------
# Monthly: mean ± SD over time, faceted by SPEI window
year_bands <- summary_monthly |>
  dplyr::distinct(year) |>
  dplyr::transmute(
    year, xmin = as.Date(sprintf("%04d-01-01", year)),
    xmax = as.Date(sprintf("%04d-12-31", year)),
    is_even = (year %% 2) == 0
  ) |>
  dplyr::filter(is_even)
 
p_monthly <- summary_monthly %>% 
  dplyr::filter(grepl("(?i)harg\\s*?1([^0-9]|$)", scene, perl = TRUE)) %>% 
  ggplot(aes(x = date, y = mean,  group = scene)) +
  # geom_rect(
  #   data = year_bands,
  #   aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
  #   inherit.aes = FALSE, fill = "grey85", alpha = 0.25
  # ) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, alpha = 0.7) +
  geom_line() +
  geom_point(size = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ spei, ncol = 1, scales = "free_y") +
  labs(x = "Date", y = "SPEI (monthly mean ± SD)",# color = "Scene", fill = "Scene",
       title = "Monthly SPEI by scene") +
  theme_bw(base_size = 9)

p_monthly


p_monthly_cols <- summary_monthly %>% 
  dplyr::filter(grepl("(?i)harg\\s*?1([^0-9]|$)", scene, perl = TRUE)) %>% 
  dplyr::mutate(sign = ifelse(mean > 0, "Positive (> 0)", "Negative (< 0)")) %>%
  ggplot(aes(x = date, y = mean, group = scene, color = sign)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, alpha = 0.7) +
  geom_line() +
  geom_point(size = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ spei, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Negative (< 0)" = "red", "Positive (> 0)" = "blue")) +
  labs(
    x = "Date",
    y = "SPEI (monthly mean ± SD)",
    title = "Monthly SPEI by scene",
    color = NULL
  ) +
  theme_bw(base_size = 9)

p_monthly_cols


p_monthly_cols_sub <- summary_monthly %>% 
  dplyr::filter(grepl("(?i)harg\\s*?1([^0-9]|$)", scene, perl = TRUE)) %>% 
  dplyr::filter(year > 1999) %>% 
  
  dplyr::mutate(sign = ifelse(mean > 0, "Positive (> 0)", "Negative (< 0)")) %>%
  ggplot(aes(x = date, y = mean, group = scene, color = sign)) +
  geom_rect(
    data = filter(year_bands,year>1999),
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "grey85", alpha = 0.25
  ) +
  
  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, alpha = 0.7) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ spei, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Negative (< 0)" = "red", 
                                "Positive (> 0)" = "blue")) +
  labs(
    x = "Date",
    y = "SPEI (monthly mean ± SD)",
    title = "Monthly SPEI by scene",
    color = NULL
  ) +
  theme_bw(base_size = 9)

p_monthly_cols_sub




# Yearly: mean ± SD per year, faceted by SPEI window
p_yearly <- summary_monthly %>% 
  dplyr::filter(grepl("(?i)harg\\s*?1([^0-9]|$)", scene, perl = TRUE)) %>%  
  ggplot(aes(x = year, y = mean, color = scene, fill = scene, group = scene)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25, alpha = 0.7) +
  geom_line() +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_wrap(~ spei, ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "SPEI (annual mean ± SD)", color = "Scene", fill = "Scene",
       title = "Annual SPEI by scene", subtitle = "Facets are SPEI windows") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

p_monthly
p_yearly
