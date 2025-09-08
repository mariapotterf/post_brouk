

# Analysize SPEI12 value

library(terra)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)

# --- auto-detect candidate subfolders under raw_SPEI12 ----------------------
base_dir <- "raw/SPEI12"

# If you want to force specific names, set them here (comment out auto-detect)
# folders <- c(file.path(base_dir, "spei12_harg1_split"),
#              file.path(base_dir, "spei12_harg2_split"))

# Otherwise, detect all immediate subfolders under base_dir
folders <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
print(folders)

# --- helper to parse YYYY_MM_DD from filenames like spei12_2000_01_01.tif ---
get_year_month_day <- function(x) {
  m <- str_match(basename(x), "((19|20)\\d{2})[-_ ]?(\\d{1,2})[-_ ]?(\\d{1,2})")
  tibble(
    year  = suppressWarnings(as.integer(m[,2])),
    month = suppressWarnings(as.integer(m[,4])),
    day   = suppressWarnings(as.integer(m[,5]))
  )
}

# --- summarize function (unchanged, except recursive search & skip empty) ---
summarize_folder <- function(folder_name) {
  files <- list.files(folder_name,
                      pattern = "\\.(tif|tiff)$",
                      full.names = TRUE, ignore.case = TRUE,
                      recursive = TRUE)  # <— search nested folders too
  
  if (length(files) == 0) {
    warning("No GeoTIFFs found in: ", folder_name, " — skipping.")
    return(list(monthly = tibble(), yearly = tibble()))
  }
  
  meta <- get_year_month_day(files)
  keep <- which(!is.na(meta$year) & !is.na(meta$month))
  if (!length(keep)) {
    warning("No files with parsable year+month in: ", folder_name, " — skipping.")
    return(list(monthly = tibble(), yearly = tibble()))
  }
  
  files <- files[keep]
  meta  <- meta[keep, , drop = FALSE]
  
  # quick inventory
  inv <- meta %>% count(year, month, name = "n_files")
  message("Folder: ", folder_name, " | total monthly rasters: ", nrow(meta))
  # print(head(inv, 6))
  
  r_all <- rast(files)
  names(r_all) <- basename(files)
  
  # (1) monthly summaries (across cells)
  monthly_df <- imap_dfr(seq_len(nlyr(r_all)), function(i, idx) {
    vals <- values(r_all[[i]], mat = FALSE)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      tibble(mean = NA_real_, sd = NA_real_,
             median = NA_real_, q25 = NA_real_, q75 = NA_real_)
    } else {
      qs <- quantile(vals, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      tibble(mean = mean(vals, na.rm = TRUE),
             sd = sd(vals, na.rm = TRUE),
             median = median(vals, na.rm = TRUE),
             q25 = qs[1], q75 = qs[2])
    }
  }) %>%
    bind_cols(meta[, c("year","month")]) %>%
    mutate(folder = basename(folder_name)) %>%
    relocate(folder, year, month) %>%
    arrange(year, month)
  
  # (2) yearly summaries (across months & cells)
  idx_by_year <- split(seq_len(nlyr(r_all)), monthly_df$year)
  
  yearly_df <- purrr::imap_dfr(idx_by_year, function(layer_idx, yr) {
    r_year <- r_all[[layer_idx]]
    vals <- values(r_year, mat = TRUE)
    vals <- as.vector(vals)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      tibble(year = as.integer(yr),
             mean = NA_real_, sd = NA_real_,
             median = NA_real_, q25 = NA_real_, q75 = NA_real_)
    } else {
      qs <- quantile(vals, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      tibble(year   = as.integer(yr),
             mean   = mean(vals, na.rm = TRUE),
             sd     = sd(vals, na.rm = TRUE),
             median = median(vals, na.rm = TRUE),
             q25    = qs[1], q75 = qs[2])
    }
  }) %>%
    mutate(folder = basename(folder_name)) %>%
    relocate(folder, year) %>%
    arrange(year)
  
  list(monthly = monthly_df, yearly = yearly_df)
}

# --- run and bind -----------------------------------------------------------
res_list <- map(folders, summarize_folder)

summary_monthly <- bind_rows(lapply(res_list, `[[`, "monthly"))
summary_yearly  <- bind_rows(lapply(res_list, `[[`, "yearly"))

# Clean 'type' column from folder name (keep only harg1 or harg2)
get_type <- function(x) {
  stringr::str_extract(x, "harg[0-9]+")
}

# simplify the different type
summary_monthly <- summary_monthly %>% 
  mutate(type = get_type(folder),
         # keep proper ordering
         date = as.Date(sprintf("%04d-%02d-01", year, as.integer(month))))

summary_yearly <- summary_yearly %>% 
  mutate(type = get_type(folder))




# --- PLOT -------------------------------------------------------------------
#   - Ribbon: IQR (25–75%)
#   - Error bars: mean ± SD
#   - Line/points: annual mean
# Build "even-year" bands (or flip to odd if you prefer)
year_bands <- summary_monthly %>%
  distinct(year) %>%
  transmute(
    year,
    xmin = as.Date(sprintf("%04d-01-01", year)),
    xmax = as.Date(sprintf("%04d-12-31", year))
  ) %>%
  mutate(is_even = (year %% 2) == 0) %>%
  filter(is_even)   # keep only every other year for shading




ggplot(summary_monthly, aes(x = date, y = mean, color = type, fill = type)) +
  # background grey bands per (every other) year
  geom_rect(
    data = year_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey80", alpha = 0.25
  ) +
  #geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_line() +
  geom_point(size = 1.2) +
  geom_hline(yintercept=0, linetype = 'dashed', 'grey') +
  #facet_wrap(~ folder, ncol = 1, scales = "free_y") +
  labs(x = "Year&Month",
       y = "SPEI12 (monthly mean)",
       title = "",
       subtitle = "") +
  theme_bw(base_size = 8)


ggplot(summary_yearly, aes(x = year, y = mean, color = type, fill = type)) +
  #geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_line() +
  geom_point(size = 1.2) +
  geom_hline(yintercept=0, linetype = 'dashed', 'grey') +
  #facet_wrap(~ folder, ncol = 1, scales = "free_y") +
  labs(x = "Year",
       y = "SPEI12 (annual mean)",
       title = "",
       subtitle = "") +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom")



