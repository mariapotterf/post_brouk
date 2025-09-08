

# get SPEI
# summarize per okres level

# make gif of all rasters over time
# make gif for summarized SPEIs over okres

# ---- Minimal GIF maker for SPEI rasters (ggplot, fixed scale -2..2) --------
library(terra)
library(dplyr)
library(stringr)
library(ggplot2)
library(magick)
library(fs)
library(scales)

# INPUT: pick ONE folder
spei_folder <- "raw/SPEI12/spei12_harg1_split"
out_dir     <- "outData"
#dir_create(out_dir)



# Parse YYYY_MM_DD from filenames like spei12_2000_01_01.tif
get_meta <- function(folder) {
  files <- list.files(folder, pattern = "\\.(tif|tiff)$",
                      full.names = TRUE, ignore.case = TRUE)
  stopifnot(length(files) > 0)
  m <- str_match(basename(files), "([12][0-9]{3})[-_ ]?(0?[1-9]|1[0-2])[-_ ]?(0?[1-9]|[12][0-9]|3[01])")
  tibble(
    file  = files,
    year  = as.integer(m[,2]),
    month = as.integer(m[,3]),
    day   = as.integer(m[,4])
  ) %>%
    filter(!is.na(year), !is.na(month)) %>%
    arrange(year, month, day)
}

# Terra raster -> data.frame for ggplot
r_to_df <- function(r) {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  val_col <- names(r)[1]
  names(df)[names(df) == val_col] <- "value"
  df
}

# Shared ggplot canvas with fixed scale -2..2
plot_frame <- function(df, title_text) {
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_equal(expand = FALSE) +
    scale_fill_gradient2(
      low = "red", mid = "white", high = "blue",
      limits = c(-3, 3), midpoint = 0,
      oob = scales::squish, name = "SPEI12"
    ) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),  # white background
      panel.background = element_rect(fill = "white", colour = NA), # white panel
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 14,
        margin = margin(b = 8, t = 6)
      ),
      plot.title.position = "panel",
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}


# MONTHLY GIF (one frame per month)
make_gif_monthly_gg <- function(folder, outfile, width = 7, height = 5, dpi = 150, fps = 10) {
  meta <- get_meta(folder)
  frames <- lapply(seq_len(nrow(meta)), function(i) {
    r  <- rast(meta$file[i])
    df <- r_to_df(r)
    p  <- plot_frame(df, sprintf("SPEI12 — %04d-%02d", meta$year[i], meta$month[i]))
    tmp <- tempfile(fileext = ".png")
    ggsave(tmp, p, width = width, height = height, dpi = dpi)
    image_read(tmp)
  })
  gif <- image_animate(image_join(frames), fps = fps, loop = 0)
  image_write(gif, outfile)
  message("Saved monthly GIF: ", outfile)
}

# YEARLY GIF (mean of months per year)
make_gif_yearly_gg <- function(folder, outfile, width = 7, height = 5, dpi = 150, fps = 10) {
  meta <- get_meta(folder)
  years <- unique(meta$year)
  frames <- lapply(years, function(yr) {
    r_year <- rast(meta$file[meta$year == yr])
    r_mean <- mean(r_year, na.rm = TRUE)
    df <- r_to_df(r_mean)
    p  <- plot_frame(df, sprintf("SPEI12 — %d (yearly mean)", yr))
    tmp <- tempfile(fileext = ".png")
    ggsave(tmp, p, width = width, height = height, dpi = dpi)
    image_read(tmp)
  })
  gif <- image_animate(image_join(frames), fps = fps, loop = 0)
  image_write(gif, outfile)
  message("Saved yearly GIF: ", outfile)
}

# ---- RUN (choose what you need) --------------------------------------------
#make_gif_monthly_gg(spei_folder, file.path(out_dir, "SPEI12_monthly_fixedScale.gif"))
make_gif_yearly_gg (spei_folder, file.path(out_dir, "SPEI12_yearly_fixedScale.gif", fps = 5))


