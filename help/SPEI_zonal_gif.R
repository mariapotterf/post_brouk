

# get SPEI
# summarize per okres level

# make gif of all rasters over time
# make gif for summarized SPEIs over okres


library(terra)
library(ggplot2)
library(magick)   # for GIF creation
library(dplyr)
library(stringr)

# --- SETTINGS ---------------------------------------------------------------
folder <- "raw/SPEI12/spei12_harg1_split"   # choose one folder
out_gif_monthly <- "outData/SPEI12_monthly.gif"
out_gif_yearly  <- "outData/SPEI12_yearly.gif"

# --- HELPERS ----------------------------------------------------------------
get_year_month_day <- function(x) {
  m <- str_match(basename(x), "([12][0-9]{3})[-_ ]?(0?[1-9]|1[0-2])[-_ ]?(0?[1-9]|[12][0-9]|3[01])")
  tibble(
    year  = suppressWarnings(as.integer(m[,2])),
    month = suppressWarnings(as.integer(m[,3])),
    day   = suppressWarnings(as.integer(m[,4])),
    file  = x
  )
}

# --- READ FILES -------------------------------------------------------------
files <- list.files(folder, pattern = "\\.tif$|\\.tiff$", full.names = TRUE, ignore.case = TRUE)
stopifnot(length(files) > 0)
files


meta <- get_year_month_day(files) %>% mutate(file = files)
meta <- meta %>% arrange(year, month, day)

# --- (A) MONTHLY ANIMATION --------------------------------------------------
# Caution: heavy if many rasters!
# Use terra::plot() + magick to combine

# If you want to reuse the yearly zlim, pass it via the zlim= argument.
make_gif_monthly <- function(meta, out_file, fps = 6, zlim = NULL) {
  # ensure proper order
  meta <- meta %>% arrange(year, month, day)
  
  # 1) Determine global symmetric z-limits across ALL monthly rasters (if not provided)
  if (is.null(zlim)) {
    mm <- do.call(rbind, lapply(meta$file, function(f) {
      r  <- rast(f)
      gm <- terra::global(r, c("min","max"), na.rm = TRUE)
      as.numeric(gm[1, ])  # c(min, max)
    }))
    lim  <- max(abs(range(mm, na.rm = TRUE)))
    zlim <- c(-lim, lim)
  }
  
  # 2) Diverging palette with exact white at zero
  pal  <- colorRampPalette(c("#B2182B", "#FFFFFF", "#2166AC"))  # red - white - blue
  ncol <- 101  # odd => middle color exactly white
  
  # 3) Build frames
  frames <- lapply(seq_len(nrow(meta)), function(i) {
    r <- rast(meta$file[i])
    
    png_file <- tempfile(fileext = ".png")
    png(png_file, width = 1200, height = 900, res = 150)
    par(mar = c(3,3,3,6))
    plot(r,
         main = sprintf("SPEI12 — %04d-%02d", meta$year[i], meta$month[i]),
         col  = pal(ncol),
         zlim = zlim)
    # mtext("SPEI12 (− = dry … + = wet)", side = 4, line = 3, cex = 0.9)
    dev.off()
    
    image_read(png_file)
  })
  
  gif <- image_animate(image_join(frames), fps = fps, loop = 0)
  image_write(gif, out_file)
  message("Saved monthly GIF: ", out_file, " | zlim: ", paste0(round(zlim, 3), collapse = " to "))
}

# --- (B) YEARLY ANIMATION ---------------------------------------------------
make_gif_yearly <- function(meta, out_file, fps = 2) {
  yrs <- sort(unique(meta$year))
  
  # 1) Compute yearly mean rasters first (to reuse & to scan range)
  yearly_means <- lapply(yrs, function(yr) {
    r <- rast(meta$file[meta$year == yr])
    mean(r, na.rm = TRUE)
  })
  
  # 2) Global symmetric z-limits around 0 (consistent across time)
  all_vals <- unlist(lapply(yearly_means, values))
  all_vals <- all_vals[is.finite(all_vals)]
  lim <- max(abs(range(all_vals, na.rm = TRUE)))
  zlim <- c(-lim, lim)
  
  # 3) Diverging palette with exact white at zero
  pal <- colorRampPalette(c("#B2182B", "#FFFFFF", "#2166AC"))  # red - white - blue
  ncol <- 101  # odd number ensures center color is white
  
  # 4) Draw frames with fixed zlim & palette
  frames <- lapply(seq_along(yrs), function(i) {
    r_m <- yearly_means[[i]]
    png_file <- tempfile(fileext = ".png")
    png(png_file, width = 1200, height = 900, res = 150)
    par(mar = c(3,3,3,6))
    plot(r_m,
         main = sprintf("SPEI12 — %d (yearly mean)", yrs[i]),
         col  = pal(ncol),
         zlim = zlim)
    #mtext("SPEI12 (− = dry … + = wet)", side = 4, line = 3, cex = 0.9)
    dev.off()
    image_read(png_file)
  })
  
  gif <- image_animate(image_join(frames), fps = fps, loop = 0)
  image_write(gif, out_file)
  message("Saved yearly GIF: ", out_file, " | zlim: ", paste0(round(zlim, 3), collapse = " to "))
}


# --- RUN --------------------------------------------------------------------
# monthly version (can be heavy)
make_gif_monthly(meta, out_gif_monthly)

# safer yearly version
make_gif_yearly(meta, out_gif_yearly)
