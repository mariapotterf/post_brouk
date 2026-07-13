

# Make a barplot for Karim

library(tidyverse)

df <- read_csv("C:/Users/potterf/OneDrive - CZU v Praze/Dokumenty/2026_Karim/data/workflow_steps.csv")

df <- df %>%
  mutate(category = str_trim(iconv(category, "UTF-8", "ASCII", sub = "")))

# updated color map — A to G progressively
color_map <- c(
  "A. Platform Initialization"              = "#4A7FA0",
  "B. Input Data Retrieval and Preprocessing" = "#6BAED6",
  "C. Feature Selection and Data Assembly"  = "#AED6E8",
  "D. Patch Extraction (U-Net Only)"        = "#F2A97E",
  "E. Data Splitting and Augmentation"      = "#E8C87A",
  "F. Model Training"                       = "#89B4C9",
  "G. Evaluation and Output"                = "#C4903A")

# factor order: A at bottom, G at top of stack
df$category <- factor(df$category, levels = c(
  "G. Evaluation and Output",
  "F. Model Training",
  "E. Data Splitting and Augmentation",
  "D. Patch Extraction (U-Net Only)",
  "C. Feature Selection and Data Assembly",
  "B. Input Data Retrieval and Preprocessing",
  "A. Platform Initialization"))

# configurations ordered lowest to highest total steps
df$configuration <- factor(df$configuration, levels = c(
  "AEF+RF", "AEF+U-Net", "Sentinel+RF", "Sentinel+U-Net"))

ggplot(df, aes(x = configuration,
               y = steps,
               fill = category)) +
  geom_bar(stat = "identity", col = "black") +
  scale_fill_manual(values = color_map) +
  labs(x = NULL, y = "Number of steps", fill = NULL) +
  theme_classic() +
  theme(legend.position = "right")

