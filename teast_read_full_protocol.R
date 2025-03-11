
# test -----
# export protocol: can i see aliases? and answers? 

library(sf)        # For spatial data handling
library(terra)     # For `vect()` method
library(dplyr)
library(tidyr)
library(data.table)

# read data
#path = "C:/Users/potterf/OneDrive - CZU v Praze\Dokumenty\2025_CZU_postbrouk\example_from_Christian\full_data_for_tablet\slk""
protocol <- vect("raw/data_SVK_protocol.gpkg")
aliases <- vect("raw/data_SVK_protocol_aliases.gpkg")
#aliases <- fread("raw/test_CSV.csv")

head(protocol)
head(aliases)

ncol(protocol)  # 326
ncol(aliases)   # 175

# tehy do not fit: but the questions are enought to export
aliases_df <- as.data.frame(aliases)
class(t(aliases_df))

# Convert vector to data frame
aliases_clean <- data.frame(Alias = colnames(aliases_df))

# export
fwrite(aliases_clean, 'outTable/aliases_names.csv')


