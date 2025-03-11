
# test -----
# export protocol: can i see aliases? and answers? 

library(sf)        # For spatial data handling
library(terra)     # For `vect()` method
library(RSQLite)   # For accessing GPKG metadata

# read metadata:

# Connect to the GPKG file as a SQLite database
con <- dbConnect(RSQLite::SQLite(), "raw/data_SVK_protocol.gpkg")

# List all tables in the GPKG (find the correct metadata table)
dbListTables(con)

dbGetQuery(con, "PRAGMA table_info('rap_ass_v2_a021dd1f_4cf0_4c71_ae0a_6cf6431fda28')")

# Check for alias data
aliases_data <- dbGetQuery(con, "SELECT * FROM log_added_attrs")
print(aliases_data)

search_aliases <- function(con) {
  tables <- dbListTables(con)
  results <- list()
  
  for (table in tables) {
    query <- paste("SELECT * FROM", table)
    try({
      data <- dbGetQuery(con, query)
      if (any(grepl("alias|description|comment", names(data), ignore.case = TRUE))) {
        results[[table]] <- data
      }
    }, silent = TRUE)
  }
  
  results
}

alias_data <- search_aliases(con)
print(alias_data)


# read data
#path = "C:/Users/potterf/OneDrive - CZU v Praze\Dokumenty\2025_CZU_postbrouk\example_from_Christian\full_data_for_tablet\slk""
protocol <- vect("raw/data_SVK_protocol.gpkg")

head(protocol)
