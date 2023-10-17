library("dataverse")
require(sf)
library(dplyr)
library(tmap)

# Sys.setenv("DATAVERSE_SERVER" = "dataverse.harvard.edu")

VBG <- get_dataframe_by_name(
  filename = "VBG.gpkg",
  dataset = "10.7910/DVN/IME0M5",
  original = TRUE,
  .f = sf::read_sf)

# not straightforward to load .rda files with this method
#get_dataframe_by_name(
#  filename = "SIG.rda",
#  dataset = "10.7910/DVN/IME0M5",
#  original = TRUE,
#  .f = base::load)

gpstrack <- get_dataframe_by_name(
  filename = "05.gpx",
  dataset = "10.7910/DVN/Y1AQKS",
  original = TRUE,
  .f = sf::read_sf)
