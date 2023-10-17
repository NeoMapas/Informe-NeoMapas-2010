## change format to make this more download-friendly
## ogr2ogr VBG.gpkg VBG.shp

library(sp)
library(rgdal) ## also rgeos will be loaded when calling the function 
load("SIG.rda")
writeOGR(CNEB, dsn="CNEB.gpkg", driver="GPKG", layer="CNEB")