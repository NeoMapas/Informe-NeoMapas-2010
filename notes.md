# Quarto book

```sh
quarto create project book informe_qmd
```

## Import data

Need to change format to make this more download-friendly.

Files in shapefile format can be converted to gpkg with OGR:
```sh
 ogr2ogr VBG.gpkg VBG.shp
```

R object in the old `sp` format, can be written to spatial objects with `rgdal`:

```{r}
library(sp)
library(rgdal) ## also rgeos will be loaded when calling the function 
load("SIG.rda")
writeOGR(CNEB, dsn="CNEB.gpkg", driver="GPKG", layer="CNEB")
```