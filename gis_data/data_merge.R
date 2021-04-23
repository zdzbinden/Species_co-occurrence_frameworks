
############################
### Reading spatial data ###
############################
setwd("")
# reading HydroRIVERS shapefile into R (modified this file in ArcGIS)
library(raster)
hydroRiv <- shapefile("./rivers.shp", verbose=TRUE, warnPRJ=TRUE)
plot(hydroRiv) # check map looks OK
# reading fish collection site coordinates N=142 (this file was also projected in ArcGIS)
sites <- shapefile("./site_shapefiles/sites.shp")
points(sites) # check map looks OK
### snapping points to river lines ###
library(sf)
sites<-st_as_sf(sites)
hydroRiv <-st_as_sf(hydroRiv)
join<-st_join(sites, hydroRiv, st_nearest_feature,left=TRUE) # combined in situ ENV variables with hydroATLAS variables
join <-st_set_geometry(join, NULL) #325 vars
rm(sites)

write.csv(join, "./DATA.csv")
