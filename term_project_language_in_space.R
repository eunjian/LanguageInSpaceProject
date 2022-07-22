library(tmap)
library(terra)
library(gstat)
library(mgcv)
library(tidyverse)
library(sf)
library(spdep)
library(spatialreg)
library(dplyr)
library(modelr)

# load the data
pad_voronoi <- read_sf("pad_voronoi.shp")
germany <- pad_voronoi %>% 
  st_geometry() %>%
  st_union()

pad_mds <- read_csv("pad_mds.csv") %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE")) # 183 obs in total
st_crs(pad_mds) <- 4326

pad_grid <- germany %>% 
  st_bbox() %>%
  st_as_sfc() %>%
  st_make_grid(
    cellsize = c(.05, .05),
    what="centers"
  ) %>% 
  st_as_sf(crs=4326)
st_crs(pad_mds) <- 4326

pad_train = sample_n(pad_mds, 165) # 90% of the data for training
st_crs(pad_train) <- 4326

outside <- sapply(st_intersects(pad_mds, pad_train),function(x){length(x)==0})
pad_test = pad_mds[outside, ] # 10% of the data for testing
st_crs(pad_test) <- 4326


# Inverse distance weighted interpolation
pad.idw <- idw(r ~ 1, location = pad_train, newdata = pad_test, idp=2)

# Ordinary Kriging
pad.v <- pad_train %>% variogram(r ~ 1, ., cloud=F, cutoff=1000) 

myVariogramModel <- vgm(psill=0.28, "Sph", range=1000, nugget=0.02)

pad.vfit <- fit.variogram(pad.v, myVariogramModel, fit.ranges=F)

pad.krige <- krige(r ~ 1, pad_train, pad_test, pad.vfit)

# Generalized Additive Models
gam.fit <- pad_train %>%
  cbind(., st_coordinates(.)) %>%
  select(r, X, Y) %>%
  mgcv::gam(r ~ s(X,Y), data=.) 

gam.prediction <- predict(gam.fit, newdata = data.frame(st_coordinates(pad_test)))

pad.gam <- pad.krige %>%
  cbind(., st_coordinates(.)) %>%
  st_drop_geometry() %>%
  mutate(Z = gam.prediction) %>%
  select(X,Y,Z)
# select(X, Y, Z) %>%
# raster::rasterFromXYZ(crs=4326) %>%
# as("SpatRaster")