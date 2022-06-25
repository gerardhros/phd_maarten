#set working directory
setwd(paste0("C:/Users/", Sys.info()[["user"]], "/SPRINGG/Gerard Ros - NMI-PROJ/JustPmaps/"))

#load libraries
library(ranger)
library(terra)
library(ggplot2)
library(ggspatial)

#load models
alox <- readRDS("02 models/psc/alox/final/rf.model.RDS")
feox <- readRDS("02 models/psc/feox/final/rf.model.RDS")

#load isric rasterstack
isric <- rast("01 data/isric rasterstack/isric_stack_mean_0_30.tif")

#set names of isric to correspond to rf model
names(isric) <- gsub("phw", "phh2o", names(isric))

#make xy dataframe
isric.df <- as.data.frame(isric, na.rm = TRUE, xy = TRUE)

#remove values of 0
isric.df <- subset(isric.df, isric.df$bdod > 0)

#predict
isric.df$alox <- predict(alox, isric.df)$predictions
isric.df$feox <- predict(feox, isric.df)$predictions
isric.df$fe_alox <- isric.df$feox + isric.df$alox

#convert to raster again
isric.rast <- rast(isric.df)

#set crs to crs of raster
crs(isric.rast) <- crs(isric)

#save
terra::writeRaster(isric.rast$alox, filename = "03 output/alox.tiff")
terra::writeRaster(isric.rast$feox, filename = "03 output/feox.tiff")
terra::writeRaster(isric.rast$fe_alox, filename = "03 output/fe_alox.tiff")

