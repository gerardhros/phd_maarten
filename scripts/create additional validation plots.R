#load packages
library(terra) 
library(data.table)
library(sf)
library(ranger)
library(ggplot2)

#create ref to directory 
setwd(paste0("C:/Users/", Sys.info()[["user"]], "/SPRINGG/Gerard Ros - NMI-PROJ/JustPmaps/"))

#load maps
feox <- rast("03 output/feox.tiff")
alox <- rast("03 output/alox.tiff")

#read validation sets and models
alox.rf <- readRDS("02 models/psc/alox/final/rf.model.RDS") #model
alox.val <- fread("02 models/psc/alox/final/data/validation.csv") #validation (for ML, not map validation)

feox.rf <- readRDS("02 models/psc/feox/final/rf.model.RDS") #model
feox.val <- fread("02 models/psc/feox/final/data/validation.csv") #validation (for ML, not map validation)

#extract fe+alox from maps (0.5 x 0.5 degrees resolution)
feox.val$map_feox <- terra::extract(feox, feox.val[, c("lon", "lat")])[["feox"]]
alox.val$map_alox <- terra::extract(alox, alox.val[, c("lon", "lat")])[["alox"]]

#predict values with ML models
alox.val$ranger_alox <- predict(alox.rf, alox.val)[["predictions"]]
feox.val$ranger_feox <- predict(feox.rf, feox.val)[["predictions"]]

#create plots
plotfun <- function(dt, name){
  #set names
  setnames(dt, old = name, new = "measured")

  #pivot longer, plot
  dt.long <- tidyr::pivot_longer(dt, cols = c(paste0("map_", name), paste0("ranger_", name)))
  ggplot(dt.long, aes(x = measured, y = value)) + facet_wrap(~name, scales = "fixed") +
    geom_point(alpha = 0.1) + theme_bw() + geom_abline(slope = 1, intercept = 0) +
    theme(text = element_text(size = 15))
}
plotfun(alox.val, "alox")
plotfun(feox.val, "feox")
