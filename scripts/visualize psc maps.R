####visualization of raster images

#load packages
library(terra) 
library(ggplot2)
library(ggspatial)

#set wd
setwd(paste0("C:/Users/", Sys.info()[["user"]], "/SPRINGG/Gerard Ros - NMI-PROJ/JustPmaps/"))

#set ggplot theme
theme_set(theme_bw())

#homosline projection
homosoline <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs" #homosoline projection

#load rasters
alox <- rast("03 output/alox.tiff")
feox <- rast("03 output/feox.tiff")
both <- rast("03 output/fe_alox.tiff")

#visualize
visualize <- function(raster, name, breaks){
  #convert rast to homosline
  rasterfile <- project(raster, homosoline)
  
  #raster to xy
  df <- as.data.frame(rasterfile, xy = TRUE)
  
  #colnames
  colnames(df) <- c("x", "y", "variable")
  
  #plot
  ggplot(df, aes(x = x, y = y))  + geom_tile(fill = NA, col = "black") + geom_tile(aes(fill = cut(variable, breaks))) +
    coord_sf(crs = homosoline, lims_method = "box") +
    scale_fill_brewer(palette = "Oranges") +
    labs(fill = name) + theme(text = element_text(size = 15), legend.position = c(0.1,0.3))
}

#create plots
plot.alox <- visualize(alox, name = "alox\n(mmol/kg)", breaks = 10)
plot.feox <- visualize(feox, name = "feox\n(mmol/kg)", breaks = 10)
plot.both <- visualize(both, name = "fe+alox\n(mmol/kg)", breaks = 10)

#save to output
if(!dir.exists("03 output/plots")){
  dir.create("03 output/plots")
}
ggsave(filename = "03 output/plots/alox.png", plot = plot.alox, width = 40, height = 30, units = c("cm"), dpi = 1200)
ggsave(filename = "03 output/plots/feox.png", plot = plot.feox, width = 40, height = 30, units = c("cm"), dpi = 1200)
ggsave(filename = "03 output/plots/al_feox.png", plot = plot.both, width = 40, height = 30, units = c("cm"), dpi = 1200)
