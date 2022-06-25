####visualization of raster images

#load packages
require(terra); require(ggplot2); require(ggspatial) ; require(data.table)

  # clear environment
  rm(list= ls())

  # set wd if no link to local data is required
  # setwd(paste0("C:/Users/", Sys.info()[["user"]], "/SPRINGG/Gerard Ros - NMI-PROJ/JustPmaps/"))

  #set ggplot theme
  theme_set(theme_bw())

  #homosline projection
  homosoline <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs" 

  # load the poxptot model
  poxptotlm <- readRDS('data/poxptot_lm.RDS')
    
  # load rasters with predicted Al and Fe oxides (made by Maarten)
  alox <- rast("data/alox.tiff")
  feox <- rast("data/feox.tiff")
  both <- rast("data/fe_alox.tiff")
  
  # load crop P uptake in arable land (from IMAGE) and convert to boolean (there is arable crop or not)
  crop <- rast('data/crop_puptake.asc')
  crop[crop==0] <- NA
  crop <- terra::resample(crop,alox,method ='near')
  crop[!is.na(crop)] <- 1
  
  # update alox and feox
  alox2 <- alox * crop
  feox2 <- feox * crop
  alfeox2 <- both * crop
  
  # read in relevant data from IMAGE
  dens <- rast('D:/DATA/07 image/bulkdensity_cr.asc')
  dens <- terra::resample(dens,alox,method ='bilinear')
  ptot <- rast('D:/DATA/07 image/Ptot.asc')
  ptot <- terra::resample(ptot,alox,method ='bilinear')
  
  # stack them
  spr <- c(alox2,feox2,alfeox2,ptot,dens)
  
  # estimate just P levels
  
    # set to data.table
    dt <- as.data.table(as.data.frame(spr, xy=TRUE, na.rm=F))
  
    # estimate pox-ptot fraction given fe_alox in mmol/kg, function trained over full PSI range in 'pox-pot-regression.R'
    dt[, pox_ptot := 0.16352 + 0.24635 * log10(fe_alox)]

    # estiamate current PSI (as fraction (-)
    dt[,psi_current := (Ptot / 31) * pox_ptot / fe_alox]
      
    # set PSI target to 0.1
    dt[,psi_target := 0.1]
    
    # depth of the soil layer
    dt[, depth := 0.3]
    
    # duration to reach the target (years)
    dt[,duration := 35]
    
    # estimate additional P input given PSI target = 0.1
    dt[, preq1 := bulkdensity_cr * 1000 * depth * pmax(0,0.1 - psi_current) * fe_alox * (1/pox_ptot) * 31 * 0.01 / duration]
    
    # estimate additional P input given PSI target = 0.1
    dt[, preq2 := bulkdensity_cr * 1000 * depth * pmax(0,0.15 - psi_current) * fe_alox * (1/pox_ptot) * 31 * 0.01 / duration]
    
    
  
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
