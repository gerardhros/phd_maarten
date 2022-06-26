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
  names(alfeox2) <- 'alfeox'
  
  # read in relevant data from IMAGE
  
    # read density (ton/m3) and convert to kg/m3
    dens <- rast('D:/DATA/07 image/bulkdensity_cr.asc')
    dens <- terra::resample(dens,alox,method ='bilinear')
    dens <- dens * 1000
    names(dens) <- 'density'
  
    # read in total P levels (mg P /kg)
    ptot <- rast('D:/DATA/07 image/Ptot.asc')
    ptot <- terra::resample(ptot,alox,method ='bilinear')
    names(ptot) <- 'ptot'
    
    # read in cropland area (ha)
    croparea <-  rast('D:/DATA/07 image/croparea_2015.asc')
    croparea <- terra::resample(croparea,alox,method ='bilinear')
    names(croparea) <- 'croparea'
    
    # read in crop P uptake in kg P2O5 / year and divide by area to get kg P / ha / yr
    pup <- rast('data/crop_puptake.asc')
    pup <- terra::resample(pup,alox,method ='bilinear')
    pup <- pup * 0.4365 / croparea
    names(pup) <- 'puptake'
    
    # read in P inputs in 2015 (kg P / year) and divide by area to get kg P/ ha / year
    pinput <- rast('D:/DATA/07 image/inputs_2015.asc')
    pinput <- terra::resample(pinput,alox,method ='bilinear')
    pinput <- pinput * crop
    # pinput <- pinput / croparea => it seemed unit is kg P / ha rather than kg P / year
    names(pinput) <- 'pinput'
    
    # read in P surplus (kg P2O5 per year) and divide by cropland (ha) and P-P2O5 converstion to get kg P / ha /yr
    psurplus <- rast('D:/DATA/07 image/balance_p_arable.asc')
    psurplus <- terra::resample(psurplus,alox,method ='bilinear')
    psurplus <- psurplus * 0.4365 / croparea
    names(psurplus) <- 'psurplus'
    
    # read in lena stack 
    r.lena1 <- rast('../02 data/Fig 4a uptake_P_arable_current_kgP_perha.tif')
    r.lena2 <- rast('../02 data/Fig 4b uptake_P_arable_target_kgP_perha.tif')
    r.lena3 <- rast('../02 data/Fig 4c uptake_P_arable_finaltarget_kgP_perha.tif')
    names(r.lena1) <- 'pup1'  
    names(r.lena2) <- 'pup2' 
    names(r.lena3) <- 'pup3' 
    r.lena <- c(r.lena1,r.lena2,r.lena3)
    r.lena <- terra::resample(r.lena,alox,method='bilinear')
    
    
  # stack them
  dt.stack <- c(croparea,alox2,feox2,alfeox2,ptot,dens,pup,pinput,psurplus,r.lena)
  
  # estimate just P levels
  
    # set to data.table
    dt <- as.data.table(as.data.frame(dt.stack, xy=TRUE, na.rm=F))
  
    # estimate pox-ptot fraction given fe_alox in mmol/kg, function trained over full PSI range in 'pox-pot-regression.R'
    dt[, pox_ptot := 0.16352 + 0.24635 * log10(alfeox)]

    # estiamate current PSI (as fraction (-)
    dt[,psi_current := (ptot / 31) * pox_ptot / alfeox]
      
    # depth of the soil layer
    dt[, depth := 0.3]
    
    # duration to reach the target (years)
    dt[,duration := 35]
    
    # estimate additional P input (kg P / ha) given PSI target = 0.1
    dt[, preq1 := density * depth * pmax(0,0.1 - psi_current) * alfeox * (1/pox_ptot) * 31 * 0.01 / duration]
    
    # estimate additional P input  (kg P / ha) given PSI target = 0.1
    dt[, preq2 := density * depth * pmax(0,0.15 - psi_current) * alfeox * (1/pox_ptot) * 31 * 0.01 / duration]
    
    # get just P level for P equilibrium fertilization (kg P / ha)
    dt[, pfert1 := puptake + preq1]
    dt[, pfert2 := puptake + preq2]
  
    # difference with current P surplus (kg P ha-1)
    dt[,pdiff1 := preq1 - pinput]
    dt[,pdiff2 := preq2 - pinput]
    
    # total P input needed
    dt[, req_ptotal1 := croparea * pfert1]
    dt[, req_ptotal2 := croparea * pfert2]
    dt[, req_ptotal0 := croparea * pinput]
    
    # long-term plus mid-term P input
    dt[, pfert1lt := pup3 + preq1]
    dt[, pfert2lt := pup3 + preq2]
    dt[, req_ptotal3 := croparea * pfert1lt]
    dt[, req_ptotal4 := croparea * pfert2lt]
    dt[, req_preq1 := croparea * preq1]
    dt[, req_preq2 := croparea * preq2]
    
  # make raster
  r.fin <- terra::rast(dt,type='xyz')
  terra::crs(r.fin) <- 'epsg:4326'
  terra::writeRaster(r.fin,'products/justplevels.tif', overwrite = TRUE)
  
  # get base world map
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # create plot for just P level for PSI target = 0.1
  p1 <- visualize(raster = r.fin, 
                 layer = 'preq1', 
                 name = "P dose\n(kg P / ha)", 
                 breaks = c(-1000,5,50,100,500),
                 labels = c('<5','5-50','50-100','>100'),
                 ftitle = 'Mid-term required P surplus to reach PSI target of 10%')
  ggsave(filename = "products/preq_target1.png", 
         plot = p1, width = 32, height = 24, units = c("cm"), dpi = 1200)
  plot.pfert1 <- visualize(raster = r.fin, 
                           layer = 'pfert1', 
                           name = "P dose\n(kg P / ha)", 
                           breaks = c(-1000,5,50,100,500),
                           labels = c('<5','5-50','50-100','>100'),
                           ftitle = 'just P level with PSI target of 10%')
  ggsave(filename = "products/pfert_target1.png", 
         plot = plot.pfert1, width = 32, height = 24, units = c("cm"), dpi = 1200)
  plot.pfert1lt <- visualize(raster = r.fin, 
                           layer = 'pfert1lt', 
                           name = "P dose\n(kg P / ha)", 
                           breaks = c(-1000,5,50,100,500),
                           labels = c('<5','5-50','50-100','>100'),
                           ftitle = 'just P level for PSI target of 10% and corr. target P uptake from yield gap')
  ggsave(filename = "products/pfert_target1lt.png", 
         plot = plot.pfert1lt, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  # create plot for just P level for PSI target = 0.15
  p2 <- visualize(raster = r.fin, 
                  layer = 'preq2', 
                  name = "P dose\n(kg P / ha)", 
                  breaks = c(-1000,5,50,100,500),
                  labels = c('<5','5-50','50-100','>100'),
                  ftitle = 'Mid-term required P surplus to reach PSI target of 15%')
  ggsave(filename = "products/preq_target2.png", 
         plot = p2, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  plot.pfert2 <- visualize(raster = r.fin, 
                           layer = 'pfert2', 
                           name = "P dose\n(kg P / ha)", 
                           breaks = c(-1000,5,50,100,1000),
                           labels = c('<5','5-50','50-100','>100'),
                           ftitle = 'just P level with PSI target of 15%')
  ggsave(filename = "products/pfert_target2.png", 
         plot = plot.pfert2, width = 32, height = 24, units = c("cm"), dpi = 1200)
  plot.pfert2lt <- visualize(raster = r.fin, 
                             layer = 'pfert2lt', 
                             name = "P dose\n(kg P / ha)", 
                             breaks = c(-1000,5,50,100,500),
                             labels = c('<5','5-50','50-100','>100'),
                             ftitle = 'just P level for PSI target of 15% and corr. target P uptake from yield gap')
  ggsave(filename = "products/pfert_target2lt.png", 
         plot = plot.pfert2lt, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  
  
  # create plot for current P surplus
  p3 <- visualize(raster = r.fin, 
                 layer = 'psurplus', 
                 name = "P surplus\n(kg P / ha)", 
                 breaks = c(-1000,-50,-5,5,50,1000),
                 labels = c('< -50','-50 to -5','-5 to +5','5 to 50', '>50'),
                 ftitle = 'P surplus on cropland')
  ggsave(filename = "products/psurplus.png", 
         plot = p3, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  
  # create additive P input for soil fertility improvement
  p4 <- visualize(raster = r.fin, 
                  layer = 'preq1', 
                  name = "P input\n(kg P / ha)", 
                  breaks = c(-1000,5,50,100,1000),
                  labels = c('<5','5-50','50-100','>100'),
                  ftitle = 'Extra P needed for soil fertility improvement (target PSI 10%)')
  ggsave(filename = "products/pinput_sf_target1.png", 
         plot = p4, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  
  # create additive P input for soil fertility improvement
  p5 <- visualize(raster = r.fin, 
                  layer = 'preq2', 
                  name = "P input\n(kg P / ha)", 
                  breaks = c(-1000,5,50,100,1000),
                  labels = c('<5','5-50','50-100','>100'),
                  ftitle = 'Extra P needed for soil fertility improvement (target PSI 15%)')
  ggsave(filename = "products/pinput_sf_target2.png", 
         plot = p5, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  
  # create map for current input 
  p6 <- visualize(raster = r.fin, 
                  layer = 'pinput', 
                  name = "P input\n(kg P / ha)", 
                  breaks = c(-1000,5,50,100,1000),
                  labels = c('<5','5-50','50-100','>100'),
                  ftitle = 'Current P inputs')
  ggsave(filename = "products/pinput_current.png", 
         plot = p6, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  # create map for difference between target and current P surplus
  p7 <- visualize(raster = r.fin, 
                  layer = 'pdiff1', 
                  name = "delta P surplus\n(kg P / ha)", 
                  breaks = c(-1000,-25,0,25,100,1000),
                  labels = c('< -25','-25 to 0','0 to 25','25 to 100','>100'),
                  ftitle = 'Change between target (10% PSI) and current P surplus')
  ggsave(filename = "products/psp_difference_target1.png", 
         plot = p7, width = 32, height = 24, units = c("cm"), dpi = 1200)
  p8 <- visualize(raster = r.fin, 
                  layer = 'pdiff2', 
                  name = "delta P surplus\n(kg P / ha)", 
                  breaks = c(-1000,-25,0,25,100,1000),
                  labels = c('< -25','-25 to 0','0 to 25','25 to 100','>100'),
                  ftitle = 'Change between target (15% PSI) and current P surplus')
  ggsave(filename = "products/psp_difference_target2.png", 
         plot = p8, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  # maps with soil data
  p9 <- visualize(raster = r.fin, 
                  layer = 'alfeox', 
                  name = "Al and Fe oxides\n(mmol / kg)", 
                  breaks = c(-1000,25,75,125,200,1000),
                  labels = c('<25','25-75','75-125','125-200','>200'),
                  ftitle = 'Aluminium and iron oxides')
  ggsave(filename = "products/alfeox.png", 
         plot = p9, width = 32, height = 24, units = c("cm"), dpi = 1200)
  p10 <- visualize(raster = r.fin, 
                  layer = 'psi_current', 
                  name = "PSI (-)", 
                  breaks = c(-1,0.05,0.10,0.15,0.5,5),
                  labels = c('<0.05','0.05-0.10','0.10-0.15','0.15-0.50','>0.5'),
                  ftitle = 'Phosphate Saturation Index (-)')
  ggsave(filename = "products/psi_current.png", 
         plot = p10, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  p11 <- visualize(raster = r.fin, 
                   layer = 'ptot', 
                   name = "total P/(mg/kg)", 
                   breaks = c(0,100,250,750,1500,10000),
                   labels = c('<100','100-250','250-750','750-1500','>1500'),
                   ftitle = 'Total Phosphorus in Arable Soils')
  ggsave(filename = "products/ptot_soil.png", 
         plot = p11, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
 

  pa1 <- visualize(raster = r.lena, 
                   layer = 'pup1', 
                   name = "P uptake\n(kg P / ha)", 
                   breaks = c(-100,5,10,20,30,40,50,1000),
                   labels = c('<5','5-10','10-20','20-30','30-40','40-50','>50'),
                   ftitle = 'Current P uptake (2015)')
  ggsave(filename = "products/pup1.png", 
         plot = pa1, width = 26, height = 14, units = c("cm"), dpi = 1200)
  pa2 <- visualize(raster = r.lena, 
                   layer = 'pup2', 
                   name = "Target P uptake\n(kg P / ha)", 
                   breaks = c(-100,5,10,20,30,40,50,1000),
                   labels = c('<5','5-10','10-20','20-30','30-40','40-50','>50'),
                   ftitle = 'Target P uptake based on yield gap')
  ggsave(filename = "products/pup2.png", 
         plot = pa2, width = 26, height = 14, units = c("cm"), dpi = 1200)
  pa3 <- visualize(raster = r.lena, 
                   layer = 'pup3', 
                   name = "Target P uptake\n(kg P / ha)", 
                   breaks = c(-100,5,10,20,30,40,50,1000),
                   labels = c('<5','5-10','10-20','20-30','30-40','40-50','>50'),
                   ftitle = 'Corrected target P uptake based on yield gap')
  ggsave(filename = "products/pup3.png", 
         plot = pa3, width = 26, height = 14, units = c("cm"), dpi = 1200)
  
      
# --- PlOT FIGURES IN HOMOSOLINE PROJECTION --------
# https://wilkelab.org/practicalgg/articles/goode.html
  
  # create plot for just P level for PSI target = 0.1
  p1 <- visualize(raster = r.fin, 
                 layer = 'pfert1', 
                 name = "P dose\n(kg P / ha)", 
                 breaks = c(-1000,5,50,100,500),
                 labels = c('<5','5-50','50-100','>100'),
                 ftitle = 'just P level with PSI target of 10%',
                 homo = TRUE)
  
  ggsave(filename = "products/pfert_target1.png", 
         plot = plot.pfert1, width = 32, height = 24, units = c("cm"), dpi = 1200)
  
  
  
  
  #visualize
  visualize <- function(raster, layer, name, breaks, labels, ftitle, homo = FALSE){
   
    # select the correct layer
    raster.int <- raster[layer]
    
    # convert rast to homosline
    if(homo == TRUE){
      raster.int <- terra::project(raster.int, homosoline)
      world2 <- sf::st_transform(world,homosoline)
      world2 <- world2[world2$region_wb != 'Antarctica',]
      world2 <- world2[world2$name_sort != 'Greenland',]
      plotcrs <- coord_sf(crs = homosoline, lims_method = "box")
    } else {
      world2 <- world
      plotcrs <- coord_sf(crs = 4326, lims_method = "box")
    }
    
    #raster to xy
    df <- as.data.frame(raster.int, xy = TRUE)
    
    #colnames
    colnames(df) <- c("x", "y", "variable")
    
    #plot
    ggplot() + 
      geom_sf(data = world2, color = "black", fill = "white",show.legend = FALSE) +
      geom_tile(data = df, aes(x = x, y = y,fill = cut(variable, breaks,labels = labels))) +
      plotcrs +
      scale_fill_viridis_d() +
      xlab("") + ylab("")+
      #xlab("Longitude") + ylab("Latitude") +
      labs(fill = name) + 
      theme(text = element_text(size = 12), 
            legend.position = c(0.1,0.4),
            legend.background = element_rect(fill = "white",color='white'),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle(ftitle)
  }
  
  
  
  
  
  
  
  
# make plots
visualize <- function(raster, layer,name, breaks, labels = NULL, homo = FALSE){
  
  # select the correct layer
  raster.int <- raster[layer]
  
  # convert rast to homosline
  if(homo == TRUE){
    raster.int <- terra::project(raster.int, homosoline)
    world <- sf::st_transform(world,homosoline)
  }
  
  # raster to xy
  df <- as.data.frame(raster.int, xy = TRUE)
  
  # colnames
  colnames(df) <- c("x", "y", "variable")
  
  # add categories for discrete color scale
  if(length(breaks)>1){
    df$pclass <- cut(df$variable,breaks = breaks, labels = labels,na.rm=TRUE)
    df$pclass <- factor(df$pclass,levels = labels)
    df$pclass2 <- as.numeric(df$pclass)
  } else {
    df$pclass <- cut(df$variable, breaks)
  }
  
  # scale_fill_brewer(palette = "Oranges") +
  # plot
  ggplot(data = world) + 
    geom_sf(color = "black", fill = "white",show.legend = FALSE) +
    geom_tile(data = df, aes(x = x, y = y, fill = NA, color = "black"),show.legend = FALSE) + 
    geom_tile(data = df, aes(x = x, y = y, fill = variable)) +
    scale_fill_viridis_c()+ theme_void()
  
    coord_sf(crs = 4326, lims_method = "box") +
    labs(fill = name) + 
    theme(text = element_text(size = 12), 
          legend.position = c(0.1,0.4),
          legend.background = element_rect(fill = "white",color='white')) +
    ggtitle(titel)
}

#create plots
plot.pfert1 <- visualize(raster = r.fin, 
                         layer = 'pfert1', 
                         name = "P dose\n(kg P/ha)", 
                         breaks = c(-1000,5,50,100,500),
                         labels = c('<5','5-50','50-100','>100'))
ggsave(filename = "03 output/plots/alox.png", plot = plot.alox, width = 40, height = 30, units = c("cm"), dpi = 1200)

plot.feox <- visualize(feox, name = "feox\n(mmol/kg)", breaks = 10)
plot.both <- visualize(both, name = "fe+alox\n(mmol/kg)", breaks = 10)

#save to output
if(!dir.exists("03 output/plots")){
  dir.create("03 output/plots")
}
ggsave(filename = "03 output/plots/alox.png", plot = plot.alox, width = 40, height = 30, units = c("cm"), dpi = 1200)
ggsave(filename = "03 output/plots/feox.png", plot = plot.feox, width = 40, height = 30, units = c("cm"), dpi = 1200)
ggsave(filename = "03 output/plots/al_feox.png", plot = plot.both, width = 40, height = 30, units = c("cm"), dpi = 1200)
