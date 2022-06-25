#####convert ISRIC rasters to one 0-30cm rasterstack

#load packages
library(terra) 
library(data.table)

#create ref to directory 
dir <- paste0("C:/Users/", Sys.info()[["user"]], "/SPRINGG/Gerard Ros - NMI-PROJ/JustPmaps/01 data")

#create folder to save rasterstack to
if(!dir.exists(paste0(dir, "/isric rasterstack"))){
  dir.create(paste0(dir, "/isric rasterstack"))
}

#create data.table of filepaths
rasterfiles <- data.table(raster = list.files(dir, pattern = ".tif$"), path = list.files(dir, full.names = T, pattern = ".tif$"))
rasterfiles <- rasterfiles[grepl("^isric", raster)]
rasterfiles[, parameter := sapply(strsplit(raster, "_"), `[`, 2)] #variable

#create 0-30cm rasterfiles per variable
rasters.depthweighed <- list()
for(i in unique(rasterfiles$parameter)){
  print(i)
  #subset rasterfiles
  tmp <- rasterfiles[parameter == i]
  
  #load as stack
  rel.rasters <- rast(tmp$path)
  
  #change names
  names(rel.rasters) <- gsub("isric_|mean_|soc|bdod|cec|clay|phw|silt|sand|ntot", "", names(rel.rasters))
  names(rel.rasters) <-  gsub("^_", "", names(rel.rasters))
  
  #create 0-30cm raster
  rel.rasters$weighed <- (5 * rel.rasters$`0_5` + 10 * rel.rasters$`5_15` + 15 * rel.rasters$`15_30`) / 30
  
  #select depth weighed raster and change names
  rel.rasters <- rel.rasters$weighed
  names(rel.rasters) <- paste0(i, "_mean_0-30")  
  
  #bind to list
  rasters.depthweighed[[i]] <- rel.rasters
}

#check extents and resolution
lapply(rasters.depthweighed, terra::ext)
lapply(rasters.depthweighed, terra::res)
  #resolution and extents of SOC and BDOD is slightly different, resampling required

#1) resample bdod and soc rasters
rasters.depthweighed$bdod <- terra::resample(rasters.depthweighed$bdod, rasters.depthweighed$cec)
rasters.depthweighed$soc <- terra::resample(rasters.depthweighed$soc, rasters.depthweighed$cec)

#2) stack rasters (from list)
stacked <- terra::rast(rasters.depthweighed)

#3) save rasterstack
terra::writeRaster(stacked, filename = paste0(dir, "/isric rasterstack/isric_stack_mean_0_30.tif"))
