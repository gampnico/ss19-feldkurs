library("raster")
library("rgdal")
raster_list <- list.files('~/PycharmProjects/bkp/ss19-feldkurs/Pycode/rasters/static', full.names = TRUE)
master_raster <- raster(raster_list[1])
for(i in 1:length(raster_list)){
  
  # get file name
  file_name <- raster_list[i]
  
  # read raster in
  road_rast_i <- projectRaster(raster(file_name), master_raster, method='ngb')
  
  
  if(i == 1){
    
    combined_raster <- road_rast_i
    
  } else {
    # merge rasters and calc overlap
    # combined_raster <- merge(combined_raster, road_rast_i, 
    #                          fun = function(x, y){sum(x@data@values, y@data@values)})
    combined_raster <- mosaic(combined_raster, road_rast_i, fun = sum)
  }
}
max.raster = combined_raster@data@max
min.raster = combined_raster@data@min

breakpoints <- c(min.raster, 0.001*max.raster, 0.4*max.raster,0.5*max.raster,
                 0.6*max.raster, 0.7*max.raster, 0.8*max.raster, max.raster)
colors <- c("white","lightgoldenrod","lightgoldenrod2","lightgoldenrod3",
            "lightgoldenrod4","darkorange4","firebrick4", "black")
plot(combined_raster,breaks=breakpoints,col=colors)
