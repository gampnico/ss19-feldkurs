mosaicList <- function(rasList){
  
  #Internal function to make a list of raster objects from list of files.
  ListRasters <- function(list_names) {
    raster_list <- list() # initialise the list of rasters
    for (i in 1:(length(list_names))){ 
      grd_name <- list_names[i] # list_names contains all the names of the images in .grd format
      raster_file <- flip(raster::raster(grd_name),2)
      # raster_file <- flip(raster_file,2)
      
      raster_file <- projectRaster(raster_file, snap, method = "ngb")
    }
    raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
  }
  
  #convert every raster path to a raster object and create list of the results
  raster.list <-sapply(rasList, FUN = ListRasters)
  
  # edit settings of the raster list for use in do.call and mosaic
  names(raster.list) <- NULL
  #####This function deals with overlapping areas
  raster.list$fun <- sum
  #raster.list$tolerance <- 0.1
  
  #run do call to implement mosaic over the list of raster objects.
  mos <- do.call(raster::mosaic, raster.list)
  
  #set crs of output
  crs(mos) <- crs(x = raster(rasList[1]))
  return(mos)
}

raster_files <- list.files(path ="~/PycharmProjects/bkp/ss19-feldkurs/Pycode/rasters/gen_tif",
                           pattern = ".tif$",full.names = TRUE )
# snap <- raster(resolution = c(5,5), xmn = 180000, xmx = 300000, ymn = 60000,
#                ymx = 100000, crs = "+init=epsg:31254 +units=m")
snap <- raster(resolution = c(10,10), xmn = (79219 - 15000), xmx = (79219 + 15000),
               ymn = (237011- 10000), ymx = (237011 + 10000), crs = "+init=epsg:31254 +units=m") 

national_layer <- mosaicList(raster_files)
max.raster = national_layer@data@max
min.raster = national_layer@data@min
layer.mean = mean(national_layer@data@values)
# breakpoints <- c(0.01*min.raster, 0.01*max.raster, 0.3*max.raster, 0.4*max.raster,0.5*max.raster,
#                  0.6*max.raster, 0.7*max.raster, 0.8*max.raster, max.raster)
# colors <- c("white","lightgoldenrod","lightgoldenrod2","lightgoldenrod3",
#             "lightgoldenrod4","darkorange4","firebrick4", "black", "red")


plot(national_layer)
plot(national_layer, col=terrain.colors(10))
levels <- c(0.1*max.raster, 0.2*max.raster, 0.3*max.raster,
            0.4*max.raster, 0.5 * max.raster, 0.6*max.raster,
            0.7*max.raster, 0.8*max.raster, 0.9*max.raster)

levels_2 <- c(0.02*max.raster, 0.04*max.raster, 0.06*max.raster,
              0.05* max.raster, 0.2*max.raster, 0.4*max.raster,
              0.6*max.raster, 0.9 * max.raster)
# levels_3 <- c(layer.mean * 0.4, layer.mean * 0.6, layer.mean * 0.8,
#              layer.mean,layer.mean * 1.2,layer.mean*1.4,layer.mean*1.6)

# contour plot ------------------------------------------------------------
# cont <- rasterToContour(national_layer, maxpixels=1000000000, levels=levels)
cont <- rasterToContour(national_layer, maxpixels=100000000, levels=levels_2)

plot(cont)

