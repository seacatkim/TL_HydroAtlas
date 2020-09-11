######################################
#### Code written by Felipe Suarez ###
#### Version 3:  14-01-2020       ####
######################################


library(tidyverse)
library(dplyr)
library(rgdal)
library(raster)
library(gdistance)
library(reshape2)
library(sf)
library(rgrass7)
library(parallel)
library(tmap)

# To find nearest non-NA pixels (when pour points are not along the coast):
#install.packages("remotes")
#remotes::install_github("SEEG-Oxford/seegSDM")
library(seegSDM)


#IGNORE LINES 25 - 99 (unless you want to modify the function) 

### change SUM_Risk_1 at line 38,101,119 to SUM_Risk_1 or SUM_Risk_2 or SUM_RISK_0 - SUM_sed_BL, SUM_1_LULC

########################################################################


sed_flow<-function(pp_sf,#list of pour points
                   ocean_raster,
                   buffer_value,#buffer from pour point
                   resolution_raster,#in meters
                   searchdistance #max searching distance to snap the points to the coast
                   
){
  #set initial sediment value
  initial.value<-pp_sf@data$dsprs
  #create buffer around point and clip ocean raster
  point<-st_as_sf(pp_sf)#first need to convert to sf object
  point_buff<-st_buffer(point, buffer_value)#
  r1<-crop(ocean_raster,extent(point_buff),snap = "near")
  
  #create raster to match point location to ocean raster
  #for a new version we can use the function nearest.raster.point, but I had already solved this issue 
  #rasterizing the point
  
  point_1000<-st_buffer(point, 500)#
  #create raster of pour point
  r2  <- raster(ncol=1, nrow=1, xmn=extent(point_1000)[1], xmx=extent(point_1000)[2], ymn=extent(point_1000)[3], ymx=extent(point_1000)[4])
  crs(r2)<-crs(ocean_raster)
  r2[]<-initial.value
  
  #now sum to match rasters of different extent
  
  extend_all =
    function(rasters){
      extent(Reduce(extend,rasters))
    }
  
  sum_all =
    function(rasters, extent){
      re = lapply(rasters, function(r){extend(r,extent)})
      Reduce("+", re)
    }
  
  rasterOptions(tolerance = 1)
  r3 = sum_all(list(r1,r2), extend_all(list(r1,r2)))
  
  #identify if points are on the coast
  
  if(length(which(r3[]>0))< 1){
    snap_point<-nearestLand(coordinates(r2)[],r1,searchdistance)
    pp_sf@coords<-snap_point
    point<-st_as_sf(pp_sf,coords=c("longitude.y","latitude.y"))
    point_1000<-st_buffer(point, 500)#
    #create raster of pour point
    r2  <- raster(ncol=1, nrow=1, xmn=extent(point_1000)[1], xmx=extent(point_1000)[2], ymn=extent(point_1000)[3], ymx=extent(point_1000)[4])
    crs(r2)<-crs(ocean_raster)
    r2[]<-initial.value
    r3 = sum_all(list(r1,r2), extend_all(list(r1,r2)))
  }
  
  #find exact pour point position
  pour_point <- which(r3[] >0)
  
  fromCoords<-pp_sf@coords
  
  raster_temp<-r3
  r<-raster_temp
  
  #create transition layer to be used in cost layer
  transition_layer<-transition(r1,transitionFunction=mean, directions=8)
  
  transition_layer<-geoCorrection(transition_layer,scl=FALSE)
  
  cost_layer<-accCost(transition_layer, fromCoords)
  #calculate distance frequencies
  
  
  initial.value<-pp_sf@data$dsprs
  
  values<-rasterToPoints(cost_layer)[,3]
  frequency_distances<-as.data.frame(table(round(values,1)))
  frequency_distances<-subset(frequency_distances,frequency_distances[1] != Inf)
  frequency_distances[,1]<-as.character(frequency_distances[,1])
  frequency_distances[,1]<-as.numeric(frequency_distances[,1])
  new.values<-initial.value
  raster_temp[pour_point]<-initial.value*0.005  ##Decay parameter
  
  for(i in 2:nrow(frequency_distances))
  {
    cellstoreplace<-which(round(cost_layer[],1)==frequency_distances[i,1])
    percell = initial.value * 0.005  ## decay parameter
    sumcells = percell * length(cellstoreplace)
    remain = initial.value - sumcells
    initial.value = remain
    raster_temp[cellstoreplace]<-percell
    remain < pp_sf@data$dsprs*0.0005
    if(remain < pp_sf@data$dsprs*0.0005) {
      break
    }
  }
  raster_temp[which(raster_temp[]==0)]<-NA
  if(length(unique(raster_temp[]))>1){
    raster_temp<-trim(raster_temp)
  }
  
  return(raster_temp)
}

###########################################################
#APPLY YOU FUNCTION TO POINT LIST########################
#READ OCEANS RASTER
ocean_raster<-raster("data/CT_ocean.tif")  # TLocean_raster.tif
#plot(ocean_raster)

#READ POUR POINTS SHAPEFILE
pour_points<-st_read("data/NTimor_pour_points.shp")   #50R_pp 
plot(pour_points)

pour_points<-as(pour_points,"Spatial")

#Assign unique ID to pour point
pour_points@data$ID<-1:nrow(pour_points@data)
#create list of points
point_list<-split(pour_points,as.character(pour_points@data$ID))

pb <- txtProgressBar(min = 0, max = length(pour_points), style = 3) #To check the progress

test<-list()
###To split process, change "length(however many points)"for i in XX:XX to the range and then change 1 to whatever the starting point is
system.time(for (i in 1:length(point_list)){
  setTxtProgressBar(pb, i) 
  test[[i]]<-  sed_flow(  
    point_list[[i]],#list of pour points
    ocean_raster,
    80000,#buffer from pour point (in meters)-change value to make nutrients go farther 
    1,#in km
    90000#search distance to snap points to the coast
  )
})

newextent<-extent(Reduce(extend,test))

re = lapply(test, function(r){extend(r,newextent)})

for (i in 1:length(re)){
  r<-re[[i]]
 r[is.na(r[])] <- 0 # try commenting out line to set no data to NULL
  re[[i]]<-r
}

sed_plume_merge<-Reduce("+", re)


writeRaster(sed_plume_merge, "output/sed_plume_CT1_Timor.tif", overwrite = TRUE) 
  



