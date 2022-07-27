create_raster_from_raster = function(oriRas,  sampleR){
  library(sp)
  library(raster)
  
  # raste to points
  oriRasXY = rasterToPoints(oriRas)
  
  # assign numeric values for NA
  oriRasXY[is.na(oriRasXY)] = -9999
  
  # first two column and lon(x) and lat (y) coordinates
  oriRasXY = data.frame(oriRasXY)
  colnames(oriRasXY)[1:3]  = c('X','Y','Val')
  coordinates(oriRasXY) =  ~ X +Y
  
  r = rasterize(oriRasXY,sampleR, field= 'Val')
  r[r == -9999] = NA
  
  return(r)
}

calc_sat_stastics_metrics_single = function(obsdat, satdat){
  # remove na
  dt = data.frame(obsdat, satdat)
  dt = na.omit(dt)
  obsdat = dt[,1]; satdat = dt[,2]
  
  evalMetric = rep(NA,3) # CC, RB, RMSE
  # statistical metrics/ quantitative metrics
  evalMetric[1] = round(cor(satdat, obsdat),3)  # CC
  evalMetric[2] = round(mean(satdat)/mean(obsdat)-1,2)  #RB
  evalMetric[3] = round(sqrt(sum(1/(length(satdat))*(satdat - obsdat)^2)),1)  #RMSE
  
  # evalMetric - CC, BIAS, RMSE
  return(evalMetric) 
}

create_raster_from_datXY = function(datXY, fieldName, sampleR){
  library(sp)
  library(raster)
  
  # assign numeric values for NA
  datXY[is.na(datXY)] = -9999
  
  # first two column and lon(x) and lat (y) coordinates
  colnames(datXY)[1:2]  = c('X','Y')
  coordinates(datXY) =  ~ X +Y
  
  r = rasterize(datXY,sampleR, field= fieldName)
  r[r == -9999] = NA
  
  return(r)
}


t_col <- function(color, percent = 70, name = NULL) {
  # create transparent colors
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255,
               alpha = (100 - percent) * 255 / 100,names = name)
  invisible(t.col)
}

deg2rad = function(deg){
  ### Convert degrees to radians
  ### source: http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
  return(deg*pi/180)
}

gcd = function(long1, lat1, long2, lat2) {
  ### Calculates the geodesic distance between two points specified by
  ### degree latitude/longitude using the Spherical Law of Cosines (
  ### based on : http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
  long1 = deg2rad(long1)
  long2 = deg2rad(long2)
  lat1 = deg2rad(lat1)
  lat2 = deg2rad(lat2)    
  R = 6371 # Earth mean radius [km]
  d = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)
  if (d-1 == .Machine$double.eps) d = 1
  d = acos(d) * R
  return(d) # Distance [km]
}

gcd_matrix = function(x.lat,x.lon,y.lat,y.lon){
  ### Calculates the geodesic distance matrix
  if(length(x.lat)!=length(x.lon) | length(y.lat)!=length(y.lon))
    stop("x or y vectors do not match\n")    
  station.dist = matrix(NA,nrow=length(x.lat),ncol=length(y.lat))
  ## station.dist = matrix(NA,nrow=length(y.lat),ncol=length(x.lat))
  for(xi in 1:nrow(station.dist)){
    for(yi in 1:ncol(station.dist)){
      station.dist[xi,yi] = gcd(long1=y.lon[yi],
                                lat1=y.lat[yi],
                                long2=x.lon[xi],
                                lat2=x.lat[xi])
    }
  }
  return(station.dist)
}
