#Functions to use when analyzing hyena data

library(circular)
library(dismo)
library(lubridate)
library(gplots)
library(viridis)
library(raster)
library(rgdal)
library(zoo)


#LAT/LON TO UTM CONVERSIONS (AND VICE VERSA)
#Converts a matrix of lons and lats (lons first column, lats second column) to UTM
#Inputs:
#	LonsLats: [N x 2 matrix] of lons (col 1) and lats (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	EastingsCol1: whether eastings should be given in first column of output (default) or not
#Outputs:
#	EastNorths or NorthEasts: [N x 2 matrix] of Eastings and Northings - eastings are first column by default
latlon.to.utm <- function(LonsLats,EastingsCol1 = TRUE,utm.zone='34',southern_hemisphere=TRUE){
	latlons <- data.frame(X=LonsLats[,2],Y=LonsLats[,1])
	non.na.idxs <- which(!is.na(latlons$X) & !is.na(latlons$Y))
	len <- nrow(latlons)
	non.na.latlons <- latlons[non.na.idxs,]
	coordinates(non.na.latlons) <- ~Y + X
	proj4string(non.na.latlons) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')
	if(southern_hemisphere){
		projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
	} else{
		projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
	}
	utm <- spTransform(non.na.latlons,CRS(projection.string))
	EastNorths <- matrix(NA,nrow=len,ncol=2)
	EastNorths[non.na.idxs,] <- utm@coords
	if(!EastingsCol1){
		NorthEasts <- EastNorths[,c(2,1)]
		return(NorthEasts)
	} else{
		return(EastNorths)
	}
}

#Converts a matrix of eastings and northings (eastings first column, northings second column) to UTM
#Inputs:
#	EastNorths: [N x 2 matrix] of eastings (col 1) and northings (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	LonsCol1: whether lons should be given in first column of output (default) or not
#Outputs:
#	LonLats or LatLons: [N x 2 matrix] of longitudes and latitudes - lons are first column by default 
utm.to.latlon <- function(EastNorths,LonsCol1=TRUE,utm.zone = '34',southern_hemisphere=TRUE){
	utms <- data.frame(X=EastNorths[,1],Y=EastNorths[,2])
	non.na.idxs <- which(!is.na(utms$X) & !is.na(utms$Y))
	len <- nrow(utms)
	if(southern_hemisphere){
		projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
	} else{
		projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
	}
	non.na.utms <- SpatialPoints(utms[non.na.idxs,],proj4string=CRS(projection.string))
	lonlat <- spTransform(non.na.utms,CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
	LonLats <- matrix(NA,nrow=len,ncol=2)
	LonLats[non.na.idxs,] <- lonlat@coords
	if(!LonsCol1){
		LatLons <- LonLats[,c(2,1)]
		return(LatLons)
	} else{
		return(LonLats)
	}	
}

#Get state matrix from data frame
#INPUTS:
# gps.vedba.all: data frame of gps and vedba data, created from link_gps_and_vedba.R
# state.to.use: 1 if using single threshold defined states ($state), 2 if using double threshold defined states ($state2)
#OUTPUTS:
# out, an object containing
#   $state.mat: [n.inds x n.times] matrix of states where state.mat[i,t] = state of individual i at time t
#   $times: vector of times corresponding to the indices of state.mat
#   state.used: 1 or 2, directly taken from state.to.use
get.state.mat <- function(gps.vedba.all,state.to.use = 2){
  #create a matrix N x T whose element gives the state of a given hyena
  n <- length(unique(gps.vedba.all$hyena.id))
  times <- seq(min(gps.vedba.all$time),max(gps.vedba.all$time),by=1)
  state.mat <- matrix(NA,nrow=n,ncol=length(times))
  for(i in 1:n){
    dat.i <- gps.vedba.all[which(gps.vedba.all$hyena.id==i),]
    tmin.idx <- which(times == dat.i$time[1])
    if(state.to.use==1){
      state.mat[i,seq(tmin.idx,tmin.idx+nrow(dat.i)-1,1)] <- dat.i$state
    } else{
      state.mat[i,seq(tmin.idx,tmin.idx+nrow(dat.i)-1,1)] <- dat.i$state2
    }
  }
 
  out <- list()
  out$state.mat <- state.mat
  out$times <- times
  out$state.used <- state.to.use
  return(out)
   
}

#Get activity state over a time window by computing how many high activity seconds within that window and testing whether it is over a threshold
#INPUTS:
# state.vec: vector of states over time (0 or 1)
# t.win: time window over which to compute aggregated states (must be odd)
# thresh: threshold number of 1's in order to conside a hyena "active" during a set of seconds
#OUTPUTS:
# agg.state.vec: aggregated state vector, thresholded in the manner described over the time window given
get.aggregated.activity.state <- function(state.vec, t.win, thresh){
  return(rollapply(data=state.vec,width=t.win,FUN=function(x){return(agg.state.helper(x,thresh))},align='center',fill=NA))
}

#helper function for get.aggregated.activity.state - for a vector, returns 1 if mean(vec) > thresh and 0 otherwise
agg.state.helper <- function(vec,thresh){
  return(as.numeric(mean(vec,na.rm=T) > thresh))
}
