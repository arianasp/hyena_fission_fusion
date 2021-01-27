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

#FUNCTIONS
#' Computes the first passage time out of a circle of radius R based on x and y vectors
#' Assumes that the initial point from which to compute is (x[1],y[1])
# @param x numeric vector of the x coordinates of the trajectory
# @param y numeric vector of the y coordinates of the trajectory
# @param R numeric value specifying the radius to use
# @returns the index of the first element in the (x,y) vectors at which R is exceeded
first.passage.time <- function(x, y, R){
  
  #get initial point
  x0 <- x[1]
  y0 <- y[1]
  
  #if initial point is NA, return NA
  if(is.na(x0) | is.na(y0)){
    return(NA)
  }
  
  #throw error if x and y are not the same length
  if(length(x) != length(y)){
    stop('x and y not the same length')
  }
  
  #get distance to initial point over time
  dist <- sqrt((x - x0)^2 + (y - y0)^2)
  
  #get the distances that are greater than R from the starting point
  passed <- which(dist >= R)
  
  #if no passages found, return NA, otherwise return the index of the first passage
  if(length(passed)==0){
    return(NA)
  } else{
    return(min(passed)) 
  }
  
}


plot_dyad_dists <- function(indices){
  for(r in indices){
    plot(y = dyad.dists[together.seqs$i[r], together.seqs$j[r], together.seqs$t.start[r]:together.seqs$t.end[r]], type = 'l',
         x = together.seqs$t.start[r]:together.seqs$t.end[r],
         xlab = 'Time (s)', ylab = 'Distance between individuals (m)')
  }
}

plot_canonical_shape <- function(rows, dyad.dists, together.seqs){
  for(r in rows){
    
    y <- dyad.dists[together.seqs$i[r], together.seqs$j[r], together.seqs$t.start[r]:together.seqs$t.end[r]]
    x <- together.seqs$t.start[r]:together.seqs$t.end[r]
    plot(y = y, type = 'l',
         x = x,
         xlab = 'Time (s)', ylab = 'Distance between individuals (m)')
    
   
  
    fp <- list(x0 = x[1], 
               y0 = y[1],
               xf = x[length(x)],
               yf = y[length(y)]) 
    
    # lines(x =x, y = fission_fusion_function(x, p$par[1], p$par[2], p$par[3], fp), col = 'red')
    y.fit <- fission_fusion_function(x = x, b1 = together.seqs$b1[r], b2 = together.seqs$b2[r],
                                     b.y.intercept = together.seqs$y.intercept[r], fixed.parameters = fp)
    lines(x =x, y = y.fit, col = 'red')
    
  }
}


plot_events <- function(indices, events, xs, ys, phase.col = TRUE){
  for(r in indices){
    x.i <- xs[events$i[r], events$t.start[r]:events$t.end[r]]
    x.i <- x.i - mean(x.i, na.rm = TRUE)
    y.i <- ys[events$i[r], events$t.start[r]:events$t.end[r]]
    y.i <- y.i - mean(y.i, na.rm = TRUE)
    x.j <- xs[events$j[r], events$t.start[r]:events$t.end[r]]
    x.j <- x.j - mean(x.j, na.rm = TRUE)
    y.j <- ys[events$j[r], events$t.start[r]:events$t.end[r]]
    y.j <- y.j - mean(y.j, na.rm = TRUE)
    
    b1.idx <- events$b1[r] - events$t.start[r] + 1
    b2.idx <- events$b2[r] - events$t.start[r] + 1
    col1 <- rep('darkblue', length(x.i))
    if(phase.col){
      col1[b1.idx:b2.idx] <- 'blue'
      col1[b2.idx:length(x.i)] <- 'cyan'
    }
    col2 <- rep('darkred', length(x.j))
    if(phase.col){
      col2[b1.idx:b2.idx] <- 'red'
      col2[b2.idx:length(x.j)] <- 'orange'
    }
    plot(x.i, y.i,asp=1, ylim = c(min(c(y.i, y.j), na.rm =TRUE), max(c(y.i, y.j), na.rm = TRUE)),
         xlim = c(min(c(x.i, x.j), na.rm =TRUE), max(c(x.i, x.j), na.rm = TRUE)), col = scales::alpha(col1, 0.5),
         cex = seq(0.1, 1, length.out = length(x.i)), main = r)
    points(x.j, y.j, col = scales::alpha(col2, 0.5), cex = seq(0.1, 1, length.out = length(x.i)))
    x.lims <- par()$usr[1:2]
    y.lims <- par()$usr[3:4]
  }
}




