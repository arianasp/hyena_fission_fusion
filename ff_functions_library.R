#################################### LIBRARIES ####################################
library(lubridate)
library(utf8)
library(rgdal)
library(pBrackets)

#################################### HELPER FUNCTIONS ####################################

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
  coordinates(non.na.latlons) <- ~Y + X ##HERE
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

### This function outputs the "canonical shape" of a fission-fusion event based on a set of input parameters
# The canonical shape is a u-shaped function defined by the parameters described below
# The input parameters include a set of parameters that will be fitted (x, b1, and b2) as well as a 
# set that remains fixed based on the starting and ending times and distances of the two individuals (fixed.parameters)
#Inputs:
# x: [vector] of time points where the function will be evaluated
# b1: [numeric] time point associated with the first "break point" in the u shape
# b2: [numeric] time point associated with the second "break point" in the u shape
# b.y.intercept: [numeric] height of the flat part of the function (from b1 to b2)
# fixed.parameters: [named list] rest of the parameters that are needed to describe the u-shaped function, but which are held fixed and not fitted
#  fixed.parameters$x0: initial time point (start of the fission-fusion event)
#  fixed.parameters$y0: initial distance apart (at the start of the fission-fusion event)
#  fixed.parameters$xf: final time point (end of the fission-fusion event)
#  fixed.parameters$yf: final distance apart (at the end of the fission-fusion event)
#Outputs:
# val: [vector] of y values, created by applying the function to the vector of time points
fission_fusion_function <- function(x, b1, b2, b.y.intercept, fixed.parameters){
  
  #extract the fixed parameters
  x0 <- fixed.parameters$x0
  y0 <- fixed.parameters$y0
  xf <- fixed.parameters$xf
  yf <- fixed.parameters$yf
  
  #val will hold the y values of the function
  val <- rep(NA, length(x))
  
  ## Aproach segment
  m1 <- (b.y.intercept - y0)/(b1 - x0)
  val[which(x < b1)] <- y0 + m1*(x[which(x < b1)]-x0)
  
  ## Together segment (flat part)
  val[which(x < b2 &  x >= b1)] <- b.y.intercept
  
  ## Depart segment
  m2 <- (yf-b.y.intercept)/(xf - b2)
  val[which(x >= b2)] <- b.y.intercept + m2*(x[which(x >= b2)] - b2)
  
  return(val)
}


### Calculate least squared error for fission_fusion_function. This funciton is used for optimizing the
# parameters b1, b2, and b.y.intercept in the function fission.fusion.function for a given fission-fusion event
#Inputs:
# params: [vector] containing the parameters eventually being fitted (b1, b2, and b.y.intercept)
# fixed.parameters: [named list] containing the fixed parameters associated with a given event (to be passed in to fission.fusion.function)
# x: [vector] of time points from the beginning to the end of the fission-fusion event
# y: [vector] of distances between the two individuals during the fission-fusion event
#Outputs:
# error: [numeric] sum of squared error between data (x,y) and model (fission.fusion.function with the specified parameter set)
ls_error <- function(params, fixed.parameters, x, y){
  y.fit <- fission_fusion_function(x = x, b1 = params[1], b2 = params[2],
                                   b.y.intercept = params[3], fixed.parameters = fixed.parameters)
  
  error <- sum((y.fit - y)^2, na.rm = TRUE)
  return(error)
}

#Calculate the angle between two vectors (i and j) using law of cosines
# The vectors are defined by two end points each (1 and 2)
#Inputs:
# x1.i: [numeric] x value of the first point (1) of the vector i
# x2.i: [numeric] x value of the second point (2) of the vector i
# y1.i: [numeric] y value of the first point (1) of the vector i
# y2.i: [numeric] y value of the second point (2) of the vector i
# x1.j: [numeric] x value of the first point (1) of the vector j
# x2.j: [numeric] x value of the second point (2) of the vector j
# y1.j: [numeric] y value of the first point (1) of the vector j
# y2.j: [numeric] y value of the second point (2) of the vector j
#Outputs:
# angle: angle between the two vectors (in radians)
get_angle_between_vectors <- function(x1.i, x2.i, y1.i, y2.i, 
                                      x1.j, x2.j, y1.j, y2.j){
  dx.i <- x2.i - x1.i
  dy.i <- y2.i - y1.i
  
  dx.j <- x2.j - x1.j
  dy.j <- y2.j - y1.j
  
  s.i <- sqrt((x2.i - x1.i)^2 + (y2.i - y1.i)^2)
  s.j <- sqrt((x2.j - x1.j)^2 + (y2.j - y1.j)^2)
  
  cos.a <- ((dx.i * dx.j) + (dy.i * dy.j))/(s.i * s.j)
  angle <- acos(cos.a)
  return(angle)
}

#Read in den locations from a file
#Inputs:
# den.file.path: [string] full path to the file containing the den locations (as well as other locations)
# den.names: [vector of strings] the names of the dens to use (defaults to the dens used for the current study, i.e. 2017 pilot field season)
#Outputs:
# den.locs: data frame containing information about the dens, including
#  $name: name of the den
#  $lat: latitude coordinate of the den
#  $lon: longitude coordinate of the den
#  $east: easting coordinate of the den
#  $north: northing coordinate of the den
get_dens <- function(den.file.path,
                     den.names = c('DAVE D','RBEND D','RES M D1','DICK D')){
  known.locs <- read.csv(den.file.path,stringsAsFactors=F)
  eastsNorths <- latlon.to.utm(cbind(known.locs$lon,known.locs$lat),southern_hemisphere=T,utm.zone=36)
  known.locs$east <- eastsNorths[,1]
  known.locs$north <- eastsNorths[,2]
  den.locs <- known.locs[which(known.locs$name %in% den.names),]
  den.locs$name <- as.character(den.locs$name)
  return(den.locs)
}

#Calculate spatially discretized heading over time for an individual's trajectory
#The spatially discretized heading of an individual at a given time point is defined as the vector 
#pointing from its current location to its location after it has moved a fixed distance R
#Spatially discretized headings are useful to use when individuals often stop moving, 
#as the headings of stationary individuals (if calculated using a time window) will bounce around randomly with the GPS noise
#Inputs:
# x: [vector] x coordinates of an individual trajectory
# y: [vector] y coordinates of an individual trajectory
# R: [numeric] spatial radius to use for spatial discretization of the headings
# t.idxs: [vector] time indices at which to calculate the heading
# backward: [boolean] whether to calculate the heading at a given time point based on where it was in the past (backward = T)
#  or where it is heading in the future (backward = F). Defaults to F.
# fpt.thresh: [numeric] threshold for the first passage time, below which the heading will not be calculated and will be filled in with NA. This is to prevent individuals from getting 'headings' when they are stationary for too long.
# subsample: [numeric] rate to subsample data at (to save time). Defaults to 1 (use all data, no subsampling)
#Outputs:
# spat.heads: [vector] giving the headings (in radians) of the individual over time
spatial.headings <- function(x,y,R,t.idxs=1:length(x),backward=F, fpt.thresh, subsample = 1){
  
  #initialize
  n.times <- length(x)
  spat.heads <- rep(NA,n.times)
  fpt <- rep(NA,n.times)
  
  #go backwards for backwards vectors
  if(backward){
    t.idxs <- rev(t.idxs)
  }
  
  #loop over all times
  for(t in t.idxs){
    
    ## Skip time points that aren't a multiple of subsampling unit (start at 0 instead of 1)
    if((t-1) %% subsample != 0)
      next
    
    #get current location
    x0 <- x[t]
    y0 <- y[t]
    
    if(is.na(x0)){
      next
    }
    
    #move forward (or backward) until radius reached
    found <- 0
    na.found <- 0
    time.vec <- t:n.times
    if(backward){
      time.vec <- seq(t,1,-1)
    }
    for(j in 1:length(time.vec)){
      i <- time.vec[j]
      
      if(backward){
        dx <- x0 - x[i]
        dy <- y0 - y[i]
      } else{
        dx <- x[i] - x0
        dy <- y[i] - y0
      }
      dist <- sqrt(dx^2+dy^2)
      
      if(is.na(dist)){
        spat.heads[t] <- NA
        found <- 1
        na.found <- 1
        break
      }
      else{
        if(dist >= R){
          found <- 1
          break
        }
      }
    }
    
    #if you reach the end of the trajectory, return (and leave rest as NAs)
    #also make sure that first passage time is less than fpt.thresh
    if(found){
      if(!na.found & (i - t) <= fpt.thresh){
        spat.heads[t] <- atan2(dy,dx)
      }
    } else{
      return(spat.heads)
    }
  }
  return(spat.heads)
}

#Count the number of fission-fusion events per dyad
#Inputs:
# events: [data frame] containing information about the extracted fission-fusion events
# symmetrize: [boolean] indicating whether the lower part of the output matrix should be filled in or not
#Outputs:
# net: [matrix] of size n.inds x n.inds, where net[i,j] gives the number of ff events involving individuals i and j
count_events_per_dyad <- function(events, symmetrize = T){
  
  n.inds <- length(unique(c(events$i, events$j)))
  net <- matrix(NA, nrow = n.inds, ncol = n.inds)
  
  for(i in 1:(n.inds-1)){
    for(j in (i+1):n.inds){
      net[i,j] <- length(which((events$i==i & events$j==j) | (events$i==j & events$j==i)))
      if(symmetrize){
        net[j,i] <- net[i,j]
      }
    }
  }
  
  return(net)
  
}

##### Make table of contiguous runs in a vector v
# Helper function for extract ff events and phases
#Inputs:
# v: [vector] Any type of vector
# Outputs:
# runs.table: data frame containing all contiguous runs in vector v with columns values (run value), t0 (first coordinate of run in v) and tf (last coordinate of run in v)
# (Note: Runs of NA's are not recognized, but returned as multiple one-element runs)
make_run_table = function(v){
  runs = rle(v)
  tf = cumsum(runs$lengths)
  t0 = tf+1
  t0 = c(1,t0[1:(length(t0)-1)])
  if (length(tf)==1){
    t0=1
  }
  runs.table = data.frame(runs$values,t0, tf)
  return(runs.table)
} 

#################################### MAIN FUNCTIONS ####################################

##### Extract ff events and phases
#This is the main function for extracting fission-fusion events from spatial data, fitting the phases for each event, and classifying them into types 
#Inputs:
# xs: [matrix] n.inds x n.times matrix of x coordinates (eastings) for all individuals. xs[i,t] gives the x coordinate of hyena i at time index t
# ys: [matrix] same as xs, but for y coordinates (northings)
# params: [named list] of parameters for extracting ff events, including
#  params$R.fusion: [numeric] inner radius for ff events - individuals must come closer to one another than this distance to for something to count as a ff event
#  params$R.fission: [numeric] outer radius for ff events - when individuals cross this radius, the ff event begins (or ends)
#  params$max.break: [numeric] maximum number of time steps between two events connected by NAs to concatenate them together into 1 event
#  params$move.thresh: [numeric] minimum amount moved (meters) to be considered 'moving' during the fusion or fission phase
#  params$together.travel.thresh: [numeric] minimum amount moved (meters) to be considered 'moving' during the together phase
#  params$local.time.diff: [numeric] time difference between local time and UTC (hours)
#  params$den.dist.thresh: [numeric] threshold distance to consider something to be happening 'at the den'
#  params$last.day.used: [numeric] last day in the data set to use (this should be the day before the first collar died)
#verbose: [boolean] whether to print progress updates as the function runs
#Outputs:
# together.seqs: data frame containing all the fission-fusion events as rows, with information about each event contained in the columns, including
#   together.seqs$i: [numeric] index of the first individual in the ff event
#   together.seqs$j: [numeric] index of the second individual in the ff event
#   together.seqs$t0: [numeric] time point associated with crossing the R.fusion threshold for the first time in an event (note: this is NOT the beginning of the ff event - for this we need to look back in time to when the individuals crossed the R.fission threshold!)
#   together.seqs$tf: [numeric] time point associated with crossing the R.fusion threshold to end the event (note: this is NOT the end of the ff event - for this we need to look to when the individuals cross the R.fission threshold!)
#   together.seqs$t.start: [numeric] time point associated with the start of the event (R.fission crossed)
#   together.seqs$t.end: [numeric] time point associated with the end of the event (R.fission crossed)
#   together.seqs$start.exact: [boolean] whether the start time of the event is exact (will be FALSE if the distance between the pair immediately before t.start is NA)
#   together.seqs$end.exact: [boolean] whether the end time of the event is exact (will be FALSE if the distance between the pair immediately after t.end is NA)
#   together.seqs$closest.app: [numeric] closest approach distance during the event (m)
#   together.seqs$t.closest: [numeric] time index associated with the closest approach during the event
#   together.seqs$b1: [numeric] time index associated with the first break point (beginning of the "together phase")
#   together.seqs$b2: [numeric] time index assocaited with the second break point (end of the "together phase")
#   together.seqs$y.intercept: [numeric] dyadic distance associated with the together phase (the 'height' of the bottom of the "u" of the canonical ff curve)
get_ff_events_and_phases <- function(xs, ys, params, verbose = TRUE){
  
  #Get basic parameters
  n.inds <- nrow(xs)
  n.times <- ncol(xs)
  
  if(verbose)
    print('Computing dyadic distances')
  
  #Compute dyadic distance over time
  dyad.dists <- array(NA, dim=c(n.inds, n.inds, n.times))
  for (i in 1:n.inds) {
    dyad.dists[i,,] <- sqrt( (sweep(xs,2,xs[i,],FUN='-'))^2 + (sweep(ys,2,ys[i,],FUN='-'))^2 )
    dyad.dists[i,i,] <- NA
  }
  
  ################################################################################
  #### Identify together sequences
  
  if(verbose) print('Identifying fission-fusion events')
  
  #Get contiguous sequences of "together" times for each pair and store in data frame
  #data frame has columns:
  # i, j (the two individuals)
  # t0, tf (beginning and ending time of 'bout' of being together)
  
  #- find t0, tf as start and end points of runs <100m
  together.seqs <- data.frame()
  for(i in 1:(n.inds-1)){
    for(j in (i+1):n.inds){
      dists <- dyad.dists[i,j,]
      below.Rfusion = dists<params$R.fusion
      runs.table = make_run_table(below.Rfusion)  # Table of all runs and their start and end points
      true.runs = subset(runs.table, runs.values == TRUE) #  Select only runs of TRUEs (below R.fusion)
      together.seqs = rbind(together.seqs, data.frame(i=i, j=j, t0=true.runs$t0, tf=true.runs$tf))
    }
  }
  
  
  # Aggregate events that are close together in time and separated only by a sequence of NAs of maximum length max.break
  together.seqs.orig = together.seqs
  together.seqs = data.frame()
  
  for(i in 1:(n.inds-1)){
    for(j in (i+1):n.inds){
      dists = dyad.dists[i,j,]
      sub.together = together.seqs.orig[(together.seqs.orig$i == i & together.seqs.orig$j == j), ]
      
      # identify rows that should be merged (with previous row)
      for (r in 2:(nrow(sub.together))){
        end_prev = sub.together[r-1,'tf']
        start_next = sub.together[r,'t0']
        time.between = (start_next-end_prev)-1
        
        if ((time.between <= params$max.break) & sum(is.na(dists[end_prev:start_next]))==time.between){
          sub.together[r,'merge.w.prev'] = T} else {sub.together[r,'merge.w.prev']= F
          }
      }
      
      # merge
      runs.table = make_run_table(sub.together$merge.w.prev)
      true.runs = subset(runs.table, runs.values == TRUE)
      
      if(nrow(true.runs)>0){
        # replace old end with new end
        for (l in 1:nrow(true.runs)){
          sub.together[true.runs[l,'t0']-1, 'tf'] = sub.together[true.runs[l,'tf'], 'tf']
        }
        # remove the other rows that have been merged
        k=1
        while(k<=nrow(sub.together)){
          sub.together = sub.together[!((sub.together$t0>sub.together[k,'t0'] & sub.together$t0<=sub.together[k,'tf'])),]
          k=k+1
        }
      }
      together.seqs = rbind(together.seqs, sub.together)
    }
  }
  # drop merge column
  together.seqs = together.seqs[,c("i", "j", "t0", "tf")]
  
  if(verbose) print('Identifying the correct start and end times of ff events')
  
  #- find t.start and t.end by extending t0 and tf to R.fission threshold (200m)
  together.seqs$t.start = together.seqs$t0 # default value (no elongation)
  together.seqs$t.end = together.seqs$tf # default value (no elongation)
  
  for(i in 1:nrow(together.seqs)){
    ind.i <- together.seqs[i,'i']
    ind.j <- together.seqs[i,'j']
    t0 <- together.seqs$t0[i]
    tf <- together.seqs$tf[i]
    dists =  dyad.dists[ind.i, ind.j,]
    
    # Find t.start
    ## Identify time points
    prev.dists <- dists[1:t0]
    far.times <- which(prev.dists > params$R.fission)
    
    ## Exclude events that are too early to have prior data where they are apart
    if(length(far.times) > 0){
      
      t.start = max(far.times, na.rm = TRUE)+1
      
      # check if the sequence to be added contains any long NA stretches which should break the seq
      to.add = dists[t.start:t0]
      runs.table = make_run_table(is.na(to.add))  
      true.runs = subset(runs.table, runs.values == TRUE)
      true.runs$len = true.runs$tf-true.runs$t0
      break.runs = subset(true.runs, len >= params$max.break)
      
      # if it does, start after last stretch
      if (nrow(break.runs)>0){
        last.break = break.runs[nrow(break.runs), 'tf']
        t.start = t.start + last.break
      }
      together.seqs$t.start[i] = t.start
    }
    
    # Find t.end
    ## Identify time points
    following.dists <- dists[tf:length(dists)]
    far.times <- which(following.dists > params$R.fission)
    
    ## Exclude events that are too late to have following data where they are apart
    if(length(far.times) > 0){
      
      t.end = tf+min(far.times, na.rm = TRUE)-2
      
      # check if the sequence to be added contains any long NA stretches which should break the seq
      to.add = dists[(tf+1):t.end]
      runs.table = make_run_table(is.na(to.add))  
      true.runs = subset(runs.table, runs.values == TRUE)
      true.runs$len = true.runs$tf-true.runs$t0
      break.runs = subset(true.runs, len >= params$max.break)
      
      # if it does, end before first stretch
      if (nrow(break.runs)>0){
        first.break = break.runs[1, 't0']
        t.end = tf + first.break
      }
      together.seqs$t.end[i] = t.end
    }
  }
  
  #- remove duplicates
  together.seqs = together.seqs[!duplicated(together.seqs[,c("i", "j", "t.start", "t.end")]),]
  
  
  # Remove NA's from start and ends of runs
  for (r in 1:nrow(together.seqs)){
    between.dists = dyad.dists[together.seqs$i[r], together.seqs$j[r], together.seqs$t.start[r]:together.seqs$t.end[r]]
    ori.start = together.seqs$t.start[r]
    ori.end = together.seqs$t.end[r]
    first.non.na = min(which(!is.na(between.dists)))
    last.non.na = max(which(!is.na(between.dists)))
    together.seqs$t.start[r] = ori.start + first.non.na - 1
    together.seqs$t.end[r] = ori.start + last.non.na - 1
  }
  
  # Set start.exact and end.exact
  start.exact <- end.exact <- rep(TRUE, nrow(together.seqs))
  together.seqs$start.exact = start.exact
  together.seqs$end.exact = end.exact
  
  # Get whether start and end are exact by checking whether the point right before and right after are NAs
  for(r in 1:nrow(together.seqs)){
    
    ind.i <- together.seqs[r,'i']
    ind.j <- together.seqs[r,'j']
    t.start = together.seqs[r,'t.start']
    t.end = together.seqs[r,'t.end']
    
    if (is.na(dyad.dists[ind.i, ind.j, (t.start-1)])){
      together.seqs[r, 'start.exact'] = FALSE
    }
    
    if (is.na(dyad.dists[ind.i, ind.j, (t.end+1)])){
      together.seqs[r, 'end.exact'] = FALSE
    }
  }
  
  #get closest approach distance
  together.seqs$closest.app <- NA
  for(i in 1:nrow(together.seqs)){
    dists <- dyad.dists[together.seqs$i[i], together.seqs$j[i], together.seqs$t.start[i]:together.seqs$t.end[i]] 
    together.seqs$closest.app[i] <- min(dists,na.rm=T)
  }
  
  #add column for time until next event
  together.seqs$dt <- NA
  for(i in 1:(n.inds-1)){
    for(j in (i+1):n.inds){
      rows <- which(together.seqs$i == i & together.seqs$j == j)
      dt <- together.seqs$t.start[rows[2:length(rows)]] - together.seqs$t.end[rows[1:(length(rows)-1)]] -1
      together.seqs$dt[rows[1:(length(rows)-1)]] <- dt
    }
  }
  
  #Identify phases
  if(verbose)
    print('Identifying phases')
  
  together.seqs[,c('b1', 'b2', 'y.intercept')] <- NA
  for(r in 1:nrow(together.seqs)){
    
    if(!together.seqs$start.exact[r] | !together.seqs$end.exact[r])
      next
    
    y <- dyad.dists[together.seqs$i[r], together.seqs$j[r], together.seqs$t.start[r]:together.seqs$t.end[r]]
    x <- together.seqs$t.start[r]:together.seqs$t.end[r]
    
    fp <- list(x0 = x[1], 
               y0 = y[1],
               xf = x[length(x)],
               yf = y[length(y)])
    
    
    p <- constrOptim(theta = c(fp$x0+5, fp$xf-5, .001), f = ls_error, fixed.parameters = fp, x = x, y =y,
                     ui = matrix(nrow  = 7, ncol = 3, byrow = TRUE,
                                 data = c(1,0,0,
                                          -1,0,0,
                                          0,1,0,
                                          0,-1,0,
                                          0,0,1,
                                          0,0,-1,
                                          -1,1,0)),
                     ci = c(fp$x0+1, -fp$xf-1, fp$x0+1, -fp$xf-1, 0, -fp$yf , 0),
                     grad = NULL)
    
    
    together.seqs[r,c('b1', 'b2', 'y.intercept')] <-p$par
  }
  
  together.seqs$b1 <- round(together.seqs$b1)
  together.seqs$b2 <- round(together.seqs$b2)
  
  return(together.seqs)
}


#### Extract features from fission fusion events
#Inputs:
# xs: [matrix] n.inds x n.times matrix of x coordinates (eastings) for all individuals. xs[i,t] gives the x coordinate of hyena i at time index t
# ys: [matrix] same as xs, but for y coordinates (northings)
# together.seqs: [data frame] of fission-fusion events and associated information (output from get_ff_events_and_phases)
# params: [named list] of parameters for extracting ff events (see above for what it includes)
# den.names: [vector of strings] the names of the dens to use (defaults to the dens used for the current study, i.e. 2017 pilot field season)
# vedbaas: [matrix] of vedba values, of same dimensions and structure as xs and ys matrices
# get.sync.measures: [boolean] whether to calculate the synchrony measures associated with the together phase (this step is computationally intensive)
# sync.subsample: [numeric] by what factor to downsample the data when computing sync measurse (defaults to 10, 1 means no downsampling)
#Outputs:
# together.seqs: [data frame] same as the input, but updated with additional metrics computed including 
#   displacements during each phase of each individual, duration, phase types and event type, heading correlation, activity synchrony
#   (column names should be self-explanatory - see how everything is computed inside the function for details)
get_ff_features <- function(xs, ys, together.seqs, params, den.file.path, den.names, vedbas = vedbas, get.sync.measures = FALSE, sync.subsample = 10){
  
  empty.vec <- rep(NA, nrow(together.seqs))
  positions <- data.frame(x.before.i = empty.vec,
                          x.start.i = empty.vec,
                          x.b1.i = empty.vec,
                          x.b2.i = empty.vec,
                          x.end.i = empty.vec,
                          x.after.i = empty.vec,
                          x.closest.i = empty.vec,
                          y.before.i = empty.vec,
                          y.start.i = empty.vec,
                          y.b1.i = empty.vec,
                          y.b2.i = empty.vec,
                          y.end.i = empty.vec,
                          y.after.i = empty.vec,
                          y.closest.i = empty.vec,
                          x.before.j = empty.vec,
                          x.start.j = empty.vec,
                          x.b1.j = empty.vec,
                          x.b2.j = empty.vec,
                          x.end.j = empty.vec,
                          x.after.j = empty.vec,
                          x.closest.j = empty.vec,
                          y.before.j = empty.vec,
                          y.start.j = empty.vec,
                          y.b1.j = empty.vec,
                          y.b2.j = empty.vec,
                          y.end.j = empty.vec,
                          y.after.j = empty.vec,
                          y.closest.j = empty.vec)
  
  ### Get locations for relevant time points
  
  idxs.start.end.exact <- which(together.seqs$start.exact & together.seqs$end.exact)
  idxs.start.exact <- which(together.seqs$start.exact)
  idxs.end.exact <- which(together.seqs$end.exact)
  
  #####i
  # start exact
  positions$x.start.i[idxs.start.exact] <- xs[cbind(together.seqs$i, together.seqs$t.start)[idxs.start.exact,]]
  positions$y.start.i[idxs.start.exact] <- ys[cbind(together.seqs$i, together.seqs$t.start)[idxs.start.exact,]]
  
  #both exact
  positions$x.b1.i[idxs.start.end.exact] <- xs[cbind(together.seqs$i, together.seqs$b1)[idxs.start.end.exact,]]
  positions$y.b1.i[idxs.start.end.exact] <- ys[cbind(together.seqs$i, together.seqs$b1)[idxs.start.end.exact,]]
  positions$x.b2.i[idxs.start.end.exact] <- xs[cbind(together.seqs$i, together.seqs$b2)[idxs.start.end.exact,]]
  positions$y.b2.i[idxs.start.end.exact] <- ys[cbind(together.seqs$i, together.seqs$b2)[idxs.start.end.exact,]]
  
  #end exact
  positions$x.end.i[idxs.end.exact] <- xs[cbind(together.seqs$i, together.seqs$t.end)[idxs.end.exact,]]
  positions$y.end.i[idxs.end.exact] <- ys[cbind(together.seqs$i, together.seqs$t.end)[idxs.end.exact,]]
  
  positions$x.closest.i <- xs[cbind(together.seqs$i, together.seqs$t.closest)]
  positions$y.closest.i <- ys[cbind(together.seqs$i, together.seqs$t.closest)]
  
  #####j
  # start exact
  positions$x.start.j[idxs.start.exact] <- xs[cbind(together.seqs$j, together.seqs$t.start)[idxs.start.exact,]]
  positions$y.start.j[idxs.start.exact] <- ys[cbind(together.seqs$j, together.seqs$t.start)[idxs.start.exact,]]
  
  #both exact
  positions$x.b1.j[idxs.start.end.exact] <- xs[cbind(together.seqs$j, together.seqs$b1)[idxs.start.end.exact,]]
  positions$y.b1.j[idxs.start.end.exact] <- ys[cbind(together.seqs$j, together.seqs$b1)[idxs.start.end.exact,]]
  positions$x.b2.j[idxs.start.end.exact] <- xs[cbind(together.seqs$j, together.seqs$b2)[idxs.start.end.exact,]]
  positions$y.b2.j[idxs.start.end.exact] <- ys[cbind(together.seqs$j, together.seqs$b2)[idxs.start.end.exact,]]
  
  #end exact
  positions$x.end.j[idxs.end.exact] <- xs[cbind(together.seqs$j, together.seqs$t.end)[idxs.end.exact,]]
  positions$y.end.j[idxs.end.exact] <- ys[cbind(together.seqs$j, together.seqs$t.end)[idxs.end.exact,]]
  
  positions$x.closest.j <- xs[cbind(together.seqs$j, together.seqs$t.closest)]
  positions$y.closest.j <- ys[cbind(together.seqs$j, together.seqs$t.closest)]
  
  
  ##### Extract features
  ### Displacement - no before or after because they are defined a priori as 100m
  
  together.seqs$disp.fusion.i <- sqrt((positions$x.start.i-positions$x.b1.i)^2 + (positions$y.start.i-positions$y.b1.i)^2)
  together.seqs$disp.together.i <- sqrt((positions$x.b2.i-positions$x.b1.i)^2 + (positions$y.b2.i-positions$y.b1.i)^2)
  together.seqs$disp.fission.i <- sqrt((positions$x.end.i-positions$x.b2.i)^2 + (positions$y.end.i-positions$y.b2.i)^2)
  
  together.seqs$disp.fusion.j <- sqrt((positions$x.start.j-positions$x.b1.j)^2 + (positions$y.start.j-positions$y.b1.j)^2)
  together.seqs$disp.together.j <- sqrt((positions$x.b2.j-positions$x.b1.j)^2 + (positions$y.b2.j-positions$y.b1.j)^2)
  together.seqs$disp.fission.j <- sqrt((positions$x.end.j-positions$x.b2.j)^2 + (positions$y.end.j-positions$y.b2.j)^2)
  
  
  ### Time
  together.seqs$duration.fusion <- together.seqs$b1 - together.seqs$t.start
  together.seqs$duration.together <- together.seqs$b2 - together.seqs$b1
  together.seqs$duration.fission <- together.seqs$t.end - together.seqs$b2
  
  ### Angle
  together.seqs$angle.fusion <- get_angle_between_vectors(x1.i = positions$x.start.i, x2.i = positions$x.b1.i, y1.i = positions$y.start.i, y2.i = positions$y.b1.i,
                                                          x1.j = positions$x.start.j, x2.j = positions$x.b1.j, y1.j = positions$y.start.j, y2.j = positions$y.b1.j)
  
  together.seqs$angle.together <- get_angle_between_vectors(x1.i = positions$x.b1.i, x2.i = positions$x.b2.i, y1.i = positions$y.b1.i, y2.i = positions$y.b2.i,
                                                            x1.j = positions$x.b1.j, x2.j = positions$x.b2.j, y1.j = positions$y.b1.j, y2.j = positions$y.b2.j)
  
  together.seqs$angle.fission <- get_angle_between_vectors(x1.i = positions$x.b2.i, x2.i = positions$x.end.i, y1.i = positions$y.b2.i, y2.i = positions$y.end.i,
                                                           x1.j = positions$x.b2.j, x2.j = positions$x.end.j, y1.j = positions$y.b2.j, y2.j = positions$y.end.j)
  
  #### Categorize into types
  
  #defining clusters - fusion
  fusion.stay.move.idxs <- which(together.seqs$disp.fusion.i <= params$move.thresh)
  fusion.move.stay.idxs <- which(together.seqs$disp.fusion.j <= params$move.thresh)
  fusion.move.move.idxs <- which(together.seqs$disp.fusion.i > params$move.thresh & together.seqs$disp.fusion.j > params$move.thresh)
  
  #defining clusters - together
  together.local.idxs <- which(together.seqs$disp.together.i <= params$together.travel.thresh | together.seqs$disp.together.j <= params$together.travel.thresh)
  together.travel.idxs <- which(together.seqs$disp.together.i > params$together.travel.thresh & together.seqs$disp.together.j > params$together.travel.thresh)
  
  #defining clusters - fission
  fission.stay.move.idxs <- which(together.seqs$disp.fission.i <= params$move.thresh)
  fission.move.stay.idxs <- which(together.seqs$disp.fission.j <= params$move.thresh)
  fission.move.move.idxs <- which(together.seqs$disp.fission.i > params$move.thresh & together.seqs$disp.fission.j > params$move.thresh)
  
  #store clusters in events data frame
  together.seqs$fusion.type <- together.seqs$together.type <- together.seqs$fission.type <- NA
  together.seqs$fusion.type[fusion.stay.move.idxs] <- 'fusion.stay.move'
  together.seqs$fusion.type[fusion.move.stay.idxs] <- 'fusion.move.stay'
  together.seqs$fusion.type[fusion.move.move.idxs] <- 'fusion.move.move'
  together.seqs$together.type[together.local.idxs] <- 'together.local'
  together.seqs$together.type[together.travel.idxs] <- 'together.travel'
  together.seqs$fission.type[fission.stay.move.idxs] <- 'fission.stay.move'
  together.seqs$fission.type[fission.move.stay.idxs] <- 'fission.move.stay'
  together.seqs$fission.type[fission.move.move.idxs] <- 'fission.move.move'
  
  #combine into event types
  together.seqs$event.type <- paste(together.seqs$fusion.type, together.seqs$together.type, together.seqs$fission.type, sep = '__')
  
  #create a column for the collapsed event type (collapsing the symmetry of move.stay and stay.move)
  together.seqs$event.type.sym <- together.seqs$event.type
  
  #fusion move.stay --> fusion.stay.move in cases where fission has a mover and a stayer
  fusion.move.stay.idxs <- which(together.seqs$fusion.type == 'fusion.move.stay') #get indexes where fusion type is move.stay
  fission.asymmetric.idxs <- which(grepl('move',together.seqs$fission.type) & grepl('stay', together.seqs$fission.type)) #get indexes where fusion type is stay.move or move.stay
  idxs.to.swap <- intersect(fusion.move.stay.idxs, fission.asymmetric.idxs)
  together.seqs$event.type.sym[idxs.to.swap] <- textclean::swap(together.seqs$event.type[idxs.to.swap], 'move','stay')
  
  #fusion.move.stay --> fusion.stay.move in cases where fission has only movers
  fission.move.move.idxs <- which(!grepl('stay', together.seqs$fission.type))
  idxs.to.change <- intersect(fusion.move.stay.idxs, fission.move.move.idxs)
  together.seqs$event.type.sym[idxs.to.change] <- sub('fusion.move.stay', 'fusion.stay.move', together.seqs$event.type[idxs.to.change])
  
  #fission.move.stay --> fission.stay.move in cases where fusion has only movers
  fission.move.stay.idxs <- which(together.seqs$fission.type == 'fission.move.stay') #get indexes where fission type is move.stay
  fusion.move.move.idxs <- which(!grepl('stay', together.seqs$fusion.type))
  idxs.to.change <- intersect(fission.move.stay.idxs, fusion.move.move.idxs)
  together.seqs$event.type.sym[idxs.to.change] <- sub('fission.move.stay', 'fission.stay.move', together.seqs$event.type[idxs.to.change])
  
  #distance from nearest den
  den.locs <- get_dens(den.file.path, den.names)
  dist.to.den.start <- dist.to.den.end <- matrix(nrow = nrow(together.seqs), ncol = nrow(den.locs))
  for(i in 1:nrow(den.locs)){
    x.d <- den.locs$east[i]
    y.d <- den.locs$north[i]
    
    #start of event
    dxi <- xs[cbind(together.seqs$i, together.seqs$t.start)] - x.d
    dyi <- ys[cbind(together.seqs$i, together.seqs$t.start)] - y.d
    ddi <- sqrt(dxi^2 + dyi^2)
    dxj <- xs[cbind(together.seqs$j, together.seqs$t.start)] - x.d
    dyj <- ys[cbind(together.seqs$j, together.seqs$t.start)] - y.d
    ddj <- sqrt(dxj^2 + dyj^2)
    dist.to.den.start[,i] <- apply(X = cbind(ddi,ddj), MARGIN = 1, FUN = min, na.rm=T)
    
    #end of event
    dxi <- xs[cbind(together.seqs$i, together.seqs$t.end)] - x.d
    dyi <- ys[cbind(together.seqs$i, together.seqs$t.end)] - y.d
    ddi <- sqrt(dxi^2 + dyi^2)
    dxj <- xs[cbind(together.seqs$j, together.seqs$t.end)] - x.d
    dyj <- ys[cbind(together.seqs$j, together.seqs$t.end)] - y.d
    ddj <- sqrt(dxj^2 + dyj^2)
    dist.to.den.end[,i] <- apply(X = cbind(ddi,ddj), MARGIN = 1, FUN = min, na.rm=T)
    
  }
  
  #min dist to den at start and end of events
  together.seqs$dist.den.start <- apply(dist.to.den.start, 1, min, na.rm=T)
  together.seqs$dist.den.end <- apply(dist.to.den.end, 1, min, na.rm=T)
  
  ## Calculate sync measures
  if(get.sync.measures == TRUE){
    
    ## Heading correlation
    together.seqs$together.heading.similarity <- NA
    together.seqs$together.heading.samples <- NA
    together.seqs$vedba.similarity <- NA
    together.seqs$vedba.samples <- NA
    pb <- txtProgressBar(min = 0, max = nrow(together.seqs), char = '..ᕕ(ᐛ)ᕗ..', style = 3)
    for(r in 1:nrow(together.seqs)){
      if(is.na(together.seqs$b1[r]) | is.na(together.seqs$b2[r]))
        next
      
      ## Headings of individual i
      i.heads <- spatial.headings(x = xs[together.seqs$i[r], together.seqs$b1[r]:together.seqs$b2[r]],
                                  y = ys[together.seqs$i[r], together.seqs$b1[r]:together.seqs$b2[r]],
                                  R = 5,
                                  fpt.thresh = 10,
                                  subsample = sync.subsample)
      
      ## Headings of individual j
      j.heads <- spatial.headings(x = xs[together.seqs$j[r], together.seqs$b1[r]:together.seqs$b2[r]],
                                  y = ys[together.seqs$j[r], together.seqs$b1[r]:together.seqs$b2[r]],
                                  R = 5,
                                  fpt.thresh = 10,
                                  subsample = sync.subsample)
      
      
      ## xs and ys of unit vectors of headings
      x.i <- cos(i.heads)
      y.i <- sin(i.heads)
      x.j <- cos(j.heads)
      y.j <- sin(j.heads)
      
      
      ## Dot product of these two vectors
      heading.cors <- (x.i*x.j)+(y.i*y.j)
      
      
      ## Save mean of dot products and number of non NA dot products
      together.seqs$together.heading.similarity[r] <- mean(heading.cors, na.rm = TRUE)
      together.seqs$together.heading.samples[r] <- sum(!is.na(heading.cors))
      
      
      ### Vedba correlation
      i.vedba <- vedbas[together.seqs$i[r], together.seqs$b1[r]:together.seqs$b2[r]]
      j.vedba <- vedbas[together.seqs$j[r], together.seqs$b1[r]:together.seqs$b2[r]]
      if(sum(!is.na(i.vedba + j.vedba)) >= 2){
        together.seqs$vedba.similarity[r] <- cor(i.vedba, j.vedba, use = 'complete')
        together.seqs$vedba.samples[r] <- sum(!is.na(i.vedba + j.vedba))
      }
      
      setTxtProgressBar(pb, r)
    }
    cat('\n')
  }
  
  return(together.seqs)
}

########################## RANDOMIZATION ######################################

#Generate a randomization plan where hyena trajectories from each day will be shuffled
#If rand.params$ensure.no.day.matches == T, make sure that the trajectories don't 'align' on any day 
#(i.e. the hyenas are actually randomized to have the same day represented at the same time in the randomized data)
#Inputs:
# rand.params: [named list] of parameters associated with the randomization, including
#   rand.params$break.hour: [numeric] which hour to "break" at when randomizing days (0 = midnight, 12 = noon), defaults to 12
#   rand.params$last.day.used: [numeric] last day to use in the randomizations (and real data)
#   rand.params$den.blocks: [list of vectors of length n.inds] blocks to keep together for each individual (e.g. to keep den attendance roughly constant)
#   rand.params$ensure.no.day.matches: [boolean] whether to ensure that no pair of individuals is randomized to the same day
#   rand.params$n.rands: [numeric] how many randomizations to do
# max.tries: maximum number of tries at generating the randomization plan before giving up and trying again (only relevant if ensure.no.day.matches == T)
#Outputs:
# rand.plan: [matrix] of dimension n.inds x n.days x n.rands specifying indexes for day swaps. rand.plan[i,j,k] gives the day that should be swapped in for individual i on day j in randomization k
generate_randomization_plan <- function(rand.params, n.inds, max.tries = 10000){
  
  #Set seed for reproducibility
  set.seed(177823)
  
  #initialize array to hold randomization plan
  n.rands <- rand.params$n.rands
  rand.plan <- array(NA, dim = c(n.inds, rand.params$last.day.used, n.rands))
  
  #get blocks
  if(is.null(rand.params$blocks)){
    
    #if no blocks specified, all days are counted as one big block
    blocks <- list()
    for(i in 1:n.inds){
      blocks[[i]] <- c(1, rand.params$last.day.used+1)
    }
    
  } else{
    
    #if blocks were specified, use those and add end points at 1 and last.day.used
    blocks <- rand.params$blocks
    for(i in 1:n.inds){
      blocks[[i]] <- rand.params$blocks[[i]]
      blocks[[i]] <- blocks[[i]][which(blocks[[i]] < rand.params$last.day.used)] 
      blocks[[i]] <- c(1, blocks[[i]], rand.params$last.day.used+1)
    }
    
  }
  
  #generate permutation plans
  r <- 1
  while(r <= n.rands){
    rand.plan.curr <- array(NA, dim = c(n.inds, rand.params$last.day.used))
    
    #generate an initial possibility
    for(i in 1:n.inds){
      for(j in 1:(length(blocks[[i]])-1)){
        d0 <- blocks[[i]][j]
        df <- blocks[[i]][j+1] - 1
        rand.plan.curr[i,d0:df] <- sample(d0:df, replace=F) #shuffle within block #TODO - fix so don't need to use replace, possible to efficiently find possibilities satisfying constraint?
      }
    }
    
    converged <- T
    if(rand.params$ensure.no.day.matches){
      converged <- F
      tries <- 1
      while(!converged){
        uniques <- apply(rand.plan.curr, 2, FUN = function(x){return(length(unique(x)))})
        match.cols <- which(uniques < n.inds)
        if(length(match.cols)==0){
          converged <- T
          break
        }
        if(length(match.cols)>1){
          c <- sample(match.cols,1) #get a column at random from the ones that have duplicates
        } else{
          c <- match.cols
        }
        dup.inds <- which(duplicated(rand.plan.curr[,c])) #for now this prioritizes keeping the first individual in place and shuffling the next individual
        dup.ind <- dup.inds[1]
        block <- max(which(blocks[[dup.ind]] <= c))
        d0 <- blocks[[dup.ind]][block]
        df <- blocks[[dup.ind]][block+1] - 1
        rand.plan.curr[dup.ind,d0:df] <- sample(d0:df, replace = F) #try again
        tries <- tries + 1
        
        #ensure no infinite looping by setting a maximum number of tries before the whole thing is reinitialized
        if(tries > max.tries){
          warning('reached maximum number of tries when generating randomization plans - tried again')
          break
        }
      }
    }
    
    #add to the big array storing all plans (unless not converged, then try again for the same randomization index)
    if(converged){
      rand.plan[,,r] <- rand.plan.curr
      r <- r +1
    } 
  }
  
  return(rand.plan)
}



########################## ANALYSIS + PLOTTING #################################


#get distributions of phase types
#Inputs:
# together.seqs: [data frame] of information about all ff events (output from get_ff_features)
#Outputs:
# out: [vector] number of events of each type
get_event_type_distributions <- function(together.seqs){
  
  event.types.all <- c('fusion.stay.move__together.local__fission.stay.move',
                       'fusion.stay.move__together.local__fission.move.stay',
                       'fusion.stay.move__together.local__fission.move.move',
                       'fusion.move.move__together.local__fission.stay.move',
                       'fusion.move.move__together.local__fission.move.move',
                       'fusion.stay.move__together.travel__fission.stay.move',
                       'fusion.stay.move__together.travel__fission.move.stay',
                       'fusion.stay.move__together.travel__fission.move.move',
                       'fusion.move.move__together.travel__fission.stay.move',
                       'fusion.move.move__together.travel__fission.move.move')
  
  events.all <- together.seqs$event.type.sym
  
  out <- rep(NA, length(event.types.all))
  
  for(i in 1:length(event.types.all)){
    
    out[i] <- sum(events.all == event.types.all[i])
    
  }
  
  names(out) <- event.types.all
  
  return(out)
  
}

#remove events around day breaks
#Inputs:
# together.seqs: [data frame] of information about all ff events (output from get_ff_features)
# timestamps: [vector] of timestamps associated with the data in POSIX format
# rand.params: [named list] of randomization parameters (see above for what it includes)
#Outputs:
# together.seqs: [data frame] of information for all ff events, with events spanning day breaks removed
remove_events_around_day_breaks <- function(together.seqs, timestamps, rand.params){
  
  idx.rem <- c()
  for(i in 1:nrow(together.seqs)){
    
    t.interval <- seq.POSIXt(from = timestamps[together.seqs$t.start[i]],to = timestamps[together.seqs$t.end[i]], by = 'sec')
    break.time <- as.POSIXct(paste(format(t.interval[[1]], '%Y-%m-%d'), paste0(rand.params$break.hour, ':00:00')), tz = tz(timestamps[together.seqs$t.start[i]]))
    if(break.time %in% t.interval){
      idx.rem <- c(idx.rem, i)
    }
  }
  
  if(length(idx.rem) > 0){
    together.seqs <- together.seqs[-idx.rem,]
  }
  
  return(together.seqs)
}


visualize_event_type_distributions <- function(events, events.rand.list, rand.params, timestamps, remove.events.around.day.breaks = T, col){
  
  #quartz()
  
  n.rands <- length(events.rand.list)
  event.type.symbols <- get_event_type_symbols()
  distribs.rand.list <- list()
  distribs.dat <- rep(NA, 10)
  
  #remove events surrounding the 'day break'
  if(remove.events.around.day.breaks){
    events <- remove_events_around_day_breaks(events, timestamps, rand.params)
  }
  
  distribs.dat <- get_event_type_distributions(events)
  for(i in 1:n.rands){
    
    #get events associated with that randomization
    events.rand <- events.rand.list[[i]]
    
    #limit to only events with exact start and end times
    events.rand <- events.rand[events.rand$start.exact & events.rand$end.exact,]
    
    #remove events surrounding the 'day break'
    if(remove.events.around.day.breaks){
      events.rand <- remove_events_around_day_breaks(events.rand, timestamps, rand.params)
    }
    distribs.rand.list[[i]] <- get_event_type_distributions(events.rand)
  }
  
  plot.df <- data.frame(freq = do.call(c, distribs.rand.list),
                        condition = names(do.call(c, distribs.rand.list)))
  
  plot.df$condition <- factor(plot.df$condition, levels = c('fusion.stay.move__together.local__fission.stay.move',
                                                            'fusion.stay.move__together.local__fission.move.stay',
                                                            'fusion.stay.move__together.local__fission.move.move',
                                                            'fusion.move.move__together.local__fission.stay.move',
                                                            'fusion.move.move__together.local__fission.move.move',
                                                            'fusion.stay.move__together.travel__fission.stay.move',
                                                            'fusion.stay.move__together.travel__fission.move.stay',
                                                            'fusion.stay.move__together.travel__fission.move.move',
                                                            'fusion.move.move__together.travel__fission.stay.move',
                                                            'fusion.move.move__together.travel__fission.move.move'),
                              labels = event.type.symbols)
  
  
  plot.df.obs <- data.frame(condition = names(distribs.dat), freq = distribs.dat)
  plot.df.obs$condition <- factor(plot.df.obs$condition, levels = c('fusion.stay.move__together.local__fission.stay.move',
                                                                    'fusion.stay.move__together.local__fission.move.stay',
                                                                    'fusion.stay.move__together.local__fission.move.move',
                                                                    'fusion.move.move__together.local__fission.stay.move',
                                                                    'fusion.move.move__together.local__fission.move.move',
                                                                    'fusion.stay.move__together.travel__fission.stay.move',
                                                                    'fusion.stay.move__together.travel__fission.move.stay',
                                                                    'fusion.stay.move__together.travel__fission.move.move',
                                                                    'fusion.move.move__together.travel__fission.stay.move',
                                                                    'fusion.move.move__together.travel__fission.move.move'),
                                  labels = event.type.symbols)
  
  
  ggplot(data = plot.df, aes(x = condition, y = freq))+
    geom_violin(scale = 'width', fill = col, col = col)+
    geom_point(data = plot.df.obs, shape = '|', size = 7)+
    theme_classic(base_size = 12)+
    geom_vline(aes(xintercept = 5.5))+
    ylab('Frequency')+
    theme(axis.text = element_text(color = 'black'), axis.title.y = element_blank())+
    scale_x_discrete(limits = levels(plot.df$condition)[c(6:10, 1:5)])+
    geom_line(data = data.frame(x = get_event_type_symbols()[c('fusion.stay.move__together.local__fission.move.move',
                                                               'fusion.move.move__together.local__fission.stay.move',
                                                               'fusion.stay.move__together.travel__fission.move.move',
                                                               'fusion.move.move__together.travel__fission.stay.move')],
                                y = plot.df.obs$freq[match(get_event_type_symbols()[c('fusion.stay.move__together.local__fission.move.move',
                                                                                      'fusion.move.move__together.local__fission.stay.move',
                                                                                      'fusion.stay.move__together.travel__fission.move.move',
                                                                                      'fusion.move.move__together.travel__fission.stay.move')],
                                                           plot.df.obs$condition)],
                                group=c(1,1,2,2)),
              aes(x = x, y = y, group = group), lty = 3)+
    coord_flip()
}


#get ascii plotting symbols for each event type
get_event_type_symbols <- function(){
  
  event.types.all <- c('fusion.stay.move__together.local__fission.stay.move',
                       'fusion.stay.move__together.local__fission.move.stay',
                       'fusion.stay.move__together.local__fission.move.move',
                       'fusion.move.move__together.local__fission.stay.move',
                       'fusion.move.move__together.local__fission.move.move',
                       'fusion.stay.move__together.travel__fission.stay.move',
                       'fusion.stay.move__together.travel__fission.move.stay',
                       'fusion.stay.move__together.travel__fission.move.move',
                       'fusion.move.move__together.travel__fission.stay.move',
                       'fusion.move.move__together.travel__fission.move.move')
  
  # arrow <- intToUtf8(0x2191)
  # dot <- intToUtf8(0x2022) 
  
  # travel <- intToUtf8(0x27F9)
  # travel <- intToUtf8(0x21CA)
  #travel <- paste0(arrow, arrow)
  #travel <-intToUtf8(0x21D1)
  #local <- intToUtf8(0x235F)
  
  arrow <- '\u2191'
  dot <- '\u2022'
  local <- '\u2295'
  travel <- '\u21D1'
  
  
  # travel <- paste0(intToUtf8(8593),intToUtf8(8593))
  # local <- paste0(intToUtf8(8625),intToUtf8(8626))
  fission.move.move <- paste0(intToUtf8(8598),intToUtf8(8599))
  fusion.move.move <- paste0(intToUtf8(8599),intToUtf8(8598))
  fission.stay.move <- paste0(intToUtf8(8226), intToUtf8(8599))
  fission.move.stay <- paste0(intToUtf8(8598), intToUtf8(8226))
  fusion.stay.move <- paste0(intToUtf8(8226),intToUtf8(8598))
  
  
  
  event.type.symbols <- rep('', length(event.types.all))
  names(event.type.symbols) <- event.types.all
  event.type.symbols[1] <- paste0(dot, arrow, ' - ', local, ' - ', dot, arrow)
  event.type.symbols[2] <- paste0(dot, arrow, ' - ', local, ' - ', arrow, dot)
  event.type.symbols[3] <- paste0(dot, arrow, ' - ', local, ' - ', arrow, arrow)
  event.type.symbols[4] <- paste0(arrow, arrow, ' - ', local, ' - ', dot, arrow)
  event.type.symbols[5] <- paste0(arrow, arrow, ' - ', local, ' - ', arrow, arrow)
  event.type.symbols[6] <- paste0(dot, arrow, ' - ', travel, ' - ', dot, arrow)
  event.type.symbols[7] <- paste0(dot, arrow, ' - ', travel, ' - ', arrow, dot)
  event.type.symbols[8] <- paste0(dot, arrow, ' - ', travel, ' - ', arrow, arrow)
  event.type.symbols[9] <- paste0(arrow, arrow, ' - ', travel, ' - ', dot, arrow)
  event.type.symbols[10] <- paste0(arrow, arrow, ' - ', travel, ' - ', arrow, arrow)
  
  return(event.type.symbols)
}

visualize_compare_event_properties <- function(events, events.rand.list, params, rand.params, timestamps, remove.events.around.day.breaks = T, cols){
  
  #concatenate randomized events list of tables into one giant table containing data from all randomizations
  events.rand.all <- list()
  for(i in 1:length(events.rand.list)){
    tmp <- events.rand.list[[i]]
    tmp$rand <- i
    #limit to only events with exact start and end times
    tmp <- tmp[tmp$start.exact & tmp$end.exact,]
    if(remove.events.around.day.breaks){
      tmp <- remove_events_around_day_breaks(tmp, timestamps, rand.params)
    }
    events.rand.all[[i]] <- tmp
  }
  events.rand.all <- do.call(rbind, events.rand.all)
  
  #remove events around day breaks if needed (also done above for randomized data)
  if(remove.events.around.day.breaks){
    events <- remove_events_around_day_breaks(events, timestamps, rand.params)
  }
  
  good.idxs.data <- which(events$start.exact & events$end.exact)
  good.idxs.rand <- which(events.rand.all$start.exact & events.rand.all$end.exact)
  
  ## Den events based on distance to den during fusion only (change during revision)
  den.cat.data <- events$dist.den.start <= params$den.dist.thresh # | events$dist.den.end <= params$den.dist.thresh
  den.cat.rand <- events.rand.all$dist.den.start <= params$den.dist.thresh # | events.rand.all$dist.den.end <= params$den.dist.thresh
  
  #quartz(height = 12, width = 10)
  par(mfrow = c(4,2), mar = c(4,4,1,1), oma = c(0,0,2,0))
  #duration
  events.rand.all$duration <- events.rand.all$t.end - events.rand.all$t.start
  events$duration <- events$t.end - events$t.start
  compare_histograms(events$duration[good.idxs.data]/60, events.rand.all$duration[good.idxs.rand]/60, events.rand.all$rand[good.idxs.rand],
                     xlab = 'Duration (min)', logaxes = 'x', n.breaks = 100, custom.breaks=NULL, cumulative=T, 
                     categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand], cols)
  mtext(text = 'A', side = 3, line = 0, adj = 0, font = 2)
  
  #displacement during together phase (minimum of the 2 individuals used)
  events$disp.together <- suppressWarnings(apply(cbind(events$disp.together.i, events$disp.together.j), FUN = function(x){return(min(x,na.rm=T))}, 1))
  events$disp.together[which(is.infinite(events$disp.together))] <- NA
  events.rand.all$disp.together <- suppressWarnings(apply(cbind(events.rand.all$disp.together.i, events.rand.all$disp.together.j), FUN = function(x){return(min(x,na.rm=T))}, 1))
  events.rand.all$disp.together[which(is.infinite(events.rand.all$disp.together))] <- NA
  compare_histograms(events$disp.together[good.idxs.data], events.rand.all$disp.together[good.idxs.rand], events.rand.all$rand[good.idxs.rand],  n.breaks = 500, xlab = 'Displacement together (m)', logaxes = 'x',cumulative=T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand], cols = cols)
  mtext(text = 'B', side = 3, line = 0, adj = 0, font = 2)
  
  #TODO - look at cases of very large displacements together in randomized data - is this real or some weird artifact?
  
  #closest approach
  compare_histograms(events$closest.app[good.idxs.data], events.rand.all$closest.app[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 100, xlab = 'Closest approach (m)', logaxes = 'x', cumulative=T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand], cols = cols)
  mtext(text = 'C', side = 3, line = 0, adj = 0, font = 2)
  
  #whether at the den or not
  events$dist.den.min <- suppressWarnings(apply(cbind(events$dist.den.start, events$dist.den.end), FUN = function(x){return(min(x,na.rm=T))}, 1))
  events$dist.den.min[which(is.infinite(events$dist.den.min))] <- NA
  events.rand.all$dist.den.min <- suppressWarnings(apply(cbind(events.rand.all$dist.den.start, events.rand.all$dist.den.end), FUN = function(x){return(min(x,na.rm=T))}, 1))
  events.rand.all$dist.den.min[which(is.infinite(events.rand.all$dist.den.min))] <- NA
  #compare_histograms(events$dist.den.min[good.idxs.data], events.rand.all$dist.den.min[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 100, xlab = 'Distance from den (m)', logaxes='x', cumulative = T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand])
  
  #time of day (use midpoint of event)
  events$hour <- hour(timestamps[(events$t.start + events$t.end)/2] + params$local.time.diff)
  events.rand.all$hour <- hour(timestamps[(events.rand.all$t.start + events.rand.all$t.end)/2] + params$local.time.diff)
  compare_histograms(events$hour[good.idxs.data], events.rand.all$hour[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 23, xlab = 'Hour of day', logaxes='', custom.breaks = seq(0,23,1), categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand], cols = cols)
  legend('topright',legend = c('Den','Non-Den'), col = c('blue','magenta'), lwd = c(1.5,1.5), lty = c(1,1),cex = 1, bty='n')
  mtext(text = 'D', side = 3, line = 0, adj = 0, font = 2)
  
  #heading similarity
  events$heading.similarity.filt <- events$together.heading.similarity
  events$heading.similarity.filt[which(events$together.heading.samples==0)] <- NA
  events.rand.all$heading.similarity.filt <- events.rand.all$together.heading.similarity
  events.rand.all$heading.similarity.filt[which(events.rand.all$together.heading.samples==0)] <- NA
  compare_histograms(events$heading.similarity.filt[good.idxs.data], events.rand.all$heading.similarity.filt[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 100, xlab = 'Heading similarity', logaxes = '', cumulative = T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand], cols = cols)
  mtext(text = 'E', side = 3, line = 0, adj = 0, font = 2)
  
  #vedba similarity
  compare_histograms(events$vedba.similarity[good.idxs.data], events.rand.all$vedba.similarity[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 100, xlab = 'Activity similarity', logaxes = '', cumulative = T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand], cols = cols)
  mtext(text = 'F', side = 3, line = 0, adj = 0, font= 2)
  
  #distance from den at start
  compare_histograms(events$dist.den.start[good.idxs.data], events.rand.all$dist.den.start[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 100, xlab = 'Starting distance from den (m)', logaxes='x', cumulative = T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand], cols = cols)
  mtext(text = 'G', side = 3, line = 0, adj = 0, font = 2)
  
  #distance from den at end
  compare_histograms(events$dist.den.end[good.idxs.data], events.rand.all$dist.den.end[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 100, xlab = 'Ending distance from den (m)', logaxes='x', cumulative = T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand], cols = cols)
  mtext(text = 'H', side = 3, line = 0, adj = 0, font = 2)
  
}

#make a plot of randomized (black) vs data (red) histograms of a given features
compare_histograms <- function(values.dat, values.rand, randomization.idxs, n.breaks = 100, xlab = '', logaxes = 'x', custom.breaks = NULL, cumulative=F, categories.data = NULL, categories.rand = NULL, cols = cols){
  
  n.rands <- length(unique(randomization.idxs))
  breaks <- seq(min(c(values.dat, values.rand),na.rm=T),max(c(values.dat, values.rand),na.rm=T), length.out = n.breaks)
  if(!is.null(custom.breaks)){
    breaks <- custom.breaks
  }
  
  #if specified, break into two categories based on category idxs (T or F, e.g. den or non-den)
  if(!is.null(categories.data)){
    hist.data.T <- hist(values.dat[categories.data], breaks = breaks, plot = F)
    hist.data.F <- hist(values.dat[!categories.data], breaks = breaks, plot = F)
    mids <- hist.data.T$mids
  } else{
    hist.data <- hist(values.dat, breaks = breaks, plot = F)
    mids <- hist.data$mids
  }
  
  if(!is.null(categories.data)){
    if(cumulative){
      y.data.T <- cumsum(hist.data.T$counts) / sum(hist.data.T$counts)
      y.data.F <- cumsum(hist.data.F$counts) / sum(hist.data.F$counts)
    } else{
      y.data.T <- hist.data.T$density
      y.data.F <- hist.data.F$density
    }
  } else{
    if(cumulative){
      y.data <- cumsum(hist.data$counts) / sum(hist.data$counts)
    } else{
      y.data <- hist.data$density
    }
  }
  
  if(!is.null(categories.data)){
    hists.rand.T <- hists.rand.F <- matrix(NA, nrow = n.rands, ncol = length(breaks)-1)
  } else{
    hists.rand <- matrix(NA, nrow = n.rands, ncol = length(breaks)-1)
  }
  
  for(i in 1:n.rands){
    idxs.rand.i <- which(randomization.idxs == i)
    values.rand.i <- values.rand[idxs.rand.i]
    if(!is.null(categories.data)){
      categories.rand.i <- categories.rand[idxs.rand.i]
      if(cumulative){
        tmp.T <- hist(values.rand.i[which(categories.rand.i)], breaks = breaks, plot = F)
        tmp.F <- hist(values.rand.i[which(!categories.rand.i)], breaks = breaks, plot = F)
        hists.rand.T[i,] <- cumsum(tmp.T$counts) / sum(tmp.T$counts)
        hists.rand.F[i,] <- cumsum(tmp.F$counts) / sum(tmp.F$counts)
      } else{
        hists.rand.T[i,] <- hist(values.rand.i[which(categories.rand.i)], breaks = breaks, plot = F)$density
        hists.rand.F[i,] <- hist(values.rand.i[which(!categories.rand.i)], breaks = breaks, plot = F)$density
      }
      
    } else{
      if(cumulative){
        tmp <- hist(values.rand.i, breaks = breaks, plot = F)
        hists.rand[i,] <- cumsum(tmp$counts) / sum(tmp$counts)
      } else{
        hists.rand[i,] <- hist(values.rand.i, breaks = breaks, plot = F)$density
      }
    }
  }
  
  #plotting 
  if(!is.null(categories.data)){
    ymin <- min(min(cbind(hists.rand.T, hists.rand.F),na.rm=T),min(cbind(y.data.T, y.data.F),na.rm=T))
    ymax <- max(max(cbind(hists.rand.T, hists.rand.F),na.rm=T),max(cbind(y.data.T, y.data.F),na.rm=T))
    if(cumulative){
      ylab <- 'Cumulative probability'
    } else{
      ylab <- 'Density'
    }
    plot(NULL,  xlab = '', ylab = '', log = logaxes, xlim = range(mids), ylim = c(ymin, ymax), cex.lab = 1, cex.axis = 1)
    for(i in 1:n.rands){
      lines(hist.data.T$mids, hists.rand.T[i,], lwd = 0.4, col = alpha(cols[2], 0.2))
      lines(hist.data.F$mids, hists.rand.F[i,], lwd = 0.4, col = alpha(cols[1], 0.2))
    }
    lines(mids, y.data.T, col = cols[2], lwd = 1.5)
    lines(mids, y.data.F, col = cols[1], lwd = 1.5)
    title(ylab = ylab, line = 2)
    title(xlab = xlab, line = 2)
    
    
  } else{
    ymin <- min(min(hists.rand,na.rm=T),min(y.data,na.rm=T))
    ymax <- max(max(hists.rand,na.rm=T),max(y.data,na.rm=T))
    if(cumulative){
      ylab <- 'Cumulative probability'
    } else{
      ylab <- 'Density'
    }
    plot(NULL,  xlab = '', ylab = '', log = logaxes, xlim = range(hist.data$mids), ylim = c(ymin, ymax), cex.lab = 1.5, cex.axis = 1.5)
    for(i in 1:n.rands){
      lines(hist.data$mids, hists.rand[i,], lwd = 0.4, col = alpha(cols[2], 0.2))
    }
    lines(hist.data$mids, y.data, col = 'red', lwd = 1.5)
    title(ylab = ylab, line = 2)
    title(xlab = xlab, line = 2)
  }
}

#make a visualization of ff event of index r
plot_events <- function(r, events, xs, ys, cols, ...){
  i <- events$i[r]
  j <- events$j[r]
  x.i <- xs[events$i[r], events$t.start[r]:events$t.end[r]]
  x.i <- x.i - mean(range(x.i, na.rm = TRUE))
  y.i <- ys[events$i[r], events$t.start[r]:events$t.end[r]]
  y.i <- y.i - mean(range(y.i, na.rm = TRUE))
  x.j <- xs[events$j[r], events$t.start[r]:events$t.end[r]]
  x.j <- x.j - mean(range(x.j, na.rm = TRUE))
  y.j <- ys[events$j[r], events$t.start[r]:events$t.end[r]]
  y.j <- y.j - mean(range(y.j, na.rm = TRUE))
  
  b1.idx <- events$b1[r] - events$t.start[r] + 1
  b2.idx <- events$b2[r] - events$t.start[r] + 1
  
  plot.df <- data.frame(x = c(x.i, x.j), y = c(y.i, y.j), id = factor(c(rep(i, length(x.i)), rep(j, length(x.j)))), time = rep(1:length(x.i), 2), 
                        phase = NA)
  plot.df[plot.df$time < b1.idx,'phase'] <- 'fusion'
  plot.df[plot.df$time >= b1.idx & plot.df$time <= b2.idx,'phase'] <- 'together'
  plot.df[plot.df$time > b2.idx,'phase'] <- 'fission'
  
  x.range = max(plot.df$x, na.rm = TRUE) - min(plot.df$x, na.rm = TRUE)
  y.range = max(plot.df$y, na.rm = TRUE) - min(plot.df$y, na.rm = TRUE)
  x.most.extreme = range(plot.df$x, na.rm  = TRUE)[which.max(abs(range(plot.df$x, na.rm = TRUE)))]
  y.most.extreme = range(plot.df$y, na.rm  = TRUE)[which.max(abs(range(plot.df$y, na.rm = TRUE)))]
  max.range = ceiling(max(x.range, y.range) / 50) * 50
  scalebar = data.frame(x = c(max.range/2, max.range/2 - R.fusion),
                        y = rep(max.range/2 - 0.1*max.range/2, 2))
  scalebar.label = data.frame(x = mean(scalebar$x), y = mean(scalebar$y) + 0.07 * mean(scalebar$y),
                              label = paste0(R.fusion, 'm'))
  arrow.df = data.frame(x = c(-max.range/2 + 0.1*max.range/2, -max.range/2 + 0.1*max.range/2),
                     y = c(max.range/2 - 0.2*max.range/2, max.range/2- 0.07*max.range/2))
  arrow.label = data.frame(x = -max.range/2 + 0.1*max.range/2,
                           y = max.range/2- 0.1*max.range/2 + 0.07 * max.range/2,
                           label = 'N')
  
  ggplot(plot.df, aes(x =x ,y=y, col = id))+
    geom_path(size = 1)+
    geom_line(aes(x = x, y= y, group = time), inherit.aes = F, lty = 1, size = 0.5,
              data = plot.df[plot.df$phase == 'together',], alpha = 0.2, col = 'gray30')+
    geom_line(data = scalebar, aes(x = x, y = y), inherit.aes = F)+
    geom_line(arrow = arrow(type = 'closed', length = unit(0.01, units = 'npc')), data = arrow.df, aes(x = x, y = y), inherit.aes = F)+
    geom_text(data = arrow.label, aes(x = x, y = y, label = label), inherit.aes = F, size = 3)+
    geom_text(data = scalebar.label, aes(x = x, y = y, label = label), inherit.aes = F, size = 3)+
    geom_point(data = plot.df[plot.df$time == 1,], shape = 'S', size = 2)+
    geom_point(data = plot.df[plot.df$time == max(plot.df$time),], shape = 'E', size = 2)+
    theme_classic()+
    theme(legend.position = 'none', axis.title = element_blank(), aspect.ratio = 1,
          axis.line = element_blank(), panel.border = element_rect(fill = NA, color = 'black', size = 2),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.margin=grid::unit(c(0,0,0,0), "mm"))+
    geom_line(data = data.frame(x.bar = max(plot.df$x) - 100,  max(plot.df$x),
                                y.bar = rep(min(plot.df$y) + 10, 2)),
              aes(x = x.bar, y = y.bar), inherit.aes = FALSE)+
    scale_color_manual(values = cols)+
    xlim(-max.range/2, max.range/2)+
    ylim(-max.range/2, max.range/2)
}

animate_events <- function(r, events, xs, ys, cols, ...){
  i <- events$i[r]
  j <- events$j[r]
  x.i <- xs[events$i[r], events$t.start[r]:events$t.end[r]]
  x.i <- x.i - mean(range(x.i, na.rm = TRUE))
  y.i <- ys[events$i[r], events$t.start[r]:events$t.end[r]]
  y.i <- y.i - mean(range(y.i, na.rm = TRUE))
  x.j <- xs[events$j[r], events$t.start[r]:events$t.end[r]]
  x.j <- x.j - mean(range(x.j, na.rm = TRUE))
  y.j <- ys[events$j[r], events$t.start[r]:events$t.end[r]]
  y.j <- y.j - mean(range(y.j, na.rm = TRUE))
  
  b1.idx <- events$b1[r] - events$t.start[r] + 1
  b2.idx <- events$b2[r] - events$t.start[r] + 1
  
  plot.df <- data.frame(x = c(x.i, x.j), y = c(y.i, y.j), id = factor(c(rep(i, length(x.i)), rep(j, length(x.j)))), time = rep(1:length(x.i), 2), 
                        phase = NA)
  plot.df[plot.df$time < b1.idx,'phase'] <- 'fusion'
  plot.df[plot.df$time >= b1.idx & plot.df$time <= b2.idx,'phase'] <- 'together'
  plot.df[plot.df$time > b2.idx,'phase'] <- 'fission'
  
  x.range = max(plot.df$x, na.rm = TRUE) - min(plot.df$x, na.rm = TRUE)
  y.range = max(plot.df$y, na.rm = TRUE) - min(plot.df$y, na.rm = TRUE)
  x.most.extreme = range(plot.df$x, na.rm  = TRUE)[which.max(abs(range(plot.df$x, na.rm = TRUE)))]
  y.most.extreme = range(plot.df$y, na.rm  = TRUE)[which.max(abs(range(plot.df$y, na.rm = TRUE)))]
  max.range = ceiling(max(x.range, y.range) / 50) * 50
  scalebar = data.frame(x = c(max.range/2, max.range/2 - R.fusion),
                        y = rep(max.range/2 - 0.1*max.range/2, 2))
  scalebar.label = data.frame(x = mean(scalebar$x), y = mean(scalebar$y) + 0.07 * mean(scalebar$y),
                              label = paste0(R.fusion, 'm'))
  arrow.df = data.frame(x = c(-max.range/2 + 0.1*max.range/2, -max.range/2 + 0.1*max.range/2),
                        y = c(max.range/2 - 0.2*max.range/2, max.range/2- 0.07*max.range/2))
  arrow.label = data.frame(x = -max.range/2 + 0.1*max.range/2,
                           y = max.range/2- 0.1*max.range/2 + 0.07 * max.range/2,
                           label = 'N')
  
  ggplot(plot.df, aes(x =x ,y=y, col = id))+
    geom_point(size = 3, alpha = 0.5)+
    geom_point(size = 2, alpha = 0.5, data = filter(plot.df, phase == 'together'))+
    # geom_line(aes(x = x, y= y, group = time), inherit.aes = F, lty = 1, size = 0.5,
    #           data = plot.df[plot.df$phase == 'together',], alpha = 0.2, col = 'gray30')+
    geom_line(data = scalebar, aes(x = x, y = y), inherit.aes = F)+
    geom_line(arrow = arrow(type = 'closed', length = unit(0.01, units = 'npc')), data = arrow.df, aes(x = x, y = y), inherit.aes = F)+
    geom_text(data = arrow.label, aes(x = x, y = y, label = label), inherit.aes = F, size = 3)+
    geom_text(data = scalebar.label, aes(x = x, y = y, label = label), inherit.aes = F, size = 3)+
    #geom_point(data = plot.df[plot.df$time == 1,], shape = 'S', size = 2)+
    #geom_point(data = plot.df[plot.df$time == max(plot.df$time),], shape = 'E', size = 2)+
    theme_classic()+
    theme(legend.position = 'none', axis.title = element_blank(), aspect.ratio = 1,
          axis.line = element_blank(), panel.border = element_rect(fill = NA, color = 'black', size = 2),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.margin=grid::unit(c(0,0,0,0), "mm"), plot.title = element_text(size = 11, hjust = 0.5, vjust = -0.01))+
    # geom_line(data = data.frame(x.bar = max(plot.df$x) - 100,  max(plot.df$x),
    #                             y.bar = rep(min(plot.df$y) + 10, 2)),
    #           aes(x = x.bar, y = y.bar), inherit.aes = FALSE)+
    scale_color_manual(values = cols)+
    xlim(-max.range/2, max.range/2)+
    ylim(-max.range/2, max.range/2)+
    transition_components(time)+
    shadow_wake(wake = 0.3, exclude_layer = 1)+
    ggtitle(gsub(x = gsub(x = events$event.type.sym[r], '__', ' | ', fixed = TRUE), '.', ' ', fixed = TRUE))
}

#plot the fitted fission-fusion function associated with a given event r
plot_canonical_shape <- function(r, together.seqs, xs, ys){
  
  y <- sqrt( (xs[together.seqs$i[r],together.seqs$t.start[r]:together.seqs$t.end[r]] - xs[together.seqs$j[r],together.seqs$t.start[r]:together.seqs$t.end[r]])^2 +
               (ys[together.seqs$i[r],together.seqs$t.start[r]:together.seqs$t.end[r]] - ys[together.seqs$j[r],together.seqs$t.start[r]:together.seqs$t.end[r]])^2)
  
  x <- together.seqs$t.start[r]:together.seqs$t.end[r]
  
  
  fp <- list(x0 = x[1],
             y0 = y[1],
             xf = x[length(x)],
             yf = y[length(y)])
  
  y.fit <- fission_fusion_function(x = x, b1 = together.seqs$b1[r], b2 = together.seqs$b2[r],
                                   b.y.intercept = together.seqs$y.intercept[r], fixed.parameters = fp)
  
  ggplot(data = data.frame(x = x-together.seqs$t.start[r], y, y.fit), aes(x,y))+
    geom_line(size = 0.5)+
    geom_line(aes(y = y.fit), col = '#ba0c2f', size = 0.5)+
    theme_classic(base_size = 12)+
    xlab('Time (s)')+
    ylab('Distance between individuals (m)')
  
}

#generate all the figures and save output
generate_figures <- function(data.outdir, plots.outdir, code.directory){
  print('Generating figures')
  source(paste0(code.directory, 'ff_summary_figures.R'), local = TRUE, print.eval = TRUE)
}
