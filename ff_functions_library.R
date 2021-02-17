#################################### LIBRARIES ####################################
library(lubridate)
library(utf8)
library(rgdal)

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

### Find parameters defining the 'canonical shape' that best matches data for a fission-fusion event
fission_fusion_function <- function(x, b1, b2, b.y.intercept, fixed.parameters){
  x0 <- fixed.parameters$x0
  y0 <- fixed.parameters$y0
  xf <- fixed.parameters$xf
  yf <- fixed.parameters$yf
  
  val <- rep(NA, length(x))
  ## Aproach
  m1 <- (b.y.intercept - y0)/(b1 - x0)
  val[which(x < b1)] <- y0 + m1*(x[which(x < b1)]-x0)
  
  ## Together
  val[which(x < b2 &  x >= b1)] <- b.y.intercept
  
  ## Depart
  m2 <- (yf-b.y.intercept)/(xf - b2)
  val[which(x >= b2)] <- b.y.intercept + m2*(x[which(x >= b2)] - b2)
  return(val)
}


### Calculate least squared error for fission_fusion_function. For optimization. 
ls_error <- function(params, fixed.parameters, x, y){
  y.fit <- fission_fusion_function(x = x, b1 = params[1], b2 = params[2],
                                   b.y.intercept = params[3], fixed.parameters = fixed.parameters)
  
  error <- sum((y.fit - y)^2, na.rm = TRUE)
  return(error)
}

#Calculate the angle between two vectors using law of cosines
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

# #get day start indexes
# get_day_start_idxs <- function(timestamps, local.time.diff = 3){
#   
#   ts.local <- timestamps + local.time.diff * 60 * 60
#   dates <- as.Date(ts.local)
#   dates.uniq <- sort(unique(dates))
#   day.start.idxs <- rep(NA, length(dates.uniq) + 1)
#   for(i in 1:length(dates.uniq)){
#     day.start.idxs[i] <- min(which(dates == dates.uniq[i]))
#   }
#   day.start.idxs[length(day.start.idxs)] <- length(dates)
#   
#   return(day.start.idxs)
#   
# }

#read in den locations
get_dens <- function(den.file.path = '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/rawdata/metadata/hyena_isolate_dens.csv',
                     den.names = c('DAVE D','RBEND D','RES M D1','DICK D')){
  known.locs <- read.csv(den.file.path,stringsAsFactors=F)
  eastsNorths <- latlon.to.utm(cbind(known.locs$lon,known.locs$lat),southern_hemisphere=T,utm.zone=36)
  known.locs$east <- eastsNorths[,1]
  known.locs$north <- eastsNorths[,2]
  den.locs <- known.locs[which(known.locs$name %in% den.names),]
  den.locs$name <- as.character(den.locs$name)
  return(den.locs)
}

#Calculate spatially discretized heading over time using sliding window
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

#################################### MAIN FUNCTIONS ####################################

##### Extract ff events and phases
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
  
  if(verbose)
    print('Identifying fission-fusion events')
  
  #Get contiguous sequences of "together" times for each pair and store in data frame
  #data frame has columns:
  # i, j (the two individuals)
  # t0, tf (beginning and ending time of 'bout' of being together)
  # t0.exact, tf.exact (indicates whether the time step immediately before and after the 'bout' was present (TRUE) or missing/NA (FALSE))
  together.seqs <- data.frame()
  for(i in 1:(n.inds-1)){
    for(j in (i+1):n.inds){
      dists <- dyad.dists[i,j,]
      prev.together <- NA
      for(k in 1:n.times){
        
        #get current distance apart
        curr.dist <- dists[k]
        
        #if coming from an NA sequence...
        if(is.na(prev.together)){
          #if still NA, pass
          if(is.na(curr.dist)){
            next
            
            #if not NA...
          } else{
            
            #if within R.fusion, start a together sequence
            if(curr.dist < params$R.fusion){
              t0 <- k
              prev.together <- TRUE
              next
              
              #otherwise, do not start a sequence but set prev.together to FALSE
            } else{
              prev.together <- FALSE
            }
          }
        }
        
        #if previously not together
        if(prev.together == FALSE){
          
          #if you run into an NA, set prev.together to NA
          if(is.na(curr.dist)){
            prev.together <- NA
            next
          }
          
          #if curr.dist is less than R.fusion, start a sequence
          if(curr.dist < params$R.fusion){
            t0 <- k
            prev.together <- TRUE
            next
            
            #otherwise pass
          } else{
            next
          }
          
        }
        
        #if previously together...
        if(prev.together == TRUE){
          
          #if you run into an NA, end the sequence, add to table, set prev.together to NA
          if(is.na(curr.dist)){
            tf <- k - 1
            together.seqs <- rbind(together.seqs, data.frame(i = i, j = j, t0 = t0, tf = tf))
            prev.together <- NA
            next
          }
          
          #if moved away to > R.fission distance, end sequence, set prev.together to FALSE
          if(curr.dist > params$R.fission){
            tf <- k -1
            together.seqs <- rbind(together.seqs, data.frame(i = i, j = j, t0 = t0, tf = tf))
            prev.together <- FALSE
            next
          }
        }
      }
    }
  }
  together.seqs.orig <- together.seqs
  
  #Aggregate events that are close together in time and separated only by a sequence of NAs of maximum length max.break
  together.seqs <- together.seqs.orig
  converged <- F
  while(!converged){
    new.seqs <- data.frame()
    for(i in 1:(n.inds-1)){
      for(j in (i+1):n.inds){
        rows <- which(together.seqs$i == i & together.seqs$j == j)
        k <- 1
        while(k < length(rows)){
          break.start <- together.seqs$tf[rows[k]] + 1
          break.end <- together.seqs$t0[rows[k+1]] - 1
          
          #if two together periods are separated only by NAs and less than max.break seconds, aggregate into one event
          if((break.end - break.start <= params$max.break) & 
             (sum(!is.na(dyad.dists[i,j,break.start:break.end])) == 0)){
            new.seqs <- rbind(new.seqs, data.frame(i = i, j = j, t0 = together.seqs$t0[rows[k]], tf = together.seqs$tf[rows[k+1]]))
            k <- k + 2
          } else{
            new.seqs <- rbind(new.seqs, together.seqs[rows[k],])
            k <- k + 1
            
            #add the last row if needed
            if(k == length(rows)){
              new.seqs <- rbind(new.seqs, together.seqs[rows[k],])
            }
          }
        }
      }
    }
    
    #check if aggregation has converged
    if(nrow(new.seqs) == nrow(together.seqs)){
      converged <- T
    } else{
      converged <- F
    }
    
    #replace with new sequences
    together.seqs <- new.seqs
  }
  
  #add column for time until next event
  together.seqs$dt <- NA
  for(i in 1:(n.inds-1)){
    for(j in (i+1):n.inds){
      rows <- which(together.seqs$i == i & together.seqs$j == j)
      dt <- together.seqs$t0[rows[2:length(rows)]] - together.seqs$tf[rows[1:(length(rows)-1)]]
      together.seqs$dt[rows[1:(length(rows)-1)]] <- dt
    }
  }
  
  #add columns for whether start and end time are exact (whether there is an NA before or after the event)
  befores <- afters <- rep(NA, nrow(together.seqs))
  for(j in 1:nrow(together.seqs)){
    if(together.seqs$t0[j] > 1){
      befores[j] <- dyad.dists[together.seqs$i[j], together.seqs$j[j], together.seqs$t0[j]-1]
    } 
    if(together.seqs$tf[j] < ncol(xs)){
      afters[j] <- dyad.dists[together.seqs$i[j], together.seqs$j[j], together.seqs$tf[j]+1]
    } 
  }
  #befores <- dyad.dists[cbind(together.seqs$i, together.seqs$j, together.seqs$t0-1)]
  #afters <- dyad.dists[cbind(together.seqs$i, together.seqs$j, together.seqs$tf+1)]
  together.seqs$t0.exact <- !is.na(befores)
  together.seqs$tf.exact <- !is.na(afters)
  
  #get closest approach distance
  together.seqs$closest.app <- NA
  for(i in 1:nrow(together.seqs)){
    dists <- dyad.dists[together.seqs$i[i], together.seqs$j[i], together.seqs$t0[i]:together.seqs$tf[i]]
    together.seqs$closest.app[i] <- min(dists,na.rm=T)
  }
  
  
  ################################################################################
  ##### Extract start, end, and other relevant time points
  
  if(verbose)
    print('Extracting start and end time points')
  
  together.seqs[,c('t.start', 't.end',
                   't.closest', 'start.exact', 'end.exact')] <- NA
  for(i in 1:nrow(together.seqs)){
    ind.i <- together.seqs[i,'i']
    ind.j <- together.seqs[i,'j']
    t0 <- together.seqs$t0[i]
    tf <- together.seqs$tf[i]
    
    
    ## Identify time points
    prev.dists <- dyad.dists[ind.i, ind.j,1:t0]
    
    #End
    together.seqs$t.end[i] <- tf
    t.end <- tf
    together.seqs$end.exact[i] <- together.seqs[i,'tf.exact']
    
    #Start
    far.times <- which(prev.dists > params$R.fission)
    
    ## Exclude events that are too early to have prior data where they are apart
    if(!(length(far.times) > 0)){
      together.seqs$start.exact[i] <- FALSE
      next
    }
    
    together.seqs$t.start[i] <- max(far.times, na.rm = TRUE)
    t.start <- together.seqs$t.start[i]
    together.seqs$start.exact[i] <- !is.na(dyad.dists[ind.i, ind.j, t.start+1])
    
    #Closest approach - (first) time at t.end individuals are closest
    together.seqs$t.closest[i] <- t.start + base::which.min(dyad.dists[ind.i, ind.j, t.start:t.end]) - 1
  }
  
  ### Remove events (only 1) with undefined start and end times
  together.seqs <- together.seqs[!is.na(together.seqs$t.start) & !is.na(together.seqs$t.end),]
  
  ################################################################################
  ### Loop over all events and extract time points that define phases
  
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
get_ff_features <- function(xs, ys, together.seqs, params, den.file.path, den.names, get.sync.measures = FALSE, sync.subsample = 10){
  
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
    for(r in 1:nrow(together.seqs)){
      pb <- txtProgressBar(min = 0, max = nrow(together.seqs), char = '|', style = 3)
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
      
      setTxtProgressBar(pb, i)
    }
  }
  
  return(together.seqs)
}

#Compute the datetimes of den attendance at each den by each individual
#Creates a plot of this if specified
plot_den_attendance_by_ind <- function(xs, ys, den.file.path, den.names, params, den.blocks = NULL){

  #get number of inds and times
  n.inds <- dim(xs)[1]
  n.times <- dim(xs)[2]
  n.dens <- length(den.names)
  
  #hard coded boundaries
  if(is.null(den.blocks)){
    den.blocks <- list()
    den.blocks[[1]] <- c(12, 22, 30, 33)
    den.blocks[[2]] <- c(13, 26)
    den.blocks[[3]] <- c(14, 35)
    den.blocks[[4]] <- c(15)
    den.blocks[[5]] <- c(12, 31)
  }
  
  #for now set a single threshold or arrival and leaving (= den.dist.thresh) - consider changing this to two threhsolds for less 'flicker'
  thresh.arrive <- params$den.dist.thresh
  thresh.leave <- params$den.dist.thresh

  #Get data frame with den locations and names
  den.locs <- get_dens(den.file.path = den.file.path, den.names = den.names)
  dens <- den.locs[,c('name','lat','lon','east','north')]
  dens$col <- c('green','magenta','blue','red')

  #Get distance of each hyena to each den over time
  dist.dens <- array(NA,dim=c(dim(xs),nrow(den.locs)))
  for(i in 1:nrow(den.locs)){
    dist.dens[,,i] <- sqrt((xs - den.locs$east[i])^2 + (ys - den.locs$north[i])^2)
  }
  
  #shift time by 12 hours (or other break.hour time) - shift to the night basically
  nights.all <- date(timestamps + params$local.time.diff * 60 * 60 + rand.params$break.hour * 60 *60)
  nights <- unique(nights.all)
  
  #for each night compute secs spent at each den
  den.nights <- array(dim = c(n.inds, length(nights), n.dens))
  for(i in 1:n.inds){
    for(j in 1:length(nights)){
      idxs <- which(nights.all == nights[j])
      for(k in 1:n.dens){
        den.nights[i,j,k] <- sum(dist.dens[i,idxs,k] < params$den.dist.thresh, na.rm=T)
      }
    }
  }
  
  #find the most used den of that night (or NA if no dens were used for > 1 hr)
  most.used.dens <- matrix(NA,nrow = n.inds, ncol = length(nights))
  for(i in 1:n.inds){
    for(j in 1:length(nights)){
      if(sum(den.nights[i,j,]) > 60*60){
        most.used.dens[i,j] <- which.max(den.nights[i,j,])
      }
    }
  }
  

  #get days and minutes into day for all ts
  step <- 60 #subsample to one sample per minute
  dist.dens.sub <- dist.dens[,seq(1,length(timestamps),step),]
  ts.sub <- timestamps[seq(1,length(timestamps),step)]
  days <- as.Date(ts.sub)
  uniq.days <- unique(days)
  mins <- minute(ts.sub) + hour(ts.sub)*60
  
  dist.dens.thresh <- dist.dens.sub < thresh.arrive
  quartz(width=12,height=6)
  par(mfrow=c(1,5),mar=c(4,4,1,0))
  cols <- dens$col
  for(i in 1:n.inds){
    plot(NULL,xlim=c(0,24),ylim=c(1,length(uniq.days)),main=hyena.ids$name[i],xlab='hour of day',ylab='day')
    for(d in 1:dim(dist.dens.thresh)[3]){
      for(day in 1:length(uniq.days)){
        idxs.curr <- which(days == uniq.days[day])
        mins.curr <- mins[idxs.curr]
        den.presence.curr <- which(dist.dens.thresh[i,idxs.curr,d]==TRUE)
        mins.at.den <- mins.curr[den.presence.curr]
        points(mins.at.den/60,rep(day,length(mins.at.den)),col=cols[d],pch=19,cex=0.3)
        
      }
    }
    
    #plot boundary marks
    for(j in 1:length(den.blocks[[i]])){
      b <- den.blocks[[i]][j]
      lines(c(0,12,12,24),c(b,b,b-1,b-1)+.5, lty = 2)
    }
    if(i==4){
      legend('topright',legend=c('Rbend Den','Dave Den','Res Den','Dick Den'),fill=c('blue','green','red','magenta'))
    }
  }
  
  
}

########################## RANDOMIZATION ######################################

#Generate a randomization plan where hyena trajectories from each day will be shuffled
#If ensure.no.day.matches == T, make sure that the trajectories don't 'align' on any day 
#(i.e. the hyenas are actually randomized to have the same day represented at the same time in the randomized data)
generate_randomization_plan <- function(rand.params, n.inds, ensure.no.day.matches = F, max.tries = 10000){

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
    if(ensure.no.day.matches){
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
get_phase_type_distributions <- function(together.seqs){
  
  fusion.types <- c('fusion.stay.move','fusion.move.move')
  together.types <- c('together.local','together.travel')
  fission.types <- c('fission.stay.move','fission.move.move')
  
  complete.events <- which(!grepl('NA',together.seqs$event.type))
  
  events.fusion <- together.seqs$fusion.type[complete.events]
  events.fusion[which(events.fusion=='fusion.move.stay')] <- 'fusion.stay.move'
  events.together <- together.seqs$together.type[complete.events]
  events.fission <- together.seqs$fission.type[complete.events]
  events.fission[which(events.fission=='fission.move.stay')] <- 'fission.stay.move'
  
  events.all <- c(events.fusion, events.together, events.fission)
  event.types.all <- c(fusion.types, together.types, fission.types)
  
  out <- rep(NA, length(event.types.all))
  
  for(i in 1:length(event.types.all)){
    
    out[i] <- sum(events.all == event.types.all[i])
    
  }
  
  names(out) <- event.types.all
  
  return(out)
  
}

#get distributions of phase types
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
remove_events_around_day_breaks <- function(together.seqs, timestamps, rand.params){
  
  idx.rem <- c()
  for(i in 1:nrow(together.seqs)){

    t.interval <- seq.POSIXt(from = timestamps[together.seqs$t0[i]],to = timestamps[together.seqs$tf[i]], by = 'sec')
    hrs <- hour(t.interval)
    if(rand.params$break.hour %in% t.interval){
      idx.rem <- c(idx.rem, i)
    }
  }
  
  if(length(idx.rem) > 0){
    together.seqs <- together.seqs[-idx.rem,]
  }
  
  return(together.seqs)
}

visualize_event_type_distributions <- function(events, events.rand.list, rand.params, timestamps, remove.events.around.day.breaks = T){
  
  quartz()
  
  n.rands <- length(events.rand.list)
  event.type.symbols <- get_event_type_symbols()
  distribs.rand <- matrix(NA, nrow = 10, ncol = n.rands)
  distribs.dat <- rep(NA, 10)
  
  #remove events surrounding the 'day break'
  if(remove.events.around.day.breaks){
    events <- remove_events_around_day_breaks(events, timestamps, rand.params)
  }
  
  distribs.dat <- get_event_type_distributions(events)
  for(i in 1:n.rands){
    
    #get events associated with that randomization
    events.rand <- events.rand.list[[i]]
    
    #remove events surrounding the 'day break'
    if(remove.events.around.day.breaks){
      events.rand <- remove_events_around_day_breaks(events.rand, timestamps, rand.params)
    }
    distribs.rand[,i] <- get_event_type_distributions(events.rand)
  }
  
  visualize_yvals_vs_event_type(distribs.dat, distribs.rand, 'Frequency')
  
}

#general plotting function for plotting event types on x axis and some value on y axis
#yvals.dat is a vector of length 10 (n event types)
#yvals.rand is a matrix of 10 x n.rands containing data from randomizations
visualize_yvals_vs_event_type <- function(yvals.dat, yvals.rand, ylab){
  
  n.rands <- ncol(yvals.rand)
  event.type.symbols <- get_event_type_symbols()
  
  par(mar = c(6,5,1,1))
  max.freq <- max(max(yvals.rand,na.rm=T), max(yvals.dat,na.rm=T))
  plot(NULL, xlim = c(.8, 10.2), ylim = c(0, max.freq*1.1), xlab = '', ylab = ylab, xaxt='n', cex.lab=1.5, cex.axis = 1.5)
  axis(side = 1, at = seq(1,10), labels = event.type.symbols, line = 2.3, tick = F, srt = 180, cex.axis = 1.5)
  abline(v=5.5, lty = 1)
  polygon(c(2.5,2.5,4.5,4.5), c(-max.freq*.03,max.freq*1.1, max.freq*1.1,-max.freq*.03), border='blue', lty = 2)
  polygon(c(7.5,7.5,9.5,9.5), c(-max.freq*.03,max.freq*1.1, max.freq*1.1,-max.freq*.03), border='blue', lty = 2)
  
  for(i in 1:n.rands){
    jitter <- rnorm(10, mean=0, sd = .1)
    points(seq(1,10) + jitter, yvals.rand[,i], col = '#00000044', cex = .8)
  }
  points(yvals.dat, col = 'red', cex = 1.5, pch = 8)
  
}

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
  
  arrow <- intToUtf8(0x2191)
  dot <- intToUtf8(0x2022) 
  travel <- intToUtf8(0x27F9)
  travel <- intToUtf8(0x21CA)
  travel <- paste0(arrow, arrow)
  #travel <-intToUtf8(0x21D1)
  #local <- intToUtf8(0x235F)
  local <- paste0(dot, dot)
  
  event.type.symbols <- rep('', length(event.types.all))
  names(event.type.symbols) <- event.types.all
  event.type.symbols[1] <- paste0(dot, arrow, '\n', local, '\n', dot, arrow)
  event.type.symbols[2] <- paste0(arrow, dot, '\n', local, '\n', dot, arrow)
  event.type.symbols[3] <- paste0(arrow, arrow, '\n', local, '\n', dot, arrow)
  event.type.symbols[4] <- paste0(dot, arrow, '\n', local, '\n', arrow, arrow)
  event.type.symbols[5] <- paste0(arrow, arrow, '\n', local, '\n', arrow, arrow)
  event.type.symbols[6] <- paste0(dot, arrow, '\n', travel, '\n', dot, arrow)
  event.type.symbols[7] <- paste0(arrow, dot, '\n', travel, '\n', dot, arrow)
  event.type.symbols[8] <- paste0(arrow, arrow, '\n', travel, '\n', dot, arrow)
  event.type.symbols[9] <- paste0(dot, arrow, '\n', travel, '\n', arrow, arrow)
  event.type.symbols[10] <- paste0(arrow, arrow, '\n', travel, '\n', arrow, arrow)
  
  return(event.type.symbols)
}

#plot fraction of time each event type occurs at den in beginning or end, compare to randomized
visualize_den_association_by_event_type <- function(events, events.rand.list, timestamps, rand.params, params, assess.at.start = T, remove.events.around.day.breaks=T){
  
  quartz()
  
  n.rands <- length(events.rand.list)
  event.type.symbols <- get_event_type_symbols()
  event.types.all <- names(event.type.symbols)
  den.frac.dat <- rep(NA, length(event.types.all))
  den.frac.rand <- matrix(NA, nrow = length(event.types.all), ncol = n.rands)
  
  if(remove.events.around.day.breaks){
    events <- remove_events_around_day_breaks(events, timestamps, rand.params)
  }
  
  for(i in 1:length(event.types.all)){
    if(assess.at.start){
      den.frac.dat[i] <- mean(events$dist.den.start[which(events$event.type.sym==event.types.all[i])] <= params$den.dist.thresh)
    } else{
      den.frac.dat[i] <- mean(events$dist.den.end[which(events$event.type.sym==event.types.all[i])] <= params$den.dist.thresh)
    }
    for(j in 1:n.rands){
      events.rand <- events.rand.list[[j]]
      #remove events surrounding the 'day break'
      if(remove.events.around.day.breaks){
        events.rand <- remove_events_around_day_breaks(events.rand, timestamps, rand.params)
      }
      if(assess.at.start){
        den.frac.rand[i,j] <- mean(events.rand$dist.den.start[which(events.rand$event.type.sym==event.types.all[i])] <= params$den.dist.thresh)
      } else{
        den.frac.rand[i,j] <- mean(events.rand$dist.den.end[which(events.rand$event.type.sym==event.types.all[i])] <= params$den.dist.thresh)
      }
    }
  }
  
  if(assess.at.start){
    ylab <- 'Fraction start at den'
  } else{
    ylab <- 'Fraction end at den'
  }
  visualize_yvals_vs_event_type(den.frac.dat, den.frac.rand, ylab = ylab)
  
  
}

#plot event type frequencies in real data vs permuted data - separate out den vs non-den events
visualize_event_type_frequency_den_vs_nonden <- function(events, events.rand.list, params, rand.params, timestamps, remove.events.around.day.breaks = T){
  
  quartz()
  
  #setup
  n.rands <- length(events.rand.list)
  event.type.symbols <- get_event_type_symbols()
  distribs.rand.den <- distribs.rand.nonden <- matrix(NA, nrow = 10, ncol = n.rands)
  distribs.dat.den <- distribs.dat.nonden <- rep(NA, 10)
  
  #remove events surrounding the 'day break'
  if(remove.events.around.day.breaks){
    events <- remove_events_around_day_breaks(events, timestamps, rand.params)
  }
  
  #get indexes to events associated with dens (at either beginning or end) and ones not assoc with dens
  den.idxs <- which(events$dist.den.start <= params$den.dist.thresh | events$dist.den.end <= params$den.dist.thresh)
  nonden.idxs <- which(events$dist.den.start > params$den.dist.thresh & events$dist.den.end > params$den.dist.thresh)
  
  distribs.dat.den <- get_event_type_distributions(events[den.idxs,])
  distribs.dat.nonden <- get_event_type_distributions(events[nonden.idxs,])
  
  for(i in 1:n.rands){
    
    #get events associated with that randomization
    events.rand <- events.rand.list[[i]]
    
    #remove events surrounding the 'day break'
    if(remove.events.around.day.breaks){
      events.rand <- remove_events_around_day_breaks(events.rand, timestamps, rand.params)
    }
    
    den.idxs <- which(events.rand$dist.den.start <= params$den.dist.thresh | events.rand$dist.den.end <= params$den.dist.thresh)
    nonden.idxs <- which(events.rand$dist.den.start > params$den.dist.thresh & events.rand$dist.den.end > params$den.dist.thresh)
    
    distribs.rand.den[,i] <- get_event_type_distributions(events.rand[den.idxs,])
    distribs.rand.nonden[,i] <- get_event_type_distributions(events.rand[nonden.idxs,])
  }
  
  par(mfrow=c(2,1))
  visualize_yvals_vs_event_type(distribs.dat.den, distribs.rand.den, 'Frequency (at den)')
  visualize_yvals_vs_event_type(distribs.dat.nonden, distribs.rand.nonden, 'Frequency (not at den)')
  
}

visualize_compare_event_properties <- function(events, events.rand.list, params, rand.params, timestamps, property.colname, remove.events.around.day.breaks = T){
  
  #concatenate randomized events list of tables into one giant table containing data from all randomizations
  events.rand.all <- events.rand.list[[1]]
  events.rand.all$rand <- 1
  for(i in 2:length(events.rand.list)){
    tmp <- events.rand.list[[i]]
    tmp$rand <- i
    if(remove.events.around.day.breaks){
      tmp <- remove_events_around_day_breaks(tmp, timestamps, rand.params)
    }
    events.rand.all <- rbind(events.rand.all, tmp)
  }
  
  #remove events around day breaks if needed (also done above for randomized data)
  if(remove.events.around.day.breaks){
    events <- remove_events_around_day_breaks(events, timestamps, rand.params)
  }
  
  good.idxs.data <- which(events$start.exact & events$end.exact)
  good.idxs.rand <- which(events.rand.all$start.exact & events.rand.all$end.exact)
  
  den.cat.data <- events$dist.den.start <= params$den.dist.thresh | events$dist.den.end <= params$den.dist.thresh
  den.cat.rand <- events.rand.all$dist.den.start <= params$den.dist.thresh | events.rand.all$dist.den.end <= params$den.dist.thresh
  
  quartz(height = 6, width = 20)
  par(mfrow = c(1,5), mar = c(5,6,1,1))
  #duration
  events.rand.all$duration <- events.rand.all$t.end - events.rand.all$t.start
  events$duration <- events$t.end - events$t.start
  compare_histograms(events$duration[good.idxs.data]/60, events.rand.all$duration[good.idxs.rand]/60, events.rand.all$rand[good.idxs.rand], xlab = 'Duration (min)', logaxes = 'x', n.breaks = 100, custom.breaks=NULL, cumulative=T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand])
  
  #displacement during together phase (minimum of the 2 individuals used)
  events$disp.together <- suppressWarnings(apply(cbind(events$disp.together.i, events$disp.together.j), FUN = function(x){return(min(x,na.rm=T))}, 1))
  events$disp.together[which(is.infinite(events$disp.together))] <- NA
  events.rand.all$disp.together <- suppressWarnings(apply(cbind(events.rand.all$disp.together.i, events.rand.all$disp.together.j), FUN = function(x){return(min(x,na.rm=T))}, 1))
  events.rand.all$disp.together[which(is.infinite(events.rand.all$disp.together))] <- NA
  compare_histograms(events$disp.together[good.idxs.data], events.rand.all$disp.together[good.idxs.rand], events.rand.all$rand[good.idxs.rand],  n.breaks = 500, xlab = 'Displacement together (m)', logaxes = 'x',cumulative=T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand])
  
  #TODO - look at cases of very large displacements together in randomized data - is this real or some weird artifact?
  
  #closest approach
  compare_histograms(events$closest.app[good.idxs.data], events.rand.all$closest.app[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 100, xlab = 'Closest approach (m)', logaxes = 'x', cumulative=T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand])
  
  #whether at the den or not
  events$dist.den.min <- suppressWarnings(apply(cbind(events$dist.den.start, events$dist.den.end), FUN = function(x){return(min(x,na.rm=T))}, 1))
  events$dist.den.min[which(is.infinite(events$dist.den.min))] <- NA
  events.rand.all$dist.den.min <- suppressWarnings(apply(cbind(events.rand.all$dist.den.start, events.rand.all$dist.den.end), FUN = function(x){return(min(x,na.rm=T))}, 1))
  events.rand.all$dist.den.min[which(is.infinite(events.rand.all$dist.den.min))] <- NA
  compare_histograms(events$dist.den.min[good.idxs.data], events.rand.all$dist.den.min[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 100, xlab = 'Distance from den (m)', logaxes='x', cumulative = T, categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand])
  
  #time of day (use midpoint of event)
  events$hour <- hour(timestamps[(events$t.start + events$t.end)/2] + params$local.time.diff)
  events.rand.all$hour <- hour(timestamps[(events.rand.all$t.start + events.rand.all$t.end)/2] + params$local.time.diff)
  compare_histograms(events$hour[good.idxs.data], events.rand.all$hour[good.idxs.rand], events.rand.all$rand[good.idxs.rand], n.breaks = 23, xlab = 'Hour of day', logaxes='', custom.breaks = seq(0,23,1), categories.data = den.cat.data[good.idxs.data], categories.rand = den.cat.rand[good.idxs.rand])
  legend('topright',legend = c('Den','Non-Den'), col = c('blue','magenta'), lwd = c(2,2), lty = c(1,1),cex = 2)
  
}

#make a plot of randomized (black) vs data (red) histograms of a given features
compare_histograms <- function(values.dat, values.rand, randomization.idxs, n.breaks = 100, xlab = '', logaxes = 'x', custom.breaks = NULL, cumulative=F, categories.data = NULL, categories.rand = NULL){
  
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
    plot(NULL,  xlab = xlab, ylab = ylab, log = logaxes, xlim = range(mids), ylim = c(ymin, ymax), cex.lab = 1.5, cex.axis = 1.5)
    for(i in 1:n.rands){
      lines(hist.data.T$mids, hists.rand.T[i,], lwd = 0.2, col = '#0000FF33')
      lines(hist.data.F$mids, hists.rand.F[i,], lwd = 0.2, col = '#FF00FF33')
    }
    lines(mids, y.data.T, col = 'blue', lwd = 3)
    lines(mids, y.data.F, col = 'magenta', lwd = 3)
    
    
  } else{
    ymin <- min(min(hists.rand,na.rm=T),min(y.data,na.rm=T))
    ymax <- max(max(hists.rand,na.rm=T),max(y.data,na.rm=T))
    if(cumulative){
      ylab <- 'Cumulative probability'
    } else{
      ylab <- 'Density'
    }
    plot(NULL,  xlab = xlab, ylab = ylab, log = logaxes, xlim = range(hist.data$mids), ylim = c(ymin, ymax), cex.lab = 1.5, cex.axis = 1.5)
    for(i in 1:n.rands){
      lines(hist.data$mids, hists.rand[i,], lwd = 0.2, col = '#00000033')
    }
    lines(hist.data$mids, y.data, col = 'red', lwd = 3)
  }
}

