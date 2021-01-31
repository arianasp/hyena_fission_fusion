#################################### LIBRARIES ####################################
library(lubridate)

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
  befores <- dyad.dists[cbind(together.seqs$i, together.seqs$j, together.seqs$t0-1)]
  afters <- dyad.dists[cbind(together.seqs$i, together.seqs$j, together.seqs$tf+1)]
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
    
    y = dyad.dists[together.seqs$i[r], together.seqs$j[r], together.seqs$t.start[r]:together.seqs$t.end[r]]
    x = together.seqs$t.start[r]:together.seqs$t.end[r]
    fp <- list(x0 = x[1], 
               y0 = y[1],
               xf = x[length(x)],
               yf = y[length(y)])
    
    
    p <- constrOptim(theta = c(fp$x0+5, fp$xf-5, 1), f = ls_error, fixed.parameters = fp, x = x, y =y,
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
get_ff_features <- function(xs, ys, together.seqs, params, den.file.path, den.names){
  
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
  
  return(together.seqs)
}

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
  
  out <- list()
  
  fusion.table <- table(events.fusion)
  out$fusion <- fusion.table[match(names(fusion.table), fusion.types)]
  together.table <- table(events.together)
  out$together <- together.table[match(names(together.table), together.types)]
  fission.table <- table(events.fission)
  out$fission <- fission.table[match(names(fission.table), fission.types)]
  
  return(out)
  
}

visualize_phase_type_distributions <- function(events, events.rand.list){
  
  n.rands <- length(events.rand.list)
  distribs.rand <- matrix(NA, nrow = 6, ncol = n.rands)
  distribs.dat <- rep(NA, 6)
  
  for(i in 1:n.rands){
    
    freqs.rand <- get_phase_type_distributions(events.rand.list[[i]])
    distribs.rand[1,i] <- freqs.rand$fusion[1] 
    distribs.rand[2,i] <- freqs.rand$fusion[2]
    distribs.rand[3,i] <- freqs.rand$together[1] 
    distribs.rand[4,i] <- freqs.rand$together[2] 
    distribs.rand[5,i] <- freqs.rand$fission[1] 
    distribs.rand[6,i] <- freqs.rand$fission[2] 
    
  }
  
  freqs.dat <- get_phase_type_distributions(events)
  distribs.dat[1] <- freqs.dat$fusion[1] 
  distribs.dat[2] <- freqs.dat$fusion[2] 
  distribs.dat[3] <- freqs.dat$together[1]
  distribs.dat[4] <- freqs.dat$together[2] 
  distribs.dat[5] <- freqs.dat$fission[1] 
  distribs.dat[6] <- freqs.dat$fission[2] 
  
  par(mar = c(6,4,1,4))
  max.freq <- max(max(distribs.rand), max(distribs.dat))
  plot(NULL, xlim = c(.8, 6.2), ylim = c(0, max.freq+50), xlab = '', ylab = 'Frequency', xaxt='n')
  abline(v=4.5)
  abline(v=2.5)
  for(i in 1:n.rands){
    points(c(1,2), distribs.rand[1:2,i], col = 'black', cex = .8)
    points(c(3,4), distribs.rand[3:4,i], col = 'black', cex = .8)
    points(c(5,6), distribs.rand[5:6,i], col = 'black', cex = .8)
  }
  points(c(1,2), distribs.dat[1:2], col = 'black',pch = 19, cex = 2)
  points(c(3,4), distribs.dat[3:4], col = 'black',pch = 19, cex = 2)
  points(c(5,6), distribs.dat[5:6], col = 'black',pch = 19, cex = 2)
  
  par(new=T)
  plot(NULL, ylim = c(0, 1.2), xlim = c(.8,6.2), ylab = '', xlab = '', xaxt='n', yaxt = 'n', las =2)
  text(1.5, 1.2, labels='Fusion', adj = c(.5,.5))
  text(3.5, 1.2, labels='Together', adj = c(.5,.5))
  text(5.5, 1.2, labels='Fission', adj = c(.5,.5))
  for(i in 1:n.rands){
    lines(c(1,2), distribs.rand[1:2,i]/sum(distribs.rand[1:2,i]), col = 'red', lty = 2)
    lines(c(3,4), distribs.rand[3:4,i]/sum(distribs.rand[3:4,i]), col = 'red', lty = 2)
    lines(c(5,6), distribs.rand[5:6,i]/sum(distribs.rand[5:6,i]), col = 'red', lty = 2)
  }
  lines(c(1,2), distribs.dat[1:2]/sum(distribs.dat[1:2]), col = 'red',lwd=3)
  lines(c(3,4), distribs.dat[3:4]/sum(distribs.dat[3:4]), col = 'red',lwd=3)
  lines(c(5,6), distribs.dat[5:6]/sum(distribs.dat[5:6]), col = 'red',lwd=3)
  axis(side = 1, at = seq(1,6), labels = c('stay/move','move/move','local','travel','stay/move','move/move'), las =2)
  axis(side = 4, at = seq(0,1,.2), labels = seq(0,1,.2), col.axis = 'red', col = 'red')
  mtext(text = 'Fraction', side = 4, line=2.5, col = 'red')
  
}
