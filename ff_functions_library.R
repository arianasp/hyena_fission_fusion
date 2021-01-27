
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


##### Extract ff events and phases
get_ff_events_and_phases <- function(xs, ys, R.fusion = 100, R.fission = 200, max.break = 60*30, verbose = TRUE){
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
            if(curr.dist < R.fusion){
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
          if(curr.dist < R.fusion){
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
          if(curr.dist > R.fission){
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
          if((break.end - break.start <= max.break) & 
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
    far.times <- which(prev.dists > R.fission)
    
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

