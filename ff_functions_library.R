#DIRECTORIES - MODIFY THIS TO YOUR OWN COMPUTER SETUP
name <- readline(prompt = "Enter name:")

if(grepl('a|A', name)){
  basedir <- '~/Dropbox/hyenas/hyena_fission_fusion/' # <-- CHANGE TO THE MAIN PROJECT FOLDER
  codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/' # <-- CHANGE TO THE DIRECTORY WHERE CODE REPOSITORY IS
}else if(grepl('e|E', name)){
  basedir <- 'C:/Users/strau/Dropbox/hyena_fission_fusion/' # <-- CHANGE TO THE MAIN PROJECT FOLDER
  codedir <- 'C:/Users/strau/Dropbox/Documents/Research/Partial_projects/hyena_fission_fusion/' # <-- CHANGE TO THE DIRECTORY WHERE CODE REPOSITORY IS
}else{
  cat('Name not recognized')
}

#SUBDIRECTORIES
datadir <- paste(basedir, 'data/raw', sep = '')
processeddir <- paste(basedir, 'data/processed', sep = '') #this is where the output will be saved (the fission fusion events RData file)

#FUNCTIONS
setwd(codedir)
source('hyena_functions.R')

#PARAMETERS
R.before.after <- 100

#----------- MAIN -----------------
#Set working directory for main analysis
setwd(datadir)

#Load data
load('hyena_xy_level1.RData')
load('hyena_timestamps.Rdata')
load('hyena_ids.RData')
load('hyena_vedba.RData')

#Read in den data
den.file <- 'hyena_isolate_dens.csv'
den.names <- c('DAVE D','RBEND D','RES M D1','DICK D')
known.locs <- read.csv(den.file,stringsAsFactors=F)
eastsNorths <- latlon.to.utm(cbind(known.locs$lon,known.locs$lat),southern_hemisphere=T,utm.zone=36)
known.locs$east <- eastsNorths[,1]
known.locs$north <- eastsNorths[,2]
den.locs <- known.locs[which(known.locs$name %in% den.names),]
den.locs$name <- as.character(den.locs$name)

setwd(processeddir)
load('fission_fusion_events.RData')

together.seqs[,c('t.start', 't.end', 't.before.i', 
                 't.after.i', 't.before.j', 't.after.j',
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
  
  #Before - different for each individual
  before.steps <- seq(t.start, 1, -1)
  together.seqs$t.before.i[i] <- t.start - first.passage.time(x = xs[ind.i,before.steps],
                     y = ys[ind.i,before.steps],
                     R = R.before.after) + 1
  
  together.seqs$t.before.j[i] <- t.start - first.passage.time(x = xs[ind.j,before.steps],
                                             y = ys[ind.j,before.steps],
                                             R = R.before.after) + 1
  
  #After - different for each individual
  after.steps <- seq(t.end, ncol(xs), 1)
  together.seqs$t.after.i[i] <- t.end + first.passage.time(x = xs[ind.i,after.steps],
                                             y = ys[ind.i,after.steps],
                                             R = R.before.after) - 1 
  
  together.seqs$t.after.j[i] <- t.end + first.passage.time(x = xs[ind.j,after.steps],
                                             y = ys[ind.j,after.steps],
                                             R = R.before.after) - 1 
  
  #Closest approach - (first) time at t.end individuals are closest
  together.seqs$t.closest[i] <- t.start + base::which.min(dyad.dists[ind.i, ind.j, t.start:t.end]) - 1
  
  print(i)
}

### Remove events (only 1) with undefined start and end times
together.seqs <- together.seqs[!is.na(together.seqs$t.start) & !is.na(together.seqs$t.end),]




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
positions$x.before.i[idxs.start.exact] <- xs[cbind(together.seqs$i, together.seqs$t.before.i)[idxs.start.exact,]]
positions$y.before.i[idxs.start.exact] <- ys[cbind(together.seqs$i, together.seqs$t.before.i)[idxs.start.exact,]]
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
positions$x.after.i[idxs.end.exact] <- xs[cbind(together.seqs$i, together.seqs$t.after.i)[idxs.end.exact,]]
positions$y.after.i[idxs.end.exact] <- ys[cbind(together.seqs$i, together.seqs$t.after.i)[idxs.end.exact,]]

positions$x.closest.i <- xs[cbind(together.seqs$i, together.seqs$t.closest)]
positions$y.closest.i <- ys[cbind(together.seqs$i, together.seqs$t.closest)]

#####j
# start exact
positions$x.before.j[idxs.start.exact] <- xs[cbind(together.seqs$j, together.seqs$t.before.j)[idxs.start.exact,]]
positions$y.before.j[idxs.start.exact] <- ys[cbind(together.seqs$j, together.seqs$t.before.j)[idxs.start.exact,]]
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
positions$x.after.j[idxs.end.exact] <- xs[cbind(together.seqs$j, together.seqs$t.after.j)[idxs.end.exact,]]
positions$y.after.j[idxs.end.exact] <- ys[cbind(together.seqs$j, together.seqs$t.after.j)[idxs.end.exact,]]

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
together.seqs$duration.before.i <- together.seqs$t.start - together.seqs$t.before.i
together.seqs$duration.after.i <- together.seqs$t.after.i - together.seqs$t.end

together.seqs$duration.fusion <- together.seqs$b1 - together.seqs$t.start
together.seqs$duration.together <- together.seqs$b2 - together.seqs$b1
together.seqs$duration.fission <- together.seqs$t.end - together.seqs$b2

together.seqs$duration.before.j <- together.seqs$t.start - together.seqs$t.before.j
together.seqs$duration.after.j <- together.seqs$t.after.j - together.seqs$t.end

# ### Speed
# together.seqs.speed.before.i <- together.seqs$disp.before.i / together.seqs$duration.before.i
# together.seqs$speed.fusion.i <- together.seqs$disp.fusion.i / together.seqs$duration.fusion
# together.seqs$speed.together.i <- together.seqs$disp.together.i / together.seqs$duration.together
# together.seqs$speed.fission.i <- together.seqs$disp.fission.i / together.seqs$duration.fission
# together.seqs$speed.after.i <- together.seqs$disp.after.i / together.seqs$duration.after.i
# 
# together.seqs.speed.before.j <- together.seqs$disp.before.j / together.seqs$duration.before.j
# together.seqs$speed.fusion.j <- together.seqs$disp.fusion.j / together.seqs$duration.fusion
# together.seqs$speed.together.j <- together.seqs$disp.together.j / together.seqs$duration.together
# together.seqs$speed.fission.j <- together.seqs$disp.fission.j / together.seqs$duration.fission
# together.seqs$speed.after.j <- together.seqs$disp.after.j / together.seqs$duration.after.j

### Angle
together.seqs$angle.before <- get_angle(x1.i = positions$x.before.i, x2.i = positions$x.start.i, y1.i = positions$y.before.i, y2.i = positions$y.start.i,
                                        x1.j = positions$x.before.j, x2.j = positions$x.start.j, y1.j = positions$y.before.j, y2.j = positions$y.start.j)

together.seqs$angle.fusion <- get_angle(x1.i = positions$x.start.i, x2.i = positions$x.b1.i, y1.i = positions$y.start.i, y2.i = positions$y.b1.i,
                                        x1.j = positions$x.start.j, x2.j = positions$x.b1.j, y1.j = positions$y.start.j, y2.j = positions$y.b1.j)

together.seqs$angle.together <- get_angle(x1.i = positions$x.b1.i, x2.i = positions$x.b2.i, y1.i = positions$y.b1.i, y2.i = positions$y.b2.i,
                                          x1.j = positions$x.b1.j, x2.j = positions$x.b2.j, y1.j = positions$y.b1.j, y2.j = positions$y.b2.j)

together.seqs$angle.fission <- get_angle(x1.i = positions$x.b2.i, x2.i = positions$x.end.i, y1.i = positions$y.b2.i, y2.i = positions$y.end.i,
                                         x1.j = positions$x.b2.j, x2.j = positions$x.end.j, y1.j = positions$y.b2.j, y2.j = positions$y.end.j)

together.seqs$angle.after <- get_angle(x1.i = positions$x.end.i, x2.i = positions$x.after.i, y1.i = positions$y.end.i, y2.i = positions$y.after.i,
                                       x1.j = positions$x.end.j, x2.j = positions$x.after.j, y1.j = positions$y.end.j, y2.j = positions$y.after.j)

### VEDBA

for(i in 1:nrow(together.seqs)){
  together.seqs$vedba.before.i <- mean(vedbas[together.seqs$i[i], together.seqs$t.before.i[i]:together.seqs$t.start[i]], na.rm = TRUE)
}





### Save features object
events <- together.seqs
setwd(processeddir)
save(list=c('events'),file='fission_fusion_features.RData')
