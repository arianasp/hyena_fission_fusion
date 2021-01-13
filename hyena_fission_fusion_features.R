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

### Get locations for relevant time points
x.start.i <- xs[cbind(together.seqs$i, together.seqs$t.start)]
x.before.i <- xs[cbind(together.seqs$i, together.seqs$t.before.i)]
x.end.i <- xs[cbind(together.seqs$i, together.seqs$t.end)]
x.after.i <- xs[cbind(together.seqs$i, together.seqs$t.after.i)]
x.closest.i <- xs[cbind(together.seqs$i, together.seqs$t.closest)]

y.start.i <- ys[cbind(together.seqs$i, together.seqs$t.start)]
y.before.i <- ys[cbind(together.seqs$i, together.seqs$t.before.i)]
y.end.i <- ys[cbind(together.seqs$i, together.seqs$t.end)]
y.after.i <- ys[cbind(together.seqs$i, together.seqs$t.after.i)]
y.closest.i <- ys[cbind(together.seqs$i, together.seqs$t.closest)]

x.start.j <- xs[cbind(together.seqs$j, together.seqs$t.start)]
x.before.j <- xs[cbind(together.seqs$j, together.seqs$t.before.j)]
x.end.j <- xs[cbind(together.seqs$j, together.seqs$t.end)]
x.after.j <- xs[cbind(together.seqs$j, together.seqs$t.after.j)]
x.closest.j <- xs[cbind(together.seqs$j, together.seqs$t.closest)]

y.start.j <- ys[cbind(together.seqs$j, together.seqs$t.start)]
y.before.j <- ys[cbind(together.seqs$j, together.seqs$t.before.j)]
y.end.j <- ys[cbind(together.seqs$j, together.seqs$t.end)]
y.after.j <- ys[cbind(together.seqs$j, together.seqs$t.after.j)]
y.closest.j <- ys[cbind(together.seqs$j, together.seqs$t.closest)]


### Displacement
together.seqs$disp.during.i <- sqrt((x.start.i-x.end.i)^2 + (y.start.i-y.end.i)^2)
together.seqs$disp.before.i <- sqrt((x.before.i-x.start.i)^2 + (y.before.i-y.start.i)^2)
together.seqs$disp.after.i <- sqrt((x.end.i-x.after.i)^2 + (y.end.i-y.after.i)^2)

together.seqs$disp.during.j <- sqrt((x.start.j-x.end.j)^2 + (y.start.j-y.end.j)^2)
together.seqs$disp.before.j <- sqrt((x.before.j-x.start.j)^2 + (y.before.j-y.start.j)^2)
together.seqs$disp.after.j <- sqrt((x.end.j-x.after.j)^2 + (y.end.j-y.after.j)^2)

### Time
together.seqs$duration.before.i <- together.seqs$t.start - together.seqs$t.before.i
together.seqs$duration.after.i <- together.seqs$t.after.i - together.seqs$t.end

together.seqs$duration.during <- together.seqs$t.end - together.seqs$t.start

together.seqs$duration.before.j <- together.seqs$t.start - together.seqs$t.before.j
together.seqs$duration.after.j <- together.seqs$t.after.j - together.seqs$t.end

### Speed
together.seqs.speed.before.i <- together.seqs$disp.before.i / together.seqs$duration.before.i
together.seqs$speed.during.i <- together.seqs$disp.during.i / together.seqs$duration.during
together.seqs$speed.after.i <- together.seqs$disp.after.i / together.seqs$duration.after.i

together.seqs.speed.before.j <- together.seqs$disp.before.j / together.seqs$duration.before.j
together.seqs$speed.during.j <- together.seqs$disp.during.j / together.seqs$duration.during
together.seqs$speed.after.j <- together.seqs$disp.after.j / together.seqs$duration.after.j






r = 1

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

# p <- optim(par = c(fp$x0+5, fp$xf-5, 0), fn = ls_error, fixed.parameters = fp, x = x, y = y, method = "L-BFGS-B",
#            lower = c(fp$x0+1, fp$x0+1, 0), upper = c(fp$xf-1, fp$xf-1, fp$yf))


plot_dyad_dists(r)
lines(x =x, y = fission_fusion_function(x, p$par[1], p$par[2], p$par[3], fp), col = 'red')


together.seqs.exact <- together.seqs[together.seqs$start.exact & together.seqs$end.exact,]
together.seqs.exact[,c('b1', 'b2', 'y.intercept')] <- NA
for(r in 1:nrow(together.seqs.exact)){
  
  ### Remove one case at the very start of the data
  if(is.na(together.seqs.exact$t.start[r]))
    next
  
  y = dyad.dists[together.seqs.exact$i[r], together.seqs.exact$j[r], together.seqs.exact$t.start[r]:together.seqs.exact$t.end[r]]
  x = together.seqs.exact$t.start[r]:together.seqs.exact$t.end[r]
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
  
  
  together.seqs.exact[r,c('b1', 'b2', 'y.intercept')] <-p$par
}











#### Random plots for discussion 2021-01-11
hist(together.seqs$duration, breaks = 1000,
     main = 'Distribution of duration of fusion events', xlab = 'Duration')

hist(together.seqs$disp.start.end.i, breaks = 1000,
     main = 'Distribution of distances travelled', xlab = 'Distance travelled')
hist(together.seqs$disp.start.end.i, xlim = c(0,1000), breaks = 1000,
     main = 'Distribution of distances travelled (zoomed in)', xlab = 'Distance travelled')

hist(together.seqs.exact$y.intercept, breaks = 100)
hist(together.seqs.exact$b2 - together.seqs.exact$b1, breaks = 100)


disp.200 <- which(together.seqs$disp.during.i >= 195 & together.seqs$disp.during.i <= 205)
plot_events(disp.200)

disp.400 <- which(together.seqs$disp.during.i >= 395, together.seqs$disp.during.i <= 405)
plot_events(disp.400)

random.50 <- sample(1:nrow(together.seqs), 50)
plot_events(random.20)

random.10.exact <- sample(1:nrow(together.seqs.exact), 10)
plot_events(10)
plot_canonical_shape(10)

disp.over.500 <- which(together.seqs$disp.during.i >= 500)
plot_events(disp.over.500)



r=200
plot_events(r)

plot_canonical_shape(rows = c(200, 100), dyad.dists = dyad.dists, together.seqs.exact = together.seqs.exact)



y = dyad.dists[together.seqs.exact$i[r], together.seqs.exact$j[r], together.seqs.exact$t.start[r]:together.seqs.exact$t.end[r]]
x = together.seqs.exact$t.start[r]:together.seqs.exact$t.end[r]
plot(y = y, type = 'l',
     x = x,
     xlab = 'Time (s)', ylab = 'Distance between individuals (m)')



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

lines(x =x, y = fission_fusion_function(x, p$par[1], p$par[2], p$par[3], fp), col = 'red')



