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
  #Start
  far.times <- which(prev.dists > R.fission)
  
  ## Exclude events that are too early to have prior data where they are apart
  if(!(length(far.times) > 0))
    next
  together.seqs$t.start[i] <- max(far.times, na.rm = TRUE)
  t.start <- together.seqs$t.start[i]
  together.seqs$start.exact[i] <- !is.na(dyad.dists[ind.i, ind.j, t.start+1])
  #End
  together.seqs$t.end[i] <- tf
  t.end <- tf
  together.seqs$end.exact[i] <- together.seqs[i,'tf.exact']
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
together.seqs$disp.start.end.i <- sqrt((x.start.i-x.end.i)^2 + (y.start.i-y.end.i)^2)
hist(log(together.seqs$disp.start.end.i))
hist(together.seqs$disp.start.end.i, xlim = c(0,1000), breaks = 1000)




#### Random plots for discussion 2021-01-11
hist(together.seqs$t.end - together.seqs$t.start, breaks = 1000,
     main = 'Distribution of duration of fusion events', xlab = 'Duration')

hist(together.seqs$disp.start.end.i, breaks = 1000,
     main = 'Distribution of distances travelled', xlab = 'Distance travelled')
hist(together.seqs$disp.start.end.i, xlim = c(0,1000), breaks = 1000,
     main = 'Distribution of distances travelled (zoomed in)', xlab = 'Distance travelled')

disp.200 <- which(together.seqs$disp.start.end.i >= 195 & together.seqs$disp.start.end.i <= 205)

for(r in disp.200){
  x.i <- xs[together.seqs$i[r], together.seqs$t.start[r]:together.seqs$t.end[r]]
  y.i <- ys[together.seqs$i[r], together.seqs$t.start[r]:together.seqs$t.end[r]]
  x.j <- xs[together.seqs$j[r], together.seqs$t.start[r]:together.seqs$t.end[r]]
  y.j <- ys[together.seqs$j[r], together.seqs$t.start[r]:together.seqs$t.end[r]]
  
  
  plot(x.i, y.i,asp=1, ylim = c(min(c(y.i, y.j), na.rm =TRUE), max(c(y.i, y.j), na.rm = TRUE)),
       xlim = c(min(c(x.i, x.j), na.rm =TRUE), max(c(x.i, x.j), na.rm = TRUE)), col = 'blue')
  points(x.j, y.j, col = 'red')
}



