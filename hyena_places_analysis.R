#This code is an analysis of hyena repeated use sites (i.e. locations that they stay at for extended periods on multiple days)
#vs sites that are only used once (but for an extended period). We are interested in looking at:
# - the birth and death of these "places" over time
# - whether movement patterns when going to / leaving from these places differ from typical movement patterns (e.g. are more directed)

#FILES LOADED
#'hyena_xy_level1.RData'
#'hyena_timestamps.RData'
#'hyena_ids.RData'

#FILES GENERATED
#'first_passage_times.RData' #first passage times
#'hyena_visits_10minthresh_R200.RData' #visits data (also contains place ids and info about their locs)

#PARAMRETERS
#analysis parameters
downsamp.rate <- 60 #downsample to 1 point per downsamp.rate seconds (to save processing time)
Rs <- c(10,50, 100, 200, 500, 1000) #radius to use for first passage time (in meters)
fpt.thresh <- 10*60 #threshold of first passage time to consider something a 'place' (default 10 min = 600 sec)
R.place.idx <- 4 #index to the R value in the vector Rs that is used for identifying a 'place' (default = 4, ie 200 m to be consistent with field definitions)
min.separation.distance <- 100 #places less than this distance apart will be aggregated into one place using the chain rule / DBSCAN

#what to compute
recompute.first.passage.times <- FALSE #if true, recompute fpts, if false, load from file
save.first.passage.times <- FALSE #whether to overwrite fpt saved file
save.visits <- TRUE #whether to overwrite visits data frame
make.plots <- FALSE #whether to show the plots of places and histograms of path directedness when entering or leaving

#Directories and files
basedir <- '~/Dropbox/hyenas/hyena_fission_fusion/'
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'

#File save names
fpt.savename <- 'first_passage_times.RData' #where to store fpt data
visits.savename <- 'hyena_visits_10minthresh_R200.RData' #where to save visits data
den.file <- 'hyena_isolate_dens.csv'
den.names <- c('DAVE D','RBEND D','RES M D1','DICK D')

#--------------LIBRARIES-------------
library(viridis)
library(dbscan)

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

#Identify subgroups using chain rule

#cluster using chain rule and radius R
chain.cluster.groups <- function(x,y,R){
  
  n.inds <- length(x)
  
  groups <- rep(NA,n.inds)
  non.nas <- which(!is.na(x))
  
  r <- cbind(x[non.nas],y[non.nas])
  
  #dyadic dists
  dyad.dists <- dist(r)
  
  #hierarchical clustering to get group identities
  clusts <- hclust(dyad.dists,method='single')
  groups.non.na <- cutree(clusts,h=R)
  groups[non.nas] <- groups.non.na
  
  return(groups)
  
}

#compute directedness of a trajectory, defined as displacement over path length
path.directedness <- function(x, y){
  
  displacement <- sqrt((x[length(x)] - x[1])^2 + (y[length(y)] - y[1])^2)
  path.len <- sum(diff(x)^2 + diff(y)^2)
  
  pd <- displacement / path.len
  
  return(pd)
  
}

#----------------------- MAIN -----------------------------

#SOURCE FUNCTIONS
setwd(codedir)
source('hyena_functions.R')

#LOAD DATA
datadir <- paste(basedir, 'data/raw',sep='')
setwd(datadir)
load('hyena_xy_level1.RData')
load('hyena_timestamps.RData')
load('hyena_ids.RData')

#Read in den data
known.locs <- read.csv(den.file,stringsAsFactors=F)
eastsNorths <- latlon.to.utm(cbind(known.locs$lon,known.locs$lat),southern_hemisphere=T,utm.zone=36)
known.locs$east <- eastsNorths[,1]
known.locs$north <- eastsNorths[,2]
den.locs <- known.locs[which(known.locs$name %in% den.names),]
den.locs$name <- as.character(den.locs$name)

if(recompute.first.passage.times){
  #save fpt parameters in object
  fpt.params <- list()
  fpt.params$downsamp.rate <- downsamp.rate
  fpt.params$Rs <- Rs
  
  #Evaluate first passage time at only 1 point per 10 sec to reduce processing time (don't need that much resolution, and reduce autocorrelation anyway)
  t.idxs.samp <- seq(1,ncol(xs),downsamp.rate)
  
  #numer of individuals and number of times
  xs.sub <- xs[,t.idxs.samp]
  ys.sub <- ys[,t.idxs.samp]
  n.inds <- nrow(xs.sub)
  n.times <- ncol(xs.sub)
  n.Rs <- length(Rs)
  
  #Get all first passage times for the different radii
  fp.times <- array(NA, dim=c(n.inds,length(t.idxs.samp),n.Rs))
  for(r in 1:n.Rs){
    print('R')
    print(r)
    for(i in 1:n.inds){
      print(i)
      for(j in 1:length(t.idxs.samp)){
        fp.times[i,j,r] <- first.passage.time(x = xs.sub[i, j:n.times], y = ys.sub[i, j:n.times], R = Rs[r]) * downsamp.rate
      }
    }
  }
  
  #Make histograms of the first passage times for different scales (R)
  for(ind in 1:5){
  time.bins.log <- seq(0,5,.1) #time bins
  hists <- matrix(NA,nrow=n.Rs,ncol=length(time.bins.log)-1)
  for(i in 1:n.Rs){
    histo <- hist(log(fp.times[ind,,i],base=10),breaks=time.bins.log, plot = F)
    mids <- histo$mids #get midpoints of time bins (will be same for all Rs)
    hists[i,] <- histo$density #get probability densities for each bin (for each R)
  }
  
  quartz(width=10,height = 6)
  plot(NULL,xlim=range(time.bins.log),ylim=c(0,max(hists)),xlab='First passage time (log scale)', ylab = 'Probability density',xaxt='n', cex.lab=1.5,cex.axis=1.5)
  xlabels.sec <- c(10,60,10*60,60*60,12*60*60)
  xlabels.text <- c('10 sec','1 min','10 min','1 hr','10 hr')
  xlocs <- log(xlabels.sec,base=10)
  axis(side=1,at=xlocs, labels =xlabels.text,cex.axis=1.5)
  cols <- viridis(n.Rs)
  for(i in 1:n.Rs){
    lines(mids,hists[i,],col=cols[i],lwd=3)
  }
  legend('topleft',legend = paste(Rs, 'm'),lty=1,col=cols, lwd = 3,cex=1.2)
  }
  
  if(save.first.passage.times){
    setwd(paste(basedir,'data/processed',sep=''))
    save(list=c('fpt.params','fp.times','t.idxs.samp'), file = fpt.savename)
  }
  
} else{
  setwd(paste(basedir,'data/processed',sep=''))
  load(fpt.savename)
  
}


#Identify times when hyena first passage times were > a time threshold
n.inds <- nrow(xs)

places <- data.frame()
for(i in 1:n.inds){
  fp.times.ind <- fp.times[i,,R.place.idx]
  place.time.idxs <- which(fp.times.ind > fpt.thresh)
  tmp.places <- data.frame(ind = rep(i, length(place.time.idxs)), t = t.idxs.samp[place.time.idxs])
  places <- rbind(places, tmp.places)
}

#find contiguous sequences within each individual's data, call them "visits"
visits <- data.frame()
for(i in 1:n.inds){
  places.ind <- places[which(places$ind==i),]
  places.ind <- places.ind[order(places.ind$t),]
  
  t0 <- t.curr <- places.ind$t[1] #start at beginning
  for(j in 1:nrow(places.ind)){
    dt <- places.ind$t[j] - t.curr
    if(dt <= downsamp.rate){
      t.curr <- places.ind$t[j]
    } else{
      visits <- rbind(visits, data.frame(i = i, t0 = t0, tf = t.curr))
      t0 <- t.curr <- places.ind$t[j]
    }
  }
  
}

visits$dur <- visits$tf - visits$t0

#get "place" locations (mean location for the time they were in the place)
visits$east <- visits$north <- NA
for(i in 1:nrow(visits)){
  xp <- xs[visits$i[i], visits$t0[i]:visits$tf[i]]
  yp <- ys[visits$i[i], visits$t0[i]:visits$tf[i]]
  visits$east[i] <- mean(xp, na.rm=T)
  visits$north[i] <- mean(yp, na.rm=T)
  
}


# #NEW: Find clusters using DBSCAN (testing out)
# 
# cl <- dbscan(cbind(visits$east,visits$north),eps=min.separation.distance)
# cols <- rainbow(max(cl$cluster))
# cols <- c('black',cols)
# #plot(visits$east,visits$north,col=cols[cl$cluster+1],pch=19,cex=0.1,asp=1, )
# 
# #get x and y coordinates associated with each place
# visits$place.id <- cl$cluster
# 
# #cluster 0 is the 'noise' cluster but in this case these should be considered as separate places
# max.cluster.label <- max(visits$place.id)
# single.point.clusts <- which(visits$place.id==0)
# visits$place.id[single.point.clusts] <- seq(max.cluster.label+1,max.cluster.label + length(single.point.clusts))
# places <- data.frame(place.id = unique(visits$place.id))
# places$hull <- NA
# places$area <- NA
# for(p in 1:nrow(places)){
#   place <- places$place.id[p]
#   place.visit.idxs <- which(visits$place.id==place)
#   if(length(place.visit.idxs) >= 3){
#     place.points <- as.matrix(visits[place.visit.idxs,c('east','north')])
#     shape <- ashape(place.points, alpha = min.separation.distance)
#     edges <- shape$edges
#     places$hull[p] <- list(shape)
#   }
# }
# 
# #PLOT PLACES
# quartz(width = 26, height = 10)
# par(mfrow=c(2,3))
# for(i in 1:nrow(hyena.ids)){
# plot(NULL, asp=1, xlim = c(742000,758000),ylim=c(9830000,9843000))
# cols <- sample(rainbow(nrow(places)))
# for(p in 1:nrow(places)){
#   if(!is.na(places$hull[p])){
#     plot(places$hull[p][[1]], add = T, col='black',pch=19,cex=0.01, lwd = 2)
#   }
# }
# ind.idxs <- which(visits$i==i)
# points(visits$east[ind.idxs],visits$north[ind.idxs],col=hyena.ids$color[i],pch=19,cex=0.1)
# }
# 

#cluster places with x distance of one another
visits$place.id <- NA
chain.clusts <- chain.cluster.groups(x = visits$east, y = visits$north, R = min.separation.distance)
visits$place.id <- chain.clusts

#for each cluster of "sub-places" find the mean location, define this as the center of the overall "place"
#also get the minimum and maximum of the x and y values to potentially create a bounding box around the place later
place.ids <- unique(visits$place.id)
visits$place.mean.east <- visits$place.mean.north <- NA
visits$place.min.east <- visits$place.min.north <- NA
visits$place.max.east <- visits$place.max.north <- NA
for(i in place.ids){
  rows <- which(visits$place.id == i)
  xp <- visits$east[rows]
  yp <- visits$north[rows]

  visits$place.mean.east[rows] <- mean(xp, na.rm=T)
  visits$place.mean.north[rows] <- mean(yp, na.rm=T)

  visits$place.min.east[rows] <- min(xp, na.rm=T)
  visits$place.min.north[rows] <- min(yp, na.rm=T)

  visits$place.max.east[rows] <- max(xp, na.rm=T)
  visits$place.max.north[rows] <- max(yp, na.rm=T)

}

#plot all the places on one plot
#size = number of times place was visited
#color = number of unique hyenas visiting the place
if(make.plots){
  quartz(width = 12, height = 12)
  plot(NULL, xlim = range(visits$place.mean.east,na.rm=T), ylim = range(visits$place.mean.north,na.rm=T), asp = 1, xlab = 'Easting (m)', ylab = 'Northing (m)')
  cols <- c('black','blue','dark green','orange','red')
  for(i in place.ids){
    rows <- which(visits$place.id == i)
    x <- visits$place.mean.east[rows[1]]
    y <- visits$place.mean.north[rows[1]]
    xleft <- visits$place.min.east[rows[1]]
    xright <- visits$place.max.east[rows[1]]
    ybottom <- visits$place.min.north[rows[1]]
    ytop <- visits$place.max.north[rows[1]]
    #rect(xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, col = '#00000022', border=NA)
    col <- cols[length(unique(visits$i[rows]))]
    n.stays <- length(rows)
    size <- log(n.stays + 1)/3
    points(x, y, cex = size, col = col, pch = 19)
  }
  for(d in 1:nrow(den.locs)){
    points(den.locs$east[d], den.locs$north[d], cex = 1, col = 'black')
    text(x = den.locs$east[d], y = den.locs$north[d], den.locs$name[d], col = 'black', pos = 4)
  }
  legend('topleft',fill=cols,legend=c('1 hyena', '2 hyenas','3 hyenas','4 hyenas','5 hyenas'))
  legend('topright',legend=c('1 visit','3 visits','10 visits','30 visits','100 visits','300 visits'), col=rep('black',6),pt.cex = log(c(1,3,10,30,100,300)+1)/3, pch = rep(19,6))
  
  #add previous hour trajectories to the map
  # nsamps <- 200
  # for(i in sample(1:nrow(visits),nsamps)){
  #   ind <- visits$i[i]
  #   tprev <- seq(visits$t0[i] - 60*60, visits$t0[i])
  #   taft <- seq(visits$tf[i], visits$tf[i] + 10*60)
  #   if(sum(tprev<1)==0){
  #     lines(xs[ind,tprev], ys[ind,tprev],col='gray',lwd=.5)
  #   }
  # }
}

#calculate some metrics for the 10 min before getting to the place
visits$pd.prev10 <- visits$pd.prev60 <- NA
for(i in 1:nrow(visits)){
  ind <- visits$i[i]
  t.arrive <- visits$t0[i]
  t.leave <- visits$tf[i] + 60*60
  
  if(t.arrive > 60*10){
    visits$pd.prev10[i] <- path.directedness(xs[ind, seq(t.arrive - 60*10, t.arrive)], ys[ind, seq(t.arrive - 60*10, t.arrive)])
  }
  if(t.arrive > 60*60){
    visits$pd.prev60[i] <- path.directedness(xs[ind, seq(t.arrive - 60*60, t.arrive)], ys[ind, seq(t.arrive - 60*60, t.arrive)])
  }
  
  if(t.arrive > 60*10){
    visits$pd.next10[i] <- path.directedness(xs[ind, seq(t.leave, t.leave + 60*10)], ys[ind, seq(t.leave, t.leave + 60*10)])
  }
  if(t.arrive > 60*60){
    visits$pd.next60[i] <- path.directedness(xs[ind, seq(t.leave, t.leave + 60*60)], ys[ind, seq(t.leave, t.leave + 60*60)])
  }
}

#get number of visits each place has had
place.n.visits <- table(visits$place.id)
places.sorted <- sort(place.n.visits,decreasing=T)

#NOTE: ASSUMING DENS ARE THE TOP 3 MOST VISITED PLACES - this assumption valid for the current parameters, however it should be checked if anything is changed!
popular.nondens <- names(places.sorted)[places.sorted >= 10]
dens <- names(places.sorted)[1:3]
rest <- names(places.sorted)[places.sorted < 10]

visits$type <- NA
visits$type[which(visits$place.id %in% popular.nondens)] <- 'popular.nonden'
visits$type[which(visits$place.id %in% dens)] <- 'den'
visits$type[which(visits$place.id %in% rest)] <- 'nonden'

#make a histogram of path directedness for dens, popular non-dens (>10 visits), and other places (< 10 visits)
if(make.plots){
  quartz(width=8,height=8)
  bandwidth <- .1
  plot(NULL,xlim=c(0,1),ylim=c(0,2.5),xlab = 'Path directedness',ylab='Probability density',cex.lab=1.5, cex.axis=1.5)
  kd.dens <- density(visits$pd.next10[which(visits$type=='den')],na.rm=T,from=0,to=1, bw=bandwidth)
  kd.popular.nondens <- density(visits$pd.next10[which(visits$type=='popular.nonden')],na.rm=T,from=0,to=1, bw=bandwidth)
  kd.nondens <- density(visits$pd.next10[which(visits$type=='nonden')],na.rm=T,from=0,to=1, bw=bandwidth)
  lines(kd.dens,col='red',lwd=3)
  lines(kd.popular.nondens, col='blue',lwd=3)
  lines(kd.nondens,col='black',lwd=3)
  legend('topleft',legend=c('Dens','> 10 visits (non-den)','< 10 visits'), col=c('red','blue','black'),lwd=3,cex=1.5)
}

#save hyena places
if(save.visits){
  save(file=visits.savename, list = c('visits'))
}
