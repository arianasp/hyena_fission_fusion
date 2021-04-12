#This script generates a figure showing which den was attended by which individual over time

#------PARAMETERS------
datadir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps/'
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'

den.file.path <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/rawdata/metadata/hyena_isolate_dens.csv'
last.day <- 35
den.dist.thresh <- 200 #how close an individual needs to be to be considered "at the den"
local.time.diff <- 3 #offset from UTC in hrs

#-----LIBRARIES-------
library(rhdf5)
library(lubridate)

#-----SOURCE FUNCTIONS------
setwd(codedir)
source('ff_functions_library.R')

#-----DIRECTORY----

#set working directory
setwd(datadir)

#-----LOAD DATA-----

#gps data
load('hyena_xy_level1.RData')
load('hyena_timestamps.RData')
load('hyena_day_start_idxs.RData')
load('hyena_ids.RData')

#dens
den.locs <- get_dens(den.file.path = den.file.path)

#----PROCESS----
#remove data from after day 35
t.idxs <- seq(1, day.start.idxs[last.day+1]-1)
xs <- xs[,t.idxs]
ys <- ys[,t.idxs]
timestamps <- timestamps[t.idxs]

#local time
ts.local <- timestamps + local.time.diff * 60 * 60

#Get distance of each hyena to each den over time
dist.dens <- array(NA,dim=c(dim(xs),nrow(den.locs)))
for(i in 1:nrow(den.locs)){
  dist.dens[,,i] <- sqrt((xs - den.locs$east[i])^2 + (ys - den.locs$north[i])^2)
}

#Den blocks to show on plot
den.blocks <- list()
den.blocks[[1]] <- c(11, 21, 29, 33)
den.blocks[[2]] <- c(12, 25)
den.blocks[[3]] <- c(13)
den.blocks[[4]] <- c(15)
den.blocks[[5]] <- c(11, 30)

#-----PLOT------
#get days and minutes into day for all ts
step <- 60 #subsample to one sample per minute
dist.dens.sub <- dist.dens[,seq(1,length(ts.local),step),]
ts.sub <- ts.local[seq(1,length(ts.local),step)]
days <- as.Date(ts.sub)
uniq.days <- unique(days)
mins <- minute(ts.sub) + hour(ts.sub)*60

dist.dens.thresh <- dist.dens.sub < den.dist.thresh
quartz(width=12,height=6)
par(mfrow=c(1,5),mar=c(5,5,1,0))
cols <- c('green','magenta','blue','red')
for(i in 1:nrow(hyena.ids)){
  plot(NULL,xlim=c(0,24),ylim=c(1,length(uniq.days)),main=hyena.ids$name[i],xlab='hour of day',ylab='day', cex.lab = 2, cex.axis = 1.5)
  for(d in 1:dim(dist.dens.thresh)[3]){
    for(day in 1:length(uniq.days)){
      idxs.curr <- which(days == uniq.days[day])
      mins.curr <- mins[idxs.curr]
      den.presence.curr <- which(dist.dens.thresh[i,idxs.curr,d]==TRUE)
      mins.at.den <- mins.curr[den.presence.curr]
      points(mins.at.den/60,rep(day,length(mins.at.den)),col=cols[d],pch=19,cex=0.3)
    }
    den.blocks.i <- den.blocks[[i]]
    for(j in 1:(length(den.blocks.i))){
      lines(c(0, 12, 12, 24), c(den.blocks.i[j]+1, den.blocks.i[j]+1, den.blocks.i[j], den.blocks.i[j])-.5, lty = 2)
    }
  }
  if(i==4){
    legend('topright',legend=c('Rbend Den','Dave Den','Res Den','Dick Den'),fill=cols)
  }
}
