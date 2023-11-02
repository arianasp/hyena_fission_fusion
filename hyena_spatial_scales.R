#What are the relevant spatial scales associated with hyena interactions?

#This script produces the following plots:

#1. Histogram of log dyadic distance between pairs of hyenas (relevant spatial scales)
#2. Mean heading correlation as a function of dyadic distance between hyenas (directional coordination)
#3. Mean difference in speed between pairs of hyenas as a function of distance between them (speed coordination)

#PARAMS
subsamp <- 60 #subsample rate for headings - to speed up computation
speed_dt <- 60 #time step for computing individual speed
fpt_thresh <- 600 #threshold of first passage time over which a heading is not computed (because the hyena is not moving so has undefined heading)
dist_bins <- c(0, 10^seq(0,4,.2)) #distance bins to use (log spaced)

#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'
dir <- '~/Dropbox/hyenas/hyena_fission_fusion/data/raw/' #directory where all data is stored

#FUNCTIONS
#read in functions
setwd(codedir)
source('ff_functions_library.R')

#LOAD DATA
setwd(dir)
load('hyena_xy_level1.RData')
load('hyena_timestamps.Rdata')

#CALCULATE METRICS 
n_inds <- nrow(xs)
n_times <- ncol(xs)
heads <- speeds <- matrix(NA, nrow = n_inds, ncol = n_times)
ts.sub <- seq(1, length(timestamps), subsamp) #time indexes for subsampled data (headings)

#spatially discretized headings of individuals and speeds of individuals
for(i in 1:n_inds){
  print(i)
  heads[i,] <- spatial.headings(x=xs[i,],y=ys[i,],R=10,t.idxs=1:length(timestamps),backward = F, fpt.thresh = fpt_thresh, subsample=subsamp)
  idxs_fut <- seq(speed_dt+1,length(timestamps),1)
  idxs_now <- seq(1,length(timestamps)-speed_dt,1)
  speeds[i,idxs_now] <- sqrt((xs[i, idxs_fut] - xs[i, idxs_now])^2 + (ys[i, idxs_fut] - ys[i, idxs_now])^2)
}

#dyadic distances and heading correlations
dyad_dists <- head_corrs <- speed_diffs <- log_speed_diffs <- array(NA, dim = c(n_inds,n_inds,n_times))
for(i in 1:(n_inds-1)){
  for(j in (i+1):n_inds){
    print(paste(i,j,sep=','))
    dyad_dists[i,j,] <- sqrt((xs[i,]-xs[j,])^2 + (ys[i,]-ys[j,])^2)
    head_corrs[i,j,] <- cos(heads[i,])*cos(heads[j,]) + sin(heads[i,])*sin(heads[j,])
    log_speed_diffs[i,j,] <- abs(log(speeds[i,]) - log(speeds[j,]))
    speed_diffs[i,j,] <- abs(speeds[i,] - speeds[j,])
    idxs_fut <- seq(speed_dt+1,length(timestamps),1)
    idxs_now <- seq(1,length(timestamps)-speed_dt,1)
  }
}

#getting means for headings and speed diffs as a function of distance 
freq_dyad_dists <- heads_mean <- speeds_mean <- log_speeds_mean <- rep(NA, length(dist_bins)-1)
for(i in 1:(length(dist_bins)-1)){
  idxs <- which(dyad_dists >= dist_bins[i] & dyad_dists < dist_bins[i+1])
  heads_mean[i] <- mean(head_corrs[idxs], na.rm=T)
  speeds_mean[i] <- mean(speed_diffs[idxs], na.rm=T)
  log_speeds_mean[i] <- mean(log_speed_diffs[idxs], na.rm=T)
  freq_dyad_dists[i] <- length(idxs)
}

#PLOTS

#middle of the bins for plotting
mids <- (dist_bins[1:(length(dist_bins)-1)] + dist_bins[2:length(dist_bins)])/2

#set up plot
quartz(width=12, height = 5)
par(mfrow=c(1,3), mar = c(5,5,1,1))

#Plot 1: Distribution of log(dyadic distances) between hyenas
frac_dyad_dists <- freq_dyad_dists / sum(freq_dyad_dists)
plot(mids,frac_dyad_dists, xlab = 'Distance apart (m)', ylab = 'Probability', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,.23),cex.lab=2, cex.axis=1.5, log='x')
abline(v=100, col ='blue', lwd = 3, lty=2)
abline(v=200, col='darkgreen', lwd = 3, lty = 2)

#Plot 2: Mean heading correlation as a function of distance 
plot(mids,heads_mean, xlab = 'Distance apart (m)', ylab = 'Mean heading correlation', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,1),cex.lab=2, cex.axis=1.5, log='x')
abline(h=0, lty = 2)
abline(v=100, col = 'blue', lwd = 3, lty = 2)
abline(v=200, col = 'darkgreen', lwd = 3, lty = 2)

#Plot 3: Mean difference in speed as a function of distance 
plot(mids,speeds_mean, xlab = 'Distance apart (m)', ylab = 'Mean absolute speed difference (m/min)', pch = 19, col = '#00000066', cex = 2, ylim = c(0,max(speeds_mean)),cex.lab=2, cex.axis=1.5, log='x')
abline(h=0, lty = 2)
abline(v=100, col = 'blue', lwd = 3, lty = 2)
abline(v=200, col = 'darkgreen', lwd = 3, lty = 2)
