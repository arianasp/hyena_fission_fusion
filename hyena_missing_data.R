#CHARACTERIZING HYENA MISSING DATA

#This script generates some plots describing the patterns of missing GPS data in our hyena collaring dataset from 2017

#-----DIRECTORIES------
gps.dir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps'
ff.dir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/data/main_output/'
plot.dir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/results/missing_data_plots/'

#-----LIBRARIES----
library(lubridate)
library(viridis)

#------LOAD---------
setwd(gps.dir)
load('hyena_xy_level1.RData')
load('hyena_timestamps.Rdata')
load('hyena_ids.Rdata')

setwd(ff.dir)
load('fission_fusion_events_features.RData')

#-----PARAMETERS----
local.time.diff <- 3 #difference of local time from UTC (3 hours)

#-----PROCESS----
#Create a data frame to hold missing data information
missing.df <- data.frame()
for(i in 1:nrow(xs)){
  missing.df <- rbind(missing.df, data.frame(id = i, time.utc = timestamps, missing = is.na(xs[i,])))
}
missing.df$time.local <- missing.df$time.utc + local.time.diff*60*60
missing.df$date <- as.Date(missing.df$time.local)
missing.df$hour <- hour(missing.df$time.local)



#-----PLOTS-----

#---Set plotting directory---
setwd(plot.dir)

#---Plot 0: GPS coverage over time per individual---

#subsample for the purposes of visualization
step <- 60
cols <- viridis(nrow(hyena.ids))
xs.sub <- xs[,seq(1, ncol(xs), step)]
tracked <- !is.na(xs.sub)
mins.into.day <- hour(timestamps[seq(1, ncol(xs), step)] + local.time.diff*60*60) * 60 + minute(timestamps[seq(1, ncol(xs), step)] + local.time.diff*60*60)
day <- as.Date(timestamps[seq(1, ncol(xs), step)] + local.time.diff*60*60)
day <- day - min(day)

png(filename = 'missingdata_0_allcoverage.png',width = 1400, height = 600, units = 'px')
par(mfrow=c(1,nrow(hyena.ids)), mar = c(5,5,4,1))
for(i in 1:nrow(hyena.ids)){
  plot(mins.into.day/60, day*tracked[i,], pch = 19, cex = 0.1, col = cols[i], xlab = 'Hour of day', ylab = 'Day', cex.lab = 2, cex.axis = 1.5, main = hyena.ids$name[i])
}
dev.off()

#---Plot 1: How much missing data per individual?---
png(filename = 'missingdata_1_missingperind.png',width = 600, height = 600, units = 'px')
par(mar=c(5,5,1,1))
missing.per.ind <- rowMeans(is.na(xs)) * 100
plot(1:nrow(hyena.ids), missing.per.ind, ylim = c(0,50), xaxt = 'n', pch = 19, xlab = "Individual", ylab = '% Missing data', cex = 1.5, cex.lab = 2, cex.axis=1.5)
axis(1, at = 1:nrow(hyena.ids), labels = hyena.ids$name, cex.axis = 1.5)
dev.off()

#---Plot 2: Missing data by day per individual---
missing.by.date <- aggregate(x = missing.df$missing, by = list(missing.df$id, missing.df$date), FUN = mean, na.action = na.omit)
png(filename = 'missingdata_2_missingperind_byday.png',width = 1400, height = 600, units = 'px')
par(mar=c(8,5,1,1))
plot(NULL, xlim = range(as.Date(missing.by.date$Group.2)), ylim = c(0,100), xaxt = 'n', pch = 19, xlab = "", ylab = '% Missing data', cex = 1.5, cex.lab = 2, cex.axis=1.5)
axis(side = 1, at = unique(missing.by.date$Group.2), labels = unique(missing.by.date$Group.2), las =2 )
legend('topleft', legend = hyena.ids$name, col = cols, pch = rep(19, nrow(hyena.ids)))
for(i in 1:nrow(hyena.ids)){
  curr <- missing.by.date[which(missing.by.date$Group.1==i),]
  lines(as.Date(curr$Group.2), curr$x*100, col = cols[i], lwd = 2)
  points(as.Date(curr$Group.2), curr$x*100, col = cols[i], pch = 19)
}
dev.off()

#---Plot 3: How much missing data vs hour of the day?---
missing.by.hour <- aggregate(x = missing.df$missing, by = list(missing.df$id, missing.df$hour), FUN = mean, na.action = na.omit)
png(filename = 'missingdata_3_missingperhour.png',width = 1400, height = 600, units = 'px')
par(mar=c(8,5,1,1))
plot(NULL, xlim = range(missing.by.hour$Group.2), ylim = c(0,100), xaxt = 'n', pch = 19, xlab = "Hour of day", ylab = '% Missing data', cex = 1.5, cex.lab = 2, cex.axis=1.5)
axis(side = 1, at = unique(missing.by.hour$Group.2), labels = unique(missing.by.hour$Group.2),  cex.axis=1.5)
cols <- viridis(nrow(hyena.ids))
legend('topleft', legend = hyena.ids$name, col = cols, pch = rep(19, nrow(hyena.ids)))
for(i in 1:nrow(hyena.ids)){
  curr <- missing.by.hour[which(missing.by.date$Group.1==i),]
  lines(curr$Group.2, curr$x*100, col = cols[i], lwd = 2)
  points(curr$Group.2, curr$x*100, col = cols[i], pch = 19)
}
dev.off()

#---Plot 4: How much missing dyadic data vs hour of the day?---
missing.mat <- is.na(xs)
missing.dyads <- rep(0, ncol(xs))
for(i in 1:(nrow(xs)-1)){
  for(j in (i+1):nrow(xs)){
    missing.dyads <- missing.dyads + (missing.mat[i,] | missing.mat[j,])
  }
}
hours <- hour(timestamps + local.time.diff*60*60)
missing.dyad.df <- data.frame(time.local = timestamps+local.time.diff*60*60, hour = hour(timestamps + local.time.diff*60*60), missing.dyads = missing.dyads)
missing.dyad.df$dyads.present <- nrow(hyena.ids)*(nrow(hyena.ids)-1)/2 - missing.dyad.df$missing.dyads
missing.dyad.df$frac.dyads.present <- missing.dyad.df$dyads.present / (nrow(hyena.ids)*(nrow(hyena.ids)-1)/2)

dyad.present.by.hour <- aggregate(missing.dyad.df$frac.dyads.present, by = list(missing.dyad.df$hour), FUN = mean, na.action = na.omit)

png(filename = 'missingdata_4_missingperhour_dyadic.png',width = 600, height = 600, units = 'px')
par(mar=c(5,5,1,1))
plot(dyad.present.by.hour$Group.1, dyad.present.by.hour$x*100, xlim = range(dyad.present.by.hour$Group.1), ylim = c(0,100), pch = 19, xlab = "Hour of day", ylab = '% Dyads tracked', cex = 1.5, cex.lab = 2, cex.axis=1.5)
dev.off()

#---Plot 5: How much time in a ff event as a function % time dyads tracked---

#get start and end times for each ff event
events$start.time.local <- timestamps[events$t.start] + local.time.diff*60*60
events$end.time.local <- timestamps[events$t.end] + local.time.diff*60*60

#loop over events to compute how many seconds spent in ff events per hour, across all dyads
hours <- seq(0,23)
ff.time <- rep(0, length(hours))
for(i in 1:nrow(events)){
  event.start <- events$start.time.local[i]
  event.end <- events$end.time.local[i]
  tseq <- seq(event.start,event.end,1)
  event.hrs <- hour(tseq)
  secs.within.hr <- table(event.hrs)
  hrs.in.event <- as.numeric(names(secs.within.hr))
  ff.time[hrs.in.event+1] <- ff.time[hrs.in.event+1] + secs.within.hr
}

#how many dyad-seconds are there in each hour window?
dyad.present.by.hour.sec <- aggregate(missing.dyad.df$dyads.present, by = list(missing.dyad.df$hour), FUN = function(x){return(sum(x))})

#of the time dyads were tracked, what fraction of that time was spent in a ff event?
frac.time.ff <- ff.time / dyad.present.by.hour.sec$x

#make a plot
png(filename = 'missingdata_5_ff_time_over_tracked_time.png',width = 600, height = 600, units = 'px')
par(mar=c(5,5,1,1))
plot(hours, frac.time.ff*100, xlim = range(hours), ylim = c(0,40), pch = 19, xlab = "Hour of day", ylab = 'Time in FF event / Total time dyads tracked', cex = 1.5, cex.lab = 1.5, cex.axis=1.5)
dev.off()

#---Plot 6: Are GPS outages correlated across individuals within the same day?---

missing.by.date <- missing.by.date[which(!as.character(missing.by.date$Group.2) %in% c('2017-02-05','2017-02-06','2017-02-06','2017-02-07','2017-02-08','2017-02-09','2017-02-10')),]

#sum of variances within a day
within.day.variance <- sum(aggregate(missing.by.date$x, by = list(missing.by.date$Group.2), FUN = var)$x)

#sum of variance within a day for shuffled days
n.rands <- 1000
within.day.variance.shuff <- rep(NA,n.rands)

for(r in 1:n.rands){
  missing.by.date.shuff <- missing.by.date
  for(i in 1:nrow(hyena.ids)){
    idxs <- which(missing.by.date$Group.1==i)
    shuff.idxs <- sample(idxs)
    missing.by.date.shuff$x[idxs] <- missing.by.date$x[shuff.idxs]
  }
  within.day.variance.shuff[r] <- sum(aggregate(missing.by.date.shuff$x, by = list(missing.by.date.shuff$Group.2), FUN = var)$x)
}

#compare variance in fraction of time missing data within days vs for randomly-shuffled days
png(filename = 'missingdata_6_withindaycorr.png',width = 600, height = 600, units = 'px')
hist(within.day.variance.shuff, xlab = 'Sum of within-day variance in fraction of missing data across individuals', ylab = 'Frequency', main = '')
abline(v = within.day.variance, col = 'red', lwd = 2)
dev.off()

#---Plot 7: Are GPS outages correlated from one day to the next?---
png(filename = 'missingdata_7_autocorr.png',width = 600, height = 600, units = 'px')
par(mar=c(5,5,1,1))
plot(NULL, xlim = c(0,1), ylim = c(0,1), xlab = 'Fraction missing on day t', ylab = 'Fraction missing on day t + 1', cex.lab = 2, cex.axis = 1.5)
for(i in 1:nrow(hyena.ids)){
  idxs <- which(missing.by.date$Group.1 == i)
  curr <- missing.by.date[idxs,]
  points(curr$x[1:(length(curr$x)-1)], curr$x[2:length(curr$x)],pch = 19, col = cols[i])
}
legend('topleft', legend = hyena.ids$name, col = cols, pch = rep(19, nrow(hyena.ids)))
dev.off()
