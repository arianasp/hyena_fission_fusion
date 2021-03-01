#Summary plots that integrate information across multiple permutations

#-----------LIBRARIES-------------
library(ggplot2)
library(gridExtra)
library(RgoogleMaps)

#------------DIRECTORIES------------

#directory of where the original (processed) movement + vedba data is stored (+ metadata on IDs and den locations)
indir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps'

#directory of where to store extracted data for fission-fusion project
outdir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/data/no_synchrony_measures'

#directory to put plots
plotdir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/results'

#directory where code for fission-fusion project is stored
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'
#codedir <- '~/Dropbox/Documents/Research/Partial_projects/hyena_fission_fusion/'

#-----------FILENAMES---------------
data_output_file <- 'fission_fusion_events_features.RData'
nightperm_output_file <- 'hyena_day_randomization_events_features_nightperm_avoidmatchTRUE.RData'
denblock_output_file <- 'hyena_day_randomization_events_features_denblock_avoidmatchTRUE.RData'

#----------FUNCTIONS---------------
setwd(codedir)
source('ff_functions_library.R')

#----------LOAD DATA-------------

print("Loading data")
setwd(indir)
load('hyena_timestamps.RData')
load('hyena_ids.RData')
load('hyena_xy_level1.RData')

setwd(outdir)

#load events from real data
load(data_output_file)
events.data <- events
rm('events')

#load events from nightperm randomization
load(nightperm_output_file)
events.rand.list.nightperm <- events.rand.list
rm('events.rand.list')

#load events from denblock randomization
load(denblock_output_file)
events.rand.list.denblock <- events.rand.list
rm('events.rand.list')

#------------PROCESS-----------------
print("Preprocessing")
n.rands <- length(events.rand.list.denblock)
n.inds <- nrow(hyena.ids)
  
#remove events surrounding the 'day break'
events.data <- remove_events_around_day_breaks(events.data, timestamps, rand.params)
  
#get total number of events in data (total, den, non-den)
events.tot.data <- nrow(events.data)
events.tot.den.data <- sum(events.data$dist.den.start <= params$den.dist.thresh | events.data$dist.den.end <= params$den.dist.thresh)
events.tot.nonden.data <- sum(events.data$dist.den.start > params$den.dist.thresh & events.data$dist.den.end > params$den.dist.thresh)

#get number of events per dyad in data
events.net.data <- count_events_per_dyad(events.data)

#get number of events per dyad in data
den.idxs <- which(events.data$dist.den.start <= params$den.dist.thresh | events.data$dist.den.end <= params$den.dist.thresh)
nonden.idxs <- which(events.data$dist.den.start > params$den.dist.thresh & events.data$dist.den.end > params$den.dist.thresh)
events.net.data.den <- count_events_per_dyad(events.data[den.idxs,])
events.net.data.nonden <- count_events_per_dyad(events.data[nonden.idxs,])
  
#get total number of events in each randomization, and events per dyad
events.tot.nightperm <- events.tot.denblock <- rep(NA, n.rands)
events.tot.den.nightperm <- events.tot.nonden.nightperm <- events.tot.den.denblock <- events.tot.nonden.denblock <- rep(NA, n.rands)
events.net.nightperm <- events.net.denblock <- array(NA, dim = c(n.inds, n.inds, n.rands))
events.net.nightperm.den <- events.net.denblock.den <- array(NA, dim = c(n.inds, n.inds, n.rands))
events.net.nightperm.nonden <- events.net.denblock.nonden <- array(NA, dim = c(n.inds, n.inds, n.rands))

for(r in 1:n.rands){
  
  #get events associated with that randomization
  events.rand.nightperm <- events.rand.list.nightperm[[r]]
  events.rand.denblock <- events.rand.list.denblock[[r]]
  
  #remove events surrounding the 'day break'
  events.rand.nightperm <- remove_events_around_day_breaks(events.rand.nightperm, timestamps, rand.params)
  events.rand.denblock <- remove_events_around_day_breaks(events.rand.denblock, timestamps, rand.params)
  
  #total number of events
  events.tot.nightperm[r] <- nrow(events.rand.nightperm)
  events.tot.denblock[r] <- nrow(events.rand.denblock)
  
  #den vs nonden
  events.tot.den.nightperm[r] <- sum(events.rand.nightperm$dist.den.start <= params$den.dist.thresh | events.rand.nightperm$dist.den.end <= params$den.dist.thresh)
  events.tot.nonden.nightperm[r] <- sum(events.rand.nightperm$dist.den.start > params$den.dist.thresh & events.rand.nightperm$dist.den.end > params$den.dist.thresh)
  events.tot.den.denblock[r] <- sum(events.rand.denblock$dist.den.start <= params$den.dist.thresh | events.rand.denblock$dist.den.end <= params$den.dist.thresh)
  events.tot.nonden.denblock[r] <- sum(events.rand.denblock$dist.den.start > params$den.dist.thresh & events.rand.denblock$dist.den.end > params$den.dist.thresh)
  
  #den indexes
  nightperm.den.idxs <- which(events.rand.nightperm$dist.den.start <= params$den.dist.thresh | events.rand.nightperm$dist.den.end <= params$den.dist.thresh)
  denblock.den.idxs <- which(events.rand.denblock$dist.den.start <= params$den.dist.thresh | events.rand.denblock$dist.den.end <= params$den.dist.thresh)
  nightperm.nonden.idxs <- which(events.rand.nightperm$dist.den.start > params$den.dist.thresh & events.rand.nightperm$dist.den.end > params$den.dist.thresh)
  denblock.nonden.idxs <- which(events.rand.denblock$dist.den.start > params$den.dist.thresh & events.rand.denblock$dist.den.end > params$den.dist.thresh)
  
  #networks - all
  events.net.nightperm[,,r] <- count_events_per_dyad(events.rand.nightperm)
  events.net.denblock[,,r] <- count_events_per_dyad(events.rand.denblock)
  
  #networks - den vs nonden
  events.net.nightperm.den[,,r] <- count_events_per_dyad(events.rand.nightperm[nightperm.den.idxs,])
  events.net.nightperm.nonden[,,r] <- count_events_per_dyad(events.rand.nightperm[nightperm.nonden.idxs,])
  events.net.denblock.den[,,r] <- count_events_per_dyad(events.rand.denblock[denblock.den.idxs,])
  events.net.denblock.nonden[,,r] <- count_events_per_dyad(events.rand.denblock[denblock.nonden.idxs,])
  
}

#--------------PLOTS----------------

print("Creating plots")
print('# Events plot')
#prep for ggplot
plotdat <- data.frame(type = c(rep('nightperm',n.rands),rep('denblock',n.rands)), nevents = c(events.tot.nightperm, events.tot.denblock), denevents = c(events.tot.den.nightperm, events.tot.den.denblock), nondenevents = c(events.tot.nonden.nightperm, events.tot.nonden.denblock))
plotdat$type <- factor(plotdat$type,
                       levels = c('nightperm','denblock'),ordered = TRUE)

#---PLOT 1 - Number of events in real data vs nightperm and denblock randomization
#Subplot 1: All events
ptot <- ggplot(plotdat, aes(x = type, y = nevents)) + 
  geom_violin(aes(fill = type)) + 
  labs(title="",x="", y = "Number of events") +
  ylim(0, 800) + 
  theme_minimal(base_size = 18) + 
  geom_hline(yintercept = events.tot.data, color = 'black', size = 1) + 
  theme(legend.position="none") + 
  scale_x_discrete(breaks=c("nightperm","denblock"), labels=c("Complete shuffle", "Den-blocked shuffle")) + 
  scale_fill_manual(values=c("gray", "gray"))

#Subplot 2: Den events
pden <- ggplot(plotdat, aes(x = type, y = denevents)) + 
  geom_violin(aes(fill = type)) + 
  labs(title="",x="", y = "Number of events") +
  ylim(0, 800) + 
  theme_minimal(base_size = 18) + 
  geom_hline(yintercept = events.tot.den.data, color = 'blue', size = 1) + 
  theme(legend.position="none") + 
  scale_x_discrete(breaks=c("nightperm","denblock"), labels=c("Complete shuffle", "Den-blocked shuffle")) + 
  scale_fill_manual(values=c("blue", "blue"))

#Subplot 3: Non den events
pnonden <- ggplot(plotdat, aes(x = type, y = nondenevents)) + 
  geom_violin(aes(fill = type)) + 
  labs(title="",x="", y = "Number of events") +
  ylim(0, 800) + 
  theme_minimal(base_size = 18) + 
  geom_hline(yintercept = events.tot.nonden.data, color = 'magenta', size = 1) + 
  theme(legend.position="none") + 
  scale_x_discrete(breaks=c("nightperm","denblock"), labels=c("Complete shuffle", "Den-blocked shuffle")) + 
  scale_fill_manual(values=c("magenta", "magenta"))

quartz(width = 16, height = 8)
grid.arrange(ptot, pden, pnonden, nrow = 1)

#---PLOT 2 Comparing social networks generated from the data to the equivalent ones in permuted data
#get netork data into ggplot form (also use SRI here)

print('Edge weights plot')

#data frame of dyads
dyads <- data.frame()
for(i in 1:(n.inds-1)){
  for(j in (i+1):n.inds){
    dyads <- rbind(dyads, c(i,j))
  }
}
colnames(dyads) <- c('i','j')

#data frame of dyads + permutation type
networkdat <- data.frame()
for(r in 1:n.rands){
  tmp <- rbind(dyads, dyads)
  tmp$type <- c(rep('nightperm',nrow(dyads)), rep('denblock',nrow(dyads)))
  tmp$rand <- r
  networkdat <- rbind(networkdat, tmp)
}

tmp <- dyads
tmp$type <- 'real'
tmp$rand <- NA
networkdat <- rbind(networkdat, tmp)

#store dyad data in this data frame for ggplot plotting
networkdat$count <- networkdat$sri <- NA
networkdat$count.den <- networkdat$sri.den <- NA
networkdat$count.nonden <- networkdat$sri.nonden <- NA
for(i in 1:(n.inds-1)){
  for(j in (i+1):n.inds){
    for(r in 1:n.rands){
      #nighpterm
      idx <- which(networkdat$i == i & networkdat$j==j & networkdat$rand == r & networkdat$type == 'nightperm')
      networkdat$count[idx] <- events.net.nightperm[i,j,r]
      networkdat$sri[idx] <- events.net.nightperm[i,j,r] / (sum(events.net.nightperm[i,,r], na.rm=T) + sum(events.net.nightperm[,j,r], na.rm=T))
      networkdat$count.den[idx] <- events.net.nightperm.den[i,j,r]
      networkdat$sri.den[idx] <- events.net.nightperm.den[i,j,r] / (sum(events.net.nightperm.den[i,,r], na.rm=T) + sum(events.net.nightperm.den[,j,r], na.rm=T))
      networkdat$count.nonden[idx] <- events.net.nightperm.nonden[i,j,r]
      networkdat$sri.nonden[idx] <- events.net.nightperm.nonden[i,j,r] / (sum(events.net.nightperm.nonden[i,,r], na.rm=T) + sum(events.net.nightperm.nonden[,j,r], na.rm=T))
      
      #denblock
      idx <- which(networkdat$i == i & networkdat$j==j & networkdat$rand == r & networkdat$type == 'denblock')
      networkdat$count[idx] <- events.net.denblock[i,j,r]
      networkdat$sri[idx] <- events.net.denblock[i,j,r] / (sum(events.net.denblock[i,,r], na.rm=T) + sum(events.net.denblock[,j,r], na.rm=T))
      networkdat$count.den[idx] <- events.net.denblock.den[i,j,r]
      networkdat$sri.den[idx] <- events.net.denblock.den[i,j,r] / (sum(events.net.denblock.den[i,,r], na.rm=T) + sum(events.net.denblock.den[,j,r], na.rm=T))
      networkdat$count.nonden[idx] <- events.net.denblock.nonden[i,j,r]
      networkdat$sri.nonden[idx] <- events.net.denblock.nonden[i,j,r] / (sum(events.net.denblock.nonden[i,,r], na.rm=T) + sum(events.net.denblock.nonden[,j,r], na.rm=T))
      
    }
    
    #real data
    idx <- which(networkdat$i == i & networkdat$j==j & networkdat$type == 'real')
    networkdat$count[idx] <- events.net.data[i,j]
    networkdat$sri[idx] <- events.net.data[i,j] / (sum(events.net.data[i,], na.rm=T) + sum(events.net.data[,j], na.rm=T))
    networkdat$count.den[idx] <- events.net.data.den[i,j]
    networkdat$sri.den[idx] <- events.net.data.den[i,j] / (sum(events.net.data.den[i,], na.rm=T) + sum(events.net.data.den[,j], na.rm=T))
    networkdat$count.nonden[idx] <- events.net.data.nonden[i,j]
    networkdat$sri.nonden[idx] <- events.net.data.nonden[i,j] / (sum(events.net.data.nonden[i,], na.rm=T) + sum(events.net.data.nonden[,j], na.rm=T))
    
  }
}

#get dyad names
networkdat$iname <- networkdat$jname <- NA
for(i in 1:n.inds){
  networkdat$iname[which(networkdat$i == hyena.ids$id[i])] <- hyena.ids$name[i]
  networkdat$jname[which(networkdat$j == hyena.ids$id[i])] <- hyena.ids$name[i]
}
networkdat$dyad <- paste(networkdat$iname, '/', networkdat$jname)

#order
networkdat$type <- factor(networkdat$type,
                       levels = c('nightperm','denblock','real'),ordered = TRUE)

#create the plot - all
psri <- ggplot(networkdat, aes(x = dyad, y = sri)) + 
  geom_violin(aes(fill = type)) + 
  geom_point(aes(x = dyad, y = sri), data = networkdat[which(networkdat$type=='real'),], shape = 23, fill = 'black', size = 4) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(legend.position="top") +
  scale_fill_manual(values=c("yellow", "light green")) +
  ylab('Edge weight') + 
  xlab('')

#create the plot - den
psriden <- ggplot(networkdat, aes(x = dyad, y = sri.den)) + 
  geom_violin(aes(fill = type)) + 
  geom_point(aes(x = dyad, y = sri.den), data = networkdat[which(networkdat$type=='real'),], shape = 23, fill = 'black', size = 4) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(legend.position="top") +
  scale_fill_manual(values=c("yellow", "light green")) +
  ylab('Edge weight') + 
  xlab('')

#create the plot - den
psrinonden <- ggplot(networkdat, aes(x = dyad, y = sri.nonden)) + 
  geom_violin(aes(fill = type)) + 
  geom_point(aes(x = dyad, y = sri.nonden), data = networkdat[which(networkdat$type=='real'),], shape = 23, fill = 'black', size = 4) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(legend.position="top") +
  scale_fill_manual(values=c("yellow", "light green")) +
  ylab('Edge weight') + 
  xlab('')

quartz(width = 16, height = 8)
grid.arrange(psri, psriden, psrinonden, nrow = 1)


#---------PLOT 3------- FF locations
#Where do ff events take place?
good.idxs <- which(events.data$start.exact & events.data$end.exact)
events.data.exact <- events.data[good.idxs,]
events.data.exact$x.fusion <- (xs[cbind(events.data.exact$i, events.data.exact$t.start)] + xs[cbind(events.data.exact$j, events.data.exact$t.start)]) / 2
events.data.exact$y.fusion <- (ys[cbind(events.data.exact$i, events.data.exact$t.start)] + ys[cbind(events.data.exact$j, events.data.exact$t.start)]) / 2
events.data.exact$x.fission <- (xs[cbind(events.data.exact$i, events.data.exact$t.end)] + xs[cbind(events.data.exact$j, events.data.exact$t.end)]) / 2
events.data.exact$y.fission <- (ys[cbind(events.data.exact$i, events.data.exact$t.end)] + ys[cbind(events.data.exact$j, events.data.exact$t.end)]) / 2

lonlat.fusion <- utm.to.latlon(cbind(events.data.exact$x.fusion, events.data.exact$y.fusion), southern_hemisphere = T, utm.zone = '36')
lonlat.fission <- utm.to.latlon(cbind(events.data.exact$x.fission, events.data.exact$y.fission), southern_hemisphere = T, utm.zone = '36')
events.data.exact$lon.fusion <- lonlat.fusion[,1]
events.data.exact$lat.fusion <- lonlat.fusion[,2]

#get den locations
den.locs <- get_dens()

#colors
cols <- rep('#FF00FF66', nrow(events.data.exact))
cols[which(events.data.exact$dist.den.start <= params$den.dist.thresh)] <- '#0000FF66'

#get google map
map <- GetMap(center=c(-1.469269, 35.22431), zoom = 13, maptype = 'satellite')

#plot it
quartz()
PlotOnStaticMap(lat = events.data.exact$lat.fusion, lon = events.data.exact$lon.fusion, MyMap = map, col = cols, pch = 3, add = F)
PlotOnStaticMap(lat = den.locs$lat, lon = den.locs$lon, MyMap = map, col = 'blue', pch = 1, add = T, cex = 4, lwd = 3)



