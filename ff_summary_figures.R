#Summary plots that integrate information across multiple permutations

#-----------LIBRARIES-------------
library(ggplot2)
library(gridExtra)
library(ggmap)
library(RColorBrewer)
library(lubridate)
library(viridis)
library(alluvial)
library(dplyr)
library(patchwork)
library(magick)


#------------DIRECTORIES------------

user <- Sys.info()['user']
if(user == 'strau'){
  remote.stem <- 'Z:\\'
  code.stem <- '~/code/'
  den.file.path <- 'Z:\\hyena/archive/hyena_pilot_2017/rawdata/metadata/hyena_isolate_dens.csv' 
}else if(user == 'straussed'){
  remote.stem <- '/Volumes/'
  code.stem <- '~/../Dropbox/Documents/Research/Partial_projects/'
  den.file.path <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/rawdata/metadata/hyena_isolate_dens.csv' 
}else{
  remote.stem <- '/Volumes/EAS_shared/'
  code.stem <- '~/Dropbox/code_ari/'
  den.file.path <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/rawdata/metadata/hyena_isolate_dens.csv' 
}

#directory of where the original (processed) movement + vedba data is stored (+ metadata on IDs and den locations)
indir <- paste0(remote.stem, 'hyena/archive/hyena_pilot_2017/processed/gps')

#directory of where to store extracted data for fission-fusion project
outdir <- paste0(remote.stem, 'hyena/working/hyena_fission_fusion/data/main_output')

#directory where the satelite map is stored
#Load pre-downloaded map
mapdir <- paste0(remote.stem, 'hyena/working/hyena_fission_fusion/data/')

#directory to put plots
plotdir <- paste0(remote.stem, 'hyena/working/hyena_fission_fusion/results')

#directory where code for fission-fusion project is stored
codedir <- paste0(code.stem, 'hyena_fission_fusion/')

#-----------FILENAMES---------------
data_output_file <- 'fission_fusion_events_features.RData'
nightperm_output_file <- 'hyena_day_randomization_events_features_nightperm_avoidmatchTRUE.RData'
denblock_output_file <- 'hyena_day_randomization_events_features_denblock_avoidmatchTRUE.RData'
map_file <- 'hyena_satellite_map.RData'

#----------FUNCTIONS---------------
setwd(codedir)
source('ff_functions_library.R')

#----------LOAD DATA-------------

print("Loading data")
setwd(indir)
load('hyena_timestamps.RData')
load('hyena_ids.RData')
load('hyena_xy_level1.RData')
load('hyena_latlon_level0.RData')

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

#load the map
setwd(mapdir)
load(map_file)

# Load hyena image
hyena <- image_read(path = paste0(remote.stem, 'hyena/working/hyena_fission_fusion/collared_hyena.jpg'))

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
  
  #get events associated with that randomization\
  if(exists('events.rand.list.nightperm')){
    events.rand.nightperm <- events.rand.list.nightperm[[r]]
  }
  events.rand.denblock <- events.rand.list.denblock[[r]]
  
  #remove events surrounding the 'day break'
  if(exists('events.rand.list.nightperm')){
    events.rand.nightperm <- remove_events_around_day_breaks(events.rand.nightperm, timestamps, rand.params)
  }
  events.rand.denblock <- remove_events_around_day_breaks(events.rand.denblock, timestamps, rand.params)
  
  #total number of events
  if(exists('events.rand.list.nightperm')){
    events.tot.nightperm[r] <- nrow(events.rand.nightperm)
  }
  events.tot.denblock[r] <- nrow(events.rand.denblock)
  
  #den vs nonden 
  if(exists('events.rand.list.nightperm')){
    events.tot.den.nightperm[r] <- sum(events.rand.nightperm$dist.den.start <= params$den.dist.thresh | events.rand.nightperm$dist.den.end <= params$den.dist.thresh)
    events.tot.nonden.nightperm[r] <- sum(events.rand.nightperm$dist.den.start > params$den.dist.thresh & events.rand.nightperm$dist.den.end > params$den.dist.thresh)
  }
  events.tot.den.denblock[r] <- sum(events.rand.denblock$dist.den.start <= params$den.dist.thresh | events.rand.denblock$dist.den.end <= params$den.dist.thresh)
  events.tot.nonden.denblock[r] <- sum(events.rand.denblock$dist.den.start > params$den.dist.thresh & events.rand.denblock$dist.den.end > params$den.dist.thresh)
  
  #den indexes
  if(exists('events.rand.list.nightperm')){
    nightperm.den.idxs <- which(events.rand.nightperm$dist.den.start <= params$den.dist.thresh | events.rand.nightperm$dist.den.end <= params$den.dist.thresh)
    nightperm.nonden.idxs <- which(events.rand.nightperm$dist.den.start > params$den.dist.thresh & events.rand.nightperm$dist.den.end > params$den.dist.thresh)
  }  
  denblock.den.idxs <- which(events.rand.denblock$dist.den.start <= params$den.dist.thresh | events.rand.denblock$dist.den.end <= params$den.dist.thresh)
  denblock.nonden.idxs <- which(events.rand.denblock$dist.den.start > params$den.dist.thresh & events.rand.denblock$dist.den.end > params$den.dist.thresh)
  
  #networks - all
  if(exists('events.rand.list.nightperm')){
    events.net.nightperm[,,r] <- count_events_per_dyad(events.rand.nightperm)
  }
  events.net.denblock[,,r] <- count_events_per_dyad(events.rand.denblock)
  
  #networks - den vs nonden
  if(exists('events.rand.list.nightperm')){
    events.net.nightperm.den[,,r] <- count_events_per_dyad(events.rand.nightperm[nightperm.den.idxs,])
    events.net.nightperm.nonden[,,r] <- count_events_per_dyad(events.rand.nightperm[nightperm.nonden.idxs,])
  }
  events.net.denblock.den[,,r] <- count_events_per_dyad(events.rand.denblock[denblock.den.idxs,])
  events.net.denblock.nonden[,,r] <- count_events_per_dyad(events.rand.denblock[denblock.nonden.idxs,])
  
}

# Pre-processing for PLOT 1 and FIGURE 3a ------
#prep for ggplot
plotdat <- data.frame(type = c(rep('nightperm',n.rands),rep('denblock',n.rands)), nevents = c(events.tot.nightperm, events.tot.denblock), denevents = c(events.tot.den.nightperm, events.tot.den.denblock), nondenevents = c(events.tot.nonden.nightperm, events.tot.nonden.denblock))
plotdat$type <- factor(plotdat$type,
                       levels = c('nightperm','denblock'),ordered = TRUE)

# Pre-processing for FIGURE 5 ------
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

#------ Pre-processing for FIGURE 1, FIGURE 2
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
events.data.exact$lon.fission <- lonlat.fission[,1]
events.data.exact$lat.fission <- lonlat.fission[,2]
events.data.exact$at.den.start <- events.data.exact$dist.den.start <= params$den.dist.thresh
events.data.exact$at.end <- events.data.exact$dist.den.end <= params$den.dist.thresh

#get den locations
den.locs <- get_dens(den.file.path)

#-------------DESCRIPTIVE STATS---------
#Compute some basic stats mentioned in the paper and print the results
print(paste0('The total number of events is ', nrow(events.data)))
print(paste0('The total number of events with exact start and end times is ', sum(events.data$start.exact==T & events.data$end.exact==T)))

print(paste0('The number of events starting at the den is ', sum(events.data$dist.den.start <= params$den.dist.thresh)))
print(paste0('The number of events ending at the den is ', sum(events.data$dist.den.end <= params$den.dist.thresh)))
print(paste0('The number of events either starting or ending at a den is ', sum(events.data$dist.den.start <= params$den.dist.thresh | events.data$dist.den.end <= params$den.dist.thresh)))


print(paste0('The mean and 95% range of predicted number of events in denblock model is ', median(events.tot.denblock), ' (', quantile(events.tot.denblock, 0.025), ' to ', quantile(events.tot.denblock, 0.975), ')'))
print(paste0('The mean and 95% range of predicted number of den events in denblock model is ', median(events.tot.den.denblock), ' (', quantile(events.tot.den.denblock, 0.025), ' to ', quantile(events.tot.den.denblock, 0.975), ')'))
print(paste0('The mean and 95% range of predicted number of nonden events in denblock model is ', median(events.tot.nonden.denblock), ' (', quantile(events.tot.nonden.denblock, 0.025), ' to ', quantile(events.tot.nonden.denblock, 0.975), ')'))

print(paste0('The number of traveling events during the together phase is ', sum(events.data$together.type=='together.travel', na.rm=T)))
print(paste0('The number of stationary events during the together phase is ', sum(events.data$together.type=='together.local', na.rm=T)))

#--------------PLOTS----------------

print("Creating plots")
print('# Events plot')

colors <- c("#F72585", "#3f37c9", "#21054C", "#4895EF", "#ba0c2f",'gray30')

#-------FIGURE 1b time of fusions ----------
events.data.exact$hour.start <- hour(timestamps[events.data.exact$t.start] + params$local.time.diff*60*60) 
events.data.exact$hour.end <- hour(timestamps[events.data.exact$t.end] + params$local.time.diff*60*60)
pal <- colorRampPalette(colors = c('#292965','#6696ba','#e2e38b','#e7a553','#7e4b68','#292965'))
breaks <- seq(0,24,1)
timeplot <- ggplot(aes(x = hour.start, fill = as.logical(at.den.start-1)), data = events.data.exact) +
  geom_bar(position = 'stack', na.rm=T) + 
  theme_classic(base_size = 12) + 
  ylab('Number of fusions') + 
  xlab('Hour of day')+
  scale_fill_manual(values = c(colors[2], colors[1]), labels = c('Den', 'Non-Den'))+
  theme(legend.position = c(.75,.75), legend.title = element_blank())+
  labs(tag='B')

timeplot

################################################################################
##### I feel like this belongs in a different script but I'm leaving it for now
#check amount of time two inds tracked (to make sure this isn't driving time of day pattern - it is not)
hours.ts <- hour(timestamps + params$local.time.diff*60*60)
all.hrs <- c()
for(i in 1:(n.inds-1)){
  for(j in (i:n.inds)){
    both.tracked <- which(!is.na(xs[i,] &!is.na(xs[j,])))
    all.hrs <- c(all.hrs, hours.ts[both.tracked])
  }
}

#quartz()
hist(all.hrs)
################################################################################


#----------------------------------FIGURE 1-------------------------------------
#Where do ff events take place?

#colors
cols <- rep(alpha(colors[1], 0.5), nrow(events.data.exact))
cols[which(events.data.exact$dist.den.start <= params$den.dist.thresh)] <- alpha(colors[2], 0.5)

#make background slightly transparent
map_attr <- attributes(hyena_map13)
hyena_map_transparent <- matrix(adjustcolor(hyena_map13, 
                                                alpha.f = 0.6), 
                                    nrow = nrow(hyena_map13))
attributes(hyena_map_transparent) <- map_attr

#create data for a scale bar in lower left corner
low <- attributes(hyena_map13)$bb$ur.lat
left <- attributes(hyena_map13)$bb$ur.lon
ur.utm <- latlon.to.utm(cbind(left, low), utm.zone = '36', southern_hemisphere = T)
scale.min.x <- ur.utm[1,1] - 2000
scale.min.y <- ur.utm[1,2] - 1000
scale.max.x <- scale.min.x + 1000
scale.max.y <- scale.min.y
scalemin.lonlat <- utm.to.latlon(cbind(scale.min.x, scale.min.y), utm.zone = '36', southern_hemisphere = T)
scalemax.lonlat <- utm.to.latlon(cbind(scale.max.x, scale.max.y), utm.zone = '36', southern_hemisphere = T)
scale.lons <- c(scalemin.lonlat[,1], scalemax.lonlat[,1])
scale.lats <- c(scalemin.lonlat[,2], scalemax.lonlat[,2])
scalebar <- data.frame(lon = scale.lons, lat = scale.lats)
scalebar2 <- data.frame(lon = mean(scale.lons), lat = mean(scale.lats))

#create the plot (for fusions, in this case)
#quartz(height = 8, width = 8)
mapplot <- ggmap(hyena_map_transparent, alpha = 0.5) + 
  geom_point(aes(x = lon.fusion, y = lat.fusion, color = at.den.start), data = events.data.exact, size = 1, shape = 3, alpha = 0.8) +
  scale_color_manual(values=c(colors[1], colors[2])) + 
  theme(legend.position="none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_point(aes(x = lon, y = lat), data = den.locs, size = 7, shape = 21, color = colors[2], stroke = 1) +
  geom_line(aes(x = lon, y = lat), data = scalebar, color = 'black') +
  geom_text(aes(x = lon, y = lat), data = scalebar2, label = '1 km', nudge_y = .003)+
  labs(tag = 'C')

hyena.plot <- image_ggplot(hyena)+labs(tag = 'A')+theme(axis.text = element_blank(), plot.margin = unit(c(0,0,0,0), units = 'cm'))

layout <- "
AACCCC
AACCCC
BBCCCC
BBCCCC
"

png(paste0(plotdir, '/FIG1.png'), width = 6, height =4, units = 'in', res = 500)
hyena.plot + timeplot + mapplot + plot_layout(design = layout)
dev.off()


#-----------FIGURE 2 example of fission-fusion events + alluvial plot ----------


ev <- plot_events(r = 393, events = events.data.exact, xs = xs, ys = ys,
                  phase.col = FALSE, axes = FALSE, xlab = '', ylab = '', cols = colors[c(3,5)])+
  labs(tag= 'A')

canonical.shape <- plot_canonical_shape(393, together.seqs = events.data.exact, xs = xs, ys = ys) +
  labs(tag = 'B')


# --alluvial plot phase transitions --
alluv.data <- data.frame(full.type = events.data.exact[,c('event.type.sym')])
splittypes <- strsplit(as.character(alluv.data$full.type), split = '__')
alluv.data$fusion.word <- alluv.data$together.word <- alluv.data$fission.word <- alluv.data$Fusion <- alluv.data$Together <- alluv.data$Fission <- NA
symbols <- get_event_type_symbols()
for(i in 1:nrow(alluv.data)){
  sym <- symbols[match(alluv.data$full.type[i], names(symbols))]
  symsplit <- strsplit(sym, split = '\n')
  alluv.data$fusion.word[i] <- splittypes[[i]][1]
  alluv.data$together.word[i] <- splittypes[[i]][2]
  alluv.data$fission.word[i] <- splittypes[[i]][3]
  
  alluv.data$Fusion[i] <- symsplit[[1]][3] #note: this is because the symbols are listed backwards to read from bottom to top
  alluv.data$Together[i] <- symsplit[[1]][2]
  alluv.data$Fission[i] <- symsplit[[1]][1]
  
}
alluv.data <- na.omit(alluv.data)

alluv.plot.data <- alluv.data %>% 
  group_by(Fusion, Together, Fission) %>%
  summarize(count = length(full.type))

fusion.symbols <- unique(alluv.plot.data$Fusion)
#alluvp <- alluvial(alluv.plot.data[,c('Fusion','Together','Fission')], freq = alluv.plot.data$count, col = ifelse(alluv.plot.data$Fusion == fusion.symbols[1], colors[7], colors[5]), blocks = T)


ap <- ggplot(alluv.plot.data, aes(y = count, axis1 = Fusion, axis2 = Together, axis3 = Fission))+
  scale_x_continuous(breaks = c(1,2,3), labels = c('Fusion', 'Together', 'Fission'), expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.y = element_blank(), axis.line = element_blank(), rect = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), legend.position = 'none', axis.text.x = element_text(size = 12))+
  geom_alluvium(aes(fill = Fusion), col = 'black', curve_type = 'sine', width = 1/12, reverse = FALSE) + 
  geom_stratum(width = 1/8, alpha = 1, size = 0.5, reverse = FALSE)+
  geom_text(stat='stratum', aes(label = after_stat(stratum)), reverse = FALSE)+
  scale_fill_manual(values = colors[c(6,4)])+
  labs(tag = 'C')

ap.blank <- ggplot(alluv.plot.data, aes(y = count, axis1 = Fusion, axis2 = Together, axis3 = Fission))+
  scale_x_continuous(breaks = c(1,2,3), labels = c('Fusion', 'Together', 'Fission'), expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.y = element_blank(), axis.line = element_blank(), rect = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), legend.position = 'none', axis.text.x = element_text(size = 12))+
 # geom_alluvium(aes(fill = Fusion), col = 'black', curve_type = 'sine', width = 1/12, reverse = FALSE) + 
  geom_stratum(width = 1/8, alpha = 1, size = 0.5, reverse = FALSE)+
  geom_text(stat='stratum', aes(label = after_stat(stratum)), reverse = FALSE)+
  scale_fill_manual(values = colors[c(6,4)])+
  labs(tag = 'C')



png(filename = paste0(plotdir, '/FIG2.png'), width = 6.5, height = 7.5, units = 'in', res = 300)

layout <- "
AABB
CCCC
CCCC
"
ev + canonical.shape + ap + plot_layout(design = layout)
dev.off()
                    
#----------------------------------FIGURE 3-------------------------------------
#--FIGURE 3a--
plotdat.denblock <- plotdat[which(plotdat$type=='denblock'),]
plotdat.denblock$lab.all <- 'All'
plotdat.denblock$lab.den <- 'Den'
plotdat.denblock$lab.nonden <- 'Non-Den'
p_nevents <- ggplot(data = plotdat.denblock) + 
  geom_violin(mapping = aes(x = lab.all, y = nevents), fill = 'gray30', color = 'gray30') + 
  geom_point(x = 1, y = events.tot.data, shape = '_', size = 7)+
  #geom_line(aes(x = x, y = y), data = data.frame(x = c(0.5,1.5), y = rep(events.tot.data, 2)), col = 'gray', size = 2) +
  geom_violin(mapping = aes(x = lab.den, y = denevents), fill = colors[2], color = colors[2]) + 
  geom_point(x = 2, y = events.tot.den.data, shape = '_', size = 7)+
  #geom_line(aes(x = x, y = y), data = data.frame(x = c(1.5,2.5), y = rep(events.tot.den.data, 2)), col = 'blue', size = 2) +
  geom_violin(mapping = aes(x = lab.nonden, y = nondenevents), fill = colors[1], color = colors[1]) + 
  #geom_line(aes(x = x, y = y), data = data.frame(x = c(2.5,3.5), y = rep(events.tot.nonden.data, 2)), col = 'magenta', size = 2) +
  geom_point(x = 3, y = events.tot.nonden.data, shape = '_', size = 7)+
  labs(title="",x="", y = "Number of fission-fusion events") +
  ylim(0, events.tot.data) + 
  theme_classic(base_size = 12) + 
  theme(legend.position="none")+
  labs(tag = 'A')
p_nevents

#-FIGURE 3b--
png(paste0(plotdir, '/FIG3.png'), width = 6, height = 4, units = 'in', res = 500)
p_nevents + visualize_event_type_distributions(events.data, events.rand.list.denblock, 
                                               rand.params, timestamps, 
                                               remove.events.around.day.breaks = T,
                                               col = colors[6])+labs(tag = 'B')
dev.off()

visualize_symmetrical_event_type_comparison(events.data, events.rand.list.denblock, 
                                            rand.params, timestamps, 
                                            remove.events.around.day.breaks = T,
                                            col = colors[6])+labs(tag = 'C')

#-----------------------------------------FIGURE 4------------------------------
png(paste0(plotdir, '/FIG4.png'), width = 5, height = 7, units = 'in', res = 500)
visualize_compare_event_properties(events.data, events.rand.list.denblock,
                                   params, rand.params, timestamps, cols = colors[1:2])
dev.off()

#---------------------------------FIGURE 5--------------------------------------
#create the plot - all
networkdat.denblock <- networkdat[which(networkdat$type %in% c('denblock','real')),]
psri2 <- ggplot(networkdat.denblock, aes(x = dyad, y = sri)) + 
  geom_violin(fill = 'gray30', color = 'grey30', size = 1) + 
  geom_point(aes(x = dyad, y = sri), data = networkdat.denblock[which(networkdat.denblock$type=='real'),], shape = '|', fill = 'black', size = 7) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(legend.position = 'none') +
  ylab('Edge weight') + 
  xlab('')

png(paste0(plotdir, '/FIG5.png'), width = 6.5, height = 4, units = 'in', res = 500)
psri2
dev.off()






###############################################################################
##### Things not in the paper -- determine which are supplemental figures and 
##### delete the rest. 

# #---PLOT 1 - NOT IN PAPER - Number of events in real data vs nightperm and denblock randomization
# #Subplot 1: All events
# ptot <- ggplot(plotdat, aes(x = type, y = nevents)) + 
#   geom_violin(aes(fill = type)) + 
#   labs(title="",x="", y = "Number of events") +
#   ylim(0, events.tot.data + 50) + 
#   theme_minimal(base_size = 18) + 
#   geom_hline(yintercept = events.tot.data, color = 'black', size = 1) + 
#   theme(legend.position="none") + 
#   scale_x_discrete(breaks=c("nightperm","denblock"), labels=c("Complete shuffle", "Den-blocked shuffle")) + 
#   scale_fill_manual(values=c("gray", "gray"))+
#   theme_classic(base_size = 14)
# 
# #Subplot 2: Den events
# pden <- ggplot(plotdat, aes(x = type, y = denevents)) + 
#   geom_violin(aes(fill = type)) + 
#   labs(title="",x="", y = "Number of events") +
#   ylim(0, events.tot.data + 50) + 
#   theme_minimal(base_size = 18) + 
#   geom_hline(yintercept = events.tot.den.data, color = 'blue', size = 1) + 
#   theme(legend.position="none") + 
#   scale_x_discrete(breaks=c("nightperm","denblock"), labels=c("Complete shuffle", "Den-blocked shuffle")) + 
#   scale_fill_manual(values=c("blue", "blue"))+
#   theme_classic(base_size = 14)
# 
# #Subplot 3: Non den events
# pnonden <- ggplot(plotdat, aes(x = type, y = nondenevents)) + 
#   geom_violin(aes(fill = type)) + 
#   labs(title="",x="", y = "Number of events") +
#   ylim(0, events.tot.data + 50) + 
#   theme_minimal(base_size = 18) + 
#   geom_hline(yintercept = events.tot.nonden.data, color = 'magenta', size = 1) + 
#   theme(legend.position="none") + 
#   scale_x_discrete(breaks=c("nightperm","denblock"), labels=c("Complete shuffle", "Den-blocked shuffle")) + 
#   scale_fill_manual(values=c("magenta", "magenta"))+
#   theme_classic(base_size = 14)
# 
# quartz(width = 16, height = 8)
# grid.arrange(ptot, pden, pnonden, nrow = 1)
# 
# 
# #---PLOT 2 - NOT IN PAPER - Comparing social networks generated from the data to the equivalent ones in permuted data
# #get network data into ggplot form (also use SRI here)
# 
# print('Edge weights plot')
# 
# #create the plot - all
# psri <- ggplot(networkdat, aes(x = dyad, y = sri)) + 
#   geom_violin(aes(fill = type)) + 
#   geom_point(aes(x = dyad, y = sri), data = networkdat[which(networkdat$type=='real'),], shape = 23, fill = 'black', size = 4) +
#   coord_flip() +
#   theme_minimal(base_size = 18) +
#   theme(legend.position="top") +
#   scale_fill_manual(values=c("yellow", "light green")) +
#   ylab('Edge weight') + 
#   xlab('')
# 
# #create the plot - den
# psriden <- ggplot(networkdat, aes(x = dyad, y = sri.den)) + 
#   geom_violin(aes(fill = type)) + 
#   geom_point(aes(x = dyad, y = sri.den), data = networkdat[which(networkdat$type=='real'),], shape = 23, fill = 'black', size = 4) +
#   coord_flip() +
#   theme_minimal(base_size = 18) +
#   theme(legend.position="top") +
#   scale_fill_manual(values=c("yellow", "light green")) +
#   ylab('Edge weight') + 
#   xlab('')
# 
# #create the plot - den
# psrinonden <- ggplot(networkdat, aes(x = dyad, y = sri.nonden)) + 
#   geom_violin(aes(fill = type)) + 
#   geom_point(aes(x = dyad, y = sri.nonden), data = networkdat[which(networkdat$type=='real'),], shape = 23, fill = 'black', size = 4) +
#   coord_flip() +
#   theme_minimal(base_size = 18) +
#   theme(legend.position="top") +
#   scale_fill_manual(values=c("yellow", "light green")) +
#   ylab('Edge weight') + 
#   xlab('')
# 
# quartz(width = 16, height = 8)
# grid.arrange(psri, psriden, psrinonden, nrow = 1)
# 
# 
# #------- NOT USED - REMOVE? -------
# #create the plot - den
# psriden <- ggplot(networkdat, aes(x = dyad, y = sri.den)) +
#   geom_violin(aes(fill = type)) +
#   geom_point(aes(x = dyad, y = sri.den), data = networkdat[which(networkdat$type=='real'),], shape = 23, fill = 'black', size = 4) +
#   coord_flip() +
#   theme_minimal(base_size = 18) +
#   theme(legend.position="top") +
#   scale_fill_manual(values=c("yellow", "light green")) +
#   ylab('Edge weight') +
#   xlab('')
# 
# #create the plot - den
# psrinonden <- ggplot(networkdat, aes(x = dyad, y = sri.nonden)) +
#   geom_violin(aes(fill = type)) +
#   geom_point(aes(x = dyad, y = sri.nonden), data = networkdat[which(networkdat$type=='real'),], shape = 23, fill = 'black', size = 4) +
#   coord_flip() +
#   theme_minimal(base_size = 18) +
#   theme(legend.position="top") +
#   scale_fill_manual(values=c("yellow", "light green")) +
#   ylab('Edge weight') +
#   xlab('')
# 
# quartz(width = 16, height = 8)
# grid.arrange(psri, psriden, psrinonden, nrow = 1)


# 
# #-------PLOT NOT USED overall ranging patterns -------
# step <- 300
# n.inds <- dim(xs)[1]
# n.times <- dim(xs)[2]
# movedat <- data.frame()
# for(i in 1:n.inds){
#   t.sub <- seq(from = 1, to = n.times, by = step)
#   movedat <- rbind(movedat, data.frame(id = rep(hyena.ids$name[i], length(t.sub)), lon = lons[i, t.sub], lat = lats[i, t.sub]))
# }
# quartz(height = 8, width = 16)
# mapplot2 <- ggmap(hyena_map_transparent, alpha = 0.5) + 
#   geom_point(aes(x = lon, y = lat, color= id), data = movedat, size = 1, shape = 19, alpha = 0.3) +
#   scale_color_brewer(palette="Set1") + 
#   geom_point(aes(x = lon, y = lat), data = den.locs, size = 10, shape = 21, color = 'blue', stroke = 1.2)  +
#   geom_line(aes(x = lon, y = lat), data = scalebar, color = 'black', lwd = 1.2) +
#   geom_text(aes(x = lon, y = lat), data = scalebar2, label = '1 km', nudge_y = .003) +
#   theme(legend.position = 'bottom')
# mapplot2
# 
# grid.arrange(mapplot2, mapplot, nrow = 1)

