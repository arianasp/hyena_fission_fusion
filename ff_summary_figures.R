#Summary plots that integrate information across multiple permutations

#-----------LIBRARIES-------------
library(ggplot2)
library(gridExtra)
library(ggmap)
library(RColorBrewer)
library(lubridate)
library(viridis)
library(ggalluvial)
library(dplyr)
library(patchwork)
library(magick)
library(ggnetwork)
library(sna)


#-----------FILENAMES---------------
data_output_file <- 'fission_fusion_events_features.RData'
denblock_output_file <- 'hyena_day_randomization_events_features_denblock_avoidmatchTRUE.RData'
map_file <- 'hyena_satellite_map.RData'

#----------LOAD DATA-------------
load(paste0(processed.data.directory,'hyena_timestamps.RData'))
load(paste0(processed.data.directory,'hyena_ids.RData'))
load(paste0(processed.data.directory,'hyena_xy_level1.RData'))
load(paste0(processed.data.directory,'hyena_latlon_level0.RData'))

#load events from real data
load(paste0(data.outdir, data_output_file))
events.data <- events
rm('events')

#load events from denblock randomization
load(paste0(data.outdir, denblock_output_file))
events.rand.list.denblock <- events.rand.list
rm('events.rand.list')

load(paste0(raw.data.directory, map_file))

# Load hyena image
hyena <- image_read(path = paste0(raw.data.directory, "collared_hyena.jpg"))

#------------PROCESS-----------------
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
events.tot.denblock <- rep(NA, n.rands)
events.tot.den.denblock <- events.tot.nonden.denblock <- rep(NA, n.rands)
events.net.denblock <- array(NA, dim = c(n.inds, n.inds, n.rands))
events.net.denblock.den <- array(NA, dim = c(n.inds, n.inds, n.rands))
events.net.denblock.nonden <- array(NA, dim = c(n.inds, n.inds, n.rands))
  
for(r in 1:n.rands){
  
  #get events associated with that randomization\
  events.rand.denblock <- events.rand.list.denblock[[r]]
  
  #remove events surrounding the 'day break'
  events.rand.denblock <- remove_events_around_day_breaks(events.rand.denblock, timestamps, rand.params)
  
  #total number of events
  events.tot.denblock[r] <- nrow(events.rand.denblock)
  
  #den vs nonden 
  events.tot.den.denblock[r] <- sum(events.rand.denblock$dist.den.start <= params$den.dist.thresh | events.rand.denblock$dist.den.end <= params$den.dist.thresh)
  events.tot.nonden.denblock[r] <- sum(events.rand.denblock$dist.den.start > params$den.dist.thresh & events.rand.denblock$dist.den.end > params$den.dist.thresh)
  
  #den indexes
  denblock.den.idxs <- which(events.rand.denblock$dist.den.start <= params$den.dist.thresh | events.rand.denblock$dist.den.end <= params$den.dist.thresh)
  denblock.nonden.idxs <- which(events.rand.denblock$dist.den.start > params$den.dist.thresh & events.rand.denblock$dist.den.end > params$den.dist.thresh)
  
  #networks - all
  events.net.denblock[,,r] <- count_events_per_dyad(events.rand.denblock)
  
  #networks - den vs nonden
  events.net.denblock.den[,,r] <- count_events_per_dyad(events.rand.denblock[denblock.den.idxs,])
  events.net.denblock.nonden[,,r] <- count_events_per_dyad(events.rand.denblock[denblock.nonden.idxs,])
  
}

# Pre-processing for PLOT 1 and FIGURE 3a ------
#prep for ggplot
plotdat <- data.frame(type = rep('denblock',n.rands), nevents = events.tot.denblock, denevents = events.tot.den.denblock, nondenevents = events.tot.nonden.denblock)
# plotdat$type <- factor(plotdat$type,
#                        levels = c('denblock', 'real'),ordered = TRUE)

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
  tmp$type <- rep('denblock',nrow(dyads))
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
                          levels = c('denblock','real'),ordered = TRUE)

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
den.locs <- get_dens(paste0(raw.data.directory, 'metadata/hyena_isolate_dens.csv'))

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
print(paste0('The number of together-traveling events starting at the den is ', sum(events.data$dist.den.start <= params$den.dist.thresh & events.data$together.type=='together.travel', na.rm = T)))
print(paste0('The number of together-traveling events ending at the den is ', sum(events.data$dist.den.end <= params$den.dist.thresh & events.data$together.type=='together.travel', na.rm = T)))

#--------------PLOTS----------------

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
hist(all.hrs, main = 'GPS data available by hour', xlab = 'Hour')
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

png(paste0(plots.outdir, 'FIG1.png'), width = 6, height =4, units = 'in', res = 500)
hyena.plot + timeplot + mapplot + plot_layout(design = layout)
dev.off()


#-----------FIGURE 2 example of fission-fusion events + alluvial plot ----------


ev <- plot_events(r = 394, events = events.data.exact, xs = xs, ys = ys,
                  phase.col = FALSE, axes = FALSE, xlab = '', ylab = '', cols = colors[c(3,5)]) +
  labs(tag= 'A')

canonical.shape <- plot_canonical_shape(394, together.seqs = events.data.exact, xs = xs, ys = ys) +
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
  summarize(count = length(full.type)) %>%
  ungroup()

fusion.symbols <- unique(alluv.plot.data$Fusion)
#alluvp <- alluvial(alluv.plot.data[,c('Fusion','Together','Fission')], freq = alluv.plot.data$count, col = ifelse(alluv.plot.data$Fusion == fusion.symbols[1], colors[7], colors[5]), blocks = T)

alluv.plot.data$Event_Type = letters[1:nrow(alluv.plot.data)]
ap <- ggplot(alluv.plot.data, aes(y = count, axis1 = Fusion, axis2 = Together, axis3 = Fission))+
  scale_x_continuous(breaks = c(1,2,3), labels = c('Fusion', 'Together', 'Fission'), expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.y = element_blank(), axis.line = element_blank(), rect = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), legend.position = 'none', axis.text.x = element_text(size = 12))+
  geom_alluvium(aes(fill = Fusion), col = 'black', curve_type = 'sine', width = 1/12, reverse = FALSE) + 
  #geom_alluvium(aes(fill = Event_Type), col = 'black', curve_type = 'sine', width = 1/12, reverse = FALSE) + 
  geom_stratum(width = 1/8, alpha = 1, size = 0.5, reverse = FALSE)+
  geom_text(stat='stratum', aes(label = after_stat(stratum)), reverse = FALSE)+
  scale_fill_manual(values = colors[c(6,4)])+
  #scale_fill_manual(values = plasma(10))+
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



png(filename = paste0(plots.outdir, 'FIG2.png'), width = 6.5, height = 7.5, units = 'in', res = 300)

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

#-FIGURE 3b--
png(paste0(plots.outdir, 'FIG3.png'), width = 6, height = 4, units = 'in', res = 500)

layout <- '
AAABBBB
'
p_nevents + visualize_event_type_distributions(events.data, events.rand.list.denblock, 
                                               rand.params, timestamps, 
                                               remove.events.around.day.breaks = T,
                                               col = colors[6])+labs(tag = 'B') +
  plot_layout(design = layout)
dev.off()

# visualize_symmetrical_event_type_comparison(events.data, events.rand.list.denblock, 
#                                             rand.params, timestamps, 
#                                             remove.events.around.day.breaks = T,
#                                             col = colors[6])+labs(tag = 'C')

#-----------------------------------------FIGURE 4------------------------------
png(paste0(plots.outdir, 'FIG4.png'), width = 5, height = 7, units = 'in', res = 500)
visualize_compare_event_properties(events.data, events.rand.list.denblock,
                                   params, rand.params, timestamps, cols = colors[1:2])
dev.off()

#---------------------------------FIGURE 5--------------------------------------
netplotdat <- networkdat

netplotdat$i <- LETTERS[netplotdat$i]
netplotdat$j <- LETTERS[netplotdat$j]
netplotdat$dyad <- paste(netplotdat$i, netplotdat$j, sep = '--')


networkdat.denblock <- netplotdat[which(netplotdat$type %in% c('denblock','real')),]
psri <- ggplot(networkdat.denblock, aes(x = dyad, y = sri, fill = dyad, color = dyad)) + 
  geom_violin(size = 1) + 
  geom_point(aes(x = dyad, y = sri), data = networkdat.denblock[which(networkdat.denblock$type=='real'),], shape = '|', color = 'black', size = 7) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(legend.position = 'none') +
  ylab('Edge weight') + 
  xlab('Edge')+
  scale_color_manual(values = plasma(10))+
  scale_fill_manual(values = plasma(10))



obs.edges <- netplotdat[netplotdat$type == 'real',][1:10,]
obs.net <- network(obs.edges[c('i', 'j')], directed = F, loops = F, multiple = F)
obs.net %e% 'sri' <- obs.edges$sri^2 ## Because line thickness is scaled to sqrt of the value (see scale_size_area documentation)
obs.net %e% 'name' <- obs.edges$dyad
obs.net <- ggnetwork(obs.net, weights = 'sri', layout = 'circle')

obs.plot <- ggplot(obs.net, aes(x = x, y = y, xend = xend, yend = yend, col = name))+
  geom_edges(aes(size = sri), curvature = 0.2)+
  geom_nodes(size = 2, color = 'black')+
  geom_nodetext(aes(label = vertex.names), color = 'white', size = 2)+
  theme_blank()+
  ggtitle(label = 'Observed')+
  theme(legend.position = 'none',
        plot.title = element_text(size = 8, hjust = 0.5))+
  scale_color_manual(values = plasma(10))+
  scale_size_area(max_size = 2)

rand.counter <- 1
for(i in sample(1:n.rands, 3, replace = F)){
  edges <-  netplotdat[!is.na(netplotdat$rand) & netplotdat$rand == i,][1:10,]
  rand.net <- obs.net
  edges$dyad[match(rand.net$name[1:10], edges$dyad)] == rand.net$name[1:10]
  rand.net$sri[1:10] <- edges$sri[match(rand.net$name[1:10], edges$dyad)]^2 ## Because line thickness is scaled to sqrt of the value (see scale_size_area documentation)
  p <- ggplot(rand.net, aes(x = x, y = y, xend = xend, yend = yend, col = name))+
    geom_edges(aes(size = sri), curvature = 0.2)+
    geom_nodes(size = 2, color = 'black')+
    geom_nodetext(aes(label = vertex.names), color = 'white', size = 2)+
    theme_blank()+
    ggtitle(label = ifelse(rand.counter == 2, 'Reference model', ''))+
    theme(legend.position = 'none', plot.title = element_text(size = 8, hjust = 0.5))+
    scale_color_manual(values = plasma(10))+
    scale_size_area(max_size = 2)
  
  assign(paste0('rand.plot.', rand.counter), p)
  rand.counter <- rand.counter + 1
}

layout <- '
BCDE
AAAA
AAAA
AAAA
'


png(paste0(plots.outdir, 'FIG5.png'), width = 6.5, height = 5, units = 'in', res = 500)
psri + obs.plot + rand.plot.1 + rand.plot.2 + rand.plot.3 + plot_layout(design = layout)
grid.rect(x = 0.655, y = 0.835, width = 0.64, height = 0.2, gp = gpar(lwd = 1, col = 'black', fill = NA, lty = 2))
grid.rect(x = 0.2055, y = 0.835, width = 0.19, height = 0.2, gp = gpar(lwd = 1, col = 'black', fill = NA))
dev.off()
