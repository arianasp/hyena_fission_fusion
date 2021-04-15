#This script explores how many of the fission-fusion events we observe are actually polyadic

#----DIRECTORIES---
indir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/data/main_output/'
plotdir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/results'
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'

#-----LIBRARIES----
library(ggplot2)
library(gridExtra)

#----FUNCTIONS----
setwd(codedir)
source('ff_functions_library.R')

#----LOAD DATA-----
setwd(indir)
load('fission_fusion_events_features.RData')

#----AGGREGATE EVENTS THAT OVERLAP-----

#get number of individuals
n.inds <- max(c(events$i, events$j))

#sort events by start time, because why not
events <- events[order(events$t0),]

#create a column for whether the event has already been "dealt with"
events$accounted.for <- FALSE #flag to indicate if the event has already been placed in an aggregated event

agg.events <- data.frame()
inds.involved.list <- list()
event.ids.list <- list()
event.types.list <- list()
agg.event.idx <- 1
for(i in 1:nrow(events)){
  
  #if the event has not already been aggregated into another event...
  if(!events$accounted.for[i]){
    
    #get information about the current event
    t0 <- events$t0[i]
    tf <- events$tf[i]
    inds.involved <- c(events$i[i], events$j[i])
    events.involved <- c(i)
    event.types.involved <- c(events$event.type.sym[i])
    
    #count it as "accounted for"
    events$accounted.for[i] <- TRUE
  
    #find other overlapping events
    converged <- F
    while(!converged){
      
      #get the indexes of events that overlap (in time) the current event and involve at least one of the inds.involved
      max.t0 <- sapply(events$t0, FUN = function(x){return(max(x,t0))})
      min.tf <- sapply(events$tf, FUN = function(x){return(min(x,tf))})
      overlap.idxs <- which((max.t0 <= min.tf) & ((events$i %in% inds.involved) | (events$j %in% inds.involved)) & (!events$accounted.for))
      
      #if aggregated event has converged (no more events to add), add it to the data frame + lists, and set converged to TRUE
      if(length(overlap.idxs)==0){
        event.ids.list[[agg.event.idx]] <- events.involved
        inds.involved.list[[agg.event.idx]] <- inds.involved
        event.types.list[[agg.event.idx]] <- event.types.involved
        agg.row <- data.frame(t0 = min(events$t0[events.involved]), tf = max(events$tf[events.involved]))
        agg.events <- rbind(agg.events, agg.row)
        converged <- T
        agg.event.idx <- agg.event.idx + 1
      } else{
        #otherwise, add the events that overlap to the list and the indivduals that overlap to the list
        events.involved <- unique(c(events.involved, overlap.idxs))
        inds.involved <- unique(c(events$i[events.involved], events$j[events.involved]))
        event.types.involved <- c(event.types.involved, events$event.type.sym[overlap.idxs])
        events$accounted.for[events.involved] <- TRUE
        t0 <- min(events$t0[events.involved])
        tf <- max(events$tf[events.involved])
      }
    }
  }
}

#add the individuals involved and the events involved to the agg.events data frame
agg.events$inds.involved <- inds.involved.list
agg.events$event.ids.involved <- event.ids.list
agg.events$event.types.involved <- event.types.list

#get number of individuals involved
agg.events$n.inds <- sapply(agg.events$inds.involved, length)
agg.events$n.events <- sapply(agg.events$event.ids.involved, length)

#get the minimum distance from the den across all events making up the aggregated event (start and end)
agg.events$min.dist.den <- NA
for(i in 1:nrow(agg.events)){
  events.involved.curr <- agg.events$event.ids.involved[i][[1]]
  dists.from.den <- c(events$dist.den.start[events.involved.curr], events$dist.den.end[events.involved.curr])
  agg.events$min.dist.den[i] <- min(dists.from.den, na.rm=T)
}

#boolean column, whether at the den during the event
agg.events$at.den <- agg.events$min.dist.den <= params$den.dist.thresh

#------PLOTTING----

#---Plot 1: get distribution of number of individuals involved in aggregated events---
n.inds.involved <- agg.events$n.inds
n.inds.involved.tab <- table(n.inds.involved)

setwd(plotdir)
png(file = 'polyadic_events.png', width = 6, height = 6, unit = 'in', res = 300)
par(mar=c(5,5,1,1))
barplot(n.inds.involved.tab, names.arg = names(n.inds.involved.tab), xlab = '# of individuals involved', ylab = 'Event frequency', cex.lab = 2, cex.axis = 1.5)
dev.off()

#---Plot 2: number of dyadic events accounted for in events of each # of inds---
dyad.events.per.agg.event.type <- rep(NA, n.inds-1)
for(i in 1:(n.inds-1)){
  idxs <- which(agg.events$n.inds==(i+1))
  dyad.events.per.agg.event.type[i] <- sum(agg.events$n.events[idxs])
}

png(file = 'dyadic_events_making_up_polyadic_events.png', width = 6, height = 6, unit = 'in', res = 300)
par(mar=c(5,5,1,1))
barplot(dyad.events.per.agg.event.type, names.arg = 2:n.inds, xlab = '# of individuals involved', ylab = '# of dyadic events accounted for', cex.lab = 2, cex.axis = 1.5)
dev.off()

# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Stacked
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity")

#---Plot 3: get distribution of number of individuals involved in aggregated events - den vs. non-den---
plot.data <- aggregate(agg.events$n.events, by = list(agg.events$n.inds, agg.events$at.den), FUN = length)
colnames(plot.data) <- c('n.inds','at.den','n.agg.events')
plot.data <- plot.data[order(plot.data$at.den, decreasing = T),]
plot.data$percent <- NA
for(i in 2:n.inds){
  idxs <- which(plot.data$n.inds == i)
  plot.data$percent[idxs] <- plot.data$n.agg.events[idxs] / sum(plot.data$n.agg.events[idxs]) * 100
}

poly_den_events <- ggplot(plot.data, aes(fill = !at.den, y = n.agg.events, x = n.inds)) + 
  geom_bar(position = 'stack', stat = 'identity') + 
  theme_classic(base_size = 24) + 
  ylab('Number of aggregated events') + 
  xlab('Number of hyenas involved') +
  scale_fill_manual(values = c('blue','magenta'), labels = c('Den', 'Non-Den'))+
  theme(legend.position = c(.75,.85), legend.title = element_blank()) + 
  geom_label(aes(label = paste0(round(percent, digits = 1),"%")), position = position_stack(vjust = 0.5), colour = 'white', show.legend = F)

poly_den_events

#Plot 4: number of dyadic events accounted for in events of each # of inds and den status
plot.data <- aggregate(agg.events$n.events, by = list(agg.events$n.inds, agg.events$at.den), FUN = sum)

colnames(plot.data) <- c('n.inds','at.den','n.agg.events')
plot.data <- plot.data[order(plot.data$at.den, decreasing = T),]
plot.data$percent <- NA
for(i in 2:n.inds){
  idxs <- which(plot.data$n.inds == i)
  plot.data$percent[idxs] <- plot.data$n.agg.events[idxs] / sum(plot.data$n.agg.events[idxs]) * 100
}

poly_den_events_bydyad <- ggplot(plot.data, aes(fill = !at.den, y = n.agg.events, x = n.inds)) + 
  geom_bar(position = 'stack', stat = 'identity') + 
  theme_classic(base_size = 24) + 
  ylab('Number of dyadic events') + 
  xlab('Number of hyenas involved') +
  scale_fill_manual(values = c('blue','magenta'), labels = c('Den', 'Non-Den'))+
  theme(legend.position = c(.75,.85), legend.title = element_blank()) + 
  geom_label(aes(label = paste0(round(percent, digits = 1),"%")), position = position_stack(vjust = 0.5), colour = 'white', show.legend = F)

quartz(height = 10, width = 12)
grid.arrange(poly_den_events, poly_den_events_bydyad, nrow = 1)

#Plot 5: Event types making up each category

#calculate for each event, what type of (polyadic or dyadic) aggregated event it is in
event.types.all <- c('fusion.stay.move__together.local__fission.stay.move',
                     'fusion.stay.move__together.local__fission.move.stay',
                     'fusion.stay.move__together.local__fission.move.move',
                     'fusion.move.move__together.local__fission.stay.move',
                     'fusion.move.move__together.local__fission.move.move',
                     'fusion.stay.move__together.travel__fission.stay.move',
                     'fusion.stay.move__together.travel__fission.move.stay',
                     'fusion.stay.move__together.travel__fission.move.move',
                     'fusion.move.move__together.travel__fission.stay.move',
                     'fusion.move.move__together.travel__fission.move.move')

event.type.symbols <- get_event_type_symbols()

for(i in 1:nrow(events)){
  for(j in 1:nrow(agg.events)){
    if(i %in% agg.events$event.ids.involved[j][[1]]){
      events$poly[i] <- agg.events$n.inds[j]
    }
  }
}

events.non.na <- events[which(!grepl('NA', events$event.type.sym, ignore.case=F)),]
plot.data.type <- aggregate(events.non.na$poly, by = list(events.non.na$poly, events.non.na$event.type.sym), FUN = length)
colnames(plot.data.type) <- c('poly','type','count')

#trying to get ggplot to plot the bars in the right order!
plot.data.type <- plot.data.type[order(plot.data.type$poly, match(plot.data.type$type, event.types.all)),]
plot.data.type$type <- factor(plot.data.type$type, levels = event.types.all)

#calculate percentages
plot.data.type$percent <- NA
for(i in 2:n.inds){
  idxs <- which(plot.data.type$poly==i)
  plot.data.type$percent[idxs] <- plot.data.type$count[idxs] / sum(plot.data.type$count[idxs]) * 100
}

cols.local <- colorRampPalette(c('darkred','orange'))
cols.travel <- colorRampPalette(c('blue','#00CC00'))
poly_events_type <- ggplot(plot.data.type, aes(fill = type, y = percent, x = poly)) + 
  geom_bar(position = 'stack', stat = 'identity') + 
  theme_classic(base_size = 24) + 
  ylab('Percentage of events') + 
  xlab('Number of hyenas involved') +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c(cols.local(5),cols.travel(5)))

quartz(width = 6, height = 10)
poly_events_type
