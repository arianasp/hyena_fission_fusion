#This script computes networks based on fission-fusion events and looks at differentiation of social structure 
#depending on what types of events are used

#directory of where the data of fission-fusion events is stored
ff_dir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/data/'
raw_dir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps'

#file to load
realdata_file <- 'fission_fusion_events_features.RData'
permdata_file <- 'hyena_day_randomization_plan_denblock_avoidmatchTRUE.RData'
gpsdata_file <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/'

#LOAD DATA
setwd(ff_dir)
load(realdata_file)

#setwd(raw_dir)
#load('hyena_xy_level1.RData')

#PREPROCESS

#get number of individuals
n.inds <- 5

#get number of seconds when both individuals were tracked
#both.tracked <- array(dim = c(n.inds, n.inds))
#for(i in 1:(n.inds-1)){
#  for(j in (i+1):n.inds){
#    both.tracked[i,j] <- both.tracked[j,i] <- sum(!is.na(xs[i,]) & !is.na(xs[j,]))
#  }
#}

#initialize network
net <- array(NA, dim = c(n.inds, n.inds))

#edge definition (later pass as argument to a function, maybe)
edge.def <- 'disp'

#construct network based on ff events (different edge definitions possible)
for(a in 1:(n.inds-1)){
  for(b in (a+1):n.inds){
    
    #get events involving individual a and b
    events.ab <- events[which((events$i == a & events$j == b) | (events$i == b & events$j == a)),]
    
    #events involving individual 1 (as either the 'i' or 'j' individual in the event...)
    events.a <- events[which(events$i == a | events$j == a),]
    
    #events involving individual 2 (as either the 'i' or 'j' individual in the event...)
    events.b <- events[which(events$i == b | events$j == b),]
    
    #compute edge weight (depending on definition using) using SRI index AB_together / (A_anyone + B_anyone - AB_together)
    if(edge.def == 'dur'){
      dur.ab <- sum(events.ab$tf - events.ab$t0, na.rm=T)
      dur.a <- sum(events.a$tf - events.a$t0, na.rm=T)
      dur.b <- sum(events.b$tf - events.b$t0, na.rm=T)
      net[a,b] <- dur.ab / (dur.a + dur.b - dur.ab)
    } else if (edge.def == 'disp'){
      disp.ab <- sum(events.ab$disp.together.i + events.ab$disp.together.j, na.rm=T)/2 #take avg displacement of i and j in an event as the "disp together"
      disp.a <- sum(events.a$disp.together.i + events.a$disp.together.j, na.rm=T)/2
      disp.b <- sum(events.b$disp.together.i + events.b$disp.together.j, na.rm=T)/2
      net[a,b] <- disp.ab / (disp.a + disp.b - disp.ab)
    }
  }
}

#get the variance, as a metric of relationship differentiation
var.net <- var(c(net), na.rm=T)



