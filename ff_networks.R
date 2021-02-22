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

#rewire the links for a given data frame (data stream permutation)
#this permutation keeps the total number of events each individual is involved with constant, but randomizes who interacts with whom
#because of the way events data frame was constructed, individual i's index is always less than individual j's index.
#to remove this weird artifact, we first randomize which individual (greater or lesser index is listed first)
#we then loop through all events (in random order), find a "partner event" for them to swap with, and swap the j individual between the two partner events
rewire_links <- function(events, n.reps = 100){
  
  perm.events <- events
  
  n.events <- nrow(events)
  
  #first randomize which individual is considered i vs j
  #(due to the way events were extracted from the data, previously index of i was always less than index of j)
  for(r in 1:n.events){
    curr.i <- perm.events$i[r]
    curr.j <- perm.events$j[r]
    ids <- c(curr.i, curr.j)
    
    idxs <- sample(2,2)
    perm.events$i[r] <- ids[idxs[1]]
    perm.events$j[r] <- ids[idxs[2]]
  }
  
  for(rep in 1:n.reps){
    #then go through all the events (in a random order). for each one, pick a random partner event that would not create a self-loop, and swap the j individuals of those events
    event.order <- sample(n.events)
    for(r in 1:n.events){
      
      #get the index of the original event
      idx.orig <- event.order[r]
      
      #get the participants
      i.orig <- perm.events$i[idx.orig]
      j.orig <- perm.events$j[idx.orig]
      
      #find all possible partner events (avoid self-loops)
      possible.partner.events <- which(!(perm.events$i == j.orig | perm.events$j == i.orig))
      
      #select one partner event at random
      #note: this if statement is necessary because if the list only had one element, sample call would be interpreted incorrectly (though in reality this should never happen anyway)
      if(length(possible.partner.events) > 1){
        idx.partner <- sample(possible.partner.events, 1)
      } else{
        idx.partner <- possible.partner.events[1]
      }
      
      #rewire! note: here we just rewire j since we already randomized who i and j are
      j.partner <- perm.events$j[idx.partner]
      
      perm.events$j[idx.orig] <- j.partner
      perm.events$j[idx.partner] <- j.orig
      
    }
  }
  
  return(perm.events)
}

#from a data frame of events, produce a network based on the SRI index of frequency of events (symmetric or asymmetric)
#SRI (symmetric) = default: net[A,B] =  N_AB / (N_A + N_B - N_AB)
#SRI (asymmetric) = net[A,B] = N_AB / N_A
get_SRI_net <- function(events, n.inds = 5, asymmetric = F, counts = F){
  
  #initialize net to store adjacency matrix
  net <- matrix(NA, nrow = n.inds, ncol = n.inds)
  
  #construct network based on SRI
  for(a in 1:(n.inds-1)){
    for(b in (a+1):n.inds){
      
      #number of events involving individual a and b
      n.ab <- sum((events$i == a & events$j == b) | (events$i == b & events$j == a))
      
      #events involving individual 1 (as either the 'i' or 'j' individual in the event...)
      n.a <- sum(events$i == a | events$j == a)
      
      #events involving individual 2 (as either the 'i' or 'j' individual in the event...)
      n.b <- sum(events$i == b | events$j == b)
      
      #if doing counts instead of edge weight, run that here
      if(counts == T){
        net[a,b] <- net[b,a] <- n.ab
      } else{
        #compute edge weight (if symmetric, the default)
        if(!asymmetric){
          net[a,b] <- net[b,a] <- n.ab / (n.a + n.b - n.ab)
          
          #compute edge weight (if asymmetric)
        } else{
          net[a,b] <- n.ab / n.a
          net[a,b] <- n.ab / n.b
        }
      }
    }
  }
  
  return(net)
  
}

#get network and compare to partner swapping permutations
get_SRI_nets_real_vs_permuted <- function(events, n.perms = 100, n.inds = 5, asymmetric = F, counts = F){
  
  #get real network
  net.real <- get_SRI_net(events = events, n.inds = n.inds, asymmetric = asymmetric, counts = counts)
  
  #initialize array to hold permuted networks
  nets.perm <- array(NA, dim = c(dim(net.real), n.perms))
  
  #get permuted networks
  for(r in 1:n.perms){
    perm.events <- rewire_links(events)
    nets.perm[,,r] <- get_SRI_net(perm.events, n.inds = n.inds, asymmetric = asymmetric, counts = counts)
    
  }
  
  #store
  out <- list()
  out$real <- net.real
  out$perm <- nets.perm
  
  return(out)
  
}

get_all_nets <- function(events, params, counts = F){
  
  #get number of individuals
  n.inds <- length(unique(c(events$i, events$j)))
  
  #subset to only good indexes (start and end time exact)
  good.idxs <- which(events$t0.exact & events$tf.exact)
  events.good <- events[good.idxs,]
  
  #initialize list to store output
  nets <- list()
  
  #network of all events
  nets$all <- get_SRI_nets_real_vs_permuted(events = events.good, counts = counts)
  
  #network of all 'local' events
  events.local <- events.good[which(events.good$together.type == 'together.local'),]
  nets$local <- get_SRI_nets_real_vs_permuted(events = events.local, counts = counts)
  
  #network of all 'travel' events
  events.travel <- events.good[which(events.good$together.type == 'together.travel'),]
  nets$travel <- get_SRI_nets_real_vs_permuted(events = events.travel, counts = counts)
  
  #network of all den events
  events.den <- events.good[which(events.good$dist.den.start <= params$den.dist.thresh | events.good$dist.den.end <= params$den.dist.thresh),]
  nets$den <- get_SRI_nets_real_vs_permuted(events = events.den, counts = counts)
  
  #network of all nonden events
  events.nonden <- events.good[which(events.good$dist.den.start > params$den.dist.thresh & events.good$dist.den.end > params$den.dist.thresh),]
  nets$nonden <- get_SRI_nets_real_vs_permuted(events = events.nonden, counts = counts)
  
  return(nets)
  
}

#compare variances
compare_variances <- function(nets, plot_vars = T){
  
  #number of network types
  n.net.types <- length(nets)
  
  #number of permutations done
  n.perms <- dim(nets[1][[1]]$perm)[3]
  
  #vector and matrix to hold variances for real data and permutations, respectively
  vars.real <- rep(NA, n.net.types)
  vars.perm <- matrix(NA, nrow = n.perms, ncol = n.net.types)
  
  for(i in 1:n.net.types){
    
    #get networks - real and permuted
    reali <- nets[i][[1]]$real
    permsi <- nets[i][[1]]$perm
    
    #compute variances
    vars.real[i] <- var(c(reali), na.rm=T)
    vars.perm[,i] <- apply(X = permsi, MARGIN = 3, FUN = function(x){return(var(c(x),na.rm=T))})
    
  }
  
  vars <- list()
  vars$real <- vars.real
  vars$perm <- vars.perm
  
  if(plot_vars){
    
    quartz()
    netnames <- names(nets)
    minvar <- min(min(vars$real), min(vars$perm))
    maxvar <- max(max(vars$real), max(vars$perm))
    plot(NULL, xlim = c(0.5, n.net.types + 0.5), ylim = c(minvar, maxvar), xaxt='n', xlab = '', ylab = 'Differentiation (edge weight variane)')
    for(i in 1:length(netnames)){
      jitter <- rnorm(0, sd = .1, n = n.perms)
      points(rep(i, n.perms) + jitter, vars$perm[,i], cex = 0.5, col = '#00000055')
      points(i, vars$real[i], cex = 2, pch = 8, col = 'red')
    }
    axis(side = 1, at = 1:n.net.types, label = netnames)
    
  }
  
  return(vars)
  
}

allcounts <- get_all_nets(events, params, counts = T)
vars <- compare_variances(allnets)


# 
# #plot the networks
# quartz(height = 12, width = 8)
# par(mfrow = c(3,2))
# net.types <- names(nets)
# for(i in 1:length(nets)){
#   image(nets[i][[1]], main = paste0(net.types[i], ' (var = ', round(vars[i], digits=4), ')'), xlab = '', ylab = '', xaxt='n', yaxt='n', cex.main = 2, zlim=c(0,0.3))
#   for(a in 1:n.inds){
#     for(b in 1:n.inds){
#       text(a, b, counts[i][[1]][a,b])
#     }
#   }
# }
#   
  
  
  
  
