#Run analysis of "attraction range"
#Given that a hyena starts a distance of R_start meters from another hyena, what is the probability that it approaches
#that hyena within R_close = 5 m before exiting a radius of R_exit = R_start + buffer m?

#DIRECTORIES
datadir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps/'
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'

#PARAMS
R.close <- 5
buffer <- 100
R.start.max <- 100

#LOAD DATA
setwd(datadir)
load('hyena_xy_level1.RData')

#FUNCS
setwd(codedir)
source('ff_functions_library.R')

#Set params
params <- list(R.fusion = R.start.max, 
               R.fission = R.start.max + buffer, 
               max.break = 60*30, #max time between events to merge events connected by NAs
               move.thresh = 5, #minimum amount moved to be considered 'moving' during a phase
               together.travel.thresh = 200, #minimum amount to be considered 'moving' during the together phase
               local.time.diff = 3, # difference in hours from local time
               den.dist.thresh = 200
) 

#get events
events <- get_ff_events_and_phases(xs, ys, params)

#Get good indexes with defined start and end time
good.idxs <- which(events$start.exact & events$end.exact)

#Can also get this for R.exist = 200 and R.start = [0, 100]
Rs <- seq(R.close, R.start.max, 1)
attraction.probs <- rep(NA, length(Rs))
for(i in 1:length(Rs)){
  R.start <- Rs[i]
  attraction.probs[i] <- sum(events$closest.app[good.idxs] < R.close) / sum(events$closest.app[good.idxs] < R.start)
}

plot(Rs, attraction.probs, type = 'l', ylim = c(0,1), lwd = 2, xlab = 'Initial distance (m)', ylab = 'Probability of coming to close range', xlim =c(5,R.start.max))
abline(v=50)
plot(Rs, attraction.probs, type = 'l', ylim = c(0,1), lwd = 2, xlab = 'Initial distance (m)', ylab = 'Probability of coming to close range', xlim =c(5,100))
abline(v=50)
