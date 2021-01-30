#Driver script to run fission-fusion analyses

################################ SET UP DIRECTORIES ##################################

#directory of where the original (processed) movement + vedba data is stored (+ metadata on IDs and den locations)
indir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps'

#directory of where to store extracted data for fission-fusion project
outdir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/data/'

#directory where code for fission-fusion project is stored
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'

################################ CHOOSE ANALYSES TO RUN ##################################

run_extract_ff_events <- T
overwrite_extract_ff_events <- T
run_get_ff_features <- T
overwrite_extract_ff_features <- T
generate_day_randomization_plan <- T
overwrite_day_randomization_plan <- T
execute_day_randomization_plan <- T
overwrite_day_randomization_output <- T

################################ PARAMETERS ##########################################

params <- list(R.fusion = 100, 
               R.fission = 200, 
               max.break = 60*30, #max time between events to merge events connected by NAs
               move.thresh = 5, #minimum amount moved to be considered 'moving' during a phase
               together.travel.thresh = 200, #minimum amount to be considered 'moving' during the together phase
               last.day.used = 35) #last day to use in the randomizations (and real data)

verbose <- TRUE
n.rands <- 3
local.time.diff <- 3 #difference in hours from local time

events_filename <- 'fission_fusion_events.RData'
events_features_filename <- 'fission_fusion_events_features.RData'
day_start_idxs_filename <- 'hyena_day_start_idxs.RData'
day_randomization_plan_filename <- 'hyena_day_randomization_plan.RData'
day_randomization_output_filename <- 'hyena_day_randomization_events_features.RData'

################################# SOURCE FUNCTIONS #######################################

setwd(codedir)
source('ff_functions_library.R')

################################# MAIN #######################################

if(run_extract_ff_events){
  
  #Load data
  setwd(indir)
  if(verbose){
    print('Loading xy data')
  }
  load('hyena_xy_level1.RData')
  
  #extract events
  events <- get_ff_events_and_phases(xs = xs, 
                                     ys = ys, 
                                     params = params)
  
  if(overwrite_extract_ff_events){
    setwd(outdir)
    if(verbose){
      print(paste0('Saving events to ', outdir, '/', events_filename))
    }
    save(list = c('events','params'), file = events_filename)
  }
  
  
}

if(run_get_ff_features){
  
  #Need to load files if you didn't run the first step
  if(!run_extract_ff_events){
    
    #Load data
    setwd(indir)
    if(verbose){
      print('Loading xy data')
    }
    load('hyena_xy_level1.RData')
    
    setwd(outdir)
    load(events_filename)
    
  } 
  
  #Run feature extraction
  if(verbose){
    print('Extracting features')
  }
  events <- get_ff_features(xs = xs, 
                             ys = ys, 
                            together.seqs = events,
                            params = params)
  
  if(overwrite_extract_ff_features){
    
    setwd(outdir)
    if(verbose){
      print(paste0('Saving events + features to ', outdir, '/', events_features_filename))
    }
    save(list = c('events','params'), file = events_features_filename)
    
  }
  
}

if(generate_day_randomization_plan){
  
  setwd(indir)
  load('hyena_ids.RData')
  load('hyena_day_start_idxs.RData')
  
  n.inds <- nrow(hyena.ids)
  rand.plan <- array(NA, dim = c(n.inds, params$last.day.used, n.rands))
  for(r in 1:n.rands){
    for(i in 1:n.inds){
      rand.plan[i,,r] <- sample(1:params$last.day.used)
    }
  }
  
  if(overwrite_day_randomization_plan){
    
    setwd(outdir)
    save(file = day_randomization_plan_filename, list = c('rand.plan','params','n.rands'))
    
  }
  
  
}

if(execute_day_randomization_plan){
  
  if(verbose){
    print("Executing day randomization plan")
    print('Loading data')
  }
  
  setwd(indir)
  load('hyena_xy_level1.RData')
  load('hyena_timestamps.RData')
  load('hyena_day_start_idxs.RData')
  timestamps.local <- timestamps + local.time.diff*60*60
  
  setwd(outdir)
  load(day_randomization_plan_filename)
  
  n.inds <- nrow(xs)
  
  events.rand.list <- list()
  
  for(r in 1:n.rands){
    
    if(verbose){
      print(paste0('Running randomization ', r, ' / ', n.rands))
    }
    
    if(verbose){
      print('Setting up randomized xs and ys matrices')
    }
    xs.rand <- ys.rand <- matrix(NA, nrow = nrow(xs), ncol = ncol(xs))
    for(i in 1:n.inds){
      for(d in 1:params$last.day.used){
        
        #original time indexes
        t0 <- day.start.idxs[d]
        tf <- day.start.idxs[d+1] - 1
        
        #time indexes to swap in for this randomization
        t0.swap <- day.start.idxs[rand.plan[i,d,r]]
        tf.swap <- day.start.idxs[rand.plan[i,d,r]+1] - 1
        
        ts.orig <- seq(t0, tf)
        ts.swap <- seq(t0.swap, tf.swap)
        
        xs.rand[i, ts.orig] <- xs[i, ts.swap]
        ys.rand[i, ts.orig] <- ys[i, ts.swap]
        
      }
      
    }
    
    #get events and extract features
    if(verbose){
      print('Extracting events and getting features')
    }
    events.rand <- get_ff_events_and_phases(xs = xs.rand, ys = ys.rand, params = params)
    events.rand <- get_ff_features(xs = xs.rand, ys = ys.rand, together.seqs = events.rand, params = params)
    events.rand.list[[r]] <- events.rand
    
  }
  
  if(overwrite_day_randomization_output){
    save(list = c('events.rand.list','params'), file = day_randomization_output_filename)
  }
  
}
