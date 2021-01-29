#Driver script to run fission-fusion analyses

################################ SET UP DIRECTORIES ##################################

#directory of where the original (processed) movement + vedba data is stored (+ metadata on IDs and den locations)
indir <- '/Volumes/EAS_shared/hyena/working/hyena_pilot_2017/processed/'

#directory of where to store extracted data for fission-fusion project
outdir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/data/'

#directory where code for fission-fusion project is stored
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'

################################ CHOOSE ANALYSES TO RUN ##################################

run_get_day_start_idxs <- F
overwrite_get_day_start_idxs <- F
run_extract_ff_events <- F
overwrite_extract_ff_events <- F 
run_get_ff_features <- F
overwrite_extract_ff_features <- F
generate_day_randomization_plan <- F
overwrite_day_randomization_plan <- T
execute_day_randomization_plan <- T
overwrite_day_randomization_output <- T

################################ PARAMETERS ##########################################

params <- list(R.fusion = 100, 
               R.fission = 200, 
               max.break = 60*30, #max time between events to merge events connected by NAs
               move.thresh = 5, #minimum amount moved to be considered 'moving' during a phase
               together.travel.thresh = 200, #minimum amount to be considered 'moving' during the together phase
               last.day.used = 35) #last day to use in the randomizations

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

if(run_get_day_start_idxs){
  
  setwd(indir)
  load('hyena_timestamps.RData')
  
  day.start.idxs <- get_day_start_idxs(timestamps, local.time.diff = local.time.diff)
  
  if(overwrite_get_day_start_idxs){
    
    setwd(outdir)
    save(file = day_start_idxs_filename, list = c('day.start.idxs'))
    
  }
  
} else{
  
  load(day_start_idxs_filename)
  
}

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
  
  setwd(outdir)
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
  timestamps.local <- timestamps + local.time.diff*60*60
  
  setwd(outdir)
  load(day_randomization_plan_filename)
  load('hyena_day_start_idxs.RData')
  
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
        
        t0.dt.daystart <- as.numeric(timestamps.local[t0] - floor_date(timestamps.local[t0], 'day'), units = 'secs')
        #tf.dt.dayend <- as.numeric(ceiling_date(timestamps.local[tf], 'day') - timestamps.local[tf], units = 'secs')
        t0.swap.dt.daystart <- as.numeric(timestamps.local[t0.swap] - floor_date(timestamps.local[t0.swap], 'day'), units = 'secs')
        #tf.swap.dt.dayend <- as.numeric(ceiling_date(timestamps.local[tf.swap], 'day') - timestamps.local[tf.swap], units = 'secs')
        
        if(t0.dt.daystart > 0){
          shift <- t0.dt.daystart - t0.swap.dt.daystart
          t0.swap <- t0.swap + shift
        }
        if(t0.swap.dt.daystart > 0){
          shift <- t0.swap.dt.daystart - t0.dt.daystart
          t0 <- t0 + shift
        }
        
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
