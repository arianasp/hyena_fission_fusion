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
output_day_randomization_plots <- T

################################ PARAMETERS ##########################################

params <- list(R.fusion = 100, 
               R.fission = 200, 
               max.break = 60*30, #max time between events to merge events connected by NAs
               move.thresh = 5, #minimum amount moved to be considered 'moving' during a phase
               together.travel.thresh = 200, #minimum amount to be considered 'moving' during the together phase
               local.time.diff = 3, # difference in hours from local time
               den.dist.thresh = 200
                ) 

randomization.type <- 'nightperm' #options: denblock, nightperm

verbose <- TRUE

######################## FIXED PARAMETERS - DON'T CHANGE ##############################

#parameters for a regular night randomization
nightperm.rand.params <- list(break.hour = 12, #which hour to "break" at when randomizing days (0 = midnight, 12 = noon)
                    last.day.used = 35, #last day to use in the randomizations (and real data)
                    blocks = NULL, #blocks to keep together for each individual (e.g. to keep den attendance roughly constant)
                    n.rands = 500 #how many randomizations to do
                    )

#parameters for a den block permutation
#this bit is manual (hardcoded), - the boundaries are based on looking at the plot of den attendance
den.blocks <- list()
den.blocks[[1]] <- c(12, 22, 30, 33)
den.blocks[[2]] <- c(13, 26)
den.blocks[[3]] <- c(14, 35)
den.blocks[[4]] <- c(15)
den.blocks[[5]] <- c(12, 31)
denblock.rand.params <- list(break.hour = 12, #which hour to "break" at when randomizing days (0 = midnight, 12 = noon)
                              last.day.used = 35, #last day to use in the randomizations (and real data)
                              blocks = den.blocks, #blocks to keep together for each individual (e.g. to keep den attendance roughly constant)
                              n.rands = 500 #how many randomizations to do
                              )

#den info
den.file.path <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/rawdata/metadata/hyena_isolate_dens.csv'
den.names <- c('DAVE D','RBEND D','RES M D1','DICK D')


events_filename <- 'fission_fusion_events.RData'
events_features_filename <- 'fission_fusion_events_features.RData'
day_start_idxs_filename <- 'hyena_day_start_idxs.RData'
day_randomization_plan_filename <- paste0('hyena_day_randomization_plan_', randomization.type, '.RData')
day_randomization_output_filename <- paste0('hyena_day_randomization_events_features_',randomization.type,'.RData')

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
                            params = params,
                            den.file.path = den.file.path,
                            den.names = den.names)
  
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
  
  if(randomization.type == 'nightperm'){
    rand.params <- nightperm.rand.params
    
    rand.plan <- array(NA, dim = c(n.inds, rand.params$last.day.used, rand.params$n.rands))
    for(r in 1:rand.params$n.rands){
      for(i in 1:n.inds){
        rand.plan[i,,r] <- sample(1:rand.params$last.day.used)
      }
    }
  }
  
  if(randomization.type == 'denblock'){
    
    rand.params <- denblock.rand.params
    
    rand.plan <- array(NA, dim = c(n.inds, rand.params$last.day.used, rand.params$n.rands))
    for(r in 1:rand.params$n.rands){
      for(i in 1:n.inds){
        blocks <- rand.params$blocks[[i]]
        blocks <- blocks[which(blocks < rand.params$last.day.used)] 
        blocks <- c(1, blocks, rand.params$last.day.used+1)
        for(j in 1:(length(blocks)-1)){
          d0 <- blocks[j]
          df <- blocks[j+1] - 1
          rand.plan[i,d0:df,r] <- sample(d0:df) #shuffle within block
        }
        
      }
      
    }
    
  }
  
  if(overwrite_day_randomization_plan){
    
    setwd(outdir)
    save(file = day_randomization_plan_filename, list = c('rand.plan','rand.params'))
    
  }
  
  
}

if(execute_day_randomization_plan){
  
  complete <- F
  
  if(verbose){
    print("Executing day randomization plan")
    print('Loading data')
  }
  
  setwd(indir)
  load('hyena_xy_level1.RData')
  load('hyena_timestamps.RData')
  load('hyena_day_start_idxs.RData')
  timestamps.local <- timestamps + params$local.time.diff*60*60
  
  setwd(outdir)
  load(day_randomization_plan_filename)
  
  n.inds <- nrow(xs)
  
  events.rand.list <- list()
  
  for(r in 1:rand.params$n.rands){
    
    if(verbose){
      print(paste0('Running randomization ', r, ' / ', rand.params$n.rands))
    }
    
    if(verbose){
      print('Setting up randomized xs and ys matrices')
    }
    xs.rand <- ys.rand <- matrix(NA, nrow = nrow(xs), ncol = ncol(xs))
    for(i in 1:n.inds){
      for(d in 1:rand.params$last.day.used){
        
        #original time indexes
        t0 <- day.start.idxs[d] + rand.params$break.hour*60*60 #get initial time (day start, plus break.hour)
        tf <- day.start.idxs[d+1] - 1 + rand.params$break.hour*60*60 #get final time
        
        #time indexes to swap in for this randomization
        t0.swap <- day.start.idxs[rand.plan[i,d,r]] + rand.params$break.hour*60*60
        tf.swap <- day.start.idxs[rand.plan[i,d,r]+1] - 1 + rand.params$break.hour*60*60
        
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
    events.rand <- get_ff_features(xs = xs.rand, 
                                   ys = ys.rand, 
                                   together.seqs = events.rand, 
                                   params = params,
                                   den.file.path = den.file.path,
                                   den.names = den.names)
    events.rand.list[[r]] <- events.rand
    
    if(overwrite_day_randomization_output){
      save(list = c('events.rand.list','params','rand.params','complete'), file = day_randomization_output_filename)
    }
    
  }
  
  complete <- T
  
  if(overwrite_day_randomization_output){
    save(list = c('events.rand.list','params','rand.params','complete'), file = day_randomization_output_filename)
  }
  
}

if(output_day_randomization_plots){
  
  #load data
  setwd(outdir)
  load(day_randomization_output_filename)
  load(events_features_filename)
  setwd(indir)
  load('hyena_timestamps.RData')

  #visualizations
  
  #event type frequency for all events - real vs permuted
  visualize_event_type_distributions(events, events.rand.list, rand.params, timestamps, remove.events.around.day.breaks = T)

  #event type frequency for den events and non-den events separate - real vs. permuted
  visualize_event_type_frequency_den_vs_nonden(events, events.rand.list, params, rand.params, timestamps)
  
  #some other properties of events such as duration, displacement in together phase, time of day, den vs. nonden fraction
  visualize_compare_event_properties(events, events.rand.list, params, rand.params, timestamps)
  
}
