#Driver script to run fission-fusion analyses

################################ SET UP DIRECTORIES ##################################

#directory of where the original (processed) movement + vedba data is stored (+ metadata on IDs and den locations)
indir <- '/Volumes/EAS_shared/hyena/working/hyena_pilot_2017/processed/'

#directory of where to store extracted data for fission-fusion project
outdir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/data/'

#directory where code for fission-fusion project is stored
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'

################################ CHOOSE ANALYSES TO RUN ##################################

run_extract_ff_events <- T
overwrite_extract_ff_events <- T 
run_get_ff_features <- T
overwrite_extract_ff_features <- T

################################ PARAMETERS ##########################################

params <- list(R.fusion = 100, 
               R.fission = 200, 
               max.break = 60*30, #max time between events to merge events connected by NAs
               move.thresh = 5, #minimum amount moved to be considered 'moving' during a phase
               together.travel.thresh = 200) #minimum amount to be considered 'moving' during the together phase

verbose = TRUE

events_filename <- 'fission_fusion_events.RData'
events_features_filename <- 'fission_fusion_events_features.RData'

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

