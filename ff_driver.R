#Driver script to run fission-fusion analyses

################################ SET UP DIRECTORIES ##################################

#directory of where the original movement + acc data is stored
indir <- '/Volumes/EAS_shared/hyena/working/hyena_pilot_2017/processed/'

#directory of where to store extracted data
outdir <- '/Volumes/EAS_shared/hyena/working/hyena_fission_fusion/data/'

#directory where code is stored
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
  events <- get_ff_events_and_phases(xs, ys, 
                                     R.fusion = params$R.fusion, 
                                     R.fission = params$R.fission, 
                                     max.break = params$max.break, 
                                     verbose = verbose)
  if(overwrite_extract_ff_events){
    setwd(outdir)
    save(list = c('events','params'), file = 'fission_fusion_events.RData')
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
    load('fission_fusion_events.RData')
    
  } 
  
  #Run feature extraction
  events <- get_ff_features(xs, ys, 
                            together.seqs = events,
                            move.thresh = params$move.thresh,
                            together.travel.thresh = params$together.travel.thresh)
  
  if(overwrite_extract_ff_features){
    
    setwd(outdir)
    save(list = c('events','params'), file = 'fission_fusion_events_features.RData')
    
  }
  
}

