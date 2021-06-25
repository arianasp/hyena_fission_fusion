#Driver script to run fission-fusion analyses


#----------------------------- SET UP DIRECTORIES---------------------------------------

user <- Sys.info()['user']
if(user == 'strau'){
  remote.stem <- 'Z:\\'
  code.stem <- '~/../Dropbox/Documents/Research/Partial_projects/'
}else if(user == 'straussed'){
  raw.data.directory <- '~/Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/raw_data/'
  processed.data.directory <- '~/Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/processed_data/'
  results.directory <- '~/Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/results/'
  code.directory <- '~/Documents/code/hyena_fission_fusion/'
}else{
  remote.stem <- '/Volumes/EAS_shared/'
  code.stem <- '~/Dropbox/code_ari/'
}



#-------------------------- MAIN FUNCTION -------------------------------
runall <- function(
  ### Run parameters
  R.fusion = 100,  # threshold proximity for triggering a fission-fusion event
  R.fission = 200, # threshold proximity defining start and end of ff event
  n.rands = 100, # number of iterations of the reference model to run
  ensure.no.day.matches = T, # should randomizations exclude randomly reproduced original data?
  
  ### input/output directories
  raw.data.directory, # directory of raw data, containing metadata, gps data, acc data, den locations, hyena photograph, satllite map
  processed.data.directory, # directory for output of preprocessing steps
  results.directory, # directory for results; this will be created if it doesnt yet exist
  code.directory, # directory where code is located 
  
  ### Run options
  verbose = T, # Should progress be printed to console?
  preprocess = T, # Do you want to run preprocessing and overwrite output? (preprocessing steps are identical regardless of parameter values) 
  extract.ff.events = T, # Do you want to extract fission-fusion events and overwrite output?
  get.ff.features = T, # Do you want to extract features of ff events and overwrite output?
  execute.day.randomization = T, # Do you want to run the reference model and overwrite output? 
  get.sync.measures = T # Do you want to calculate synchrony measures?
){
  ################################# SOURCE FUNCTIONS #######################################
  
  source(paste0(code.directory, 'ff_functions_library.R'))
  
  
  ################################ PARAMETERS ##########################################
  
  params <- list(R.fusion = R.fusion,  #default 100
                 R.fission = R.fission, #default 200 
                 max.break = 60*30, #max time between events to merge events connected by NAs
                 move.thresh = 5, #minimum amount moved to be considered 'moving' during a phase
                 together.travel.thresh = 200, #minimum amount to be considered 'moving' during the together phase
                 local.time.diff = 3, # difference in hours from local time
                 den.dist.thresh = 200, #threshold to consider something 'at the den'
                 last.day.used = 35 #last day in the data set to use (this is the day before the first collar died)
  )
  
  run.params <- params
  
  ######################## FIXED PARAMETERS - DON'T CHANGE ##############################
  
  #parameters for a den block permutation
  #this bit is manual (hardcoded), - the boundaries are based on looking at the plot of den attendance
  den.blocks <- list()
  den.blocks[[1]] <- c(11, 21, 29, 33)
  den.blocks[[2]] <- c(12, 25)
  den.blocks[[3]] <- c(13)
  den.blocks[[4]] <- c(15)
  den.blocks[[5]] <- c(11, 30)
  denblock.rand.params <- list(break.hour = 12, #which hour to "break" at when randomizing days (0 = midnight, 12 = noon)
                               last.day.used = params$last.day.used - 1, #last day to use in the randomizations (and real data)
                               blocks = den.blocks, #blocks to keep together for each individual (e.g. to keep den attendance roughly constant)
                               ensure.no.day.matches = ensure.no.day.matches, #whether to ensure that no pair of individuals is randomized to the same day
                               n.rands = n.rands #how many randomizations to do
  )
  
  #den info
  den.names <- c('DAVE D','RBEND D','RES M D1','DICK D')
  den.file.path <- paste0(raw.data.directory, 'metadata/hyena_isolate_dens.csv')
  
  rand.params <- denblock.rand.params
  
  
  events_filename <- 'fission_fusion_events.RData'
  events_features_filename <- 'fission_fusion_events_features.RData'
  day_start_idxs_filename <- 'hyena_day_start_idxs.RData'
  day_randomization_plan_filename <- paste0('hyena_day_randomization_plan_denblock_avoidmatch', rand.params$ensure.no.day.matches, '.RData')
  day_randomization_output_filename <- paste0('hyena_day_randomization_events_features_denblock_avoidmatch', rand.params$ensure.no.day.matches, '.RData')
  
  
  
  # Adjust subdirectory to be where to store extracted data based on parameters. 

  if(R.fusion == 100 & R.fission == 200){
    results.directory.param <- paste0(results.directory, '1_main_output_R100_200/')
    data.outdir <- paste0(results.directory.param, 'data/')
    plots.outdir <- paste0(results.directory.param, 'plots/')
  } else {
    results.directory.param <- paste0(results.directory, 'stability_check_R', R.fusion, '_', R.fission, '/')
    data.outdir <- paste0(results.directory.param, 'data/')
    plots.outdir <- paste0(results.directory.param, 'plots/')
  }
  
  # Check if these output directories exist. If not, create them.
  for (x in list(data.outdir, plots.outdir)){
    if(!dir.exists(x)) dir.create(x, recursive=T)
  } 
  
  

  ################################# PREPROCESS DATA #######################################
  
  if(preprocess){
    scripts = list('0'='hyena_preprocess_0_extract_GPS_csv_to_R.R', 
                   '1'='hyena_preprocess_1_filter_gps.R',
                   '2'='hyena_preprocess_2_link_gps_and_vedba.R')
    
    if (verbose) print('Starting data preprocessing')
    for (i in seq_along(scripts)){
      if (verbose) print(paste0('Starting step ',names(scripts)[i],'...'))
      
      source(paste0(code.directory, scripts[[i]]), local = T)
      
      if (verbose) print(paste0('... step ',names(scripts)[i],' finished'))
    }
    if (verbose) print('Completed data preprocessing')
    
    params <- run.params ## params object is overwritten during preprocessing. reset to run parameters
  }
  
  
  ################################# MAIN #######################################
  
  # ********************** CLEAN GPS AND VEDBA DATA ***********************************  
  # this part removes some data from GPS and vedba  
  # we need to do this preprocessing if we plan to either extract events, extract features
  # or execute day randomization
  
  if(extract.ff.events | get.ff.features | execute.day.randomization){
    
    #Load data
    if(verbose) print('Loading xy data')
    load(paste0(processed.data.directory,'/hyena_xy_level1.RData')) # loads xs, xy
    load(paste0(processed.data.directory, '/hyena_day_start_idxs.RData')) # loads day.start.idxs
    
    #remove everything after the last day used
    t.idxs.use <- 1:(day.start.idxs[params$last.day.used+1]-1)
    xs <- xs[,t.idxs.use]
    ys <- ys[,t.idxs.use]
    
    #Set the first and last 12 hours to NAs (to match with randomizations, which do not include these hours)
    xs[,1:rand.params$break.hour * 60 * 60] <- NA
    ys[,1:rand.params$break.hour * 60 * 60] <- NA
    xs[,(ncol(xs) - rand.params$break.hour * 60 * 60 + 1):ncol(xs)] <- NA
    ys[,(ncol(xs) - rand.params$break.hour * 60 * 60 + 1):ncol(xs)] <- NA
    
    #since vedba file also has something called params, recover correct params 
    load(paste0(processed.data.directory, 'hyena_vedba.RData')) # loads vedba, params
    params <- run.params
    
    #remove everything after the last day used
    vedbas <- vedbas[,t.idxs.use]
  }
  
  
  # ********************** EXTRACT FF EVENTS ***********************************
  # this part extracts ff events ...
  # creates object saved under events_filename
  
  if(extract.ff.events){
    if(verbose) print('Extracting events')
    events <- get_ff_events_and_phases(xs = xs, 
                                       ys = ys, 
                                       params = params)
    
    if(verbose) print(paste0('Saving events to ', data.outdir, events_filename))
    save(list = c('events','params'), file = paste0(data.outdir, events_filename))
  }
  else{
    load(paste0(data.outdir, events_filename))
  }
  
  # ********************** GET FF FEATURES ***********************************
  # this part extracts ff features  ... 
  # creates object saved under events_features_filename
  
  #Run feature extraction
  if(get.ff.features){
    if(verbose) print('Extracting features')
    events <- get_ff_features(xs = xs, 
                              ys = ys, 
                              together.seqs = events,
                              params = params,
                              den.file.path = den.file.path,
                              den.names = den.names,
                              vedbas = vedbas,
                              get.sync.measures = get.sync.measures,
                              sync.subsample = 10)
    
    
    if(verbose) print(paste0('Saving features to ', data.outdir, events_features_filename))
    save(list = c('events','params'), file = paste0(data.outdir, events_features_filename))
  }
  else{
    load(paste0(data.outdir, events_features_filename))
  }
  
  # ********************** EXECUTE DAY RANDOMIZATION ***********************************
  # this part creates and then analyzes randomized versions of the original data
  
  if(execute.day.randomization){
    
    load(paste0(processed.data.directory, 'hyena_ids.RData')) # loads hyena.ids
    n.inds <- nrow(hyena.ids)
    rand.params <- denblock.rand.params 
    
    if (verbose) print('Generating randomization plan') 
    
    rand.plan <- generate_randomization_plan(rand.params, n.inds = n.inds)
    save(file = paste0(data.outdir, day_randomization_plan_filename), list = c('rand.plan','rand.params'))
    
    complete <- F
    
    if(verbose) writeLines("Executing day randomization plan\nLoading data")
    
    load(paste0(processed.data.directory,'hyena_timestamps.RData')) # loads timestamps
    timestamps.local <- timestamps + params$local.time.diff*60*60
    
    load(paste0(data.outdir, day_randomization_plan_filename)) # loads rand.plan and rand.params
    
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
      if(verbose) print('Extracting events and getting features')
      events.rand <- get_ff_events_and_phases(xs = xs.rand, ys = ys.rand, params = params)
      events.rand <- get_ff_features(xs = xs.rand, 
                                     ys = ys.rand, 
                                     together.seqs = events.rand, 
                                     params = params,
                                     den.file.path = den.file.path,
                                     vedbas = vedbas,
                                     get.sync.measures = get.sync.measures,
                                     den.names = den.names,
                                     sync.subsample = 10)
      
      events.rand.list[[r]] <- events.rand
      save(list = c('events.rand.list','params','rand.params','complete'), file = paste0(data.outdir, day_randomization_output_filename))
      
      complete <- T
      save(list = c('events.rand.list','params','rand.params','complete'), file = paste0(data.outdir, day_randomization_output_filename))
    }
    
    invisible(c(data.outdir, plots.outdir))
  }
}



#-----------------------------------RUN ME-------------------------------------------------
set.seed(43410)
print('--------------------------- MAIN RESULTS ---------------------------------')
output.dirs <- runall(ensure.no.day.matches = T, R.fusion = 100, R.fission = 200, n.rands = 4, 
                      raw.data.directory, processed.data.directory, results.directory, code.directory, preprocess = T,
                      execute.day.randomization = T, extract.ff.events = T, get.sync.measures = T, get.ff.features = T)
generate_figures(output.dirs[1], output.dirs[2], code.directory)

print('--------------------------- CHECK SMALLER THRESHOLD ---------------------------------')
output.dirs <- runall(ensure.no.day.matches = T, R.fusion = 50, R.fission = 100, n.rands = 4, 
                      raw.data.directory, processed.data.directory, results.directory, code.directory, preprocess = F,
                      execute.day.randomization = T, extract.ff.events = T, get.sync.measures = T, get.ff.features = T)
generate_figures(output.dirs[1], output.dirs[2], code.directory)
print('--------------------------- CHECK LARGER THRESHOLD ---------------------------------')
output.dirs <- runall(ensure.no.day.matches = T, R.fusion = 200, R.fission = 300, n.rands = 4, 
                      raw.data.directory, processed.data.directory, results.directory, code.directory, preprocess = F,
                      execute.day.randomization = T, extract.ff.events = T, get.sync.measures = T, get.ff.features = T)
generate_figures(output.dirs[1], output.dirs[2], code.directory)



