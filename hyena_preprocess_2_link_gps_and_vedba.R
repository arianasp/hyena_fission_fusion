#HYENA 2017 COLLAR DATA PREPROCESSING SCRIPT - STEP 2

#This script reads in ACC and GPS data, computes the VeDBA and links them together, outputting into one file for all hyenas
#It also outputs a second file with only the VeDBA information

#Note: This script used to be called link_gps_and_vedba.R - it was renamed to make the order of steps clear

##This script should be run after hyena_preprocess_0_extract_GPS_csv_to_R.R and hyena_preprocess_1_filter_gps.R

#The script reads ACC files (see below) at 25 Hz. 
#It first computes the VeDBA from the triaxial ACC (window size = 1 sec, following Qasem et al. 2012). 
#It then smooths them over a window of 1 second and subsamples to 1 value per second, to link with GPS data

#It outputs two files - one a data frame with GPS and vedba data for all individuals and 
#one a matrix-form of vedbas only (same structure as xs and ys in hyena_xy_level1)

#--------INPUTS----------
#Required data files to run:

#GPS and IDs files
#File: hyena_xy_level1.RData
#File: hyena_timestamps.RData
#File: hyena_ids.RData

#ACC files
#File: cc16_352a_A_25Hz.h5
#File: cc16_352b_A_25Hz.h5
#File: cc16_354a_A_25Hz.h5
#File: cc16_360a_A_25Hz.h5
#File: cc16_366a_A_25Hz.h5

#------OUTPUTS------

#First output file (hyena_gps_and_acc.RData) contains:
# gps.vedba.all: data frame containing all gps and vedba data for all hyenas, with columns
#   $hyena.id: id of hyena, 1 to 5, matches to hyena.ids table
#   $time: timestamp (1 per second)
#   $x: easting
#   $y: northing
#   $vedba: vedba smoothed over a window of params$vedba.win, subsampled to 1 sec rate
#   $state: state taken from ACC distribution, based on params$threshes
#   $state2: state taken using a double threshold method, using params$threshes2
# params: list of parameters for generating the main data frame, has elements
#   $vedba.win: window (in seconds) used to smooth the vedba
#   $threshes: threshes on vedba to generate state
#   $threshes2: threshes on vedba to generate state2
#   $gps.file: gps file used
#   $times.file: times file used
#   $ids.file: hyena ids file used
# hyena.ids: table of information about hyena identities, including columns
#   $name: hyena name
#   collar: collar id
#   id: numeric id
#   color: color to use for plotting

#Second output file (hyena_vedba.RData) contains:
# params: same as above
# vedbas: matrix of size n.inds x n.timesteps of vedba values for each individual over time (comparable in structure to xs and ys matrices of GPS data)

#----------------------------PARAMETERS-----------------------------
vedba.win <- 1 #window for computing vedba, in sec
vedba.smooth.win <- 1 #window for smoothing the vedba, after it is computed, in sec
threshes <- c(.05,0.8) #thresholds of vedba for activity states (not used in fission-fusion analysis)
double.threshes <- exp(c(-9,-4,-2.8,-.4,-.1,3)) #thresholds for vedba for activity states, when using the double threshold method (not used in fission-fusion analysis)

#----------------------FILES AND DIRECTORIES------------------------
#Inputs
indir_gps <- processed.data.directory
indir_acc <- paste0(raw.data.directory, 'acc/')
gps.file <- paste0(processed.data.directory, 'hyena_xy_level1.RData')
times.file <- paste0(processed.data.directory,'hyena_timestamps.RData')
ids.file <- paste0(processed.data.directory,'hyena_ids.RData')

#Outputs
savefile <- paste0(processed.data.directory,'hyena_gps_and_acc.RData')
savefile.hdf5 <- paste0(processed.data.directory,'hyena_gps_and_acc.h5')
savefile.vedba <- paste0(processed.data.directory,'hyena_vedba.RData')

#-----------------------------MAIN----------------------------------

#LIBRARIES
library(rhdf5)
library(pracma)
library(lubridate)
library(fields)
library(viridis)
library(scales)

#LOAD FILES

#gps, times, ids
load(gps.file)
load(times.file)
load(ids.file)

#vedba from each hyena
acc.files <- list.files(indir_acc)
gps.vedba.all <- data.frame()
for(file in 1:length(acc.files)){
  if(verbose)
    print(paste('loading in file:',acc.files[file]))
  
  #get ACC data, start time, sample rate
  utc <- h5read(file = paste0(indir_acc, acc.files[file]),'/UTC')
  fs <- as.numeric(h5read(file = paste0(indir_acc, acc.files[file]),'/fs'))
  A <- h5read(file = paste0(indir_acc, acc.files[file]), name = '/A')
  vedba_old <- h5read(file = paste0(indir_acc, acc.files[file]), name = '/vedba')
  
  #Estimate static acceleration using sliding window
  Af <- matrix(NA, nrow=nrow(A), ncol = ncol(A))
  for(j in 1:3){
    Af[,j] <- stats::filter(A[,j],rep(1,fs*vedba.win)/(fs*vedba.win),method='convolution',sides=2)
  }
  
  #Then calculate VeDBA as vectorial difference from static acceleration
  Adiff <- A - Af
  vedba <- sqrt(Adiff[,1]^2 + Adiff[,2]^2 + Adiff[,3]^2)
  
  #get hyena name and collar id
  collar.id <- strsplit(acc.files[file],'_')[[1]][2]
  hyena.name <- hyena.ids$name[which(hyena.ids$collar==collar.id)]
  
  #smooth vedba at smoothing window
  vedba.sm <- stats::filter(vedba,rep(1,fs*vedba.smooth.win)/(fs*vedba.smooth.win),method='convolution',sides=2)
  
  #subsample to once per second
  vedba.sub <- vedba.sm[seq(1,length(vedba.sm),fs)]
  
  start.date <- paste(utc[1:3],collapse='-')
  start.time <- paste(utc[4:6],collapse=':')
  start <- strptime(paste(start.date,start.time,sep=' '),format='%Y-%m-%d %H:%M:%S',tz='GMT')
  
  #make a times vector
  vedba.times <- seq.POSIXt(from=start,by=1,length.out=length(vedba.sub))
  
  #put into a data frame
  vedba.df <- data.frame(time=vedba.times,vedba=as.numeric(vedba.sub))
  
  #threshold to get 3 activity states
  vedba.df$state <- 0
  vedba.df$state[which(vedba.df$vedba>threshes[1])] <- 1
  vedba.df$state[which(vedba.df$vedba>threshes[2])] <- 2
  vedba.df$state[which(is.na(vedba.df$vedba))] <- NA
  
  #link gps and vedba
  all.dat <- vedba.df
  rm(vedba.df)
  
  all.dat$hyena.id <- hyena.ids$id[grep(collar.id,hyena.ids$collar)]
  
  time.idxs <- match(all.dat$time,timestamps)
  all.dat$x <- xs[cbind(all.dat$hyena.id,time.idxs)]
  all.dat$y <- ys[cbind(all.dat$hyena.id,time.idxs)]
  
  #double threshold version of state definitions
  all.dat$init.state <- 1
  for(i in 1:length(double.threshes)){
    all.dat$init.state[which(all.dat$vedba > double.threshes[i])] <- i
  }
  all.dat$init.state[which(is.na(all.dat$vedba))] <- NA
  
  all.dat$state2 <- all.dat$init.state
  all.dat$state2[which((all.dat$init.state %% 2)==0)] <- 0
  unknowns <- which(all.dat$state2==0)
  states <- all.dat$state2
  for(i in 1:length(unknowns)){
    curr <- states[unknowns[i]]
    if(!is.na(curr)){
      idx <- 0
      prev <- 0
      while(prev==0){
        idx <- idx + 1
        prev <- states[unknowns[i]-idx]
        if(is.na(prev)){
          break
        }
      }
      states[unknowns[i]] <- prev
    }
  }
  
  all.dat$state2[which(states==1)] <- 0
  all.dat$state2[which(states==3)] <- 1
  all.dat$state2[which(states==5)] <- 2
  
  #add to final data frame
  gps.vedba.all <- rbind(gps.vedba.all,all.dat)
  
}

#get rid of unneccesary columns
gps.vedba.all <- gps.vedba.all[,c('hyena.id','time','x','y','vedba','state','state2')]

#save output
params <- list()
params$veba.win <- vedba.win
params$threshes <- threshes
params$threshes2 <- double.threshes
params$gps.file <- gps.file
params$times.file <- times.file
params$ids.file <- ids.file

#matrix with only vedba (as a matrix, structured same as xs and ys)
vedbas <- matrix(NA, nrow = nrow(xs), ncol = ncol(xs))
t.idxs <- match(gps.vedba.all$time, timestamps)
non.nas <- which(!is.na(t.idxs))
vedbas[cbind(gps.vedba.all$hyena.id[non.nas], t.idxs[non.nas])] <- gps.vedba.all$vedba[non.nas]

#-----SAVE-----

#R file
save(list=c('gps.vedba.all','params','hyena.ids'),file=savefile)

#HDF5 file
h5createFile(file=savefile.hdf5)
h5write(savefile.hdf5,name='/gps_vedba_all',obj=gps.vedba.all)
h5write(savefile.hdf5,name='/params',obj=params)
h5write(savefile.hdf5,name='/hyena_ids',obj=hyena.ids)

save(list = c('vedbas', 'params'), file = savefile.vedba)
