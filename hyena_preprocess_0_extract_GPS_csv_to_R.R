#HYENA 2017 COLLAR DATA PREPROCESSING SCRIPT - STEP 0

#Extract hyena GPS data from raw Technosmart csv files and save as R files

#Note: this file was formerly named "hyena_extract_GPS_csv_to_R.R"
#It has been renamed to "hyena_preprocess_0_extract_GPS_csv_to_R.R" to make the order of steps clear

#--------INPUTS----------
#Required data files to run (one csv file with GPS output for each hyena):
#cc16_352aGPS.csv: WRTH
#cc16_352bGPS.csv: BORA
#cc16_354aGPS.csv: BYTE
#cc16_360aGPS.csv: MGTA
#cc16_366aGPS.csv: FAY 

#Required code files to run:
#'hyena_functions.R'

#------------OUTPUTS--------
#File: hyena_gps_level0.RData - this file contains all the GPS info in multiple formats and as a result is very large! other files contain subsets for easier accessibility
#Contents:
# lats, lons: [n.inds x n.timesteps numeric] matrices of the latitude or longitude of each hyena at each timestep
# xs, ys: [n.inds x n.timesteps numeric] matrices of the UTM easting (xs) and northing (ys) location of each hyena at each timestep
# timestamps [n.timesteps POSIX]: vector of timestamps associated with the columns of lats, lons, xs, and ys matrices
# hyena.ids: [data frame] containing information about each hyena, with row indexes matching the rows of lats, lons, xs, and ys matrices
#     - name: [char] name of the hyena
#     - collar: [char] collar id of the hyena
#     - color: [char] what color it is normally plotted with when we make visualizations
# hyena.gps: [data frame] containing raw technosmart output in table form, plus UTM coordinates
# dates: [date] vector of dates included in the data
# day.start.idxs [numeric]: vector of indexes to the starts of days (midnight in local time) in the xs, ys, lats, and lons matrices

#File: hyena_xy_level0.RData - this file contains only xs and ys matrices, as described above

#File: hyena_timestamps.Rdata: contains only timestamps, as described above

#File: hyena_day_start_idxs.RData: contains only day.start.idxs, as described above

#File: hyena_latlon_level0.RData: contains lats and lons matrices, as described above

#File: hyena_ids.RData: contains hyena.ids data frame, as described above

#---------PARAMETERS--------
#whether to overwrite saved files on the server
overwrite <- F

#specify time zone offset from GMT
tz_offset <- 3 #3 hours

#Libraries
library(lubridate)

#Get useful functions
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'
setwd(codedir)
source('hyena_functions.R')

#Directory where csv files are stored
indir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/rawdata/gps_technosmart'

#Directory where to store output
outdir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps'


#------MAIN--------
#Set working directory
setwd(indir)

#First create a data frame (hyena.ids) that maps hyenas to collar ids
hyena.ids <- data.frame(name = c('WRTH','BORA','BYTE','MGTA','FAY'), collar = c('352a','352b','354a','360a','366a'), id = seq(1,5,1),color=c('red','blue','green','magenta','orange'))
#Get all files in folder
files <- list.files()
files <- files[grep('.csv',files)]
hyena.ids$name <- as.character(hyena.ids$name)
hyena.ids$collar <- as.character(hyena.ids$collar)
hyena.ids$color <- as.character(hyena.ids$name)

#Read everything into a single data table (hyena.gps)
hyena.gps.list <- list()
for(i in 1:nrow(hyena.ids)){
	print(i)
	file <- paste('cc16_',hyena.ids$collar[i],'GPS.csv',sep='')
	collar <- hyena.ids$collar[i]
	dat.curr <- read.table(files[i],sep='\t',header=TRUE,stringsAsFactors=F)
	dat.curr$id <- i
	dat.curr$name <- hyena.ids$name[i]
	hyena.gps.list[[i]] <- dat.curr
}
hyena.gps <- do.call(rbind,hyena.gps.list)
rm(hyena.gps.list)

#add a datetime column
datetimes <- gsub('[.]00','',hyena.gps$timestamp)
hyena.gps$datetime <- strptime(datetimes,format='%d/%m/%Y %H:%M:%S',tz='UTC')

#Filter out data with NA timestamps
hyena.gps <- hyena.gps[which(!is.na(hyena.gps$datetime)),]

#Convert to UTM
utms <- latlon.to.utm(cbind(hyena.gps$location.long,hyena.gps$location.lat),utm.zone=36,southern_hemisphere=TRUE)
hyena.gps$easting <- utms[,1]
hyena.gps$northing <- utms[,2]

#Get list of timestamps (from first to last)
first.time <- floor_date(min(hyena.gps$datetime + tz_offset*60*60,na.rm=T), 'day') - tz_offset*60*60
last.time <- ceiling_date(max(hyena.gps$datetime + tz_offset*60*60 ,na.rm=T), 'day')  - tz_offset*60*60 - 1
timestamps <- seq(first.time,last.time,by=1,tz='UTC')

#Get matrices of lats and lons
n <- nrow(hyena.ids)
lats <- lons <- xs <- ys <- matrix(NA,nrow=n,ncol=length(timestamps))
for(i in 1:n){
	curr.dat <- hyena.gps[which(hyena.gps$id==i),]
	if(nrow(curr.dat)>0){
		idxs <- match(as.character(curr.dat$datetime),as.character(timestamps))
		non.nas <- which(!is.na(idxs))
		lats[i,idxs[non.nas]] <- curr.dat$location.lat[non.nas]
		lons[i,idxs[non.nas]] <- curr.dat$location.long[non.nas]
		xs[i,idxs[non.nas]] <- curr.dat$easting[non.nas]
		ys[i,idxs[non.nas]] <- curr.dat$northing[non.nas]
	}
}

#Get list of dates and indexes to starts of (complete) days (so skip the first Dec 31 hours)
days <- as.Date(timestamps + tz_offset*60*60)
dates <- sort(unique(days))
day.start.idxs <- c()
for(i in 1:length(dates)){
	day.start.idxs <- c(day.start.idxs, min(which(days==dates[i])))
}

if(overwrite){
  setwd(outdir)
  save(file='hyena_gps_level0.RData',list=c('xs','ys','lats','lons','timestamps','hyena.ids','hyena.gps','dates','day.start.idxs'))
  save(file='hyena_xy_level0.Rdata',list=c('xs','ys'))
  save(file='hyena_timestamps.Rdata',list=c('timestamps'))
  save(file='hyena_day_start_idxs.RData', list = c('day.start.idxs'))
  save(file='hyena_latlon_level0.RData',list=c('lats','lons'))
  save(file='hyena_ids.RData',list=c('hyena.ids'))
}
