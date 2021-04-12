#HYENA 2017 COLLAR DATA PREPROCESSING SCRIPT - STEP 1

#Pre-process hyena data
#This script:
#	Removes unrealistic speeds and replaces with NAs ( > 99.95 percentile = 14.8 m/s)
#	Fills in missing data gaps less than max.interp = 5 sec with linear interpolation
#	Finds instances where hyena did not move during an NA gap of < max.move time = 5 minutes (moved < max.move=5 m) and replaces 
#		them with the mean location of the hyena between start and end of seq
# Find extreme distances (0.01% quantile or > 99.99% quantile of xs or ys for that individual), if they are not next to other values within 1000 m of that point, replace with NA
# Find extreme xs and ys above mean + sd * 10 for each ind and remove those

#Note: This script used to be called preprocess_hyena_gps.R - it was renamed to make the order of steps clear
#This script should be run after hyena_preprocess_0_extract_GPS_csv_to_R.R

#--------INPUTS----------
#Required data files to run:

#File: hyena_xy_level0.RData
#File: hyena_timestamps.RData
#File: hyena_ids.RData

#Required code files to run:
#None

#------------OUTPUTS---------
#File: hyena_xy_level1.RData
# Contains xs and ys with filtered GPS (as described above) analogous to hyena_xy_level0.RData

#File: hyena_xy_level1.h5
# Contains the same data in hdf5 format, as well as the timestamps object

#---------------PARAMETERS-------------
max.interp <- 5 #maximum number of seconds of NAs to interpolate (default = 5)
max.speed.quantile <- 0.9995 #maximum speed quantile to allow - speeds above this will be replaced with NAs (default = 0.9995, corresponds to 14.8 m/s)
max.move <- 5 #maximum distance moved (in meters) during an NA gap of duration <  max.move.time seconds, to replace all with the mean of the start and end point (default = 5 m)
max.move.time <- 5*60 #maximum NA gap time (in seconds) to replace with mean value if the individual did not move more than max.move meters (default = 5*60 sec = 5 min)
max.dist.quantile <- .9999 #maximum distance quantile. distances beyond this quantile in both x and y direction for each individual will be replaced with NAs (defaults to 0.9999)
min.dist.quantile <- .0001 #minimum distance quantile. distances below this quantile in both x and y direction for each individual will be replaced with NAs (defaults to 0.0001)

overwrite <- F #wheter the overwrite the output file (defaults to F)

#----DIRECTORY---- 
#Directory where data are stored (and output is stored)
dir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps'

#set working directory
setwd(dir)

#--------MAIN-----------
#LIBRARIES
library(lubridate)
library(rhdf5)

#Load level 0 gps data
load('hyena_xy_level0.RData')
load('hyena_timestamps.RData')
load('hyena_ids.RData')

#number of indidivudals
n.inds <- nrow(xs)
n.times <- ncol(xs)

print('initial NA frac:')
print(sum(is.na(xs))/length(xs))

#Find extreme speeds and replace with NAs
speeds <- matrix(NA,nrow=n.inds,ncol=n.times)
for(i in 1:n.inds){
	speeds[i,seq(1,n.times-1)] <- sqrt(diff(xs[i,])^2 + diff(ys[i,])^2)
}
max.speed <- quantile(speeds,max.speed.quantile,na.rm=T)

#replace unrealistic speeds with NAs
for(i in 1:n.inds){
   unrealistic.speed.idxs <- which(speeds[i,] > max.speed)
   xs[i,unrealistic.speed.idxs] <- NA
   ys[i,unrealistic.speed.idxs] <- NA
   xs[i,unrealistic.speed.idxs + 1] <- NA
   ys[i,unrealistic.speed.idxs + 1] <- NA
}

print('after removing unrealistic speeds:')
print(sum(is.na(xs))/length(xs))

#Interpolate through seqs of NAs of length < max.interp
for(i in 1:n.inds){

	#data for a single individual
	x <- xs[i,]
	y <- ys[i,]
	
	#vectors to hold interpolated data
	x.interp <- x
	y.interp <- y
	
	#get runs of NAs
	runs <- rle(is.na(x)) #get runs of NAs
	vals <- runs$values
	lens <- runs$lengths
	idxs <- c(0,cumsum(runs$lengths))
	
	#for each run of NAs, fill in with linearly interp values if length is less than max.interp
	for(j in 1:length(vals)){
		if(vals[j]==TRUE){
			first.idx <- idxs[j]+1
			last.idx <- idxs[j+1]
			
			#If not too near the beginning or the end of the sequence...
			if((first.idx > 1) & (last.idx < length(x))){
			
				#Get values before and after NA sequence
				prev.val.x <- x[first.idx-1]
				next.val.x <- x[last.idx + 1]
				prev.val.y <- y[first.idx-1]
				next.val.y <- y[last.idx +1]
				
				#Fill in with linear interpolation if the NA sequence is short (< 5)
				if(lens[j] <= max.interp){
					interp.vals.x <- seq(prev.val.x,next.val.x,length.out=lens[j]+2)
					interp.vals.y <- seq(prev.val.y,next.val.y,length.out=lens[j]+2)
					x.interp[seq(first.idx-1,last.idx+1)] <- interp.vals.x
					y.interp[seq(first.idx-1,last.idx+1)] <- interp.vals.y
				}
			
				#Otherwise...if less than 5 minutes gap...
				if((lens[j] > max.interp) & (lens[j] < max.move.time)){
					#Fill in with mean value at start and end if they are close enough ( <= max.move)
					dist.moved <- sqrt((next.val.x - prev.val.x)^2 + (next.val.y - prev.val.y)^2)
					time.elapsed <- last.idx - first.idx
					if(dist.moved < max.move){
						mean.x <- mean(c(next.val.x,prev.val.x))
						mean.y <- mean(c(next.val.y,prev.val.y))
						x.interp[seq(first.idx,last.idx)] <- mean.x
						y.interp[seq(first.idx,last.idx)] <- mean.y
					}
				}
			}
		}
	}
	
	xs[i,] <- x.interp
	ys[i,] <- y.interp
	
}

print('after interpolating:')
print(sum(is.na(xs))/length(xs))

print('removing unrealistic distances...')

#Find and remove unrealistic distances (if they are not surrounded by other similar distances)
for(i in 1:n.inds){
  print('ind:')
  print(i)
  xi <- xi.new <- xs[i,]
  yi <- yi.new <- ys[i,]
  non.nas <- which(!is.na(xi))
  
  #get very large or very small values of x and y (outside their normal range)
  bigs <- unique(c(which(xi>quantile(xi,max.dist.quantile,na.rm=T)),which(yi>quantile(yi,max.dist.quantile,na.rm=T))))
  smalls <- unique(c(which(xi<quantile(xi,min.dist.quantile,na.rm=T)),which(yi<quantile(yi,min.dist.quantile,na.rm=T))))
  extremes <- unique(c(bigs,smalls))
  print('number of extremes = ')
  print(length(extremes))
  
  #for each extreme value, find the previous and next non-NA data point
  #if this previous point is more than 1 km away, just replace the unrealistic location with NA
  #otherwise, leave it
  for(j in 1:length(extremes)){
    t.idx <- extremes[j]
    prev.t <- max(non.nas[which(non.nas < t.idx)])
    next.t <- min(non.nas[which(non.nas > t.idx)])
    
    dist.prev <- sqrt((xs[i,prev.t]-xs[i,t.idx])^2 + (ys[i,prev.t]-ys[i,t.idx])^2)
    dist.next <- sqrt((xs[i,next.t]-xs[i,t.idx])^2 + (ys[i,next.t]-ys[i,t.idx])^2)
    
    if(dist.prev > 1000 & dist.next > 1000){
      print(paste('found unrealistic distance at time:',t.idx,'dist.prev = ',dist.prev,'dist.next = ',dist.next,'dt1=',(t.idx-prev.t)/60,'dt2=',(next.t-prev.t)/60))
      
      xi.new[t.idx] <- NA
      yi.new[t.idx] <- NA
    }
    
  }
  
  xs[i,] <- xi.new
  ys[i,] <- yi.new
  
}

#Find and remove super unrealistic locations (>mean + 10*sd for each ind in x or y)
for(i in 1:n.inds){
  xi <- xs[i,]
  yi <- ys[i,]
  extremes.high.x <- which(xi > mean(xi,na.rm=T)+10*sd(xi,na.rm=T))
  extremes.low.x <- which(xi < mean(xi,na.rm=T)-10*sd(xi,na.rm=T))
  extremes.high.y <- which(yi > mean(yi,na.rm=T)+10*sd(yi,na.rm=T))
  extremes.low.y <- which(yi < mean(yi,na.rm=T)-10*sd(yi,na.rm=T))
  extremes <- c(extremes.high.x,extremes.high.y,extremes.low.x,extremes.low.y)
  
  if(length(extremes)>0){
    xs[i,extremes] <- NA
    ys[i,extremes] <- NA
  }
  
}

#Save as RData and HDF5
if(overwrite){
  setwd(dir)
  save(file='hyena_xy_level1.RData',list=c('xs','ys'))
  
  h5createFile(file='hyena_xy_level1.h5')
  h5write(file='hyena_xy_level1.h5',name='/xs',obj=xs)
  h5write(file='hyena_xy_level1.h5',name='/ys',obj=ys)
  h5write(file='hyena_xy_level1.h5',name='/timestamps',obj=as.character(timestamps))
  H5close()
}



