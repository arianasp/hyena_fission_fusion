#Identify events where two hyenas are together within 200 m
#Beginning is fusion, end is fission

#First just use one threshold to identify start and end times (but maybe double threshold later to get rid of 'flicker')

#------------- SETUP ----------------

#PARAMETERS
R.fusion <- 200 # distance threhsold to start an 'interaction'
R.fission <- 300 #distance threshold to end an 'interaction' (should be > R.fusion)
max.break <- 60*30 #maximum length of NA sequence between events to clump them together
hyena.cols <- c('red','orange','blue','magenta','green')
grid.res <- 100

#DIRECTORIES - MODIFY THIS TO YOUR OWN COMPUTER SETUP
name <- readline(prompt = "Enter name:")

if(grepl('a|A', name)){
  basedir <- '~/Dropbox/hyenas/hyena_fission_fusion/' # <-- CHANGE TO THE MAIN PROJECT FOLDER
  codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/' # <-- CHANGE TO THE DIRECTORY WHERE CODE REPOSITORY IS
}else if(grepl('e|E', name)){
  basedir <- 'C:/Users/strau/Dropbox/hyena_fission_fusion/' # <-- CHANGE TO THE MAIN PROJECT FOLDER
  codedir <- 'C:/Users/strau/Dropbox/Documents/Research/Partial_projects/hyena_fission_fusion/' # <-- CHANGE TO THE DIRECTORY WHERE CODE REPOSITORY IS
}else{
  cat('Name not recognized')
}


#OVERWRITE PERMISSION
overwrite <- T #set to T if you want to save outputs of script, overwriting previous versions
generate_animations <- F #set to T if you want to generate animations of the fission-fusion events

#SUBDIRECTORIES
datadir <- paste(basedir, 'data/raw', sep = '')
processeddir <- paste(basedir, 'data/processed', sep = '') #this is where the output will be saved (the fission fusion events RData file)

#FUNCTIONS
setwd(codedir)
source('hyena_functions.R')

#----------- MAIN -----------------
#Set working directory for main analysis
setwd(datadir)

#Load data
load('hyena_xy_level1.RData')
load('hyena_timestamps.Rdata')
load('hyena_ids.RData')

#Read in den data
den.file <- 'hyena_isolate_dens.csv'
den.names <- c('DAVE D','RBEND D','RES M D1','DICK D')
known.locs <- read.csv(den.file,stringsAsFactors=F)
eastsNorths <- latlon.to.utm(cbind(known.locs$lon,known.locs$lat),southern_hemisphere=T,utm.zone=36)
known.locs$east <- eastsNorths[,1]
known.locs$north <- eastsNorths[,2]
den.locs <- known.locs[which(known.locs$name %in% den.names),]
den.locs$name <- as.character(den.locs$name)

#Get basic parameters
n.inds <- nrow(xs)
n.times <- ncol(xs)

#Compute dyadic distance over time
dyad.dists <- array(NA, dim=c(n.inds, n.inds, n.times))
for (i in 1:n.inds) {
  dyad.dists[i,,] <- sqrt( (sweep(xs,2,xs[i,],FUN='-'))^2 + (sweep(ys,2,ys[i,],FUN='-'))^2 )
  dyad.dists[i,i,] <- NA
}

#Get contiguous sequences of "together" times for each pair and store in data frame
#data frame has columns:
# i, j (the two individuals)
# t0, tf (beginning and ending time of 'bout' of being together)
# t0.exact, tf.exact (indicates whether the time step immediately before and after the 'bout' was present (TRUE) or missing/NA (FALSE))
together.seqs <- data.frame()
for(i in 1:(n.inds-1)){
  for(j in (i+1):n.inds){
    dists <- dyad.dists[i,j,]
    prev.together <- NA
    for(k in 1:n.times){
      
      #get current distance apart
      curr.dist <- dists[k]
      
      #if coming from an NA sequence...
      if(is.na(prev.together)){
        #if still NA, pass
        if(is.na(curr.dist)){
          next
          
        #if not NA...
        } else{
          
          #if within R.fusion, start a together sequence
          if(curr.dist < R.fusion){
            t0 <- k
            prev.together <- TRUE
            next
            
          #otherwise, do not start a sequence but set prev.together to FALSE
          } else{
            prev.together <- FALSE
          }
        }
      }
      
      #if previously not together
      if(prev.together == FALSE){
        
        #if you run into an NA, set prev.together to NA
        if(is.na(curr.dist)){
          prev.together <- NA
          next
        }
        
        #if curr.dist is less than R.fusion, start a sequence
        if(curr.dist < R.fusion){
          t0 <- k
          prev.together <- TRUE
          next
          
        #otherwise pass
        } else{
          next
        }
        
      }

      #if previously together...
      if(prev.together == TRUE){
        
        #if you run into an NA, end the sequence, add to table, set prev.together to NA
        if(is.na(curr.dist)){
          tf <- k - 1
          together.seqs <- rbind(together.seqs, data.frame(i = i, j = j, t0 = t0, tf = tf))
          prev.together <- NA
          next
        }
        
        #if moved away to > R.fission distance, end sequence, set prev.together to FALSE
        if(curr.dist > R.fission){
          tf <- k -1
          together.seqs <- rbind(together.seqs, data.frame(i = i, j = j, t0 = t0, tf = tf))
          prev.together <- FALSE
          next
        }
      }
    }
  }
}
together.seqs.orig <- together.seqs

#Aggregate events that are close together in time and separated only by a sequence of NAs of maximum length max.break
together.seqs <- together.seqs.orig
converged <- F
while(!converged){
  print(nrow(together.seqs.orig))
  new.seqs <- data.frame()
  for(i in 1:(n.inds-1)){
    for(j in (i+1):n.inds){
      rows <- which(together.seqs$i == i & together.seqs$j == j)
      k <- 1
      while(k < length(rows)){
        break.start <- together.seqs$tf[rows[k]] + 1
        break.end <- together.seqs$t0[rows[k+1]] - 1
        
        #if two together periods are separated only by NAs and less than max.break seconds, aggregate into one event
        if((break.end - break.start <= max.break) & 
           (sum(!is.na(dyad.dists[i,j,break.start:break.end])) == 0)){
            new.seqs <- rbind(new.seqs, data.frame(i = i, j = j, t0 = together.seqs$t0[rows[k]], tf = together.seqs$tf[rows[k+1]]))
            k <- k + 2
        } else{
            new.seqs <- rbind(new.seqs, together.seqs[rows[k],])
            k <- k + 1
            
            #add the last row if needed
            if(k == length(rows)){
              new.seqs <- rbind(new.seqs, together.seqs[rows[k],])
            }
        }
      }
    }
  }
  
  #check if aggregation has converged
  if(nrow(new.seqs) == nrow(together.seqs)){
    converged <- T
  } else{
    converged <- F
  }
  
  print(nrow(together.seqs))
  print(nrow(new.seqs))
  
  #replace with new sequences
  together.seqs <- new.seqs
}


#add column for time until next event
together.seqs$dt <- NA
for(i in 1:(n.inds-1)){
  for(j in (i+1):n.inds){
    rows <- which(together.seqs$i == i & together.seqs$j == j)
    dt <- together.seqs$t0[rows[2:length(rows)]] - together.seqs$tf[rows[1:(length(rows)-1)]]
    together.seqs$dt[rows[1:(length(rows)-1)]] <- dt
  }
}

#add columns for whether start and end time are exact (whether there is an NA before or after the event)
befores <- dyad.dists[cbind(together.seqs$i, together.seqs$j, together.seqs$t0-1)]
afters <- dyad.dists[cbind(together.seqs$i, together.seqs$j, together.seqs$tf+1)]
together.seqs$t0.exact <- !is.na(befores)
together.seqs$tf.exact <- !is.na(afters)

#get closest approach distance
together.seqs$closest.app <- NA
for(i in 1:nrow(together.seqs)){
  dists <- dyad.dists[together.seqs$i[i], together.seqs$j[i], together.seqs$t0[i]:together.seqs$tf[i]]
  together.seqs$closest.app[i] <- min(dists,na.rm=T)
}

#save sequences
setwd(processeddir)
if(overwrite){
  save(list=c('together.seqs','R.fusion','R.fission','max.break'),file='fission_fusion_events.RData')
}

#plot image sequences of each event to turn into movies (if specified)
if(generate_animations){
  #create a directory to hold animations if it doesn't exist yet
  img.outdir <- paste(resultdir,'/animations',sep='')
  if(!file.exists(img.outdir)){
    dir.create(img.outdir)
  }
  setwd(img.outdir)
  
  #for each sequence, create a subdirectory to hold images and then put images in it
  img.idx <- 1
  for(idx in 1:nrow(together.seqs)){
    print(paste(idx, '/', nrow(together.seqs)))
    before.time <- 1200
    after.time <- 1200
    i <- together.seqs$i[idx]
    j <- together.seqs$j[idx]
    t0 <- together.seqs$t0[idx]
    tf <- together.seqs$tf[idx]
    
    #subdirectory
    #subdir <- paste(img.outdir, '/event',idx, '_', hyena.ids$name[i], '_', hyena.ids$name[j], '_', t0,'-',tf, sep = '')
    #dir.create(subdir)
    #setwd(subdir)
    
    if((t0-before.time)>0){
      tstart <- t0 - before.time
    } else{
      tstart <- 1
    }
    
    if(tf + after.time <= n.times){
      tend <- tf + after.time
    } else{
      tend <- n.times
    }
    
    xi <- xs[i, tstart:tend]
    yi <- ys[i, tstart:tend]
    xj <- xs[j, tstart:tend]
    yj <- ys[j, tstart:tend]
    
    xt <- xs[,tstart:tend]
    yt <- ys[,tstart:tend]
    xrange <- range(c(xi,xj),na.rm=T)
    yrange <- range(c(yi,yj),na.rm=T)
    #img.idx <- 1
    for(t in seq(1,length(xi),60)){
      png(filename = paste(img.idx,'.png',sep=''), width = 800, height = 800, units = 'px')
      plot(NULL, xlim=xrange, ylim = yrange, asp = 1, xlab = 'Easting (m)',ylab = 'Northing (m)', main = paste('Event ',idx,' (',hyena.ids$name[i], ' + ',hyena.ids$name[j],')',sep=''))
      #lines(c(min(xrange), min(xrange + 200)), c(max(yrange), max(yrange)),lwd=3, col = 'purple')
      
      #grid lines every 100 m
      xgrid <- seq(floor(min(xrange)/grid.res)*grid.res, ceiling(max(xrange)/grid.res)*grid.res, grid.res)
      ygrid <- seq(floor(min(yrange)/grid.res)*grid.res, ceiling(max(yrange)/grid.res)*grid.res, grid.res)
      abline(h=ygrid,col='gray',lwd=.3)
      abline(v=xgrid,col='gray',lwd=.3)
      
      for(ind in 1:n.inds){
        lines(xt[ind,1:t],yt[ind,1:t],col=hyena.ids$color[ind])
      }
      for(ind in 1:n.inds){
        if(ind == i){
          points(xt[ind,t],yt[ind,t],pch=19,col=hyena.ids$color[ind],cex=2)
        } else{
          if(ind == j){
            points(xt[ind,t],yt[ind,t],pch=19,col=hyena.ids$color[ind],cex=2)
          } else{
              points(xt[ind,t],yt[ind,t],pch=19,col=hyena.ids$color[ind],cex=1)
          }
        }
      }
      if(t > before.time & t < (length(xi)-after.time)){
        lines(c(xi[t],xj[t]),c(yi[t],yj[t]),col='black',lwd=2)
      }
      points(den.locs$east,den.locs$north,pch=0,cex=2)
      legend('topright',legend = hyena.ids$name, col = hyena.ids$color, pch = 19)
      #print time (local time)
      mtext(text = timestamps[t + tstart - 1] + 3*60*60, side = 3, col = 'blue')
      dev.off()
      img.idx <- img.idx + 1
    }
    
    #3 blank frames between events
    blanks <- 3
    for(blank in seq(1,blanks)){
      png(filename = paste(img.idx,'.png',sep=''), width = 800, height = 800, units = 'px', bg='black')
      plot(NULL,xlim=c(0,1),ylim=c(0,1))
      dev.off()
      img.idx <- img.idx + 1
    }
    
    #convert into a movie using the command line
    #command <- paste('ffmpeg -y -r 20/1 -i %d.png -qscale 0 -vcodec mpeg4 ', basename(subdir), '.mp4', sep = '')
    #system(command)
    
    #get rid of pngs
    #system('rm *.png')
    
    #move the movie to the animations folder
    #command <- paste('mv ', basename(subdir), '.mp4', ' ../', basename(subdir), '.mp4', sep='')
    #system(command)
    
    #remove the subdirectory that held the pngs
    #system(paste('cd',img.outdir))
    #system(paste('rm -r', subdir, sep = ' '))
    
  }
  
  #convert into movies using command line
}