##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Fuunction to coordinate spatio-temporal species records and environmental information  ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#

Time_machine<-function(d.route, # Route to the directories or files with the temporal information
                       # d.format, # Format for the dates
                       d.range=c(2015,2020), # Range of dates to select in years
                       # d.aggregate=NULL, # Should the data by aggregated by year/month/week
                       d.file=".tif$", # format of the files
                       # ref.rast=NULL, # Arguments to build a reference raster to harmonize all the data
                       y_points=y_points, # Gbif or presence records to syncronize
                       bk_points=bk_points, # background points
                       ref_rast # Reference raster
                       ){
  
  # 0. Load the needed packages and functions----
  list.of.packages<-c("terra","data.table")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Check the points and format the dates
  # Transform the 
  y_points$eventDate <- y_points$eventDate %>% as.POSIXct(format="%Y-%m-%d")
  
  index.points <- data.frame(
              date=y_points$eventDate,
              year= y_points$eventDate %>% strftime(format="%Y") %>% as.numeric(),
              month=y_points$eventDate %>% strftime(format="%m") %>% as.numeric(),
              week=y_points$eventDate %>% strftime(format="%V") %>% as.numeric(),
              day=y_points$eventDate %>% strftime(format="%u") %>% as.numeric()
                )

  # Prepare the raster data
                                                              
      #if(dir.exists(d.route)){
      # Working with the files within the directory  
        t.full <- lapply(d.route,function(w) list.files(w,pattern = d.file,full.names = TRUE)) %>% unlist()
        t.x <- lapply(d.route,function(w) list.files(w,pattern = d.file,full.names = F)) %>% unlist()
      
      # Sort the files for the particular time-range
        r.time <- c(d.range[1]:d.range[2])
        index.t <- unlist(lapply(r.time,function(w) grep(pattern = w,t.x)))
    
      # Get the number of background points for each time-period 
        #n.bk.t <- (n_bk/length(r.time)) %>% round(digits=0)
        n.bk.t <- seq(from=1,to=n_bk,length.out=length(r.time)+2) %>% round(digits=0)
        
      # All layers in the directory that fall into the range  
        t.x <- t.x[index.t]        
        t.full <- t.full[index.t]
        
        if(length(t.x)==0){
          stop("No environmental data for the selected time-period")
        }
        
      # For each time period-combine the information and extract the data from the points  
        date.data <- list()
        bk.data <- list()
        
        for(i in 1:length(r.time)){
          r.r <- t.full[grep(t.x,pattern = r.time[i])]
          
          if(length(r.r)==0){
            warning(paste("No environmental information for the year",r.time[i]))
            next()
          }
          
          t.rast <- resample.rast(x=ref_rast, y = r.r, route.r = paste(getwd(),"time_machine",sep="/"), 
                                  results.r = paste(getwd(),"time_machine","resampled",sep="/"), 
                                  crs.r = NULL,rm.temp=FALSE) # change later
          
          if(is.na(t.rast)[1]){
            time.rast <- rast(r.r)
          }else{
            time.rast <-  paste(getwd(),"time_machine","resampled",sep="/") %>% list.files(pattern=d.file,full.names = TRUE) %>% rast()
            lyer_params <- t.rast
             
          }
          
          # Extract the information
            p.y <- y_points[index.points$year==r.time[i],]
          
          if(nrow(p.y)==0){
            next()
          }
          
          # Extract the data
          date.data[[i]] <- time.rast %>% terra::extract(p.y %>% vect(),ID=F)
          names(date.data)[i]<-t.x[i]
          
          # Extract information for the background
          bk.data[[i]]<-time.rast %>% terra::extract(bk_points[n.bk.t[i]:n.bk.t[i+1],] %>% vect(),ID=F)
          names(bk.data)[i]<-t.x[i]
          
          # remove the temporal folder
          unlink(paste(getwd(),"time_machine",sep="/"),recursive=TRUE)
          }
        
        pres.dat <- rbindlist(date.data,fill=TRUE)
        bk.dat <- rbindlist(bk.data,fill=TRUE)
        
        return(list(pres.dat,bk.dat))
        }

# 
# 
# # Example:
# # multiple variables environmental variables for time syncronization
# route_vars <- c("F:/OneDrive - Natural History Museum/Projects/MAXENT_disease/Data/DataGeomaTys/EbolaSDMdata_20241112/bioclim",
#                 "F:/OneDrive - Natural History Museum/Projects/MAXENT_disease/Data/DataGeomaTys/EbolaSDMdata_20241112/relative_humidity_dengue",
#                 "F:/OneDrive - Natural History Museum/Projects/MAXENT_disease/Data/DataGeomaTys/EbolaSDMdata_20241112/urbanicity")
# 
# 
# constant <- "F:/OneDrive - Natural History Museum/Projects/MAXENT_disease/Data/DataGeomaTys/EbolaSDMdata_20241112/elevation/wc2.1_2.5m_elev"
# 
# Time_machine(d.route = route_vars,d.range=c(2016,2018),y_points = y_points,bk_points=bk_points)
# 







