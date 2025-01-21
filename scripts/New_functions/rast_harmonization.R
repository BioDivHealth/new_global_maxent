# Super Raster processing code!   
resample.rast<-function(x,# Reference raster (should be a spatrast object or matrix)
                        res=NULL, # If x is missing we can built the reference raster using the resolution (res), extent(ex) and crs (crs.r)
                        ex=NULL, # if x is a matrix, the extent of the raster must be specified, by default it takes the value of the whole globe
                        crs.r=NULL, # Coordinate reference system, default WGS84 (lat~lon) # "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                        y, # filepath to the layers we want to re-size/sample
                        paralell=FALSE, # Use multiple when there is more than one object to resize? default=FALSE
                        cores=2,      # If paralell==TRUE, how many cores do you want to commit? default==2
                        route.r= NULL, # route to save the temporal re-sampled layers, by default it would create a temporal folder in the working environment with the name of r_resample
                        results.r=NULL, # Route to save the results, by default it create a ResultsRast folder in the working environment
                        mask_r=NULL, # polygon to use as mask, only if mask=TRUE,Should the final rasters need to be masked? default FALSE if TRUE a wrld_simpl dataset is usded to mask the rasters
                        rm.temp=TRUE
                        ){
  
  # Load the required libraries to run the function
  # a. Load the packages needed----
  list.of.packages<-c("terra","doParallel","dplyr","foreach","sf")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # b. Creates the Temporal file folder and the folder to allocated the results
  if(is.null(route.r)){Tdir<-paste(getwd(),"r_resample",sep="/")}else{Tdir<-route.r}
  if(!dir.exists(Tdir)){dir.create(Tdir,recursive = TRUE)}
  
  if(is.null(results.r)){Resdir<-paste(getwd(),"ResultsRast",sep="/")}else{Resdir<-results.r}
  if(!dir.exists(Resdir)){dir.create(Resdir,recursive = TRUE)}
  
  # Tmp dir---disposable
  Tmp_dir<-paste(getwd(),"junk",sep="/")
  if(!dir.exists(Tmp_dir)){dir.create(Tmp_dir)}
    
  # b. Reference raster, loaded or create it----
  if(missing(x)|is.null(x)){
    if(TRUE %in% lapply(list(crs.r,ex,crs.r),is.null)){
      print("Missing key components to build the reference raster, sampling the layers to select the best fit!")
    
    # CRS ----
     # a. Are the CRS of the layers the same?
      if(length(lapply(y,function(x) check_layers(x,return.var="crs")) %>% unlist() %>% unique())!=1){
        
        # Select the most frequent crs
        env_crs<-lapply(y,function(x) check_layers(x,return.var="crs") %>% as.character()) %>% unlist()
        new_crs<-xtabs(~env_crs) %>% as.data.frame(); new_crs<-new_crs[new_crs$Freq==max(new_crs$Freq),"env_crs"] %>% as.character()
        
        # rerpoject the rasters to the common CRS
        lyr_index <- new_crs==env_crs
        
        lyr_project <- y[!lyr_index] %>% list()
        
        lapply(lyr_project,function(w) to_project(x=w,
                                                  crs.n=new_crs,
                                                  export="./Data/temp_files")) # route to expor the re-projected layer
        
        # substitude the new adapted layer in the env_var list
        y[!lyr_index]<-list.files("./Data/temp_files",pattern=".tif$",full.names = TRUE)
        
        crs_check <- FALSE
        
      }else{
        crs_check <- TRUE
        }
      
      # Which is the CRS of the environmental layers?
       new_crs <- lapply(y,function(x) check_layers(x,return.var="crs")) %>% unlist() %>% unique()
      
    # Extension----
    # b. Check the extend and resolution of the environmental variables
    # Which is the lower extend and resolution?
    ext_env <- lapply(y,function(x) check_layers(x,return.var="ext") %>% as.vector())# %>% unlist()
    rast_exted_names <- ext_env[[1]] %>% names()
    
    # Are all extensions de same
    if(length(unique(ext_env))<=1){
      new_ext <- unique(ext_env) %>% unlist()
      ext_check <- TRUE
    
    }else{
      # Each element of ext_env represent the bounding box of the different environmental variables, we are going to choose a extension in which all variables
      # intersects
      pol_ext <- lapply(ext_env,poly_from_ext,crs_p=new_crs) # transform the extension into sf objecs
      pol_ext <- do.call(c,pol_ext) # combine the spatial features 
      
      # Select the area in which the maximun number of env_layers overlap
      sf_ext <- st_sf(pol_ext)
      i_ext <- st_intersection(sf_ext) ; i_ext <- i_ext %>% filter(n.overlaps==max(i_ext$n.overlaps)) 
      
      new_ext <- st_bbox(i_ext) #; new_ext <- new_ext[rast_exted_names]     
      # env_vars_index <- i_ext$origins
      
      ext_check <- FALSE
    }
    
    # Select the resolution
    res_env <- lapply(y,function(x) check_layers(x,return.var="res")[1] %>% as.vector()) %>% unlist()
    
    if(length(unique(res_env))<=1){
      new_res <- unique(res_env) %>% unlist()
      res_check <- TRUE
      
      }else{
      new_res <- res_env %>% unlist() %>% max(na.rm=TRUE)
      res_check <- FALSE
    }

    # Save the new_raster parameters
    env_params <- list(crs=new_crs,extent=new_ext,resolution=new_res)
    # env_vars_index
  }else{
    env_params <- list(crs=crs.r,extent=ex,resolution=res)
    }
    
    # Create a reference raster
    ref.rast <- rast(extent=env_params$extent,resolution=env_params$resolution,crs=env_params$crs)
    
    if(all(c(res_check,ext_check,crs_check))){
    print("Layers are compatible, no need for resampling!")
    return(NA)
    stop()
    
  }else{
    if(is.matrix(x)){
      
      if(TRUE %in% lapply(list(crs.r,ex),is.null)){
        print("ex or crs.r parameters missing!")
        stop()
      }
      
      ref.rast<-rast(x,crs=crs.r,extent=ex,resolution=res)
      env_params <- list(crs=terra::crs(ref.rast),extent=terra::ext(ref.rast),resolution=terra::res(ref.rast))
     }
    
    if(class(x)=="SpatRaster"){
      ref.rast<-x
      env_params <- list(crs=terra::crs(x),extent=terra::ext(x),resolution=terra::res(x))
      }
    }
  
    # d. Resample the layers----
    # If all parameters are the same, combine layers but skip resampling 
  
      
  }else{
    # d.1 If paralell processing is set to true and the number of layers is greater than 1
    if(length(y)>1 & paralell==TRUE){
      
      if(cores>detectCores(logical = FALSE)){
        cores<-detectCores(logical = FALSE)-1   
      }
      
      # Need to write the reference rast into disk to distribute it to the workers
      if(class(x)=="SpatRaster"){
        ref.rast<-x
        env_params <- list(crs=terra::crs(x),extent=terra::ext(x),resolution=terra::res(x))
      }
      
      ref.rast <- x
      
      terra::values(ref.rast)<-NA
      terra::writeRaster(ref.rast,
                         filename=paste(Tmp_dir,paste0("ref_rast.tif"),sep="/"),
                         overwrite=TRUE)
      
      # Create the cluster structure
      cl<-makeCluster(cores)
      registerDoParallel(cl)
      
      foreach(k=1:length(y)) %dopar% {
        list.of.packages<-c("terra","dplyr")
        new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
        
        if(length(new.packages)) install.packages(new.packages)
        
        lapply(list.of.packages,require,character.only=TRUE)
        rm(list.of.packages,new.packages)
        
        ref.rast2<-paste(Tmp_dir,paste0("ref_rast.tif"),sep="/") %>% rast() # Load the reference raster
        r.rast<-rast(y[k]) # load the data to resample
        
        # d.2 Configure Terra to load the processing into the disk rather than the ram memory
        terraOptions(memfrac=0.9,todisk=TRUE,tempdir = Tmp_dir) # The object y will be allocated in the temporal folder of 
        terra::resample(x=r.rast,y=ref.rast2,filename=paste(Tdir,paste0(names(r.rast),".tif"),sep="/"),overwrite=TRUE)
        gc()
        
        print("SUCCESS!!")
      }
      stopCluster(cl)
      
    }else{
      
      if(length(y)>1){
        
        for(i in 1:length(y)){
          r.rast<-y[i] %>% rast()
          terraOptions(memfrac=0.9,todisk=TRUE,tempdir = Tmp_dir) # The object y will be allocated in the temporal folder of 
          terra::resample(x=r.rast,y=ref.rast,filename=paste(Tdir,paste0(names(r.rast),".tif"),sep="/"),overwrite=TRUE)
          print("SUCCESS!!")
        }
        gc()
        
      }else{
        
        terraOptions(memfrac=0.9,todisk=TRUE,tempdir = Tdir)
        r.rast<-rast(y)
        terra::resample(x=r.rast,y=ref.rast,filename=paste(Tdir,paste0(names(r.rast),".tif"),sep="/"),overwrite=TRUE)
        print("SUCCESS!!")
      }
    }
    
    # e. Unify all the re-sample layers into a single stack
    f.r<-list.files(Tdir,pattern = ".tif$",full.names = TRUE,recursive = FALSE) # put always the full route!!!
    resampled.r<-rast(f.r)
    
    if(!is.null(mask_r)){
      x.mask <- mask_r %>% sf::st_transform(crs(ref.rast))
      x.mask <- x.mask %>% vect() # need to transform the sf object into a spatial vector object
      resampled.r<-terra::mask(resampled.r,mask=x.mask,touches=TRUE) %>% crop(x.mask)  
    }
    
    terra::writeRaster(resampled.r,
                       filename=paste(Resdir,paste0("Resample_rast.tif"),sep="/"),
                       names=basename(f.r %>% gsub(pattern = ".tif$",replacement ="")),
                       overwrite=TRUE) # Export the full results as a stack of layers
    if(rm.temp==TRUE){
      unlink(Tdir,recursive = TRUE)
      unlink(Tmp_dir,recursive = TRUE)
      }
    
    print("ALL DONE!")
    return(list(layers=basename(f.r %>% gsub(pattern = ".tif$",replacement ="")),
                parameters=env_params)) # these are the names of the layers
    }  
  }
  
  
  # A. Accesory Functions to process the information----
  # a. Extract information from the raster layers----
  check_layers<-function(x, # [Character] route to rast object
                         return.var="NULL" # [Character] y_r (resolution),y_ext (extension),y_crs (layer CRS)
  ){
    y<-rast(x)
    y_r<-res(y)[1]
    y_ext<-terra::ext(y)
    y_crs<-crs(y)
    
    if(is.null(return.var)){
      return(list(resolution=y_r,
                  extent=y_ext,
                  CRS=y_crs))
    }else{
      if(return.var=="res"){
        return(y_r)
      }
      if(return.var=="ext"){
        return(y_ext)
      }
      if(return.var=="crs"){
        return(y_crs)
      }
    }
  }
  
  # b. Function to reproject the raster layer in case different CRS are detected
  to_project<-function(x, # route to raster
                       crs.n, # crs to project
                       export=NULL # route to expor the re-projected layer
  ){
    xy <- rast(x)  
    y <- xy %>% terra::project(crs.n)
    
    if(is.null(export)){
      return(y)
      
    }else{
      ifelse(dir.exists(export),print(paste("Projected layer exported to\n",export)),
             dir.create(export,recursive = TRUE,showWarnings = FALSE))
      
      terra::writeRaster(y,
                         filename=paste(export,basename(x),sep="/"),
                         overwrite=TRUE)
      }
    }
  
  # c. Create polygon from extend object
  poly_from_ext<-function(x,crs_p){
    x2<-x[c("xmin", "ymin", "xmax", "ymax")]
    x_double <- as.double(x2)
    names(x_double) <- names(x2)
    class(x_double) <- "bbox"
    
    x_sf <- sf::st_as_sfc(x_double)
    sf::st_crs(x_sf) <- sf::st_crs(crs_p)
    
    return(x_sf)
  }
  
# End of the script