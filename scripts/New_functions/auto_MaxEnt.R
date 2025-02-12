##~~~~~~~~~~~~~~~~~~~~~##
#                       #
# Auto_MaxEnt function ##
#                       #
##~~~~~~~~~~~~~~~~~~~~~##

Auto_maxent<-function(
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Data preparation and variable selection
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  presence_dat, # [sf/data.frame] sf object containing the points or data.frame with coordinates to the presence data
  predictors, # [RAST] rasterstack object with the environmental or predictor variables to use (In all cases one object is needed). if `time_macth = TRUE` and `ignore_predictor = TRUE`, this spat_raster is used to define the study area resolution and the parameters for the time-series rasters
  coords.p = c("decimalLongitude","decimalLatitude"), # [CHARACTER] if the spp_points is not a SpatialPoints object, provide the coordinates fields e.g x,y for ~x+y
  min_obs = 20, # [NUMERICAL] Minimun number of observations to run the analysis
  rm.dp = TRUE, # [LOGICAL] Should duplicated coordinates/points be removed from the presence data
  name.mod = presence_dat$species %>% unique(), # [CHARACTER] Name to use to store the model, by default it takes the `species` from a gbif like data set
  sp_range = NULL, # [sf object] Polygon defining the range of the species, if no range is provided `NULL` a minimum complex polygon is created from the set of presence points
  crs.r = "EPSG:4326", # [CHARACTER] CRS for the spatial information `WGS84` is used by default
  buff_lim = 0, # [NUMERIC] Buffer to apply to the study area (bounding box around the sp_range). This would depend on the units used for the CRS (degrees, meters, km, etc...)
  n_bk = "AUTO", # [CHARACTER] Should the number of background points be established automatically
  prop_env = NULL, # [NUMERIC] 0-1 If n_bk==AUTO, the proportion of environmental area to sample
  type_bk = "Random", #[CHARACTER][Random,BwData,BwData_inv,EnvBK] Type of background sampling: Random = points are drawn at random from the study area; BwData = points are sampled from the study area based on a density kernell of the presence data; BwData_inv = points are sampled based on the inverse density kernell of the presence data; EnvBK = points are sampled based on a density kernel build using the first two Principal Components of a PCA build with the environmental data 
  Test_n = 30, # [NUMERIC] values between 0-100, percentage of data separated for testing purposes
  world_pol = NULL, # [sf] Sf polygon used to refine the study area. The main intention is to remove areas that might fall into the Ocean or to precisely delimit/remove other areas.
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Some data coordination and selection
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  time_macth = FALSE, # [LOGICAL] should we try to coordinate the variables time with the records dates? Still in the works, it would need more arguments in order to format all the dates/times formats
  select_var = "NUMERICAL", # [CHARACTER] If `NUMERICAL` variables are selected based on VIF scores and correlation values (VIF<5 and COR<0.7). If `DIST` is used, the variables are selected based on the overlap of the density distributions of presence and absences [not implemented]
  d.range=NULL, # [NUMERIC] A vector of dates (years at the moment) to coordinate the point and variable data
  d.file=".tif$", # [CHARACTER] Format of the layers or variables to retrieve from the folders
  d.route, # "[CHARACTER] Route to the time-series environmental information. It can be more than one folder with multiple variables or time series. Apart from these at least a spat raster with the reference study area need to be provided in the predictors field
  ignore_predictor=FALSE, # [LOGICAL] Do we want to exclude the `predictor` data from the model or include its information (When we have temporal and fix variables, e.g Elevation)
  
  #~~~~~~~~~~~~~~~~~~#
  # MaxEnt Modelling
  #~~~~~~~~~~~~~~~~~~#
  random_features = FALSE, # [LOGICAL] Should MaxEnt select its adjustment features at random? When FALSE, model fit features are clipped and all possible combinations are tested in independent models
  beta.val = 7, # [NUMERICAL] Which value or collection of b-multiplier values should we use for the models
  n.m = 10, # [NUMERIC] Number of models to run, if `random_features` set to FALSE this gets overwritten and all the possible combinations of features and beta-multiplier values are tested
  Mod.route = NULL, # [CHARACTER] route to store the models, if no route is provided the temporal folder will be used                      
  
  #~~~~~~~~~~~~~~~~~#
  # Model Selection
  #~~~~~~~~~~~~~~~~~#
  mod.select = TRUE, # [LOGICAL] Should we run a model selection procedure
  n.mods = 3, # [NUMERIC] if mod.select is TRUE, the number of models to be selected/returned
  use.boyce = NULL, # [NUMERIC] Can take values between -1 and 1. If not NULL the boice index is used to narrow the model selection according to the numeric threshold we select. It creates a selection that can pottentially return no models or a number of models < n.mods
  
  #~~~~~~~~~~~~~~~~~#
  # Return options
  #~~~~~~~~~~~~~~~~~#
  return.all = TRUE # [LOGICAL] All parameters and objects needed to replicate the models are returned alongside the models
){
  
  # 0. Load the needed packages and functions----
  list.of.packages<-c("dplyr","sf","terra","tidyr","ENMeval","parallel","dismo","data.table")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Preset some parameters
  # a. Number of models and model parameters
  if(mod.select==TRUE){
    if(n.m<=n.mods & random_features == T){
      warning("Number of models to select equal or greather than models to run, check n.m and n.mods paramters! Setting n.m to n.mods")
      n.m <- n.mods
      
      # mod.select=FALSE
    }
    
    if(time_macth == TRUE){
      # Check the requirements for the time-matching selection:
      # Environmental information with the right format
      
      # All the needed parameters for the spatio-temporal information extraction
      
      
    }
  }
  
  # 1. Format the presence data----
  if(nrow(presence_dat) < min_obs){
    print(paste("Number of points below threshold for",name.mod,sep=" "))
    stop()
    
  }else{
    if("sf" %in% class(presence_dat)){
      y_points <- presence_dat
      y_points <- y_points %>% st_transform(crs=crs.r)
      
    }else{
      print(paste("Preparing spatial information for",name.mod,sep=" "))
      y_points <- presence_dat %>% st_as_sf(coords=coords.p,crs=crs.r)
    }
  }
  
  # 1.b Remove duplicates----
  if(rm.dp){
    y_points <- y_points[!duplicated(y_points %>% st_geometry()),]
    
    if(nrow(y_points) < min_obs){
      stop(paste("Number of points below threshold for",name.mod,sep=" "))
      
    }
  }
  
  # 1.c Prepare study area----
  if(is.null(sp_range)){
    # Create mcp from points
    p.index <- chull(y_points %>% st_coordinates())
    xy.hull <- st_coordinates(y_points)[c(p.index,p.index[1]),]
    
    sp.pol <- st_sf(data.frame(ID=1,geom=st_sfc(st_polygon(list(xy.hull)))),crs=crs.r)
    sp.pol <- sp.pol %>% st_cast("POLYGON") %>% st_geometry()
    
    # Get the study area (bounding box)
    sty.a <- sp.pol %>% st_bbox() %>% poly_from_ext(crs_p = NULL)
    st_crs(sty.a) <- crs.r ; sty.a <- sty.a %>% st_transform(crs.r)
    
  }else{
    # Check intersection between points and species range/study area polygon
    sp.pol <- sp_range
    sp.pol <- sp.pol %>% st_cast("POLYGON") %>% st_geometry()
    
    xy<-sp.pol %>% st_intersects(x=y_points,sparse = TRUE) %>% unlist() %>% unique()
    
    sp.pol <- sp.pol[xy]
    
    # Get the study area (bounding box)
    sty.a <- sp.pol %>% st_bbox() %>% poly_from_ext(crs_p = crs.r)
    #st_crs(sty.a) <- crs.r ; sty.a <- sty.a %>% st_transform(crs.r)
    
  }
  
  # 1.d Add a buffer around the study area ----
  sf::sf_use_s2(FALSE)
  sty.a <- sty.a %>% st_buffer(dist=buff_lim,endCapStyle = "FLAT")
  
  # 1.d.2 Refine the study area by removing the portions that fall into the ocean
  if(!is.null(world_pol)){
    world_pol <- world_pol %>% st_transform(crs.r)
    sty.a <- sty.a %>% st_intersection(world_pol)
    
  }
  
  # Display the information ----
  sf::sf_use_s2(TRUE)
  # plot(sty.a %>% st_geometry())  
  
  # 1.e Crop the environmental information ----
  if("SpatRaster" %in% class(predictors)){
    pred.dat <- predictors %>% crop(sty.a %>% vect(),mask=TRUE)
    gc():gc()
  }else{
    print("WARNING: predictors data is not an SpatRaster object!")
    return(NULL)
    
  }
  
  # 2. Generate the background data and extract the environmental information ----
  # Get the number of background points
  if(n_bk=="AUTO"){
    # With the AUTO parameter we are going to sample the study area until we find a number of background points that represents the
    # distribution of variables of the study area (still optimizing this process!)
    
    if(is.null(prop_env)){
      warning("No prop_env provided for the selection of background sample size. Setting the sample size to 10000, the MaxEnt default")
      n_bk <- 10000 # use the MaxEnt Default method
    }
    
    if(!is.null(prop_env)){
      # Get the XX proportion of cells from the study area
      n_bk <- prod(dim(pred.dat)[1:2])
      n_bk <- (n_bk * prop_env) %>% round(digits=0)
    }
    
  }else{
    if(!is.numeric(n_bk)){
      print("WARINING: bk_n field not_valid, setting parameter to `AUTO`")
      # return(NULL)
      
      n_bk <- prod(dim(pred.dat)[1:2])
      n_bk <- (n_bk * prop_env) %>% round(digits=0)
      
    }
  }
  
  # 2.1 No time coordination needed ----
  if(time_macth == FALSE){
    # 2.1.a Create the background points ----
    if(type_bk %in% c("Random","BwData","BwData_inv")){
      bk_points <- backgroundPOINTS(presence = y_points,
                                    background_n = n_bk,
                                    TrainTest = 1,
                                    range_samp = sty.a,
                                    weights.p = type_bk)$Train
    }
    
    if(type_bk %in% "EnvBK"){
      # The selection of variables can have an impact on the shape and distribution of the PCA scores 
      if(select_var!=FALSE | time_macth == TRUE){
        Warning("When Select_var is set to `TRUE` Background points are firstly sampled at random and this data is used to run the variable selection. After this, the selected variables are used to produce the final set of Bk points.")
        bk_points <- backgroundPOINTS(presence = y_points,
                                      background_n = n_bk,
                                      TrainTest = 1,
                                      range_samp = sty.a,
                                      weights.p = "Random")$Train
        
      }else{
        bk_points <-  bk_env(p.points=y_points,
                             env_var=pred.dat,
                             n_bk=n_bk,
                             density_pc=TRUE)[["points"]] %>% st_cast("POINT")
      }
    }
    gc(); gc()
    
    # 2.1.b Prepare the presence absence vector----
    bk_points <- bk_points %>% mutate(presence=0,.before=0)
    y_points <- y_points %>% mutate(presence=1,.before=0)
    
    obs_sp <- rbind(y_points %>% dplyr::select("presence"),
                    bk_points %>% dplyr::select("presence"))
    
    obs_index <- obs_sp$presence %>% st_drop_geometry()
    
    # 2.1.c Extract the values from the environmental variables ----
    mod.dat <- pred.dat %>% terra::extract(obs_sp %>% vect(),ID=F)
    
    # 2.1.d Remove empty cases from the data----
    dat.index <- complete.cases(mod.dat)
    obs_index <- obs_index[dat.index]
    
    mod.dat <- mod.dat[dat.index,]
    }
  
  # 2.2 Time coordination/matched ----
  if(time_macth == TRUE){
    raw.time <- Time_matchine(x = pred.dat, y = y_points, y.t = "eventDate", 
                  time.i = "year", jump.i = 1, id.p = "key",continuous_data = "9999-01-01",
                  
                  # Sampling the temporal data accordingly?
                  bk_sampling = TRUE, sampling_type = type_bk, sampling_area = sty.a, 
                  number_points = n_bk # [[Numerica]] Number of ramdom points to sample at each time jump. This will also be the final number of Bk points returned by the function (a random sample of n bk_points will be extracted)
                             )
  
  # Select the time-period to run predictions
    obs_time <- xtabs(~raw.time$time-period)
    period<-
    mod.dat
    
    
  }
  # 3.b Variable selection ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Other methods are jet to be implemented :S
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if(select_var!=FALSE){
    if(select_var=="NUMERICAL"){
      vars_s <- var_select(x=mod.dat,VIF.threshold=5,cor.threshold=0.7)
    }
    
    if(select_var=="DIST"){
      print("NOT YET IMPLEMENTED")
      vars_s<-names(mod.dat)
    }   
    
    if(length(vars_s)==0){
      stop("WARNING! No variables selected for the analysis. Consider transform the variables or set select_var to FALSE")
    }
    
    if(length(vars_s)<3){
      print(paste("WARNING! Only",length(vars_s),"selected for the analysis!"))
    }else{
      print(paste(length(vars_s), "variables selected for the analysis"))
    } 
    
    mod.dat <- mod.dat[,colnames(mod.dat) %in% vars_s]
    
    if(type_bk %in% "EnvBK" & time_macth == FALSE){
      
      bk_points <-  bk_env(p.points=y_points,
                           env_var=pred.dat[[names(pred.dat) %in% vars_s]],
                           n_bk=n_bk,
                           density_pc=TRUE)[["points"]] %>% st_cast("POINT")
      
      # Prepare the presence absence vector
      bk_points <- bk_points %>% mutate(presence=0,.before=0)
      y_points <- y_points %>% mutate(presence=1,.before=0)
      
      obs_sp <- rbind(y_points %>% dplyr::select("presence"),
                      bk_points %>% dplyr::select("presence"))
      
      obs_index <- obs_sp$presence %>% st_drop_geometry()
      mod.dat <- pred.dat[[names(pred.dat) %in% vars_s]] %>% terra::extract(obs_sp %>% vect(),ID=F)
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 4. MaxEnt - Model adjustment and selection ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  models <- list()
  mod.ARGS <- list()
  
  # Set the different combinations of model adjustment paramters----
  if(random_features==FALSE){
    supress_features <- list()
    
    for(i in c(1:3)){
      supress_features[[i]]<-combinations(x=c("nolinear","noquadratic","noproduct","nothreshold"),
                                          choose = i) %>% as.data.frame()
    }
    
    supress_features[[1]]<-supress_features[[1]] %>% t() %>% as.data.frame()
    supress_features <- supress_features %>% rbindlist(fill=TRUE)
    
    supress_features <- supress_features[rep(seq_len(nrow(supress_features)), each = if(is.null(beta.val)) 15 else length(beta.val)),]
    beta.v<-rep(if(is.null(beta.val)) c(1:15) else beta.val,times=if(is.null(beta.val)) nrow(supress_features)/15 else nrow(supress_features)/length(beta.val))
    
    supress_features$beta.v<-beta.v
    n.m<-nrow(supress_features)
  }
  gc()
  
  # 4.1 Run the different models----  
  set.seed(185)
  for(w in 1:n.m){
    
    if(random_features==TRUE){ 
      supress_features<-c("nolinear","noquadratic","noproduct","nothreshold")
      xFeature <- sample(supress_features,size=sample(1:3,1),replace=F) # Select the number of predictive distributions to fit (this is random)
      
    }else{
      xFeature<-supress_features[w,-c("beta.v")] %>% as.vector() %>% unlist() %>% unname()
      xFeature <- xFeature[!is.na(xFeature)]
    }
    
    if(is.null(beta.val)){
      beta.m <- ifelse(random_features==TRUE,sample(1:15,1),supress_features[w,"beta.v"] %>% as.vector() %>% unlist() %>% unname())
    }else{
      if(random_features==TRUE){
        beta.m <- sample(beta.val,1)
      }else{
        beta.m <- supress_features[w,"beta.v"] %>% as.vector() %>% unlist() %>% unname()
      }
    }
    
    argsMX=c(
      "autofeature=false", # Prevent the creation of automatic features or relationships between predictors provided to the model and the response variable
      xFeature, # Features to be included in the model
      "defaultprevalence=1.00",
      paste0("betamultiplier=",beta.m), # random beta multiplier set between 1-15, if random_features==FALSE all possible multipliers are applied to each combination of autofeatures (beta.v is created in line 519
      "pictures=false",
      # "biasfile" <- need to try this again at some point! 
      paste0("replicates=",1), # number of replicates for each model
      "writeplotdata=true",
      paste0("threads=",detectCores(logical=TRUE)-2), # NEW
      "hinge=true",
      paste0("randomtestpoints=",100-(100-Test_n)), # NEW percentage of points to run out of sample AUC values 
      "writeplotdata=true", 
      paste0('threads=',ifelse((detectCores(logical=FALSE)-2)>1,detectCores(logical=FALSE)-2,1)), # Numeric. The number of processor threads to use. 
      # Matching this number to the number of cores on your
      # computer speeds up some operations, especially variable jackknifing.
      "responsecurves=false",
      # "jackknife=false", #Logical if TRUE measures the importance of each environmental variable by training 
      # with each environmental variable first omitted, then used in isolation
      "askoverwrite=false"
      # Check the parameters to write the data into a designated folder
    )
    
    # 4.1.b Run the model ----
    if(is.null(Mod.route)){
      rx <- paste(tempdir(),paste(name.mod,w,sep="_"),sep="\\")
      
    }else{
      rx <- Mod.route 
    }
    
    rx %>% dir.create(showWarnings = FALSE,recursive = TRUE)
    print(xFeature)
    
    # MaxEnt model
    Train.m <- dismo::maxent(x=mod.dat, p=obs_index,
                             removeDuplicates=rm.dp, 
                             args=argsMX, path=rx)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 5. Evaluate the model performance:----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #
    # Get the test and train data used by maxent----
    dat.maxent <- read.csv(paste(rx,"species_samplePredictions.csv",sep="/"))[,c("X","Test.or.train")]
    
    Train_index <- dat.maxent %>% filter(Test.or.train == "train") %>% dplyr::select("X") %>% unlist() %>% unname()
    Test_index <- dat.maxent %>% filter(Test.or.train == "test") %>% dplyr::select("X") %>% unlist() %>% unname()
    
    # Get the AUC values calculated by maxent and the data for the plots----
    test.AUC.MAXENT <- read.csv(paste(rx,"maxentResults.csv",sep="/"),header = TRUE)$Test.AUC
    training.AUC.MAXENT <- read.csv(paste(rx,"maxentResults.csv",sep="/"),header = TRUE)$Training.AUC
    
    diffAUC.MAXENT <- training.AUC.MAXENT - test.AUC.MAXENT # check model overfitting
    
    # Test the Train model against the test data----
    y.acc <- performance_model(mod.x=Train.m, dat=mod.dat[-Train_index,],
                               index_vals=obs_index[-Train_index],
                               threshold_seq = seq(0.01,0.99,by=0.01))  
    
    TSS.threshold.TEST <- seq(0.01,0.99,by=0.01)[which.max(y.acc[rownames(y.acc) %in% "TSS",])]
    
    TSS.max.TEST <- y.acc[rownames(y.acc) %in% "TSS",] %>% max(na.rm=TRUE)
    TSS.mean.TEST <- y.acc[rownames(y.acc) %in% "TSS",] %>% mean(na.rm=TRUE)
    
    if(length(TSS.threshold.TEST)==0){TSS.threshold.TEST <- "no.data"}
    
    # 5.2 Get model thresholds----
    # MaxEnt
    TSS.threshold.FULL <- read.csv(paste(rx,"maxentResults.csv",sep="/"),header = TRUE)$Maximum.training.sensitivity.plus.specificity.Cloglog.threshold
    
    # Check model performance----
    cbi <- boyce_index(pred = predict(Train.m,pred.dat),obs=y_points[Test_index,] %>% st_coordinates())
    gc() ; gc()
    
    # 5.3 Built the table to produce the AIC values----
    if(mod.select==TRUE){
      
      if(w==1){
        xraw <- predict(Train.m,mod.dat[-c(Train_index,Test_index),],args="raw")
        
        # Standarize the predictions if the sum of the background raw predictions is not close to 1
        
        if(sum(xraw)!=1){
          raw_preds <- predict(Train.m,mod.dat[c(Train_index),],args="raw")/sum(xraw)
        }
        
        # Get the number of non-cero parameters (lambda) for each model
        # Read lambdas file
        rf <- read.table(file.path(rx, 'species.lambdas'), sep=',', fill=TRUE)
        
        # record no. of params (restrict df to rows with four values and no 0 in 2nd column)
        ncoefs.maxent <- nrow(rf[!is.na(rf[3]) & rf[2] != 0,])  
        
        # Are AIC values valid?
        n.occs <- length(raw_preds) # Number of occurrence points
        AIC.valid <- ncoefs.maxent < n.occs # Number of non.cero coefficients < number of occurence points
        
        # calculate log likelihood
        LL <- sum(log(raw_preds), na.rm = TRUE)
        AICc <- (2 * ncoefs.maxent - 2 * LL) + (2 * (ncoefs.maxent) * (ncoefs.maxent + 1) / (n.occs - ncoefs.maxent - 1))
        
        # if determined invalid or if infinite, make AICc NA
        AICc <- ifelse(AIC.valid == FALSE | is.infinite(AICc), NA, AICc)
        
        # Create the output table
        AICc_mods <- data.frame(mod=ifelse(is.null(name.mod),paste0("mod.",w),paste0(name.mod,"_",w)),
                                ncoefs=ncoefs.maxent,
                                logOdds=LL,
                                AICc=AICc,
                                valid=AIC.valid)
        
        rm(xraw,rf,ncoefs.maxent,n.occs,AIC.valid,AICc)
        
      }else{
        xraw <- predict(Train.m,mod.dat[-c(Train_index,Test_index),],args="raw")
        
        # Standarize the predictions if the sum of the background raw predictinos is not close to 1
        
        if(sum(xraw)!=1){
          raw_preds <- predict(Train.m,mod.dat[c(Train_index),],args="raw")/sum(xraw)
        }
        
        # Get the number of non-cero parameters (lambda) for each model
        # Read lambdas file
        rf <- read.table(file.path(rx, 'species.lambdas'), sep=',', fill=TRUE)
        # record no. of params (restrict df to rows with four values and no 0 in 2nd column)
        ncoefs.maxent <- nrow(rf[!is.na(rf[3]) & rf[2] != 0,])  
        
        # Are AIC values valid?
        n.occs <- length(raw_preds) # Number of occurrence points
        AIC.valid <- ncoefs.maxent < n.occs # Number of non.cero coefficients < number of occurence points
        
        # calculate log likelihood
        LL <- sum(log(raw_preds), na.rm = TRUE)
        AICc <- (2 * ncoefs.maxent - 2 * LL) + (2 * (ncoefs.maxent) * (ncoefs.maxent + 1) / (n.occs - ncoefs.maxent - 1))
        
        # if determined invalid or if infinite, make AICc NA
        AICc <- ifelse(AIC.valid == FALSE | is.infinite(AICc), NA, AICc)
        
        # Create the output table
        
        AICc_mods <- rbind(AICc_mods,data.frame(mod=ifelse(is.null(name.mod),paste0("mod.",w),paste0(name.mod,"_",w)),
                                                ncoefs=ncoefs.maxent,
                                                logOdds=LL,
                                                AICc=AICc,
                                                valid=AIC.valid))
        
        rm(xraw,rf,ncoefs.maxent,n.occs,AIC.valid,AICc)
      }}
    
    #
    # 5.4 Save model parameters----
    y.data <- data.frame(mod.id = w,
                         model.id = ifelse(is.null(name.mod),paste0("mod.",w),paste0(name.mod,"_",w)),
                         n_presence = sum(Train_index %>% length(),Test_index %>% length()),
                         n_background = bk_points %>% nrow(),
                         bk_type=type_bk,
                         Test_prop = Test_n,
                         train.AUC.MAXENT = training.AUC.MAXENT,
                         test.AUC.MAXENT = test.AUC.MAXENT,
                         diff.AUC = diffAUC.MAXENT, 
                         TSS.threshold.MAXENT = TSS.threshold.FULL,
                         TSS.threshold.TEST = TSS.threshold.TEST,
                         TSS.mean.TEST = TSS.mean.TEST,
                         TSS.max.TEST = TSS.max.TEST,
                         C.Boyce = cbi$Boyce)
    
    # Create the data.frame
    if(w==1){
      mod.params <- y.data
      
    }else{
      mod.params <- rbind(y.data,mod.params)  
    }
    
    # Export the model
    models[[w]] <- Train.m
    names(models)[w] <- ifelse(is.null(name.mod),paste0("mod.",w),paste0(name.mod,"_",w))
    
    mod.ARGS[[w]] <- argsMX
    names(mod.ARGS)[w] <- ifelse(is.null(name.mod),paste0("mod.",w),paste0(name.mod,"_",w))
    gc() ; gc()
    
    # remove the temporal files
    unlink(rx,recursive = TRUE,force=TRUE)
  }
  
  # 6. Model selection ----
  if(mod.select==TRUE){  
    
    # Select the n best performing models to run and average the predictions
    # 6.1 Build a null model using a random sampling of the environment----
    rx <- paste(tempdir(),paste(name.mod,"NULL_model",sep="_"),sep="\\")
    
    rx %>% dir.create(showWarnings = FALSE,recursive = TRUE)
    obs_index_null<-rep(0,times=nrow(mod.dat))
    obs_index_null[sample(1:nrow(mod.dat),size=sum(obs_index==1))]<-1
    
    null.maxent<-dismo::maxent(x=mod.dat, p=obs_index_null,
                               removeDuplicates=rm.dp, 
                               # args=argsMX, # arguments set to automatic 
                               path=rx)
    
    # 6.1.2 Get the AICc of the random (null) model ----
    xraw <- predict(null.maxent,mod.dat[-c(obs_index_null==1),],args="raw")
    
    # Standarize the predictions if the sum of the background raw predictinos is not close to 1
    if(sum(xraw)!=1){
      raw_preds <- predict(null.maxent,mod.dat[c(obs_index_null==1),],args="raw")/sum(xraw)
    }
    
    # Read lambdas file
    rf <- read.table(file.path(rx, 'species.lambdas'), sep=',', fill=TRUE)
    # record no. of params (restrict df to rows with four values and no 0 in 2nd column)
    
    ncoefs.maxent <- nrow(rf[!is.na(rf[3]) & rf[2] != 0,])  
    
    # Are AIC values valid?
    n.occs <- length(raw_preds) # Number of occurrence points
    AIC.valid <- ncoefs.maxent < n.occs # Number of non.cero coefficients < number of occurence points
    
    # calculate log likelihood
    LL <- sum(log(raw_preds), na.rm = TRUE)
    AICc <- (2 * ncoefs.maxent - 2 * LL) + (2 * (ncoefs.maxent) * (ncoefs.maxent + 1) / (n.occs - ncoefs.maxent - 1))
    
    # if determined invalid or if infinite, make AICc NA
    AICc <- ifelse(AIC.valid == FALSE | is.infinite(AICc), NA, AICc)
    
    # Create the output table
    AICc_mods <- rbind(AICc_mods,data.frame(mod="NULL_mod",
                                            ncoefs=ncoefs.maxent,
                                            logOdds=LL,
                                            AICc=AICc,
                                            valid=AIC.valid))
    
    # 6.1.3 Remove invalid AICc, duplicated models, and cero coefficients models ----  
    AICc_mods_raw <- AICc_mods
    AICc_mods <- AICc_mods %>% filter(ncoefs!=0,valid==TRUE)
    AICc_mods <- AICc_mods[!paste(AICc_mods$AICc,AICc_mods$logOdds) %>% duplicated(),]
    
    # 6.2 Calculate delta-AIC for the different models ----
    AICc_mods$delta.AICc <- AICc_mods$AICc - min(AICc_mods$AICc, na.rm=TRUE)
    AICc_mods$w.AIC <- exp(-0.5 * AICc_mods$delta.AICc) / sum(exp(-0.5 * AICc_mods$delta.AICc), na.rm=TRUE)
    AICc_mods <- AICc_mods[order(AICc_mods$AICc),]
    
    # 6.3 Select the models based on its AICc and boice index
    if(n.mods > nrow(AICc_mods) | n.mods==nrow(AICc_mods)){
      paste0("WARNING: Number of models (",n.mods,") to select, greater than number of tested models(",nrow(AICc_mods),")! Setting n.mods to to number of tested models -1")
      n.mods <- nrow(AICc_mods)-1
    }
    
    select_models <-  AICc_mods[1:n.mods,"mod"]      
    
    if("NULL_mod" %in% select_models){
      print("WARNING: The null models is amongst the selected models! Dropping this model, number of selected models -1.")
      select_models <- select_models[!select_models %in% "NULL_mod"]
    }
    
    # 6.4 Should we consider the Boyce index of the models for the selection?
    if(!is.null(use.boyce)){
      select_models <- mod.params %>% filter(model.id %in% select_models,C.Boyce >= use.boyce) %>% 
        dplyr::select("model.id") %>% unlist() %>% unname()
    }
    
    if(length(select_models)==0){
      return(list(AICc_mods_raw,mod.params))
      stop("No models selected! Consider changing the selection parameters")
    }
    
  }else{
    select_models<-names(models)
  }
  
  # Reduce the data to macth the model selection
  models <- models[names(models) %in% select_models]
  mod.params <- mod.params %>% filter(model.id %in% select_models) 
  
  if(mod.select==TRUE){
    AICc_mods <- AICc_mods %>% filter(mod %in% select_models)
  }else{
    AICc_mods <- NULL
  }
  
  mod.ARGS <- mod.ARGS[names(mod.ARGS) %in% select_models]
  
  # 7. Prepare the return object ----
  if(return.all==TRUE){
    # a. Run model predictions
    preds<-lapply(models,predict,pred.dat) %>% rast()
    
    # b. Combine model predictions
    comb.pred <- mean(preds,na.rm=TRUE)
    
    xp<-list(mods=models,
             params=mod.params,
             AICc=AICc_mods,
             ARGS=mod.ARGS,
             variables=vars_s,
             study.area=sty.a,
             bk.oints=bk_points,
             mod.preds=preds,
             avr.preds=comb.pred)
    
    return(xp)
  }else{
    return(models)
  }
}

# End of the script
