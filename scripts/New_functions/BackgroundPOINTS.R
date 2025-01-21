#
###//\/\/\/\/\/\/\////\/\/\/\/\/\/\//\/\/\/\//\////\/\/\///////////////////\\\\\\##-#
#                 Create the background points for MAXENT algorithm
###///\/\/\/\/\/\/\////\/\/\/\/\/\/\\/\////\/\/\/\/\/\/\/\/\/\\\\\/\/\/\/\/\/\///##-#
#
#
backgroundPOINTS<-function(presence, # [sf points] distribution points of the species
                           background_n, # [Numeric] number of background points
                           TrainTest, # [Numeric] proportion of data used for the testing and training of the models
                           range_samp=NULL, # [Character] Route to the range information of species
                           buffer.dist=0, # [Numeric] buffer to apply in order to delimit the area for the random sampling of backgorund data 
                           # cut_area=range_samp, # A route to a simplified polygon of the continents shapes to cut limit the buffers
                           bias.sampling=NULL, # If a sampling bias layer is provided, this is used to run the backgorund psudoabsence sampling, otherwise the method selected in weights is applied
                           weights.p="Random"  # [Character]["Random","BwData","BwData_inv"] Method to draw the background points: BwData, background weighted data, a density kernell is first calculated and 
                                               # background data is created based on the probabilities of this kernell; Random, points are randomly generated within 
                                                # the limits of the specified range (default)
                           ){
   # 0. Load the needed packages
   pack<-c("dplyr","MASS","sf","classInt","terra","raptr")
   
      new.packages <- pack[!(pack %in% installed.packages()[,"Package"])]
      if(length(pack)) install.packages(new.packages)
      
      lapply(pack,require,character.only=TRUE)
      rm(pack,new.packages)
   
   # a.1 Create a buffer around the range or the MCP data, clip the polygon using the wrld_simple polygon
   if(!is.null(range_samp)){
      # Using the species range to delimit the background sampling area
      if(is.character(range_samp)){
             range_sp <- st_read(range_samp)
      }else{
               range_sp <-range_samp
            }
      
      range_sp<-st_transform(range_sp,crs(presence)) %>% st_make_valid()#sp::spTransform(SPP.range,crs.r)
      
      # Check the validity of range_sp, if is not valid create an MCP witht the points
      geo.check<-st_is_valid(range_sp %>% st_make_valid())
      
      if(FALSE %in% geo.check){
         range_sp<-range_sp[geo.check,]
      }
      
      if(length(range_sp)==0){
         r.buff<-st_convex_hull(st_union(presence)) %>% st_buffer(dist=buffer.dist)
         # Print a warning!
         print("Warning: The geometry of range_sp is not valid. Creating a convex_hull instead!")
      }
      
      r.buff<-st_buffer(range_sp,dist=buffer.dist) %>% st_union(by_feature = FALSE)
      
      }else{
      # Create an area using the MCP form by the distribution points
         r.buff<-st_convex_hull(st_union(presence)) %>% st_buffer(dist=buffer.dist)
      }
      
      # if(!is.null(cut_area)){
      #    # Clip the sampling area using the cut_area shape
      #    ifelse(is.character(cut_area),
      #           cut_area <- st_read(cut_area),
      #           cut_area <- cut_area)
      #    
      #    cut_area <- cut_area %>% st_transform(crs(r.buff))
      #    r.buff<-st_intersection(r.buff,cut_area) # %>% st_geometry() %>% plot()
      #    }
      
    # Creating the background train and test data
   if(is.null(bias.sampling)){   
   
   if(weights.p %in% c("BwData","BwData_inv")){
      # Create a sampling surface based on observations
      yy <- range_sp
      
     # Generate the density kernelTRUE# Generate the density kernell
      a.r<-d.fromXY(x=presence, # Points used to build the density kernell
                     y.size=c(1000,1000), # size of the density kernell matrix
                     inv=ifelse(weights.p=="BwData_inv",TRUE,FALSE), # should the final matrix need to be inverse?
                     scale01=TRUE, # should the final matrix need to be scaled between 0_1?
                     mask.p=yy, # Does the initial kernell need to be masked?
                     pol.ref=r.buff # Add a reference to create the final raster (this is the polygon used for the sampling)
                    )
      # plot(a.r)
      # plot(y_points %>% st_geometry(),add=TRUE)
      
      # If we have very few points where psudo-absences can be placed, add 0.0001 to the raster map
      if(sum(values(a.r) > 0.05,na.rm=TRUE)<background_n){
         a.r[a.r<0.05] <- 0.05
      }
      
      # Generating 10000 background points randomly considering the probability distribution of the density raster
      TrainBk <- raptr::randomPoints(mask=a.r,n=round(background_n*0.8,digits=0),prob=T)
      TtestBk <- raptr::randomPoints(mask=a.r,n=round(background_n*0.2,digits=0),prob=T)
      
      }
   
      if(weights.p=="Random"){
         # Creating the testing and training background
         TrainBk <- sf::st_sample(r.buff %>% st_make_valid(), size=round(background_n*TrainTest,digits=0),type="random") %>% as.data.frame() %>% st_as_sf(crs=st_crs(presence))
         TtestBk <- sf::st_sample(r.buff %>% st_make_valid(), size=round(background_n*c(1-TrainTest),digits=0),type="random") %>% as.data.frame() %>% st_as_sf(crs=st_crs(presence))
      }
   }else{
      
      # If we have very few points where psudo-absences can be placed, add 0.0001 to the raster map
      if(sum(values(bias.sampling) > 0,na.rm=TRUE)<background_n){
         bias.sampling <- bias.sampling %>% subst(from=NA,to=1) %>% mask(r.buff %>% vect())# %>% plot()
         }
      
      # Generating 10000 background points randomly considering the probability distribution of the density raster
      TrainBk <- raptr::randomPoints(mask=bias.sampling,n=round(background_n*0.8,digits=0),prob=T)
      TtestBk <- raptr::randomPoints(mask=bias.sampling,n=round(background_n*0.2,digits=0),prob=T)
      
   }
      
      # Correct compatibility problem
   if(is.matrix(TrainBk)){
      TrainBk <- TrainBk %>% as.data.frame() %>% st_as_sf(coords=c("x","y"),crs=st_crs(presence))
      TtestBk <- TtestBk %>% as.data.frame() %>% st_as_sf(coords=c("x","y"),crs=st_crs(presence)) 
      }
   
      return(list(Train=TrainBk,Test=TtestBk,Samp_area=r.buff)) 
  }
 

###//\/\/\/\/\/\/\////\/\/\/\/\/\/\//\/\/\/\//\////\/\/\///////////////////\\\\\\##-#
#     Accesory functions   /\/\/\/\/\/\/\/\/\\----\-\-\-\-\-\
###///\/\/\/\/\/\/\////\/\/\/\/\/\/\\/\////\/\/\/\/\/\/\/\/\/\\\\\/\/\/\/\/\/\///##-#
#
# 1.Calculate a density kernell based on geographical coordiantes
# Building a kernel density raster of the presence points used to create the background weighting data:

d.fromXY<-function(x, # Points used to build the density kernell
                   y.size=c(1000,1000), # size of the density kernell matrix
                   inv=FALSE, # should the final matrix need to be inverse?
                   scale01=TRUE, # should the final matrix need to be scaled between 0_1?
                   pol.ref=NULL, # reference polygon to create the extension of the kernell raster
                   mask.p=NULL # Does the initial kernell need to be masked?
){
   # 0. load the needed libraries----
   list.of.packages<-c("dplyr","terra","MASS","sf","KernSmooth")
   
   new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
   if(length(new.packages)) install.packages(new.packages)
   
   lapply(list.of.packages,require,character.only=TRUE)
   rm(list.of.packages,new.packages)
   
   # Run the function
   if(is.null(pol.ref)){
      coord <- x %>% st_coordinates()
      
      # Check the quantiles of our coordinates to avoid errors in kde2d
      if(bandwidth.nrd(coord[,1])==0){
         coord[,1]<-coord[,1]+rnorm(n=length(coord[,1]),mean=0.05,sd=0.0025)
         }
      
      if(bandwidth.nrd(coord[,2])==0){
         coord[,2]<-coord[,2]+rnorm(n=length(coord[,2]),mean=0.05,sd=0.0025)
         }
      
      # k = kde2d(x=coord[,"X"],y=coord[,"Y"],n=y.size)
      # k_mat = apply(t(k$z), 2, rev) # need to apply a transformation to make the indexes of kde2d compatible
      # 
      k<-KernSmooth::bkde2D(coord,bandwidth = y.size*0.01,range.x = list(range(coord[,"X"]),range(coord[,"Y"])),
                            gridsize = y.size,truncate=TRUE)
      
      k_mat <- apply(t(k$fhat), 2, rev)
      
      
      }else{
         coord <- x %>% st_coordinates()
         
         # Check the quantiles of our coordinates to avoid errors in kde2d
         if(bandwidth.nrd(coord[,1])==0){
            coord[,1]<-coord[,1]+rnorm(n=length(coord[,1]),mean=0.05,sd=0.0025)
         }
         
         if(bandwidth.nrd(coord[,2])==0){
            coord[,2]<-coord[,2]+rnorm(n=length(coord[,2]),mean=0.05,sd=0.0025)
         }
      
      # limits = pol.ref %>% st_bbox() ;  limits <- limits[c(1,3,2,4)]# as c(xl, xu, yl, yu)
         
      # k = kde2d(x=coord[,1],y=coord[,2],n=y.size,lims = c(pol.ref %>% st_bbox())[c(1,3,2,4)])# %>% as.vector())
      # k_mat = apply(t(k$z), 2, rev) # need to apply a transformation to make the indexes of kde2d compatible
      
      k<-KernSmooth::bkde2D(coord,bandwidth = y.size*0.001,
                            range.x = list(st_bbox(pol.ref)[c(1,3)],st_bbox(pol.ref)[c(2,4)]),
                            gridsize = y.size,truncate=TRUE)
      
      k_mat <- apply(t(k$fhat), 2, rev)
      
      }
   
   if(inv==TRUE){
      k_mat<-max(k_mat)-k_mat
   }
   
   # Create the reference raster
   if(is.null(pol.ref)){
      r1 = rast(k_mat,extent=c(min(k$x),max(k$x),min(k$y),max(k$y)))
         }else{
      r1 = rast(k_mat,extent=st_bbox(pol.ref),crs=crs(pol.ref %>% vect()))  
      }
   
   if(!is.null(mask.p)){
      r1 <- r1 %>% mask(mask.p %>% vect())
   }
   
   # Scale the matrix between 0-1
   scale.f<-function(x,...){
      
      min.m <- min(x,na.rm = TRUE)
      max.m <- max(x,na.rm = TRUE)
      
      rescaled.k_mat <- (x - min.m)/(max.m - min.m)
      return(rescaled.k_mat)
   }
   
   # Process the density matrix
   if(scale01==FALSE){
      print(paste0(ifelse(inv==TRUE,"Inverse ",""),"density of observations returned!"))
      r2 <- r1
   }
   
   if(scale01==TRUE){
      print(paste0(ifelse(inv==TRUE,"Inverse ",""),"density of observations, between 0-1, returned!"))
      r2 <- r1 %>% app(fun=scale.f,na.rm=TRUE)
   }
   return(r2)
   
}

  
# End of the Function
   
