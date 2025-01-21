
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Function to sample background data using the environmental hypervolume
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#

bk_env<-function(p.points, # Presence points for a given species
                 env_var=env_study, # raster stack of environmental variables
                 n_bk=10000, # Number of background points 
                 density_pc=FALSE, # We build the sampling based on the values of a density kernell over the study area
                 dim_k=NULL, # set the dimensions of the kernell matrix (by default it takes the dimensions of 1/2 of the envirnomental raster)
                 inv=FALSE, # if density_pc is true, should we use the inverse of the density kernell to draw bk points?
                 bandwidth.env=c(0.05), # if density_pc equals "TRUE", the values of the bandwidths at which to include environmental points
                 bandwidth.sp=c(0.75) # if density_pc equals "TRUE", the values of the bandwidths at which to include species presence points for the calculations of overlap
         #plot=F # should we plot a summary historgram of the results?
         ){


    # 0. Load the needed packages----
    list.of.packages<-c("sf","terra","MASS","KernSmooth")
    
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    
    lapply(list.of.packages,require,character.only=TRUE)
    rm(list.of.packages,new.packages)
    
    # 1. Prepare the spatial data----
    p_xy <- p.points %>% vect() # convert the points to vect and extract the cells in which they fall
    p_xy_index <- terra::extract(env_var,p_xy,cells=TRUE)$cell
    
    # 2. Build the environmental space using a PCA----
    # 2.a. Transform the environmental variables into vectors without NAs values
    env_x <- rast_to_vect(x=env_var)
    
    # 2.b. Perform a PCA over the environmental variables----
    env_y <- prcomp(env_x$tab[,-1],rank=4) 
    env_xy <- env_y$x[,c(1,2)] %>% as.matrix() ; rownames(env_xy) <- env_x$tab[,1]
    
    
    # 3. Calculate the areas of the environmental space occupied by the observations----
    if(density_pc==FALSE){
       # Create a convex hull for the species within the environmental space defined by the PCA  
        hull_complete<-chull(env_xy)
        X_sp <-env_xy[row.names(env_xy) %in% p_xy_index,]
        
        hull_sp<-chull(X_sp)
        
        coord_env<-env_xy[c(hull_complete,hull_complete[1]),]       
        coord_sp<-X_sp[c(hull_sp,hull_sp[1]),]
        
        env.pol<-st_polygon(list(coord_env))
        spp.pol<-st_polygon(list(coord_sp))
        
        ENV_pols <- st_sf(data.frame(ID=1:2,Type=c("Env","Spp"),level=NA,geom=st_sfc(env.pol,spp.pol)))
        
        # Some metrics----
        ENV_pols$Area <- st_area(ENV_pols) # areas with no dimensions
      
        # Sample the environmental space----
        if(n_bk>=nrow(env_xy)){
          print("WARNING! Number of background points greather than number of raster cells, setting this to cells/2!")
          n_bk.1<-nrow(env_xy)/2 %>% round(digits = 0)
        }else{
          n_bk.1<-n_bk
        }
        
        set.seed(99)
        Selected_cells <- sample(row.names(env_xy) %>% as.numeric(),size=n_bk.1,prob=NULL)
            
        # Transform the selection into points----
        bk_points <- data.frame(lon=env_var %>% terra::xFromCell(Selected_cells),
                                lat=env_var %>% terra::yFromCell(Selected_cells))
        
        bk_points <- bk_points %>% st_as_sf(coords=c("lon","lat")) %>% st_as_sfc(crs=crs(env_var))
        
        # check bk points cover over the environment----
        bk_xy <- env_xy[rownames(env_xy) %in% as.character(Selected_cells),]
        
        hull_bk <- chull(bk_xy)
        coord_bk <- bk_xy[c(hull_bk,hull_bk[1]),]       
        
        bk.pol <- st_sf(data.frame(ID=3,Type=c("Bk"),level=NA,geom=st_sfc(st_polygon(list(coord_bk)))))
        bk.pol$Area <- st_area(bk.pol)
        
        ENV_pols<-rbind(ENV_pols,bk.pol)
    
        }else{
        # Create a density kernell of the environmental space whit the 
          # Check the quantiles of our coordinates to avoid errors in kde2d
          if(bandwidth.nrd(env_xy[,1])==0){
            env_xy[,1]<-env_xy[,1]+rnorm(n=length(env_xy[,1]),mean=0.05,sd=0.0025)
          }
          
          if(bandwidth.nrd(env_xy[,2])==0){
            env_xy[,2]<-env_xy[,2]+rnorm(n=length(env_xy[,2]),mean=0.05,sd=0.0025)
          }
          
          # Configure the dimensions of the matrix
        if(is.null(dim_k)){
            dim_k <- (dim(env_var)[1:2]/2) %>% round(digits = 0)
        
        if(TRUE %in% c(dim_k<50)){
            print("WARNING! Original raster data dimensions lower than 100x100 pixels Kernel matrix size scalled at 100x100!")
            # dim_k <- ifelse(dim_k<100,100,dim_k)
            dim_k <-c(100,100)
              }
          }
        
       # Calculate the 2 dimension distance kernell----  
          # k <- kde2d(x=env_xy[,1],y=env_xy[,2],lims=c(range(env_xy[,1]),range(env_xy[,2])),n=dim_k*0.5)
          # gc()
          
          k<-KernSmooth::bkde2D(env_xy,bandwidth = dim_k*0.1,range.x = list(range(env_xy[,1]),range(env_xy[,2])),gridsize = dim_k,truncate=TRUE)
          #k_env <- apply(t(k$z), 2, rev) 
          
          k_env <- apply(t(k$fhat), 2, rev) 
          if(inv==TRUE){k_env<-max(k_env)- k_env} # 
          
          k_env_scaled<-scale.f(k_env)
         
      # Create the reference environmental space surface
          
          #r1 = rast(k_env_scaled,extent=c(min(k$x),max(k$x),min(k$y),max(k$y)))
          
          r1 = rast(k_env_scaled,extent=c(min(k$x1),max(k$x1),min(k$x2),max(k$x2)))
            
            # Extract the countour with a bandwidth of 0.05 (preserve most of the environmental space)
            env.pol<-lapply(bandwidth.env,function(x) countour_pols(r1,bandwidth = x))
            env.pol<-do.call("rbind",env.pol)
            
      # Create the species observation surface    
            # k_sp <- kde2d(x=env_xy[rownames(env_xy) %in% p_xy_index,1],
            #               y=env_xy[rownames(env_xy) %in% p_xy_index,2],n=dim_k,
            #               lims=c(range(env_xy[,1]),range(env_xy[,2])))    
            
            k_sp<-KernSmooth::bkde2D(env_xy[rownames(env_xy) %in% p_xy_index,],
                                     bandwidth = dim_k*0.1,
                                     range.x = list(range(env_xy[,1]),
                                     range(env_xy[,2])),gridsize = dim_k)
            
            
            k_sp_l <- apply(t(k_sp$fhat), 2, rev) 
            if(inv==TRUE){k_sp_l<-max(k_sp_l)- k_sp_l} #
            
            k_sp_scaled<-scale.f(k_sp_l) 
            
            r2 = rast(k_sp_scaled,extent=c(min(k_sp$x1),max(k_sp$x1),min(k_sp$x2),max(k_sp$x2)))
            
            spp.pol<-lapply(bandwidth.sp,function(x) countour_pols(r2,bandwidth = x))
            spp.pol<-do.call("rbind",spp.pol)
            
      # Add some arguments
            env.pol$Area <- env.pol %>% st_area() ; env.pol$Type<-"Env_dens"
            spp.pol$Area <- spp.pol %>% st_area() ; spp.pol$Type<-"Spp_dens"
            
            ENV_pols <- rbind(env.pol,spp.pol) ; ENV_pols$ID<-1:nrow(ENV_pols)
            
      # Create the background points based on the density kernels for the different bandwidths      
            dat_xy <- env_xy %>% as.data.frame() ; dat_xy <- dat_xy %>% mutate(cells=rownames(env_xy))
            dat_xy <- dat_xy %>% st_as_sf(coords=colnames(env_xy))      
            
            bk_area <- ENV_pols %>% filter(Type=="Env_dens")
            
        for(ww in unique(bk_area$level)){
              bk_pol<- bk_area %>% filter(level==ww)
              p_index <- dat_xy %>% st_intersects(bk_pol,sparse = FALSE)
             
              if(ncol(p_index)>1){
                p_index <- apply(p_index,1,sum)
                p_index <- ifelse(p_index>0,TRUE,FALSE)
              }
              
              dat_xy_int <- dat_xy %>% mutate(overlap=p_index,.before="geometry")
              dat_xy_prob <- dat_xy_int %>% mutate(dens = terra::extract(r1,dat_xy_int %>% vect(),ID=FALSE),
                                          .before="geometry") %>% st_drop_geometry()
              
              names(dat_xy_prob)[ncol(dat_xy_prob)] <- "dens"
              
              cells_bk <- dat_xy_prob[dat_xy_prob$overlap==TRUE,c("cells","dens")]
              
              # check if the selection of points fits the proportion of bk data
              if(n_bk > nrow(cells_bk)){
                print("WARNING! The number of background points is greather than the number of cells, setting number of Bk points equal to selected environmental space!")
                n_bk.1 <- nrow(cells_bk)
                
              }else{
                n_bk.1 <- n_bk
              
              }
              
              Selected_cells <- sample(cells_bk$cells %>% unlist() %>% as.numeric(),size= n_bk.1,
                                       prob = cells_bk$dens %>% unlist())
              
              points_bk_env <- dat_xy_int[dat_xy_int$cells %in% Selected_cells,]
              
              # Transform the selection into points----
              bk_geo <- data.frame(lon=env_var %>% terra::xFromCell(Selected_cells),
                                     lat=env_var %>% terra::yFromCell(Selected_cells))
              
              bk_geo <- bk_geo %>% st_as_sf(coords=c("lon","lat")) %>% st_as_sfc(crs=crs(env_var))
              
              if(ww == unique(bk_area$level)[1]){
                bk_geo_points <- st_sf(data.frame(ID=ww,Type=bk_pol$Type %>% unique(),
                                                    level=bk_pol$level %>% unique(),
                                                    geom=bk_geo %>% st_union()),crs=crs(env_var)
                                       )
                
                }else{
                bk_geo_points <- rbind(bk_geo_points,st_sf(data.frame(ID=ww,Type=bk_pol$Type %>% unique(),
                                                         level=bk_pol$level %>% unique(),
                                                         geom=bk_geo %>% st_union()),crs=crs(env_var))
                                       )
                }
              
              # check bk points cover over the environment----
              # # MCP bk points
              # hull_bk <- chull(points_bk_env %>% st_coordinates())
              # coord_bk <- points_bk_env[c(hull_bk,hull_bk[1]),] %>% st_coordinates()       
              # 
              # bk.mcp.pol <- st_sf(data.frame(ID=nrow(ENV_pols)+1,Type=c("Bk_mcp"),geom=st_sfc(st_polygon(list(coord_bk)))))
              # bk.mcp.pol <- bk.mcp.pol %>% mutate(level=bk_area$level[1],.before=1)
              # bk.mcp.pol$Area <- st_area(bk.mcp.pol)
              # 
              # ENV_pols<-rbind(ENV_pols,bk.mcp.pol)
              # 
              # # contour bk points
              # if(bandwidth.nrd(coord_bk[,1])==0){
              #   coord_bk[,1]<-coord_bk[,1]+rnorm(n=length(coord_bk[,1]),mean=0.05,sd=0.0025)
              # }
              # 
              # if(bandwidth.nrd(coord_bk[,2])==0){
              #   coord_bk[,2]<-coord_bk[,2]+rnorm(n=length(coord_bk[,2]),mean=0.05,sd=0.0025)
              # }
              # 
              # 
              # k <- kde2d(x=coord_bk[,1],y=coord_bk[,2],lims=c(range(env_xy[,1]),range(env_xy[,2])),n=dim_k)
              # k_bk <- apply(t(k$z), 2, rev) 
              # 
              # k_bk_scaled<-scale.f(k_env)
              # 
              # # Create the reference environmental space surface
              # r3 = rast(k_bk_scaled,extent=c(min(k$x),max(k$x),min(k$y),max(k$y)))
              # 
              # # Extract the countour with a bandwidth of 0.05 (preserve most of the environmental space)
              # bk.k.pol <- as.contour(r3,levels=bandwidth.env) # Include loose bandwidths
              # bk.k.pol <- bk.k.pol %>% st_as_sf() %>% st_cast("POLYGON")
              # 
              # bk.k.pol <- bk.k.pol %>% mutate(Area=st_area(bk.k.pol),Type="bk_k",ID=nrow(ENV_pols)+c(1,2),.after=1)
              # ENV_pols <- rbind(ENV_pols,bk.k.pol)
              # 
              # # Extract the set of cells within the species environmental space (buffer threshold)
              # r_empty<-rast(vals=NA,nrow=env_x$dim["rows"],ncols=env_x$dim["colums"],
              #             crs=env_x$crs,extent=env_x$entent)
              # r_env<-r_empty
              # 
              # # The general environmental polygon
              # r_env[c(cells_bk$cells %>% as.numeric())]<-1
            }
            
            bk_geo_points$ID<-1:length(unique(bk_area$level))
          }
            
            # Create the environmental surfaces for the
            r_empty<-rast(vals=NA,nrow=env_x$dim["rows"],ncols=env_x$dim["colums"],
                  crs=env_x$crs,extent=env_x$entent)
    
            rr <- list()  
            
            for(ll in 1:nrow(ENV_pols)){
                
                y_index <- dat_xy %>% st_intersects(ENV_pols[ll,],sparse = FALSE)
                cells <- dat_xy[y_index,"cells"] %>% st_drop_geometry() %>% unlist() %>% as.numeric()
                
                rast_x<-r_empty
                rast_x[cells]<-1
                
                # Create the geographical polygons
                rast_pol<- rast_x %>% terra::as.polygons(crs=crs(rast_x))
                rast_pol <- rast_pol %>% st_as_sf() %>% st_union()
                
                if(ll==1){
                  GEO_pols <- st_sf(data.frame(ID=ll,Type=ENV_pols[ll,]$Type,level=ENV_pols[ll,]$level,
                                               Area=st_area(st_as_sf(rast_pol)),
                                                geom=st_as_sf(rast_pol)))
                }else{
                  GEO_pols <- rbind(GEO_pols,st_sf(data.frame(ID=ll,Type=ENV_pols[ll,]$Type,level=ENV_pols[ll,]$level,
                                                    Area=st_area(st_as_sf(rast_pol)),
                                                      geom=st_as_sf(rast_pol))))
                  }
                
                names(rast_x) <- paste(ENV_pols[ll,]$Type,ENV_pols[ll,]$level,sep="-") %>% st_drop_geometry()
                rr[[ll]]<-rast_x
            }
            
          # Create the geo_clustering
          if(density_pc==TRUE){
            dens=terra::extract(r1,env_xy %>% as.data.frame() %>% 
                                    st_as_sf(coords=c("PC1","PC2")) %>% 
                                              vect(),ID=FALSE) 
            xx <-r_empty
            xx.1 <-r_empty
            xx.2 <-r_empty
            
            xx[rownames(env_xy)%>%as.numeric()]<-dens
            xx.1[rownames(env_xy)%>%as.numeric()] <- env_xy[,1]
            xx.2[rownames(env_xy)%>%as.numeric()] <- env_xy[,2]
            
            geo_K_pca <- list(kernell=xx,PC1=xx.1,PC2=xx.2) %>% rast()
            }
          
            gc()
          # Export the results
            if(density_pc==TRUE){
             results<-list(points = bk_geo_points,
                            env_pol = ENV_pols,
                              geo_pol = GEO_pols,
                                surfaces = rast(rr),
                                  env_kernell = list(Env=r1,Spp=r2) %>% rast(),
                                   geo_kernell =  geo_K_pca,
                                    PCA_info = env_y)
            }else{
              results<-list(points = bk_geo_points,
                            env_pol = ENV_pols,
                            geo_pol = GEO_pols) 
              
              }
            
             return(results)
        }
          
       # if(plot=T){
       #   
       #   # We are going to compare how the different distributions of data look, for a nice comparison of the data we are going to
       #   # standarize the data by substracting the mean and divided by the standard deviation so all the datasets have a mean of 0 and 
       #   # a standard deviation of 1
       #   
       #  st_env<-apply(env_x$tab[,-1],2,scale,center=T,scale=T)
       #  st_bk_env <- terra::extract(env_var,bk_geo_points %>% vect(),ID=FALSE,cells=FALSE) ; st_bk_env<-apply(st_bk_env,2,scale,center=T,scale=T)
       #  st_pres_env <- terra::extract(env_var,p_xy,ID=FALSE,cells=FALSE) ; st_pres_env<-apply(st_pres_env,2,scale,center=T,scale=T)
       #   
       #  par(mfrow=c(5,5))
       #  
       #   for(jj in 1:ncol(vals_bk_env)){
       #     hist(st_env[,jj],freq=F,col=NA,border=NA,main="",ylim=c(0,1))
       #     polygon(x=density(st_bk_env[,jj])$x,y=density(st_bk_env[,jj])$y/max(density(st_bk_env[,jj])$y),col="skyblue4" %>% adjustcolor(alpha.f = 0.2),border=NA)
       #     polygon(x=density(st_env[,jj])$x,y=density(st_env[,jj])$y/max(density(st_env[,jj])$y),col="darkgreen" %>% adjustcolor(alpha.f = 0.2),border=NA)
       #     polygon(x=density(st_pres_env[,jj])$x,y=density(st_pres_env[,jj])$y/max(density(st_pres_env[,jj])$y),col="gold" %>% adjustcolor(alpha.f = 0.2),border=NA)
       # 
       #     mtext(side=3,adj=0,line=1,names(env_x$tab[,-1])[jj],cex=1.2)
       #     legend("topright",legend=c("Environment","Background","Presence points"),
       #           col=c("darkgreen","skyblue4","gold") %>% adjustcolor(alpha.f = 2),
       #           pch=15,pt.cex=1.2,bty="n")  
       #     }

# Accesory functions-------------------------------------------------------------------#

# Create polygons from countour lines for the density maps
  countour_pols <- function(x,bandwidth){
                    xM<-max(minmax(x))
                    y <- x %>% terra::classify(rcl=c(bandwidth,xM,bandwidth))
                    yp <- as.polygons(y,round=FALSE) %>% terra::union()
                    
                    values(yp)<-bandwidth
                    names(yp) <-"level"
                    
                    ypz <- yp %>% st_as_sf() #%>% st_cast("POLYGON")
                    
                    return(ypz)
                  }

# Rescalled values between 0-1 by dividing by the maximun value of x
  scale.f <- function(x, ...){
    z<-(x)/max(x,...)
    print("values rescalled betweeen 0-1")
    return(z)
  }
