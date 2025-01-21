#
###//\/\/\/\/\/\/\////\/\/\/\/\/\/\//\/\/\/\//\////\/\/\///////////////////\\\\\\##-#
#                 Cleaning the points for the maxent analysis
###///\/\/\/\/\/\/\////\/\/\/\/\/\/\\/\////\/\/\/\/\/\/\/\/\/\\\\\/\/\/\/\/\/\///##-#
#
#
Prepare_points<-function(points_sp, # Gbif data output
                         range_sp=NULL,
                         xy.c=c("decimalLongitude","decimalLatitude"), # Variables describing the points coordinates (check GBIF variable naming to correctly assign thsi)
                         b.width=0, # if range_sp is != NULL, a buffer to apply to the distribution data in order to include or exclude points
                         crs.r="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # CRS for the spatial information, default = wgs84
                         # wrld_pol=NULL # optional we can add a polygon to remove the observations that fall into the ocean or water bodies - (Redundant with range_sp; comment and rewrite)
                         ){
   # 0. Load the required packages if there are not loaded ----
   list.of.packages<-c("dplyr","sf","data.table","CoordinateCleaner","terra")
   
   new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
   if(length(new.packages)) install.packages(new.packages)
   
   lapply(list.of.packages,require,character.only=TRUE)
   rm(list.of.packages,new.packages)
   
   T1<-Sys.time()
   
   # a. Clean the GBIF points using Coordinate cleaner----
   # Does the data came directly from GBIF?
   if(class(points_sp)=="gbif_data"){
         points_sp <-points_sp$data %>% as.data.frame()
      }
     
   # Process the spatial information   
   if(is.data.frame(points_sp)){
         
      y<-CoordinateCleaner::clean_coordinates(points_sp,
                                              lon = xy.c[1],
                                              lat = xy.c[2],
                                              tests = c("capitals", "centroids", "equal", 
                                                        "gbif", "institutions", "outliers",
                                                        "seas", "zeros"))$.summary
      points_sp<-points_sp[y,] %>% as.data.frame()
    
     }
   
   # b.1 Check the overlap of the data with the designated range information----
   if(!is.null(range_sp)){
      
      # configure the range data
      if(is.character(range_sp)){
      range_sp <- range_sp %>% sf::st_read() %>% st_transform(crs=crs.r)
      }
      
      range_sp <- range_sp %>% sf::st_buffer(b.width) %>% st_union()
      
      # Check the geometries
      if(st_is_valid(range_sp)==FALSE){
         range_sp<- range_sp %>% st_make_valid()
      }
         # Configure the point data      
         points_sp$ID_p<-1:nrow(points_sp)
         cords_points <- points_sp[,xy.c]
         
         points_sp<- points_sp %>% st_as_sf(coords=xy.c,crs=crs.r)
         
         # Discard non-overlapping/intersecting points
         index_intersect<-st_intersects(points_sp,range_sp,sparse=FALSE) #%>% st_geometry() %>% plot(col="blue",add=TRUE) #%>% unlist()
         
         points_sp <- points_sp[index_intersect,] ; cords_points <- cords_points[index_intersect,]
         points_sp <- points_sp %>% as.data.frame() ; points_sp<-points_sp[,!names(points_sp) %in% "geometry"]
         
         points_sp<-cbind(points_sp,cords_points)
         }

   return(points_sp)

}

# End of the function