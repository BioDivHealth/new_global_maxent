##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\><><><|><|><>||||~/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
## Function to retrieve the taxonomic information and download spatial information from Gbif
##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\--><><><><><|>|<|<|><|>|>|>|>|~}{];/|\[;][;];][;][;]#
#


Spatial_spp <- function(sci_sp, # Scientific name of the species from which we want to gather spatial information
                        p.route=paste(getwd(),"Data/sp_points",sep="/"), # Folder to store the spatial information
                        range_sp=NULL, # do we have a shapefile with the range of the species?
                        start_date=2015 # The initial date for the spatial query [the function would look from that date onwards till present date]
){
  
  # 0. Packages and dependencies:
  options(iucn_redlist_key="eb704359f6ea22d50235efebf7f4a2f5f843a66fa951cf6a376df71cf7268986")
  
  list.of.packages<-c("tidyr","rredlist","taxize","data.table","stringr",
                      "sp","rgbif","raster","data.table","dplyr","doParallel","parallel")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # 0.1 Create the exit folder 
  p.route %>% dir.create(recursive=TRUE,showWarnings = FALSE)
  
  # a. Check species names and extract synonims information
  y_sp <- retrieve_syns(spp_name=sci_sp,Gbif=FALSE)$Spp_syn
  
  # b. with the list of synomyns, dowload all the abailable spatial information
  Dowload_gbif(sp_list=y_sp, # (Character) List of species from which to dowload spatial information
               initial_date=start_date, # (Numeric/year) By default the function will dowload 500 records for each month of the year, from the year specified till present
               exit_route=p.route, # (Character) Route to store the dowloaded information
               area=range_sp, # (character) Searches for occurrences inside a polygon in Well Known Text (WKT) format. A WKT shape written as either
               gadm_codes=NULL, # (character) The gadm id of the area occurrences are desired from. https://gadm.org/.
               locality=NULL, # If iterating around different polygons, the name of the regions or polygon
               n_records=150000)
}



#
#
#################################################################################
#       dowload species GBIF point distribution for a given set of species      #-//\/\/\/\/\/\
#################################################################################
#
# This script automatically dowload the Gbif data for a list of species on an specified area
#
#
# 
Dowload_gbif<-function(sp_list, # (Character) List of species from which to dowload spatial information
                       initial_date, # (Numeric/year) By default the function will dowload 500 records for each month of the year, from the year specified till present
                       # n_cores=2, # (Numeric) Number of cores commited to the processing (only when sp_list>1)
                       exit_route, # (Character) Route to store the dowloaded information
                       area=NULL, # (character) Searches for occurrences inside a polygon in Well Known Text (WKT) format. A WKT shape written as either
                       gadm_codes=NULL, # (character) The gadm id of the area occurrences are desired from. https://gadm.org/.
                       locality=NULL, # If iterating around different polygons, the name of the regions or polygon
                       n_records=150000 # Maximun number of records to retrieve
){
  
  # Set the dates for the searching
  present_date<-Sys.time() %>% format("%Y") %>% as.numeric()
  years=c(initial_date:as.numeric(present_date))
  
  # Create the output folder for the data
  if(!dir.exists(exit_route)){dir.create(exit_route,showWarnings = FALSE)}
  
  if(length(sp_list)<2){
    
    for (k in 1:length(years)){
      print(paste0("year=",k))
      
      for (j in 1:12){
        
        points <- NULL
        t_11 <- 1
        while(is.null(points) && t_11 <= 1000) {
          
          points<-try(occ_search(scientificName =sp_list,
                                 hasCoordinate=TRUE,
                                 year=years[k],
                                 month=j,
                                 hasGeospatialIssue=FALSE,
                                 limit=n_records,
                                 geometry=area,
                                 gadmGid=gadm_codes),
                      silent = FALSE)
          t_11 <- t_11 + 1
        }
        rm(t_11)
        # 
        
        if(is.null(points$data)){
          next
        }else{
          if(!exists("y_points")){
            y_points<-points$data
          }else{
            y_points<-data.table::rbindlist(list(y_points,points$data),fill=TRUE)
          }}
        rm(points)
        print(paste(k,j,sep="-"))
      }
    }
    gc()
    
    if(exists("y_points")){
      
      if(is.null(sp_list)){
        sp_list<-paste("All_records",locality)
      }
      
      write.csv(y_points,paste(exit_route,paste0(sp_list,".csv"),sep="/"),row.names = FALSE)
      rm(y_points)
    }
    
  }else{
    # Create a temporal folder to host the information until its combined
    Temp_tax<-paste(getwd(),"Temp_tax",sep="/")
    dir.create(Temp_tax,showWarnings = FALSE)
    
    for(f in 1:length(sp_list)){
      t1<-Sys.time()
      sp<-sp_list[f]
      
      for (k in 1:length(years)){
        print(paste0("year=",k))
        
        for (j in 1:12){
          
          points <- NULL
          t_11 <- 1
          
          while(is.null(points) && t_11 <= 1000) {
            
            points<-try(occ_search(scientificName =sp,
                                   hasCoordinate=TRUE,
                                   year=years[k],
                                   month=j,
                                   hasGeospatialIssue=FALSE,
                                   limit=n_records,
                                   geometry=area,
                                   gadmGid=gadm_codes),
                        silent = FALSE)
            t_11 <- t_11 + 1
          }
          rm(t_11)
          
          if(is.null(points$data)){
            next
          }else{
            if(!exists("y_points")){
              y_points<-points$data
            }else{
              y_points<-data.table::rbindlist(list(y_points,points$data),fill=TRUE)
            }}
          rm(points)
          print(paste(k,j,sep="-"))
        }
      }
      gc()
      
      if(exists("y_points")){
        write.csv(y_points,paste(Temp_tax,paste0(sp,".csv"),sep="/"),row.names = FALSE)
        rm(y_points)
      }
    }
    
    # Combine the files from the temporal folder, export it to the relevant folder and 
    # erase the temporal folder
    points_spp<-lapply(list.files(Temp_tax,pattern = "csv",full.names = TRUE),fread)
    # points_spp<-points_spp %>% rbindlist(fill=TRUE)
    
    LL<-lapply(points_spp,function(x) apply(x,2,as.character)%>% as.data.frame()) # Data frames contain different column classes whihc cause problems when merging
    points_spp<-LL %>% rbindlist(fill=TRUE)
    
    write.csv(points_spp,paste(exit_route,paste0(sp_list[1],".csv"),sep="/"),row.names = FALSE) # Export the point information under the original name
    unlink(Temp_tax,recursive=TRUE) # Erase the temporal data
    
  }
  return("GBIF data downloaded") 
}