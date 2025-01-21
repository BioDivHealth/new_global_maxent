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



##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\><><><|><|><>||||~/\/\/\/\/\/\#
## Extract the taxonomic information and list of sinonyms for a given specie
##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\--><><><><><|>|<|<|><|>|>|>|>|~#
#

retrieve_syns<-function(spp_name,   # [Character] The species name from which to collect taxonomic information
                        n_times=100,  # [Numeric] Number of times the search is repeated until a data is found,default value = 1
                        Gbif=FALSE  # [Logical] Should we check Gbif for a taxonomic macthing of the species
)
{
  # 0. Load the packages
  list.of.packages<-c("tidyr","rredlist","taxize","data.table","stringr",
                      "rgbif","raster","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # a. Check the species names, if the name is binomial it runs the query to collect the taxonomic information
  # a.1 Removes special characters and fix the capitalization of the name
  spp.x <- stringr::str_trim(spp_name) %>% gsub(pattern = "[[:punct:]]", replacement = " ") %>% stringr::str_trim("both") %>% gsub(pattern = "  ",replacement = "") # Remove white spaces at the start and end of the string
  
  # Resolve capitalization
  CapSp <- function(x) {
    s <- strsplit(x, " ") %>% unlist()
    
    if(length(s)>1){
      paste(paste0(toupper(substring(s[1], 1,1)), substring(s[1], 2)),
            paste(tolower(s[-1]),collapse = " "), sep = " ")
      
    }else{
      paste0(toupper(substring(s[1], 1,1)), substring(s[1], 2))
    }
  }
  
  spp.x<-CapSp(spp.x)
  
  # a.2 Check if the name is related with a species/sub-specie or other taxonomic class (Class, Order, Family) 
  #     by analyzing the number of terms of the character string  
  correct_name <- NULL
  t_11 <- 1
  
  while(is.null(correct_name) && t_11 <= n_times) {
    
    try(correct_name<-gnr_resolve(sci =spp.x,
                                  resolve_once = FALSE,
                                  canonical=TRUE,
                                  best_match_only = TRUE,
                                  highestscore = TRUE),
        silent = TRUE)
    t_11 <- t_11 + 1
  }
  rm(t_11)  
  
  # Summarize the results of the name checking and correction
  if (!is.null(correct_name)|nrow(correct_name)!=0){
    y.d<-cbind(spp.x,correct_name[,colnames(correct_name)%in%
                                    c("data_source_title","score","matched_name2")]) #
    
    names(y.d)[1]<-"or_name"
    
    spp.x<-y.d$matched_name2  # Use the corrected name for the rest of taxnomic querys
    
    ifelse(y.d$matched_name2==y.d$or_name,
           y.d$Status<-"Correct",
           y.d$Status<-"Incorrect")
    
  } else {
    y.d<-data.frame(or_name=spp.x,
                    matched_name2=NA,
                    Status="Not_found",
                    data_source_title=NA,
                    score=NA)
  }
  
  # b.Use the corrected or original name to look for taxonomic data----
  #   b.1. Get the basic data from the IUCN red list----
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  g <- NULL
  t_11 <- 1
  
  while( is.null(g) && t_11 <= n_times){
    
    try(g <- rl_search(name = spp.x,silent=TRUE)$result)
    t_11 <- t_11 + 1
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  b <- NULL
  t_2 <- 1
  while( is.null(b) && t_2 <= n_times) {
    
    try(b <- rl_synonyms(name=spp.x,silent=TRUE))
    t_2 <- t_2 + 1
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if (length(g)==0 & b$count==0) {
    
    IUCN_id  <- NA
    IUCN_name<- NA
    IUCN_Phylum   <- NA
    IUCN_Class    <- NA
    IUCN_Order    <- NA
    IUCN_Family   <- NA
    
    IUCN_Category <- NA
    IUCN_Present <- "No"
    
    IUCN_N_syn <- NA
    IUCN_syn <- NA
    
  } else {
    # if the name of the species correspond to a synonim in the IUCN-red list, 
    # use the accepted name to retrieve the species information
    
    if(b$count!=0 & length(g)==0){ 
      
      g <- NULL
      t_11 <- 1
      
      while( is.null(g) && t_11 <= n_times){
        
        try(g <- rl_search(name = b$result$accepted_name,silent=TRUE)$result)
        t_11 <- t_11 + 1
      }
    }
    
    IUCN_id  <- unique(g$taxonid) 
    IUCN_name <- unique(g$scientific_name)
    
    if (is.null(g$phylum)){ IUCN_Phylum <- NA } else {IUCN_Phylum <- unique(g$phylum) }
    
    if (is.null(g$class)){ IUCN_Class <- NA }else {IUCN_Class <- unique(g$class) }
    
    if (is.null(g$order)){ IUCN_Order <- NA }else {IUCN_Order <- unique(g$order) }
    
    if (is.null(g$family)){ IUCN_Family <- NA } else {IUCN_Family <- unique(g$family)}
    
    IUCN_Category <-g$category
    IUCN_Present <-"Yes"
    
    
    if (b$count!=0){
      
      z_c<-as.character(sapply(strsplit(b$result$synonym,split = " "),length)==2)
      
      IUCN_N_syn <- length(z_c[z_c==TRUE]) # number of sp sinonims, Subsp excluded
      
      IUCN_syn<-paste(b$result$accepted_name %>% unique(),paste(b$result$synonym[sapply(strsplit(b$result$synonym,split = " "),length)==2],
                                                                collapse=";"),sep = ";") # combine the names into a single string with all the sinonims
      
    } else {
      
      IUCN_N_syn<-NA
      IUCN_syn<-NA
    }
  }
  
  # b.2 Combine the IUCN information----
  IUCN_data<-data.frame(IUCN_Present,IUCN_id,IUCN_name,
                        IUCN_Category,IUCN_N_syn,IUCN_syn,
                        IUCN_Phylum,IUCN_Class,
                        IUCN_Order,IUCN_Family)
  
  # b.3. Get the basic data from the ITIS red list----
  # b.3.1 Get TSN a reference number that we are going to need to gather information from ITIS----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  TSN <- NULL
  t_4 <- 1
  
  while( is.null(TSN) && t_4 <= n_times){
    try(TSN <- get_tsn_(spp.x,searchtype = "scientific",silent=TRUE))
    t_4 <- t_4 + 1
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if (is.null(TSN[[1]])) { 
    
    ITIS_Present<-"No"
    ITIS_id<-NA
    ITIS_name<-NA
    ITIS_is_valid<-NA
    ITIS_Phylum<-NA
    ITIS_Class<-NA
    ITIS_Order<-NA
    ITIS_Family<-NA
    
    ITIS_accept<-NA
    ITIS_syn<-NA
    ITIS_N_syn<-NA
    
    
  } else {
    
    tsn<-TSN[[1]]
    
    tsn_n<-tsn[tsn$scientificName==spp.x & tsn$nameUsage=="valid",]
    
    if(length(tsn_n$tsn)==0){ # if the original species name is not present in the returned list of synonims we igore the results
      
      ITIS_Present<-"Wrong"
      ITIS_id<-NA
      ITIS_name<-NA
      ITIS_is_valid<-NA
      ITIS_Phylum<-NA
      ITIS_Class<-NA
      ITIS_Order<-NA
      ITIS_Family<-NA
      
      ITIS_accept<-NA
      ITIS_syn<-NA
      ITIS_N_syn<-NA
      
    }else{
      
      ITIS_Present<-"Yes"
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      AceptName <- NULL
      t_4 <- 1
      
      while( is.null(AceptName) && t_4 <= n_times){
        
        try(AceptName <- itis_acceptname(tsn_n$tsn,silent=TRUE))
        t_4 <- t_4 + 1
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      if (AceptName$submittedtsn==AceptName$acceptedtsn){ # is the submitted Or_name name accepted by ITIS
        ITIS_name<-spp.x
        ITIS_accept<-"yes"
        
      }else{
        ITIS_name<-"No"                        
        ITIS_accept<-"No"
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      Syn <- NULL
      t_5 <- 1
      
      while( is.null(Syn) && t_5 <= n_times){
        
        try(Syn <- synonyms(tsn$tsn,db="itis"))
        t_5 <- t_5 + 1
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      # Syn_l is a list of possible codes with their respective names
      # We need to unify the list and extract the names
      
      Syn<-rbindlist(Syn,fill=TRUE)
      # count the number of FALSE and TRUE       #strsplit(, " "), length
      z<-as.character(is.na(Syn$syn_name))
      # if the length equals cero, there is no syn names
      
      if (length(z[z==TRUE])==0){
        ITIS_N_syn<-NA
        ITIS_syn<-NA
      } else {
        
        z_b<-as.character(sapply(strsplit(Syn$syn_name,split = " "),length)==2)
        
        ITIS_N_syn<-length(z_b[z_b==TRUE]) # number of sp sinonims, Subsp excluded
        ITIS_syn<-paste(Syn$syn_name[sapply(strsplit(Syn$syn_name,split = " "),length)==2],
                        collapse=";") # combine the names into a single string with all the sinonims
        
      }
      
      ITIS_id<-tsn_n$tsn
      ITIS_name<-tsn_n$scientificName
      ITIS_is_valid<-tsn_n$nameUsage
      
      # Get the upstream taxonomic information
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      j<-NULL
      t_13<-1
      
      while( is.null(j) && t_5 <= n_times){
        
        try(j<-itis_hierarchy(tsn_n$tsn,what="full",silent=TRUE))
        t_13 <- t_13 + 1
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      ITIS_Phylum<-toupper(as.character(j[j$rankname=="phylum",4])) # Change from lower to upper chase
      ITIS_Class<-toupper(as.character(j[j$rankname=="class",4]))
      ITIS_Order<-toupper(as.character(j[j$rankname=="order",4]))
      ITIS_Family<-toupper(as.character(j[j$rankname=="family",4]))
      
    }
  }
  
  # b.3.2 Unify the taxonomic information from ITIS---- 
  ITIS_data<-data.frame(ITIS_Present,ITIS_is_valid,ITIS_id,ITIS_name,
                        ITIS_Phylum,ITIS_N_syn,ITIS_syn,
                        ITIS_Class,ITIS_Order,ITIS_Family)
  
  # Should we retrieve synonim information from GBIF?
  if(Gbif==TRUE){
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    key_1 <- NULL
    t_6 <- 1
    
    while(is.null(key_1) && t_6 <= n_times){
      
      try(key_1 <- get_gbifid_(sci=spp.x)[[1]],silent = TRUE) # get the taxon key
      t_6 <- t_6 + 1
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    if (length(key_1)==0){
      
      GBIF_Present<-"No"
      GBIF_id<-NA
      
      GBIF_name<-NA
      
      GBIF_Phylum<-NA
      GBIF_Class<-NA
      GBIF_Order<-NA
      GBIF_Family<-NA
      
      GBIF_syn<-NA
      GBIF_N_syn<-NA
      
    } else{
      
      GBIF_id<-paste(key_1[key_1$status=="ACCEPTED" & key_1$matchtype=="EXACT",20],collapse="-")
      GBIF_name<-paste(unique(key_1[key_1$specieskey==GBIF_id,13]),collapse="-")
      
      GBIF_Phylum<-toupper(paste(unique(key_1[key_1$specieskey==GBIF_id,9]),collapse="-"))
      GBIF_Class<-toupper(paste(unique(key_1[key_1$specieskey==GBIF_id,22]),collapse="-"))
      GBIF_Order<-toupper(paste(unique(key_1[key_1$specieskey==GBIF_id,10]),collapse="-"))
      GBIF_Family<-toupper(paste(unique(key_1[key_1$specieskey==GBIF_id,11]),collapse="-"))
      
      if (is.na(GBIF_name)){
        
        GBIF_Present<-"No"
        GBIF_id<-NA
        
        GBIF_name<-NA
        
        GBIF_Phylum<-NA
        GBIF_Class<-NA
        GBIF_Order<-NA
        GBIF_Family<-NA
        
        GBIF_syn<-NA
        GBIF_N_syn<-NA
        
      } 
      
      if (GBIF_name==""){
        
        GBIF_Present<-"No"
        GBIF_id<-NA
        
        GBIF_name<-NA
        
        GBIF_Phylum<-NA
        GBIF_Class<-NA
        GBIF_Order<-NA
        GBIF_Family<-NA
        
        GBIF_syn<-NA
        GBIF_N_syn<-NA
        
      } else {
        
        GBIF_Present<-"Yes" 
        GBIF_syn<-paste(key_1[key_1$status=="SYNONYM",13],collapse=";")
        if(length(key_1[key_1$status=="SYNONYM",1])==0){
          GBIF_N_syn<-NA
        }else{
          GBIF_N_syn<-length(key_1[key_1$status=="SYNONYM",1]) 
        }
      }
    }
    
    GBif_data<-data.frame(GBIF_Present,GBIF_id,GBIF_name,
                          GBIF_N_syn,GBIF_syn,GBIF_Phylum,
                          GBIF_Class,GBIF_Order,GBIF_Family)
  }
  
  # C. Return the Taxonomic information
  if(Gbif==TRUE){
    Tax_dat<-cbind(Or_name=spp.x,IUCN_data,ITIS_data,GBif_data)
    Spp_syn<-c(spp.x,Tax_dat[,colnames(Tax_dat) %in% c("IUCN_name","IUCN_syn","ITIS_name",
                                                       "ITIS_syn","GBIF_name","GBIF_syn")]) %>% paste(collapse = ";") %>% 
      strsplit(split = ";") %>% unlist()
    
    Spp_syn<-Spp_syn[-c(Spp_syn%>%grep(pattern = "NA"))]
    Spp_syn<-Spp_syn[!duplicated(Spp_syn)]
    
    return(list(Spp_syn=Spp_syn,
                IUCN_spp=IUCN_data$IUCN_name,
                Or_name=spp.x,
                TaxDat=Tax_dat))
  }else{
    Tax_dat<-cbind(Or_name=spp.x,IUCN_data,ITIS_data)
    Spp_syn<-c(spp.x,Tax_dat[,colnames(Tax_dat) %in% c("IUCN_name","IUCN_syn","ITIS_name",
                                                       "ITIS_syn","GBIF_name","GBIF_syn")]) %>% paste(collapse = ";") %>% 
      strsplit(split = ";") %>% unlist()
    
    Spp_syn<-Spp_syn[-c(Spp_syn%>%grep(pattern = "NA"))]
    Spp_syn<-Spp_syn[!duplicated(Spp_syn)]
    
    return(list(Spp_syn=Spp_syn,
                IUCN_spp=IUCN_data$IUCN_name,
                Or_name=spp.x,
                TaxDat=Tax_dat))
    
  }
  # Message 
  print("Taxonomic search done")
}


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
