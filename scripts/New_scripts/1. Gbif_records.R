rm(list=ls())
gc()
#.rs.restartR()
options(java.parameters = "-Xmx15g") # increase the memory space for jave before loading any package
options("rgdal_show_exportToProj4_warnings"="none") # Silence packages updates warnings

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td<-tempdir()
dir.create(td,showWarnings = FALSE)
#
#
###//\/\/\/\/><\/\////\/\/\/\/\/\/\///\\\\\\//\/\/\/\////////////////////></////##-#
##                          Gbif Original records                               ##-#  
###///\/\/\/><\/\/\////\/\/\//\/\/\/\/\/\///\\\\\\\\\..\.\\\\\\\\\\\\\\\\><\\\\\##-#
#'
#'
# 0. Load the packages----
list.of.packages<-c("sf","terra","tidyverse","rJava","viridis","data.table","parallel","future","future.apply")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Some project parameters (need to check the information of the Raster Files)----
# CRS and spatial standarization
# crs_p <- "EPSG:4083" # USE the regCan UTM area

# 0.1 Load the functions----
  functions<-"C:/Users/Gonzalo/OneDrive - Natural History Museum/Projects/New_global_maxent/new_global_maxent/scripts/New_functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
  lapply(functions,function(x) source(x))

# 1. Get the species occurrence data ----
# a. GBIF data----
  # Parameters for the point-data
  route_points <- "C:/Users/Gonzalo/OneDrive - Natural History Museum/Projects/New_global_maxent/new_global_maxent/data/points"
  D.points <- route_points  %>% list.files(pattern=".r",full.names = TRUE)

  get.points <-function(y){
                            load(y) # load the data
                            if(exists("data1")){                              
                              x <- data1$data1 # Get the point data
                              rm(data1)
                              }
                            if(exists("pack1")){
                              x <- pack1$data1$data1 # Get the point data
                              rm(pack1)
                              }
                            
                              x <- sf::st_as_sf(x) # Transform from sp to sf object
                              return(x)
                              }
  
  
  # a.1 Load all the point data----
  n_cores<-detectCores(logical=F)-1
  if(n_cores<1) n_cores<-1 # most machinese would have more than one physical core, but just in case
  
  # Create the cluster              
  cl <- makeCluster(n_cores)
  clusterExport(cl, list("get.points","D.points")) # put the needed information into the "slaves"
  
  # Load the libraries
  list.of.packages<-c("sf","terra","tidyverse","rJava","viridis","data.table")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Runt the function in parallel
  p.Gbif <- parallel::parLapply(cl,D.points,get.points)
  
  stopCluster(cl) # close the cluster
  
 # Combine all the point data into a single object                         
  check <- lapply(p.Gbif,function(x) x %>% ncol()) ; check <- do.call("c",check)
    p.Gbif <- p.Gbif[c(check==7)]
    
    p.Gbif.comb <- do.call("rbind",p.Gbif)
   
  # Tidying and cleaning up  
    p.Gbif.comb <- p.Gbif.comb[!c(p.Gbif.comb$name %>% grepl(pattern="[[:digit:]]")|is.na(p.Gbif.comb$name)),]# Some species "name" are just numbers or `NA` (remove them?) 
    
  # Preview the spatial information
    plot(p.Gbif.comb["name"],pch=19,cex=0.05,main="Species records")
  
  # b. Get some species information----  
    spp_list_raw <- unique(p.Gbif.comb$name) %>% as.character() # Some names are just really long numbers!!!
    
      # Retrive synonims
        # The retrive synonims function doesn't like to run in paralell (possibly due to the API conexions to external servers)
        # we are going to try to Future package to run multiple sessions of R in the background
    tax_raw <- lapply(index2$Species, function(y) retrieve_syns(spp=y),future.seed = TRUE)
    
      
    
        
    
    point_route <- "C:/Users/Gonzalo/OneDrive - Natural History Museum/Projects/New_global_maxent/New_Gbif_points"
    
    
    
    
    
    
    lapply(spp_list,Spatial_spp,p.route=point_route,start_date=1800) # we are going to make a SDM for the red-necked wallaby
    list.files("./Data/sp_points",full.names = TRUE)

Gbif_points <- lapply(list.files("./Data/sp_points",full.names = TRUE),read.csv)
names(Gbif_points) <- basename(list.files("./Data/sp_points")) %>% gsub(pattern=".csv$",replacement="")

# Remove points with issues
Gbif_points <- lapply(Gbif_points,function(w) w %>% filter(basisOfRecord=="HUMAN_OBSERVATION", issues != "cdc,cdround,gass84,muluriiv"))

# b. IUCN data----
IUCN_pol <- spp_list %in% c("D:/Data/Spatial information/IUCN spatial data/All polygons" %>% list.files(pattern=".shp$",recursive=TRUE))

# c. Clean the points with coordinate cleaner----
Gbif_points <- lapply(Gbif_points,Prepare_points,crs.r = crs_p)

# d. Transform the points into spatial objects and project to the right CRS----
Gbif_points$`Fringilla polatzeki` <- Gbif_points$`Fringilla polatzeki` %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG:4326") %>% st_transform(crs_p)
Gbif_points$`Fringilla teydea` <- Gbif_points$`Fringilla teydea` %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG:4326") %>% st_transform(crs_p)

# e. SEO birdlife data----
Points_SEO <- readxl::read_xlsx("./Data/Records/Ficha 2024 - 15092024.xlsx",sheet=1) %>% as.data.frame()
Points_SEO <- Points_SEO %>% st_as_sf(coords=c("UTM Lon","UTM Lat"),crs=crs_p)

# Run some spatial thining for each date
xtabs(~Fecha,data = Points_SEO)
filtered_points <- lapply(unique(Points_SEO$Fecha),function(w) sp_thining(Points_SEO %>% filter(Fecha==w), dist.t = 150,units.d = "m",reps=1)[[1]])

Points_SEO.Thin <- do.call("rbind",filtered_points)

# Load the Nest information
Nests_SEO <- read.csv("./Data/Records/Nidos exitosos CUMBRE 2008-2023.csv",sep=";")
Nests_SEO <- Nests_SEO %>% st_as_sf(coords=c("lon","lat"),crs=crs_p)

xtabs(~Año,data = Nests_SEO)
Nests_SEO <- Nests_SEO %>% filter(Año>=2019, Pollos.que.salen>=1)
