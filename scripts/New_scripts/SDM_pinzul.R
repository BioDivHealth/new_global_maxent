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
###//\/\/\/\/><\/\////\/\/\/\/\/\/\///\\\\\\//\/\/\/\////////////////////////></////##-#
##                        Single species SDM pipeline                               ##-#  
###///\/\/\/><\/\/\////\/\/\//\/\/\/\/\/\///\\\\\\\\\..\.\\\\.\\\\\\\\\\\\\\\><\\\\\##-#
#'
#'
# 0. Load the packages----
list.of.packages<-c("sf","terra","tidyverse","geodata","rJava","viridis")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Some project parameters----
# CRS and spatial standarization
crs_p <- "EPSG:4083" # USE the regCan UTM area

# 0.1 Load the needed functions----
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Get the species occurrence data ----
# a. GBIF data----
spp_list <- c("Fringilla teydea polatzeki","Fringilla teydea","Fringilla polatzeki")

lapply(spp_list,Spatial_spp,start_date=2015) # we are going to make a SDM for the red-necked wallaby
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

# Points are pretty densely package, we are going to run some spatial thining to reduce the potential impacts of these concentration of observations----
xtabs(~Año,data = Nests_SEO)
filtered_Nests_SEO <- lapply(unique(Nests_SEO$Año),function(w) sp_thining(Nests_SEO %>% filter(Año==w), dist.t = 150,units.d = "m",reps=1)[[1]])
Nest_SEO.thin <- do.call("rbind",filtered_Nests_SEO)

# f. Study area polygon----
std_border <- st_read("./Data/Spatial/Vector_information/Pro_BaseLayer.shp") %>% st_transform(crs=crs_p)
co_pol <- std_border[std_border$Island %in% c("Gran_Canaria"),]

# Check the spatial data
plot(co_pol %>% st_geometry())
# plot(Gbif_points$`Fringilla polatzeki` %>% st_geometry(),add=T,col="orange",pch=19,cex=0.75)
# plot(Gbif_points$`Fringilla teydea` %>% st_geometry(),add=T,col="skyblue3",pch=19,cex=0.75)
plot(Points_SEO.Thin %>% st_geometry(),add=T,col="tomato",pch=19,cex=0.75)
plot(Nest_SEO.thin %>% st_geometry(),add=T,col="skyblue",pch=19,cex=0.75)

legend("topleft",legend = c("nidos exitosos (al menos un pollo)","Individuos detectados"),bty="n",
       col=c("skyblue","tomato"),pch=19,pt.cex = 1.2,title="Distribucion de observaciones de Fingilla polatzeki\n periodo 2018 - 2024")

# 3. Load the environmental data ---- 
# 3.1 Bioclimatic and land-cover variables
env_dat<-"./Data/Environmental_information" %>% list.files(pattern = ".tif$",full.names = TRUE) %>% rast()
crs(env_dat)

env_dat <- env_dat %>% terra::project(y=crs_p)
names(env_dat)

# 3.2 Load the satellite data for the calculation of the NDVI metric----
Sentinel_routes <- "./Data/Environmental_information/Sentinel-2" %>% list.files(pattern=".zip",recursive = TRUE,full.names = TRUE)

  for(i in 1:length(Sentinel_routes)){
          # Routes and main jp2 files
          xtemp <- paste(td,"Sentinel",sep="/") ; dir.create(xtemp,showWarnings = FALSE,recursive = TRUE)
          unzip(Sentinel_routes[i],exdir = xtemp)
          list.files(xtemp,pattern = "jp2$",recursive = TRUE)
          
          # Load and export the Sentinel-2 information
          # 
          B04 <- list.files(xtemp,pattern = "B04_10m.jp2$",recursive = TRUE,full.names = TRUE) %>% rast()
          B08 <- list.files(xtemp,pattern = "B08_10m.jp2$",recursive = TRUE,full.names = TRUE) %>% rast()
          
          names(B04) <- paste("B04",i,sep="-")
          names(B08) <- paste("B08",i,sep="-")
          
          writeRaster(c(B04,B08),paste("./Data/Environmental_information/Sentinel-2",paste0("Sentinel",i,".tif"),sep="/"),overwrite=TRUE)
          unlink(xtemp,recursive = TRUE)
          
          rm(B04,B08,xtemp)
          }

# Combine the rasters and harmonize with the rest of spatial information----
sentinel <- "./Data/Environmental_information/Sentinel-2" %>% list.files(pattern = "tif$",full.names = TRUE)
resample.rast(x=env_dat$bio_01,y=sentinel,results.r="./Data/Environmental_information/Sentinel-2/Resampled",
              paralell=TRUE,rm.temp=T,crs.r=crs_p)

sentinel <- "./Data/Environmental_information/Sentinel-2/Resampled" %>% list.files(pattern = ".tif$",full.names = TRUE)

strip1 <- sentinel %>% rast() 

# combine crop and mask
B04 <- merge(strip1$`B04-1`,strip1$`B04-2`) %>% project(crs_p) %>% crop(co_pol %>% vect) %>% mask(co_pol %>% vect)
B08 <- merge(strip1$`B08-1`,strip1$`B08-2`) %>% project(crs_p) %>% crop(co_pol %>% vect) %>% mask(co_pol %>% vect)

# Calculate the NDVI using the formula NDVI = B8-B4B8+B4
NDVI_GranCanaria <- (B08-B04)/(B08+B04) ; names(NDVI_GranCanaria) <- "NDVI"
sentinelBands <- c(B04,B08,NDVI_GranCanaria)

panel(sentinelBands[[c(1,2)]]) ; plot(NDVI_GranCanaria,col=viridis(1000),box="n",axes=F,main="NDVI Gran Canaria")

# Variable selection----
# We have 31 variables, since our points are coarsely distributed it is better to reduce the number of variables and let the Model some freedom for
# prediction
# From Carrascal et al., 2017; The important variables are height of pine trees, forest cover, 
#                              altitude, and rainfall during the driest trimester
#
vars_select <- c("bio_01","bio_04","bio_09","bio_05","bio_17","bio_15","bio_12",
                 "Grass","Shrub","Tree",
                 "PNOA_MDT25_REGCAN95_HU28_1079_LID","NDVI")

env_dat <- env_dat %>% crop(co_pol %>% vect()) %>% mask(co_pol %>% vect())                  
env_dat <- c(env_dat,NDVI_GranCanaria)

env_analysis <- env_dat[[vars_select]]
plot(env_analysis %>% crop(co_pol %>% vect()),col=viridis(1000))

# 3.3.a Display the spatial information----
#~~~~~~~~~~~~~~~~
dir.create("junk") # Terra is having a hard time right now, this solve the errors with the mask function
#~~~~~~~~~~~~~~~~
# Some parameters
# colors <- colorRampPalette(c("#fdee9dff","#99531eff","#91ad52ff","#287624ff")) ; colors<-colors(10000)
colors <- viridis(1000)

results_r <- "Results/Figures/" ; results_r %>% dir.create(recursive=TRUE,showWarnings = FALSE)

png(paste(results_r,paste0("sp_general_map",".png"),sep="/"),width=26,height = 25,res=800,units="cm")
  par(bg=NA,xpd=TRUE)
  
  s<-terra::terrain(env_analysis[["PNOA_MDT25_REGCAN95_HU28_1079_LID"]],v="slope",unit="radians")
  a<-terra::terrain(env_analysis[["PNOA_MDT25_REGCAN95_HU28_1079_LID"]],v="aspect",unit="radians")
  
  hill_shade <- terra::shade(s,a,angle = 45)
  
  hill_shade %>% writeRaster(paste(results_r,"hillshade.tif",sep="/"),overwrite=T)
  
  env_analysis[["NDVI"]] %>% terra::crop(co_pol %>% vect()) %>% terra::mask(co_pol %>% vect()) %>% 
    plot(axes=F,alpha=0,col=colors,legend=T,plg = list(loc = "rigth", size=c(0.5,1), title = "NDVI")) 
  
  hill_shade %>% terra::crop(co_pol %>% vect()) %>% terra::mask(co_pol %>% vect()) %>% 
    plot(axes=F,add=T,alpha=1,col=grey(0:100/100),legend=F)
  
  env_analysis[["NDVI"]] %>% terra::crop(co_pol %>% vect()) %>% terra::mask(co_pol %>% vect()) %>% 
    plot(axes=F,alpha=0.55,col=colors,legend=T,add=T,plg = list(loc = "rigth", size=c(0.5,1), title = "NDVI")) 
  
  plot(co_pol %>% st_geometry(),add=TRUE)
  # Add the points
    plot(Points_SEO.Thin %>% st_geometry(),add=TRUE,pch=19,cex=0.55,col="tomato")
    plot(Nest_SEO.thin %>% st_geometry(),add=TRUE,pch=19,cex=0.55,col="skyblue")
  
  legend("bottomleft",legend=c("Nidos","Observaciones"),title="Tipo de registro",pch=19,col=c("skyblue","tomato"),bty="n")
  #  mtext(side=3,adj=0,"A",font=2,cex=2,xpd=TRUE,bty="n")

dev.off()

# 4. Run MaxEnt for the species and select the best performing models----
# Run some preeliminary test-models:
# For the nest distribution
r.nests <- Auto_maxent(presence_dat=Nest_SEO.thin, 
                        predictors=env_analysis, 
                        rm.dp = TRUE,
                        crs.r = crs_p,
                        name.mod = "Pinzul_random", 
                        type_bk = "Random", #[Random,BwData,BwData_inv,EnvBK]
                        world_pol = co_pol, select_var = F, 
                        sp_range=co_pol,
                        random_features = FALSE, beta.val = c(1:15),
                        n_bk = 10000,
                        Test_n = 20,
                        # Model Selection
                        mod.select = T, n.mods = 10, use.boyce = 0.5
                        )

r.nests$avr.preds %>% plot()
r.nests$params

# For the observations distribution
r.points <- Auto_maxent(presence_dat=Points_SEO.Thin, 
                        predictors=env_analysis, 
                        rm.dp = TRUE,
                        crs.r = crs_p,
                        name.mod = "Pinzul_random", 
                        type_bk = "Random", #[Random,BwData,BwData_inv,EnvBK]
                        world_pol = co_pol, select_var = F, 
                        sp_range=co_pol,
                        random_features = TRUE, beta.val = c(1:15), n.m=1,
                        n_bk = 10000,
                        Test_n = 20,
                        # Model Selection
                        mod.select = T, n.mods = 10, use.boyce = 0.5
)

r.points$avr.preds %>% plot()
r.points$params

# 4.1 Test different parameters for the models ----
models_r <- "Results/Models/Nests" ; models_r %>% dir.create(recursive=TRUE,showWarnings = FALSE)

# 4.1.a Random background sampling----
r.maxent <- Auto_maxent(presence_dat=Nest_SEO.thin, 
                        predictors=env_analysis, 
                        rm.dp = TRUE,
                        crs.r = crs_p,
                        name.mod = "Pinzul_random", 
                        type_bk = "Random", #[Random,BwData,BwData_inv,EnvBK]
                        world_pol = co_pol, select_var = F, 
                        sp_range=co_pol,
                        random_features = FALSE, beta.val = c(1:15),
                        n_bk = 10000,
                        Test_n = 20,
                        # Model Selection
                        mod.select = T, n.mods = 10, use.boyce = 0.5
                        )
# Export the results
# Save the rast.files
route_r <- paste(models_r,"r.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(r.maxent,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(r.maxent[-c(8,9)],paste(route_r,"r.maxent.rds",sep="/"))

# 4.1.b Presence Sampling background sampling----
pres.maxent <- Auto_maxent(presence_dat=Nest_SEO.thin, 
                           predictors=env_analysis, 
                           rm.dp = TRUE,
                           crs.r = crs_p,
                           name.mod = "Pinzul_pres", 
                           type_bk = "BwData", #[Random,BwData,BwData_inv,EnvBK]
                           world_pol = co_pol, select_var = F, 
                           sp_range=co_pol,
                           random_features = FALSE, beta.val = c(1:15), 
                           n_bk = 10000,
                           Test_n = 20,
                           # Model Selection
                           mod.select = T, n.mods = 10, use.boyce = 0.5
)

# Save the rast.files
route_r <- paste(models_r,"pres.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(pres.maxent,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(pres.maxent[-c(8,9)],paste(route_r,"pres.maxent.rds",sep="/"))

# 4.1.c Presence Sampling background sampling----
pres.inv.maxent <- Auto_maxent(presence_dat=Nest_SEO.thin, 
                               predictors=env_analysis, 
                               rm.dp = TRUE,
                               crs.r = crs_p,
                               name.mod = "Pinzul_presInv", 
                               type_bk = "BwData_inv", #[Random,BwData,BwData_inv,EnvBK]
                               world_pol = co_pol, select_var = F, 
                               sp_range=co_pol,
                               random_features = FALSE, beta.val = c(1:15), 
                               n_bk = 10000,
                               Test_n = 20,
                               # Model Selection
                               mod.select = T, n.mods = 10, use.boyce = 0.5
                            )

# Save the rast.files
route_r <- paste(models_r,"pres.inv.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(pres.inv.maxent,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(pres.inv.maxent[-c(8,9)],paste(route_r,"pres.inv.maxent.rds",sep="/"))

# 4.1.d Presence Sampling background sampling----
env.maxent <- Auto_maxent(presence_dat=Nest_SEO.thin, 
                          predictors=env_analysis, 
                          rm.dp = TRUE,
                          crs.r = crs_p,
                          name.mod = "Pinzul_EnvBK", 
                          type_bk = "EnvBK", #[Random,BwData,BwData_inv,EnvBK]
                          world_pol = co_pol, select_var = F, 
                          sp_range=co_pol,
                          random_features = FALSE, beta.val = c(1:15), 
                          n_bk = 10000,
                          Test_n = 20,
                          # Model Selection
                          mod.select = T, n.mods = 10, use.boyce = 0.5
                            )

# Save the rast.files
route_r <- paste(models_r,"env.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(env.maxent,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(env.maxent[-c(8,9)],paste(route_r,"env.maxent.rds",sep="/"))

# 4.2 Models for the sampling records distribution----
models_r <- "Results/Models/Dist_points" ; models_r %>% dir.create(recursive=TRUE,showWarnings = FALSE)

# 4.2.a Random background sampling----
r.maxentP <- Auto_maxent(presence_dat=Points_SEO.Thin, 
                        predictors=env_analysis, 
                        rm.dp = TRUE,
                        crs.r = crs_p,
                        name.mod = "Pinzul_random", 
                        type_bk = "Random", #[Random,BwData,BwData_inv,EnvBK]
                        world_pol = co_pol, select_var = F, 
                        sp_range=co_pol,
                        random_features = FALSE, beta.val = c(1:15),
                        n_bk = 10000,
                        Test_n = 20,
                        # Model Selection
                        mod.select = T, n.mods = 10, use.boyce = 0.5
)
# Export the results
# Save the rast.files
route_r <- paste(models_r,"r.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(r.maxentP,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(r.maxentP[-c(8,9)],paste(route_r,"r.maxent.rds",sep="/"))

# 4.2.b Presence Sampling background sampling----
pres.maxentP <- Auto_maxent(presence_dat=Points_SEO.Thin, 
                           predictors=env_analysis, 
                           rm.dp = TRUE,
                           crs.r = crs_p,
                           name.mod = "Pinzul_pres", 
                           type_bk = "BwData", #[Random,BwData,BwData_inv,EnvBK]
                           world_pol = co_pol, select_var = F, 
                           sp_range=co_pol,
                           random_features = FALSE, beta.val = c(1:15), 
                           n_bk = 10000,
                           Test_n = 20,
                           # Model Selection
                           mod.select = T, n.mods = 10, use.boyce = 0.5
)

# Save the rast.files
route_r <- paste(models_r,"pres.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(pres.maxentP,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(pres.maxentP[-c(8,9)],paste(route_r,"pres.maxent.rds",sep="/"))

# 4.2.c Presence Sampling background sampling----
pres.inv.maxentP <- Auto_maxent(presence_dat=Points_SEO.Thin, 
                               predictors=env_analysis, 
                               rm.dp = TRUE,
                               crs.r = crs_p,
                               name.mod = "Pinzul_presInv", 
                               type_bk = "BwData_inv", #[Random,BwData,BwData_inv,EnvBK]
                               world_pol = co_pol, select_var = F, 
                               sp_range=co_pol,
                               random_features = FALSE, beta.val = c(1:15), 
                               n_bk = 10000,
                               Test_n = 20,
                               # Model Selection
                               mod.select = T, n.mods = 10, use.boyce = 0.5
)

# Save the rast.files
route_r <- paste(models_r,"pres.inv.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(pres.inv.maxentP,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(pres.inv.maxentP[-c(8,9)],paste(route_r,"pres.inv.maxent.rds",sep="/"))

# 4.2.d Presence Sampling background sampling----
env.maxentP <- Auto_maxent(presence_dat=Points_SEO.Thin, 
                          predictors=env_analysis, 
                          rm.dp = TRUE,
                          crs.r = crs_p,
                          name.mod = "Pinzul_EnvBK", 
                          type_bk = "EnvBK", #[Random,BwData,BwData_inv,EnvBK]
                          world_pol = co_pol, select_var = F, 
                          sp_range=co_pol,
                          random_features = FALSE, beta.val = c(1:15), 
                          n_bk = 10000,
                          Test_n = 20,
                          # Model Selection
                          mod.select = T, n.mods = 10, use.boyce = 0.5
)

# Save the rast.files
route_r <- paste(models_r,"env.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(env.maxentP,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(env.maxentP[-c(8,9)],paste(route_r,"env.maxent.rds",sep="/"))

# 6. ENFA analysis, Nests and Distribution points----
# 6.1 Run the ENFA analysis for the Nests----
r.dat <- rast_to_vect(env_analysis)

# Normalize the varibles
n.dat <- apply(r.dat$tab[,-1],2,function(w){
            y<-best_normalization(x=w, # data to normalize
              allow.norm=T,center = T,
                n_cores=detectCores(logical=FALSE)-1)
            return(y$t.values)}
            )

n.dat <- cbind(cell=r.dat$tab$cell,n.dat[,-c(1)])
n.row.dat <- prod(r.dat[["dim"]])

# Prepare the obs index for the analysis
pres_index<-rep(0,times=n.row.dat)

obs <- terra::extract(x=env_analysis, y=Nest_SEO.thin %>% vect(), cells=T)$cell
pres_index[obs]<-1

obs_index <- pres_index[-r.dat$index_missin]

ENFA.r <- ENFA_function(data = n.dat[,-1], # Data.frame containing the environmental information with no NAs
                        presence_index = obs_index)
ENFA.r$marginality

# Transform the ENFA_results into rasters for the export
empty_rast <- env_analysis[[1]] ; empty_rast[-is.na(empty_rast)] <- NA

# Get the different ENFA values
maha <- empty_rast ; maha[r.dat$tab$cell] <- ENFA.r$prediction
Marginality <- empty_rast ; Marginality[r.dat$tab$cell] <- ENFA.r$marginality_specificity_vals$Marginality
Specialization <- empty_rast ; Specialization[r.dat$tab$cell] <- ENFA.r$marginality_specificity_vals$Specialization1

ENFA_rast <- c(maha,Marginality,Specialization) ; names(ENFA_rast) <- c("Mahalanobis_dist","Marginality","Specificity")

# 6.1.b Get the rest of the ENFA parameters
ENFA_extra<- plot_enfa(mar=ENFA.r$marginality_specificity_vals$Marginality, # Marginality vector
                       spc=ENFA.r$marginality_specificity_vals$Specialization1, # Specialization vector
                       m=ENFA.r$niche_centroid_coordinates, # Niche centroid
                       sp_rec=obs_index, # Species records index
                       plot_sp=TRUE, # should we plot the results
                       pts=FALSE)

ENFA_extra$prop.overlap

# 6.1.c Export the results----
ENFA_route <-"./Results/ENFA" ; dir.create(ENFA_route,showWarnings = FALSE,recursive = TRUE)

saveRDS(list(ENFA.r,ENFA_extra),paste(ENFA_route,"ENFA_nest.rds",sep="/"))
writeRaster(ENFA_rast,paste(ENFA_route,paste0("ENFA.nests",".tif"),sep="/"),overwrite=T)

# 6.2 For the Observation of individuals----
# Prepare the obs index for the analysis
pres_index<-rep(0,times=n.row.dat)

obs <- terra::extract(x=env_analysis, y=Points_SEO.Thin %>% vect(), cells=T)$cell
pres_index[obs]<-1

obs_index <- pres_index[-r.dat$index_missin]

ENFA.r2 <- ENFA_function(data = n.dat[,-1], # Data.frame containing the environmental information with no NAs
                        presence_index = obs_index)
ENFA.r2$marginality

# Transform the ENFA_results into rasters for the export
# Get the different ENFA values
maha <- empty_rast ; maha[r.dat$tab$cell] <- ENFA.r2$prediction
Marginality <- empty_rast ; Marginality[r.dat$tab$cell] <- ENFA.r2$marginality_specificity_vals$Marginality
Specialization <- empty_rast ; Specialization[r.dat$tab$cell] <- ENFA.r2$marginality_specificity_vals$Specialization1

ENFA_rast2 <- c(maha,Marginality,Specialization) ; names(ENFA_rast2) <- c("Mahalanobis_dist","Marginality","Specificity")

# 6.1.b Get the rest of the ENFA parameters
ENFA_extra2<- plot_enfa(mar=ENFA.r2$marginality_specificity_vals$Marginality, # Marginality vector
                       spc=ENFA.r2$marginality_specificity_vals$Specialization1, # Specialization vector
                       m=ENFA.r2$niche_centroid_coordinates, # Niche centroid
                       sp_rec=obs_index, # Species records index
                       plot_sp=TRUE, # should we plot the results
                       pts=FALSE)

ENFA_extra2$prop.overlap

# 6.1.c Export the results
saveRDS(list(ENFA.r2,ENFA_extra2),paste(ENFA_route,"ENFA_observations.rds",sep="/"))
writeRaster(ENFA_rast2,paste(ENFA_route,paste0("ENFA.observations",".tif"),sep="/"),overwrite=T)

# End of the script
unlink("junk",recursive = TRUE,force = TRUE)
unlink(td,recursive = TRUE,force=TRUE)