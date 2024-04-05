# Install required libraries ----
library(raster)
library(rgdal)
#library(ncdf4)
#library(lubridate)
#library(grr)
library(doParallel)
library(maptools) # deprecated in 2023
library(terra)

# Load in functions and map data ---- 
# load in helper functions 
source("scripts/functions/functions6.r")

# load in wrld_simpl tif file
data(wrld_simpl)
data("wrld_simpl", package = "terra")   
tt <- raster("X:\\DRtemp\\wrld_simpl.tif")
plot(tt)
world_shp <- st_read("data/wrld_simpl/wrld_simpl3.shp")
plot(world_shp)

files2<-list.files("C:\\temp\\da2\\",pattern="_points.r",full.names=TRUE)
files2a<-list.files("C:\\temp\\da2\\",pattern="_points.r",full.names=FALSE)
files2a<-gsub("_points.r","",files2a)

# Get current bioclim data ----
spdata<-stack("C:\\temp\\predictorsX.tif")
names(spdata)<-read.csv("E:\\Dropbox\\names1.csv",stringsAsFactors=FALSE)$x
spdata<-subset(spdata,c(8,10:length(names(spdata))))
Altitude<-subset(spdata,1)
#names(Altitude)<-"Altitude"


z=426
while (z<length(files2)){
  qq<-1
  load(files2[z])
  name1<-files2a[z]	
  
  
  
  xm<-do_MAX(datax=data1$data1,n=15,predictors=spdata,background_n=10000,vary_back=TRUE,rnd_mask=tt,cluster_n=5)
  gc()
  
  if(class(xm)[1]!="list"){qq=qq+1;if(qq>2){z=z+1};next}else{qq=1;
  #save(pack1,file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),"all_points.r",sep=""))
  #write.csv(data.frame(disease=dis,type=hv2$type[z],name1=name1),file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),".csv",sep=""))
  #save(data1, file=paste("C:\\temp\\disease_analyses\\",name1,"_points.r",sep=""))
  save(xm, file=paste("/Volumes/OS/Users/xxxx/Dropbox/New_Global_MAXENT\disease_analyses2",name1,"_maxent.r",sep=""))}
  print(z);z=z+1
  gc();removeTmpFiles(5)
}

#xms2<-unlist(unlist(xm))
#xm3<-xms2[seq(from=1,to=length(xms2),by=2)]
#auc2<-unlist(xms2[seq(from=2,to=length(xms2),by=2),drop=TRUE])
#xms<-xm3[(1:length(xm3))[rank(1-auc2)<6]]