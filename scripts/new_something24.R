jj=0

library(spocc)
library(rgbif)
library(taxize)
library(raster)
library(dismo)
library(rgdal)
library(maptools)
library(doParallel)
library(rgeos)
library(rJava)
library(XML)
library(rgbif)
library(rgdal)


#library(parallel)
#cl <- makeCluster(2)
#stopCluster(cl)

######run script#############

source("E:\\Dropbox\\R_scripts\\functions6.r")

conts<-stack(list.files("E:\\Dropbox\\Public\\conts\\",pattern=".tif",full.names=T))
rand1<-calc(conts,max,na.rm=TRUE)
#template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

##models
points1a<-list.files("C:\\temp\\disease_analyses\\", pattern="all_points",full.names=TRUE)
points1<-list.files("C:\\temp\\disease_analyses\\", pattern="_points",full.names=TRUE)
points1<-points1[!points1 %in% points1a]

points2a<-list.files("C:\\temp\\disease_analyses\\", pattern="all_points",full.names=FALSE)
points2<-list.files("C:\\temp\\disease_analyses\\", pattern="_points",full.names=FALSE)
points2<-points2[!points2 %in% points2a]
points2<-gsub("_points.r","_maxent.r",points2,fixed=TRUE)

spdata<-stack("C:\\temp\\predictorsX.tif")
names(spdata)<-read.csv("E:\\Dropbox\\names1.csv",stringsAsFactors=FALSE)$x
spdata<-subset(spdata,c(8,10:length(names(spdata))))
		
############## PROCESS SPECIES #################
z=2
while(z<=length(points1)){
		
		load(points1[z])

		xm<-suppressWarnings(do_MAX(datax=data1$data1,n=36,predictors=spdata,background_n=10000,cluster_n=6,vary_back=TRUE,rnd_mask=rand1))
		
		auc1<-unlist(unlist(xm)[(seq(from=2,to=length(unlist(unlist(xm))),by=2))])
		
		xms<-unlist(unlist(xm)[(seq(from=1,to=length(unlist(unlist(xm))),by=2))])#[auc1==max(auc1,na.rm=TRUE)]]   )

		if(class(xms[[1]])!="MaxEnt"){qq=qq+1;if(qq>2){z=z+1};next}else{qq=1;
		xms<-xms[order(auc1,decreasing=T)<=5]
		save(xms, file=paste("C:\\temp\\disease_analyses3\\",points2[z],sep=""))}
		print(z);z=z+1;removeTmpFiles(1)
		gc()
		}



	
