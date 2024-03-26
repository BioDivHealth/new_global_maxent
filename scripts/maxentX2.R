
library(raster)
library(rgdal)
library(ncdf4)
library(lubridate)
library(grr)
library(doParallel)
library(maptools)

source("C:\\Users\\xxxx\\Dropbox\\R_scripts\\functions6.r")
data(wrld_simpl)

tt<-raster("X:\\DRtemp\\wrld_simpl.tif")



files2<-list.files("X:\\DRtemp\\da2\\",pattern="_points.r",full.names=TRUE)
files2a<-list.files("X:\\DRtemp\\da2\\",pattern="_points.r",full.names=FALSE)

###get current bioclim data
spdata<-stack("X:\\DRtemp\\predictorsX.tif")
names(spdata)<-read.csv("C:\\Users\\xxxx\\Dropbox\\names1.csv",stringsAsFactors=FALSE)$x
spdata<-subset(spdata,c(8,10:length(names(spdata))))
Altitude<-subset(spdata,1)
#names(Altitude)<-"Altitude"


z=693#length(files2)
while (z!=1){
	qq<-1
	load(files2[z])
	name1<-files2a[z]	



	xm<-do_MAX(datax=data1$data1,n=15,predictors=spdata,background_n=10000,vary_back=TRUE,rnd_mask=tt,cluster_n=15)
	gc()
	
	if(class(xm)[1]!="list"){qq=qq+1;if(qq>2){z=z+1};next}else{qq=1;
		#save(pack1,file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),"all_points.r",sep=""))
		#write.csv(data.frame(disease=dis,type=hv2$type[z],name1=name1),file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),".csv",sep=""))
		#save(data1, file=paste("C:\\temp\\disease_analyses\\",name1,"_points.r",sep=""))
		save(xm, file=paste("D:\\disease_analyses2\\",name1,"_maxent.r",sep=""))}
	print(z);z=z-1
	gc();removeTmpFiles(5)
	}


	