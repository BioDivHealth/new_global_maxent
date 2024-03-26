
library(raster)
library(rgdal)
library(ncdf4)
library(lubridate)
library(grr)
library(doParallel)
library(maptools)

source("E:\\Dropbox\\R_scripts\\functions6.r")
data(wrld_simpl)



###look up all diseases
link1<-list.files("C:\\temp\\disease_analyses\\", pattern=".csv",full.names=TRUE)

for (i in 1:length(link1)){

	tt<-read.csv(link1[i],stringsAsFactors=FALSE)
	l2<-gsub("C:\\temp\\disease_analyses\\","",link1[i],fixed=TRUE)
	l2<-gsub(".csv","",l2,fixed=TRUE)
	l2<-gsub("_XXX","",l2,fixed=TRUE)
	l3<-strsplit(l2,"-")[[1]]
	if(length(l3)==2){tt$id=l3[2];tt$tag=l3[1]} else {tt$id=l3[1];tt$tag=NA}
	tt<-tt[,c("disease","type","name1","id","tag")]
	if(i==1){res1<-tt}else{res1<-rbind(res1,tt)}
	print(i)
	}
unique(res1$disease)


####start predictions

##models
sdm1<-list.files("C:\\temp\\disease_analyses2\\", pattern="maxent",full.names=TRUE)

##ciimates
fut1<-list.files("C:\\FutureClimateScenarios\\",recursive=TRUE,pattern="w001001.adf",full.names=FALSE)
	fut2<-gsub("s_bio_2_5min",";",fut1,fixed=TRUE)
	fut3z<-read.table(text=fut2,sep=";",stringsAsFactors=FALSE)
	fut3z$V1<-gsub("_rcp",";",fut3z$V1,fixed=TRUE)
	fut4<-read.table(text=fut3z$V1,sep=";",stringsAsFactors=FALSE)
	fut5<-read.table(text=fut4$V2,sep="_",stringsAsFactors=FALSE)
	fut2<-gsub("s_bio_2_5min",";",fut1,fixed=TRUE)
	clims<-data.frame(RCP=paste(fut5$V1,".",fut5$V2,sep=""),year=fut5$V3,model=fut4$V1,filen=fut1)

yrs<-paste("X",unique(clims$year),sep="")
```
###land covers
fut2x<-list.files("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\",recursive=TRUE,pattern="states.nc")
	fut2b<-gsub("/states.nc","",fut2x,fixed=TRUE)
	fut3<-read.table(text=fut2b,sep="_",stringsAsFactors=FALSE)
	names(fut3)<-c("scrap","scrap2","scrap3","SSP","RCP","model")
	fut3$RCP<-gsub("RCP","",fut3$RCP,fixed=TRUE)
	fut3$filen<-fut2x
	fut3$year<-NA
	land<-fut3

var1<-c("primf","primn","secdf","secdn","urban","c3ann","c4ann","c3per","c4per","c3nfx","pastr","range","secmb","secma")

for (j in 1:nrow(land)){

	landj<-land[j,]
	
	for (i in 1:12){
	latest<-brick(paste("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\",landj$filen,sep=""),varname=var1[i])
	l1<-(subset(latest,yrs))
	names(l1)<-paste(var1[i],names(l1),landj$RCP,sep="_")	
	if(i==1){l3<-l1}else{l3<-stack(l3,l1)}
	}
	landj2<-landj[rep(1,nlayers(l3)) , ]
	landj2$layer<-names(l3)
	
	if(j==1){l4<-l3;landj3<-landj2}else{l4<-stack(l3,l4); landj3<-rbind(landj3,landj2)}
	print(j)
	}	

writeRaster(l4,format="GTiff",file="C:\\land_future.tiff")

s2<-brick("C:\\soil_layers\\soil_stack2.tif")
load(file="C:\\soil_layers\\s2_names.r"); names(s2)<-nnn2

alt<-raster("C:\\soil_layers\\alt_2-5m_esri\\alt\\alt\\w001001.adf")

###bioclim
bcp<-
names(bcp)<-n1


l4<-subset
names(l4)<

#tt<-rasterize(wrld_simpl,invar)
#writeRaster(tt,format="GTiff",file="C:\\soil_layers\\wrld_simpl.tif")
tt<-raster("C:\\soil_layers\\wrld_simpl.tif")

predictors2<-stack(subset(s2,8:14),alt,bcp,l4)
names(predictors2)<-gsub("_2015","",names(predictors2))

pred1<-predict(predictors2,xm$xm)


	if(class(xm)[1]!="list"){qq=qq+1;if(qq>2){z=z+1};next}else{qq=1;
		#save(pack1,file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),"all_points.r",sep=""))
		#write.csv(data.frame(disease=dis,type=hv2$type[z],name1=name1),file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),".csv",sep=""))
		#save(data1, file=paste("C:\\temp\\disease_analyses\\",name1,"_points.r",sep=""))
		save(xm, file=paste("C:\\temp\\disease_analyses2\\",name1,"_maxent.r",sep=""))}
	print(z);z=z+1
	gc();removeTmpFiles(5)
	}


	