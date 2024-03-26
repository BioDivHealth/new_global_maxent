
library(raster)
library(rgdal)
#library(ncdf4)
#library(lubridate)
#library(grr)
#library(doParallel)
#library(maptools)
#library(dismo)

source("E:\\Dropbox\\R_scripts\\functions6.r")
data(wrld_simpl)



###look up all diseases
#link1<-list.files("C:\\temp\\disease_analyses\\", pattern=".csv",full.names=TRUE)

#for (i in 1:length(link1)){

#	tt<-read.csv(link1[i],stringsAsFactors=FALSE)
#	l2<-gsub("C:\\temp\\disease_analyses\\","",link1[i],fixed=TRUE)
#	l2<-gsub(".csv","",l2,fixed=TRUE)
#	l2<-gsub("_XXX","",l2,fixed=TRUE)
#	l3<-strsplit(l2,"-")[[1]]
#	if(length(l3)==2){tt$id=l3[2];tt$tag=l3[1]} else {tt$id=l3[1];tt$tag=NA}
#	tt<-tt[,c("disease","type","name1","id","tag")]
#	if(i==1){res1<-tt}else{res1<-rbind(res1,tt)}
#	print(i)
#	}
#unique(res1$disease)


####start predictions

##models
sdm1<-list.files("C:\\temp\\disease_analyses3\\", pattern="maxent",full.names=TRUE)
sdm2<-list.files("C:\\temp\\disease_analyses3\\", pattern="maxent",full.names=FALSE)
sdm2<-gsub("_maxent.r","",sdm2,fixed=TRUE)


##ciimates ##CMIP5 projections
fut1<-list.files("C:\\FutureClimateScenarios\\",recursive=TRUE,pattern="w001001.adf",full.names=FALSE)
	fut2<-gsub("s_bio_2_5min",";",fut1,fixed=TRUE)
	fut3z<-read.table(text=fut2,sep=";",stringsAsFactors=FALSE)
	fut3z$V1<-gsub("_rcp",";",fut3z$V1,fixed=TRUE)
	fut4<-read.table(text=fut3z$V1,sep=";",stringsAsFactors=FALSE)
	fut5<-read.table(text=fut4$V2,sep="_",stringsAsFactors=FALSE)
	fut2<-gsub("s_bio_2_5min",";",fut1,fixed=TRUE)
	clims<-data.frame(RCP=paste(fut5$V1,".",fut5$V2,sep=""),year=fut5$V3,model=fut4$V1,filen=fut1)
	clims$rep<-paste(clims$model,clims$year,clims$RCP,sep="_")
	clims$filen<-as.vector(clims$filen)
	bios1<-read.table(text=clims$filen,sep="/",stringsAsFactors=FALSE)
	bios1$V3<-gsub("bio_","X",bios1$V3,fixed=TRUE)
	clims$name<-bios1$V3


yrs<-paste("X",unique(clims$year),sep="")

####land covers this is for CMIP6  not CMIP5
#fut2x<-list.files("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\",recursive=TRUE,pattern="states.nc")
#	fut2b<-gsub("/states.nc","",fut2x,fixed=TRUE)
#	fut3<-read.table(text=fut2b,sep="_",stringsAsFactors=FALSE)
#	names(fut3)<-c("scrap","scrap2","scrap3","SSP","RCP","model")
#	fut3$RCP<-gsub("RCP","",fut3$RCP,fixed=TRUE)
#	fut3$filen<-fut2x
#	fut3$year<-NA
#	land<-fut3
#?focal
#var1<-c("primf","primn","secdf","secdn","urban","c3ann","c4ann","c3per","c4per","c3nfx","pastr","range","secmb","secma")
#for (j in 1:nrow(land)){
#	landj<-land[j,]
#	
#	for (i in 1:12){
#	latest<-brick(paste("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\",landj$filen,sep=""),varname=var1[i])
#	l1<-(subset(latest,yrs))
#	names(l1)<-paste(var1[i],names(l1),landj$RCP,landj$SSP,sep="_")	
#	if(i==1){l3<-l1}else{l3<-stack(l3,l1)}
#	}
#	#landj2<-landj[rep(1,nlayers(l3)) , ]
#	#landj2$layer<-names(l3)
#	if(j==1){l4<-l3;landj3<-landj2}else{l4<-stack(l3,l4); landj3<-rbind(landj3,landj2)}
#	print(j)
#	}	

spdata<-stack("C:\\temp\\predictorsX.tif")
names(spdata)<-read.csv("E:\\Dropbox\\names1.csv",stringsAsFactors=FALSE)$x
spdata<-subset(spdata,c(8,10:length(names(spdata))))
Altitude<-subset(spdata,1)
#names(Altitude)<-"Altitude"

reps2<-unique(clims$rep)
reps2<-c(reps2,"present")

for (x in 1:length(reps2)){


	if(reps2[x]=="present"){predictors2<-spdata}else{
	clims2<-stack(paste("C:\\FutureClimateScenarios\\",clims[clims$rep==reps2[x],"filen"],sep=""))
	names(clims2)<-clims[clims$rep==reps2[x],"name"]
	predictors2<-stack(clims2,Altitude)
	}

	for(k in 1:length(sdm1)){

		load(sdm1[k])

		if(paste(sdm2[k],"_",reps2[x],"_XXX.tif",sep="") %in% list.files("C:\\temp\\resultsY\\",pattern="XXX.tif",full.names=FALSE)){next}
			
			for(j in 1:5){ 
		
			system.time(
			pred1<-predict(predictors2,xms[[j]])
			)

			if(j==1){pred2<-pred1} else {pred2<-stack(pred2,pred1)}
		}

	fin_pred<-calc(pred2,mean,na.rm=TRUE)
	writeRaster(fin_pred,method="GTiff",file=paste("C:\\temp\\resultsY\\",sdm2[k],"_",reps2[x],"_XXX.tif",sep=""))
	gc();removeTmpFiles(1)
	}
}

	