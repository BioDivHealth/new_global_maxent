
# Install required libraries ----
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
library(Hmisc)

# Load libraries ----
#library(parallel)
#cl <- makeCluster(2)
#stopCluster(cl)

# Read in data ----

source("scripts/functions/find_synonyms6.R")

conts <- stack(list.files("data/conts/",pattern=".tif",full.names=T))

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

spec1 = read.csv("data/disease_table28.csv",stringsAsFactors=FALSE)
spec1b<-spec1#[spec1$run==1,]

load(file="scripts/functions/wrld_simpl2.R")
load(file="scripts/functions/regions.r")
load(file="scripts/functions/subregions.r")
source(file="scripts/functions/functions7.R")

#xxx1<-list.files("E:/Documents/climate_data",full.names=T,pattern=".tif")

#xxx<-gsub("bi",";",xxx1)
#xxx<-gsub(".tif",";",xxx)
#xxx<-gsub("_data/",";",xxx)

#yyy<-read.table(text=xxx,sep=";",stringsAsFactors=FALSE)[,c(2,3)]
#yyy$model<-substr(yyy$V2, 1, 2)
#yyy$rcp<-substr(yyy$V2, 3, 4)
#yyy$year<-substr(yyy$V3, 1, 2)
#yyy$layer<-substr(yyy$V3, 3, 4)
#yyy$filen<-xxx1
#yyy$replicate<-paste(yyy$model,yyy$rcp,yyy$year)

#lcf<-stack(list.files("E:\\Documents\\landcover_data\\",full.names=T,pattern=".tif"))
#names(lcf)<- c("modis2050max","modis2050min","modis2050norm","modis2070max","modis2070min","modis2070norm")
#lc2<-stack(list.files("E:\\Documents\\landcover_pres\\",full.names=T,pattern=".tif"))
#names(lc2)<-"present"
#lc<-stack(lc2,lcf)
#lc<-mask(lc,spdata$X1)
#writeRaster(lc,format="GTiff",file="E:\\Documents\\landcover_data\\lc.tif")
#lc<-stack("E:\\Documents\\landcover_data\\lc.tif")
#names(lc)<- c("present","modis2050max","modis2050min","modis2050norm","modis2070max","modis2070min","modis2070norm")

#zzz<-data.frame(layer=names(lc),year=c(15,50,50,50,70,70,70),rcp=c(99,85,45,60,85,45,60))

#invar<-stack(list.files("E:\\Documents\\invariant_data\\",pattern=".tif",full.names=TRUE))
#names(invar)<-c("DEMSRE1a","INMSRE1a","INSSRE1a","OPISRE1a","SBDHWS1a","SCLHWS1a","SGRHWS1a","SLPSRT1a","SOCHWS1a","SSLHWS1a","SSNHWS1a","TBDHWS1a","TCLHWS1a","TGRHWS1a","TOCHWS1a","TPHHWS1a","TSLHWS1a","TSNHWS1a","TWISRE1a")
#invar<-crop(invar,lc2)
#invar2<-resample(invar,presbi)
#writeRaster(invar2,format="GTiff",file="E:\\Documents\\invariant_data\\invar2.tif")
#invar<-stack(file="C:\\Users\\User\\Documents\\invariant_data\\invar.tif")
#for (i in 1:nlayers(invar)){
#	invar2<-(subset(invar,i))
#	writeRaster(invar2,format="GTiff",paste(file="E:\\Documents\\invariant_data\\invar_",names(invar)[i],".tif",sep=""))
#	}

#invar<-stack(list.files("E:\\Dropbox\\Documents\\desktop_crap\\soil_layers\\",pattern=".tif",full.names=TRUE))
#alt<-raster("E:\\Documents\\invasion_data\\alt\\alt\\w001001.adf")
#alt<-aggregate(alt,fact=6)
#alt<-resample(alt,invar)
#invar2<-stack(subset(invar,2:8),alt)
#writeRaster(invar2,format="GTiff",file="E:\\Documents\\invariant_data\\invar2.tif")
#
##load invariant data
#invar<-stack("E:\\Documents\\invariant_data\\invar2.tif")
#names(invar)<-c("DEMSRE1a","INMSRE1a","INSSRE1a","OPISRE1a","SBDHWS1a","SCLHWS1a","SGRHWS1a","SLPSRT1a","SOCHWS1a","SSLHWS1a","SSNHWS1a","TBDHWS1a","TCLHWS1a","TGRHWS1a","TOCHWS1a","TPHHWS1a","TSLHWS1a","TSNHWS1a","TWISRE1a")
#names(invar)<-c("TBDHWS1a","TCLHWS1a","GRHWS1a","TOCHWS1a","TPHHWS1a","TSLHWS1a","TSNHWS1a","Altitude")
#invar<-resample(invar,lc$present)

###check for errors
#for(i in 1:100){
#invar<-stack("E:\\Documents\\invariant_data\\invar2.tif")
#ttt<-randomPoints(invar,1000,tryf=5)
#invar1<-extract(invar,ttt)
#print(i)}
#source("scripts/functions/Code_S7_remove_autocorrelation.R")
##remove unecessary invariants
#invar2<-autocor(invar1,threshold=0.7)
#invar<-subset(invar,colnames(invar2))
## build predictors in loop
#e <- simpleError("test error")

### bioclim variables
#presbi<-stack(yyy[yyy$year=="15","filen"])
#names(presbi)<-paste("X",yyy[yyy$year=="15","layer"],sep="")
#presbi<-resample(presbi,lc$present)

#spdata<-stack(invar,lc$present,presbi)
#writeRaster(spdata,format="GTiff",file=".//predictorsX.tif",overwrite=TRUE)
#names1<-names(spdata)
#write.csv(names1,file=".//names1.csv")
spdata<-stack("C:\\temp\\predictorsX.tif")
names(spdata)<-read.csv("data/names1.csv",stringsAsFactors=FALSE)$x
		
############## PROCESS SPECIES #################

ID<-function() return(paste(sample(c(0:9, letters, LETTERS),size=12, replace=TRUE),collapse=""))

i=1
for (i in i:(nrow(spec1b))){
#for (i in i:((jj+1)*ff)){
#for (i in ((jj*ff)+1):((jj+1)*ff)){
#for (i in ((jj*ff)+1):nrow(spec1b)){

	spec1<-spec1b[i,]
	countries<-spec1$countries
	subs2<-wrld_simpl2[wrld_simpl2$ISO2 %in% strsplit(countries,",")[[1]],]
	if(is.na(countries)){next}
	dis<-spec1$name
	geogX<-geog(countries,wrld_simpl2,regions,subregions)	

	hosts<-capitalize(strsplit(spec1$all_host,",")[[1]])
	vectors<-capitalize(strsplit(spec1$all_vectors,",")[[1]])
	second<-strsplit(spec1$SECONDARY,",")[[1]]
	hv<-data.frame(species=c(hosts,vectors,second),type=c(rep("hosts",length(hosts)),rep("vectors",length(vectors)),rep("second",length(second))),stringsAsFactors=FALSE)
	hv<-hv[!is.na(hv$species),]
	hv<-hv[!(hv$species=="NA"),]
	dom<-hv[hv$species %in% c(c("human","cattle","ducks","chickens","pigs","sheep","goats"),capitalize(c("human","cattle","ducks","chickens","pigs","sheep","goats"))),]
	hv2<-hv[!hv$species %in% c(c("human","cattle","ducks","chickens","pigs","sheep","goats"),capitalize(c("human","cattle","ducks","chickens","pigs","sheep","goats"))),]

	if(!nrow(dom)==0){for (y in 1:nrow(dom)){
			write.table(data.frame(disease=dis,type=dom$type[y],name1=tolower(dom$species[y])),file=paste("C:\\temp\\disease_analyses\\",ID(),"_XXX.csv",sep=""),sep=",")
			}}

	if(!nrow(hv2)==0){
	z=1
	while(z<=nrow(hv2)){
	#for (z in 1:nrow(hv2)){

		#spdata<-stack(".//predictorsX.tif")
		#names(spdata)<-read.csv(".//names1.csv",stringsAsFactors=FALSE)$x
	
		key1<-find_synonyms(hv2$species[z])
		data1<-process_species(key=key1$Taxa_synonyms,keyX=key1$Genus_synonyms,geog1=geogX,template)
		if(is.na(data1[1])){print(z);z=z+1;next}
		geog1<-geog2(data1,geogX,wrld_simpl2,regions,subregions)	
		rand1<-extend(geog1$r1,spdata)
		name1=tolower(paste(c(data1$key[1],data1$key[length(data1$key)],data1$sens,"SUBS",sort(geog1$sr1)),collapse="_"))
		pack1<-list(data1=data1,geog1=geog1)
		

		#xm<-suppressWarnings(do_MAX(datax=data1$data1,n=1,predictors=spdata,background_n=5000,rnd_mask=rand1))
		
		#if(class(xm)[1]!="MaxEnt"){qq=qq+1;if(qq>2){z=z+1};next}else{qq=1;
		#save(pack1,file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),"all_points.r",sep=""))
		write.csv(data.frame(disease=dis,type=hv2$type[z],name1=name1),file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),".csv",sep=""))
		if(paste("C:\\temp\\disease_analyses\\",name1,"_points.r",sep="") %in% list.files("C:\\temp\\disease_analyses\\",full.names=TRUE)){print(z);z=z+1;next}
		save(data1, file=paste("C:\\temp\\disease_analyses\\",name1,"_points.r",sep=""))
		#save(xm, file=paste("C:\\temp\\disease_analyses\\",name1,"_maxent.r",sep=""))}
		print(z);z=z+1
		gc()
		}}
print("#############################")
	print(i)
print("#############################")
}
