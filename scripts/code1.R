
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

#latest<-nc_open("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\LUH2_v2f_beta_SSP1_RCP2.6_IMAGE\\states.nc")

for (i in 1:12){
	latest<-brick("E:/Downloads/states.nc",varname=var1[i])
	l1<-(subset(latest,1166))
	names(l1)<-paste(var1[i],2015,sep="_")	
	if(i==1){l3<-l1}else{l3<-stack(l3,l1)}
	}
#writeRaster(l3,format="GTiff",file="C:\\soil_layers\\land_2015.tiff")
#nnn<-names(l3);save(nnn,file="C:\\soil_layers\\l3_names.r")
l3<-brick("C:\\soil_layers\\land_2015.tif")
load(file="C:\\soil_layers\\l3_names.r"); names(l3)<-nnn

for (i in 1:12){
	latest<-brick("E:/Downloads/states.nc",varname=var1[i])
	l1<-(subset(latest,966:1166))
	names(l1)<-paste(var1[i],1815:2015,sep="_")	
	if(i==1){l2<-l1}else{l2<-stack(l2,l1)}
	}

l2_names<-read.table(text=names(l2),sep="_",stringsAsFactors=FALSE)
l2_names$col<-names(l2)
names(l2_names)[1:2]<-c("layer","year")

#l3<-calc(l2,sum)

#load("E:\\Dropbox\\R_scripts\\functions5.txt")
alt<-raster("C:\\soil_layers\\alt_2-5m_esri\\alt\\alt\\w001001.adf")
bcp<-stack(list.files("C:\\soil_layers\\bio_2-5m_esri\\",recursive=TRUE,pattern="w001001.adf",full.names=TRUE))
n1<-list.files("C:\\soil_layers\\bio_2-5m_esri\\",recursive=TRUE,pattern="w001001.adf",full.names=TRUE)
n1<-gsub("C:\\soil_layers\\bio_2-5m_esri\\/bio/","",n1,fixed=TRUE)
n1<-gsub("/w001001.adf","",n1)
names(bcp)<-n1
#l4<-resample(l3,invar);gc()
#writeRaster(l4,format="GTiff",file="C:\\soil_layers\\land2015_resampled.tiff")
l4<-brick("C:\\soil_layers\\land2015_resampled.tif")
names(l4)<-nnn

#tt<-rasterize(wrld_simpl,invar)
#writeRaster(tt,format="GTiff",file="C:\\soil_layers\\wrld_simpl.tif")
tt<-raster("C:\\soil_layers\\wrld_simpl.tif")

#files1<-list.files("C:\\temp\\disease_analyses\\",pattern="_points.r",full.names=TRUE)
files2<-list.files("C:\\temp\\disease_analyses\\",pattern="all_points.r",full.names=TRUE)
files2a<-read.table(text=files2,sep="-",stringsAsFactors=F)
files2a$filen<-files2
files2<-files2a[ !duplicated(files2a$V1) ,"filen"]

predictors2<-stack(subset(s2,8:14),alt,bcp,l4)
names(predictors2)<-gsub("_2015","",names(predictors2))
system.time(
pred1<-predict(predictors2,xm$xm)
)

z=1
while (z<length(files2)){
	qq<-1
	load(files2[z])
	name1<-files2[z]	
	name1<-gsub("C:\\temp\\disease_analyses\\","",name1,fixed=TRUE)	
	name1<-gsub("all_points.r","",name1,fixed=TRUE)	
	points1<-pack1$data1$data1
	year1<-as.vector(year(as.Date(points1$date)))
	year1[year1<1815]<-1815
	year1[year1>2015]<-2015
	year1[is.na(year1)]<-1980
	pn<-l2_names[matches(as.numeric(year1),as.numeric(l2_names$year),all.y=F)$y,]
	pn$row<-matches(as.numeric(year1),as.numeric(l2_names$year),all.y=F)$x
	pn$col<-matches(as.numeric(year1),as.numeric(l2_names$year),all.y=F)$y
	lu1<-raster::extract(l2,points1);gc()
	pn$value<-lu1[cbind(pn$row, pn$col)]
	pn<-pn[order(pn$row),]
	mat1<-matrix(ncol=length(pn$row[pn$row==1]),nrow=nrow(points1),data=pn$value,byrow=TRUE)
	colnames(mat1)<-pn$layer[1:length(pn$row[pn$row==1])]
	rownames(mat1)<-year1
	df1<-data.frame(points1@data,as.data.frame(mat1))
	points1@data<-df1
	#(lu1[273,970,drop=FALSE])
	for (w in 1:ncol(mat1)){
		temp1<-rasterize(points1,l4,field=colnames(mat1)[w],mean,na.rm=TRUE)
		present<-subset(l4,paste(colnames(mat1)[w],2015,sep="_"))
		present[!is.na(temp1)]<-temp1[!is.na(temp1)]
		names(present)<-colnames(mat1)[w]
		if(w==1){lu2<-present}else{lu2<-stack(lu2,present)}
	}
	#luX<-raster::extract(lu2,points1);gc()
	#cor(mat1,luX); summary(luX);summary(mat1)
	predictors<-stack(subset(s2,8:14),alt,bcp,lu2)

	xm<-do_MAX(datax=points1,n=50,predictors=predictors,background_n=10000,vary_back=TRUE,rnd_mask=tt,cluster_n=5)
	gc()
	
	if(class(xm)[1]!="list"){qq=qq+1;if(qq>2){z=z+1};next}else{qq=1;
		#save(pack1,file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),"all_points.r",sep=""))
		#write.csv(data.frame(disease=dis,type=hv2$type[z],name1=name1),file=paste("C:\\temp\\disease_analyses\\",name1,"-",ID(),".csv",sep=""))
		#save(data1, file=paste("C:\\temp\\disease_analyses\\",name1,"_points.r",sep=""))
		save(xm, file=paste("C:\\temp\\disease_analyses2\\",name1,"_maxent.r",sep=""))}
	print(z);z=z+1
	gc();removeTmpFiles(5)
	}


	