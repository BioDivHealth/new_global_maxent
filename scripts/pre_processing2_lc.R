
library(raster)
library(rgdal)


##template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#removeTmpFiles(1)

## forest not forest land
#fnf<-raster("C:\\temp\\fnf_map.txt")
#fnf2<-crop(fnf,template)
#fnf3<-disaggregate(fnf2,12,method="")
#fpop<-extend(fpop,template)

###al possible masks
ff<-list.files("C:\\temp\\masks\\",full.names=TRUE)
ff2<-list.files("C:\\temp\\masks\\", pattern="tif.aux.xml",full.names=TRUE)
ff3<-ff[!ff %in% ff2]

##population
#write.csv(futpop,file="C:\\temp\\futpop.csv")
#futpop<-read.csv(file="C:\\temp\\futpop.csv",stringsAsFactors=FALSE)
#fpop2<-stack(futpop$filen)
#fpop<-disaggregate(fpop2,3,method="")/9
#fpop<-extend(fpop,template)


##land cover
write.csv(futlg,file="C:\\temp\\futlg.csv")
futlg<-read.csv(file="C:\\temp\\futlg.csv",stringsAsFactors=FALSE)
lc<-stack(futlg$filen)
lc2<-crop(lc,template)
lc2<-disaggregate(lc2,12,method="")#/144
#lc2<-resample(lc,template,method="ngb")

##poppop
#gdppop<-stack("C:\\POP_GDP_SSP13\\gdp_pop_ssp1_3.tif")
#gpt<-read.csv(file="C:\\POP_GDP_SSP13\\gdp_pop_ssp1_3.csv",stringsAsFactors=FALSE)$x
#names(gdppop)<-gpt
#gdppop2<-crop(gdppop,template)
#gdppop2<-disaggregate(gdppop2,12,method="")/144
#gpt2<-resample(gdppop,template,method="ngb")

#sum(as.numeric(values(subset(gdppop,16))),na.rm=TRUE)
#sum(as.numeric(values(subset(gdppop2,16))),na.rm=TRUE)
#sum(as.numeric(values(subset(,1))),na.rm=TRUE)
#sum(as.numeric(values(subset(,1))),na.rm=TRUE)

#rrr<-sample(1:length(ff3),length(ff3),replace=FALSE)

for (x in length(ff3):1) {

	template2<-raster(ff3[x])
	nam1<-ff3[x]	
	#nam1<-gsub(".tif","",nam1)
	nam1<-gsub("C:\\temp\\masks\\","",nam1,fixed=TRUE)

	xxx<-crop(lc2,template2)

	writeRaster(xxx,format="GTiff",file=paste("C:\\temp\\input_masked\\lcXXX",nam1,sep=""))
	
	gc()
	
	#removeTmpFiles(1)

	print(x)

	}	


