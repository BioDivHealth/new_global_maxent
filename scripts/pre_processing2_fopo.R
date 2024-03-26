
library(raster)
library(rgdal)


##template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#removeTmpFiles(1)

## forest not forest land
fnf<-raster("X:\\DRtemp\\temp\\fnf_map.txt")
fnf2<-crop(fnf,template)
fnf3<-disaggregate(fnf2,12,method="")
fpop<-extend(fnf3,template)
namx="fnfXXX"

###al possible masks
ff<-list.files("X:\\DRtemp\\temp\\masks\\",full.names=TRUE)
ff2<-list.files("X:\\DRtemp\\temp\\masks\\", pattern="tif.aux.xml",full.names=TRUE)
ff3<-ff[!ff %in% ff2]

##population
#write.csv(futpop,file="X:\\DRtemp\\temp\\futpop.csv")
#futpop<-read.csv(file="X:\\DRtemp\\temp\\futpop.csv",stringsAsFactors=FALSE)
#futpop$filen<-gsub("C:","X:/DRtemp",futpop$filen)
#fpop2<-stack(futpop$filen)
#fpop<-disaggregate(fpop2,3,method="")/9
#fpop<-extend(fpop,template)
#namx="fpopXXX"

##land cover
write.csv(futlg,file="X:\\DRtemp\\temp\\futlc.csv")
futlg<-read.csv(file="X:\\DRtemp\\temp\\futlc.csv",stringsAsFactors=FALSE)
futlg$filen<-gsub("C:","X:/DRtemp",futlg$filen)
lc<-stack(futlg$filen)
lc2<-crop(lc,template)
lc2<-disaggregate(lc2,12,method="")#/144
namx="lcXXX"
#lc2<-resample(lc,template,method="ngb")
#lc2<-resample(lc,template,method="ngb")

##poppop
#gdppop<-stack("X:\\DRtemp\\POP_GDP_SSP13\\gdp_pop_ssp1_3.tif")
#gpt<-read.csv(file="X:\\DRtemp\\POP_GDP_SSP13\\gdp_pop_ssp1_3.csv",stringsAsFactors=FALSE)$x
#gpt$filen<-gsub("C:","X:/DRtemp",gpt$filen)
#names(gdppop)<-gpt
#gdppop2<-crop(gdppop,template)
#gdppop2<-disaggregate(gdppop2,12,method="")/144
#gpt2<-resample(gdppop,template,method="ngb")
namx="gdppopXXX"

#sum(as.numeric(values(subset(gdppop,16))),na.rm=TRUE)
#sum(as.numeric(values(subset(gdppop2,16))),na.rm=TRUE)
#sum(as.numeric(values(subset(,1))),na.rm=TRUE)
#sum(as.numeric(values(subset(,1))),na.rm=TRUE)

for (x in 1:length(ff3)) {

	template2<-raster(ff3[x])
	nam1<-ff3[x]	
	#nam1<-gsub(".tif","",nam1)
	nam1<-gsub("X:\\DRtemp\\temp\\masks\\","",nam1,fixed=TRUE)

	if(paste("X:\\DRtemp\\temp\\input_masked\\",namx,nam1,sep="") %in% list.files("X:\\DRtemp\\temp\\input_masked\\",pattern=namx,full.names=TRUE)){next} 	
	xxx<-crop(fpop,template2)
	writeRaster(xxx,format="GTiff",file=paste("X:\\DRtemp\\temp\\input_masked\\",namx,nam1,sep=""))

	gc()
	
	#removeTmpFiles(1)

	print(x)

	}	


##fnf
for (x in 1:length(ff3)) {
  
  template2<-raster(ff3[x])
  nam1<-ff3[x]	
  #nam1<-gsub(".tif","",nam1)
  nam1<-gsub("X:\\DRtemp\\temp\\masks\\","",nam1,fixed=TRUE)
  
  if(paste("X:\\DRtemp\\temp\\input_masked\\",namx,nam1,sep="") %in% list.files("X:\\DRtemp\\temp\\input_masked\\",pattern=namx,full.names=TRUE)){next} 	
  xxx<-crop(fpop,template2)
  writeRaster(xxx,format="GTiff",file=paste("X:\\DRtemp\\temp\\input_masked\\",namx,nam1,sep=""))
  
  gc()
  
  removeTmpFiles(1)
  
  print(x)
  
}	

