
require(raster)
require(rgdal)
require(dismo)
library(data.table)
library(sp)
library(taxize)
library(maptools)
library(doParallel)

#setwd('/scratch/scratch/ucbtdw0/sdms')

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
values(template)<-1:ncell(template)

#ws_ras<-rasterize(wrld_simpl,template,field="UN")
#writeRaster(ws_ras,format="GTiff",file="X:\\DRtemp\\temp\\wrld_simpl_raster.tif")
##load land area mask
ws_ras<-raster("X:\\DRtemp\\temp\\wrld_simpl_raster.tif")

##mask template by area
template<-mask(template,ws_ras)

source("X:\\DRtemp\\functions6.r")
#source("E:\\Dropbox\\R_scripts\\functions6.r")
#source("C:\\Users\\xxxx\\Dropbox\\R_scripts\\functions6.r")

##data.table function
cbind2 <- function(...) (setattr(do.call(c,list(...)),"class",c("data.table","data.frame")))

years=c(2030,2050,2070,2080,2010)
RCPS<-c(2.6,8.5,6.0,4.5)
lcs<-c("gcrop","gsecd","gpast","gurbn","gothr","gfsh1","gfsh2","gfsh3")

data(wrld_simpl)
wrld_simpl[wrld_simpl$NAME=="Russia","SUBREGION"]<-161

###read disease data
d1<-read.csv("X:\\DRtemp\\disease_table30.csv",stringsAsFactors=FALSE)

###look up all diseases
link1<-list.files("X:\\DRtemp\\temp\\disease_analyses/", pattern=".csv",full.names=TRUE)

###look up all diseases
for (i in 1:length(link1)){
	tt<-read.csv(link1[i],stringsAsFactors=FALSE)
	l2<-gsub("C:/temp/disease_analyses/","",link1[i],fixed=TRUE)
	l2<-gsub(".csv","",l2,fixed=TRUE)
	l2<-gsub("_XXX","",l2,fixed=TRUE)
	l3<-strsplit(l2,"-")[[1]]
	if(length(l3)==2){tt$id=l3[2];tt$tag=l3[1]} else {tt$id=l3[1];tt$tag=NA}
	tt<-tt[,c("disease","type","name1","id","tag")]
	if(i==1){diseasesX<-tt}else{diseasesX<-rbind(diseasesX,tt)}
	print(i)
}
diseases<-diseasesX
diseases$name2<-paste(diseases$name1,";",sep="")
diseases$name2<-gsub("_subs_",";",diseases$name2)
diseases$name2<-gsub(";;",";",diseases$name2)
dd<-read.table(text=diseases$name2,sep=";",stringsAsFactors = F)
dd$V1<-gsub("_1","",dd$V1)
dd$V1<-gsub("_2","",dd$V1)
dd$V1<-gsub("_3","",dd$V1)
dd$V1<-gsub("_4","",dd$V1)
names(dd)<-c("species","subregions")
diseases<-cbind(diseases,dd)
diseases$accu<-nchar(diseases$subregions)
diseases$ID<-paste(diseases$disease,diseases$species,diseases$subregions,sep=";")
diseases$ID2<-paste(diseases$disease,diseases$species,sep=";")
rm(dd,diseasesX,tt,l2,l3,link1)
dis1<-unique(diseases$disease)


#remove duplications
#diseases[diseases$ID=="rocio;accipiter_erythrotriorchis;5_13",]
diseases<-diseases[!duplicated(diseases$ID),]
diseases<-diseases[order(rank(diseases$accu)),]
#diseases[diseases$ID2=="tanganya;crocidura theresae_crocidura theresae",]
diseases<-diseases[!duplicated(diseases$ID2),]

##get all points data for each
points1a<-list.files("X:\\DRtemp\\temp\\points\\", pattern="all_points",full.names=TRUE)
points1<-list.files("X:\\DRtemp\\temp\\points\\", pattern="_points",full.names=TRUE)
points2b<-points1[!points1 %in% points1a]
points1a<-list.files("X:\\DRtemp\\temp\\points\\", pattern="all_points",full.names=FALSE)
points1<-list.files("X:\\DRtemp\\temp\\points\\", pattern="_points",full.names=FALSE)
points2<-points1[!points1 %in% points1a]
points2<-gsub("_points.r","",points2, fixed=TRUE)
points3<-data.frame(species=points2,filen=points2b)
rm(points1a,points1,points2b,points2)

#load gridded livestock of the world
#load("C:\\temp\\livestock_future_2030_2050_2070_2080.r")
#fwrite(lt2,file="C:\\temp\\livestock_future_2030_2050_2070_2080.csv")
#lt2<-fread("X:/DRtemp/temp\\livestock_future_2030_2050_2070_2080.csv")
#names(lt2)[1:6]<-paste(names(lt2)[1:6],"2010",sep="")
#names(lt2)<-gsub("2","_2",names(lt2))
#lt3<-reshape(lt2, varying = c(1:6,8:ncol(lt2)), sep = "_", direction = 'long')
#lt3[,cell_by_year:=paste(cell.id,time,sep="_")]
#fwrite(lt3,file="X:/DRtemp/temp\\livestock_future_2030_2050_2070_2080b.csv")
lt2<-fread("X:/DRtemp/temp\\livestock_future_2030_2050_2070_2080b.csv")
setkey(lt2,cell_by_year)

###MAKE THE SAME AS ABOVE FOR GDP AND FUTPOP
## use coorindates of year and cell.id to make just a couple of columns pop1, pop2, gdp
## melt and then match by cell.id_year (BUT by SSP)
##get points
#pt1X<-xyFromCell(template,na.omit(values(template)),spatial=TRUE)

##extract people   
#futpopx<-stack(list.files("X:\\DRtemp\\Spatial_population_scenarios_GeoTIFF\\GeoTIFF\\",pattern="_2",recursive=TRUE,full.names=TRUE)[seq(2,250,by=5)])
#futpopx<-subset(futpopx,c(c(1,3,5,7,8),c(1,3,5,7,8)+10,c(1,3,5,7,8)+20,c(1,3,5,7,8)+30,c(1,3,5,7,8)+40))
#futpopx<-extend(futpopx,template)
#futpopx<-disaggregate(futpopx,3,method ='')/9
#ft2<-setDT(as.data.frame(values(futpopx)))
#rem1<-rowSums(ft2,na.rm=TRUE)
#ft2[,cell.id:=1:ncell(template)]
#head(ft2)
#ft2<-ft2[!is.na(rem1) & rem1>0,] # NOTE THIS WILL REMOVE ALL CLIMS THAT ### PROLLY should change


#ft2<-setDT(extract(futpopx,pt1X,method="bilinear",na.rm=TRUE,df=TRUE))
#ft2[,cell.id:=na.omit(values(template))]
#ft2b<-reshape(ft2, varying = 1:25, sep = "_", direction = 'long')
#ft2b[,cell_by_year:=paste(cell.id,time,sep="_")]
#fwrite(ft3,file="X:/DRtemp/temp\\future_pop_maskedb.csv")
##ft2b<-fread("X:/DRtemp/temp\\future_pop_maskedb.csv")
#ft3<-ft2b[ft2b$time==2010,c("cell.id","ssp1")]
#names(ft3)[2]<-"humans2010"
#setkey(ft2b,cell.id)
#ft4<-ft2b[ft2b$time!=2010,]
#setkey(ft4,cell.id)
#ft2<-ft4[ft3]
#ft2<-ft2[, c("cell.id","time","ssp1","ssp2","ssp3","ssp4","ssp5","cell_by_year","humans2010")]
#fwrite(ft2,file="X:/DRtemp/temp\\future_pop_maskedb.csv")
#ft2<-fread("X:/DRtemp/temp\\future_pop_maskedb.csv")
#setkey(ft2,cell_by_year) ##create cell.id by code??


##futpopgdp
#GP2<-setDT(extract(gdppop,pt1X,method="bilinear",na.rm=TRUE,df=TRUE))
#GP2[,cell.id:=na.omit(values(template))]
#names(lt2)[1:6]<-paste(names(lt2)[1:6],"2010",sep="")
#names(GP2)<-gsub("ssp","_ssp",names(GP2))
#names(GP2)<-gsub("20","_20",names(GP2))
#names(GP2)[1]<-"ID__"
#names(GP2)[32]<-"cell.id__"
#anm<-read.table(text=names(GP2),sep="_",stringsAsFactors = F)
#names(GP2)<-paste(anm$V1,"-",anm$V3,"_",anm$V2,sep="")
#names(GP2)<-gsub("-_NA","",names(GP2))
#GP3<-reshape(GP2, varying = 2:31, sep = "_", direction = 'long')
#GP3[,cell_by_year:=paste(cell.id,time,sep="_")]
#fwrite(GP3,file="X:/DRtemp/temp\\future_gdppop_maskedb.csv")
#GP2<-fread("X:/DRtemp/temp\\future_gdppop_maskedb.csv")
#GP3<-GP2[GP2$time==2010,c("cell.id","gdp-ssp1","pop-ssp1")]
#names(GP3)[2:3]<-c("gdp2010","pop2010")
#setkey(GP3,cell.id)
#GP4<-GP2[GP2$time!=2010,]
#setkey(GP4,cell.id)
#GP2<-GP4[GP3]
#GP2<-GP2[, c("cell.id","time","gdp-ssp1","gdp-ssp2","gdp-ssp3","pop-ssp1","pop-ssp2","pop-ssp3","cell_by_year","gdp2010","pop2010")]
#fwrite(GP2,file="X:/DRtemp/temp\\future_gdppop_maskedb.csv")
#GP2<-fread("X:/DRtemp/temp\\future_gdppop_maskedb.csv")
#setkey(GP2,cell_by_year) ##create cell.id by code??

####Get all future climate niches
link2b<-list.files("X:/DRtemp/resultsY", pattern=".tif",full.names=TRUE)
link2<-list.files("X:/DRtemp/resultsY", pattern=".tif",full.names=FALSE)
link2<-gsub("_ncc_",";",link2, fixed=TRUE)
link2<-gsub("_mri_",";",link2, fixed=TRUE)
link2<-gsub("_bcc_",";",link2, fixed=TRUE)
link2<-gsub("_nimr_",";",link2, fixed=TRUE)
link2<-gsub("_cccma_",";",link2, fixed=TRUE)
link2<-gsub("_cesm1_",";",link2, fixed=TRUE)
link2<-gsub("_csiro_",";",link2, fixed=TRUE)
link2<-gsub("_ec_",";",link2, fixed=TRUE)
link2<-gsub("_fio_",";",link2, fixed=TRUE)
link2<-gsub("_gfdl_",";",link2, fixed=TRUE)
link2<-gsub("_bnu_",";",link2, fixed=TRUE)
link2<-gsub("_giss_",";",link2, fixed=TRUE)
link2<-gsub("_ipsl_",";",link2, fixed=TRUE)
link2<-gsub("_inm_",";",link2, fixed=TRUE)
link2<-gsub("_miroc_",";",link2, fixed=TRUE)
link2<-gsub("_lasg_",";",link2, fixed=TRUE)
link2<-gsub("_ncar_",";",link2, fixed=TRUE)
link2<-gsub("_mpi_",";",link2, fixed=TRUE)
link2<-gsub("_mohr_",";",link2, fixed=TRUE)
link2<-gsub("_mri_",";",link2, fixed=TRUE)
link2<-gsub("_mohc_",";",link2, fixed=TRUE)
link2<-gsub("_present_",";present;2010;",link2, fixed=TRUE) ##check when in
link2<-gsub("_XXX.tif","",link2, fixed=TRUE)
link2<-gsub("_2.6",";2.6",link2, fixed=TRUE)
link2<-gsub("_4.5",";4.5",link2, fixed=TRUE)
link2<-gsub("_6.0",";6.0",link2, fixed=TRUE)
link2<-gsub("_8.5",";8.5",link2, fixed=TRUE)
link2<-gsub("_2030",";2030",link2, fixed=TRUE)
link2<-gsub("_2050",";2050",link2, fixed=TRUE)
link2<-gsub("_2070",";2070",link2, fixed=TRUE)
link2<-gsub("_2080",";2080",link2, fixed=TRUE)
link2<-gsub("-astOJIRBNlQ0","",link2, fixed=TRUE)
link2<-gsub("_points.r","",link2, fixed=TRUE)
link3<-read.table(text=link2,sep=";",stringsAsFactors=FALSE)
link3$filen<-link2b
names(link3)[1:4]<-c("species","model","year","RCP")
link3$RCP<-as.numeric(gsub("XXX.tif","999",link3$RCP))
###make other species regions available
link3$speciesX<-link3$species
link3$speciesX<-gsub("_1_subs_",";",link3$speciesX, fixed=TRUE)
link3$speciesX<-gsub("_2_subs_",";",link3$speciesX, fixed=TRUE)
link3$speciesX<-gsub("_3_subs_",";",link3$speciesX, fixed=TRUE)
link3$speciesX<-gsub("_4_subs_",";",link3$speciesX, fixed=TRUE)
others1<-read.table(text=link3$speciesX,sep=";",stringsAsFactors=FALSE)
names(others1)<-c("species2","region")
link3<-cbind(link3,others1)
rm(link2,link2b)

link3$uni1<-paste(link3$year,link3$RCP,sep="_")

sss<-unique(link3$species)

#beginCluster(8)

##randomise
sss<-sample(sss,length(sss),replace=FALSE)

for (hh in 1:sss){
  
  l3<-link3[link3$year>2010 & link3$species==sss[hh],]
  sss2<-unique(l3$uni1)
  
           for (jj in 1:length(sss2)) {
            
             l4<-l3[l3$uni1==sss2[jj],]
             
             if(paste("V:\\ResultsY2\\",l4$species[1],"_XXX_",sss2[jj],".tif",sep="") %in% list.files("V:\\ResultsY2\\",full.names=TRUE)){next}
             
              
          tt<-stack(l4[,"filen"])
          
          if(nrow(l4)==1){
            
            names(tt)<-paste(l4$species[1],"_XXX_",sss2[jj],sep="")
            writeRaster(tt,format="GTiff",file=paste("V:\\ResultsY2\\",l4$species[1],"_XXX_",sss2[jj],".tif",sep=""))
            next
           }        
          
          #ras.mean <- clusterR(tt, calc, args=list(mean, na.rm=T))
          ras.mean <- calc(tt,mean, na.rm=TRUE)
          
          names(ras.mean)<-paste(l4$species[1],"_XXX_",sss2[jj],sep="")
          
          #ras.sd <- clusterR(tt, calc, args=list(sd, na.rm=T))
  
          writeRaster(ras.mean,format="GTiff",file=paste("V:\\ResultsY2\\",l4$species[1],"_XXX_",sss2[jj],".tif",sep=""))
          rm(l4)
           }
  print(paste(hh," of ",length(sss),sep=""))
  
}

#endCluster()

