
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
d1<-read.csv("X:\\DRtemp\\disease_table28.csv",stringsAsFactors=FALSE)

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
#ft2<-setDT(extract(futpopx,pt1X,method="bilinear",na.rm=TRUE,df=TRUE))
#ft2[,cell.id:=na.omit(values(template))]
#ft3<-reshape(ft2, varying = 2:26, sep = "_", direction = 'long')
#ft3[,cell_by_year:=paste(cell.id,time,sep="_")]
#fwrite(ft3,file="X:/DRtemp/temp\\future_pop_maskedb.csv")
ft2<-fread("X:/DRtemp/temp\\future_pop_maskedb.csv")
setkey(ft2,cell_by_year) ##create cell.id by code??

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
GP2<-fread("X:/DRtemp/temp\\future_gdppop_maskedb.csv")
setkey(GP2,cell_by_year)

####Get all future climate niches
link2b<-list.files("X:/DRtemp/resultsY", pattern=".tif",full.names=TRUE)
link2<-list.files("X:/DRtemp/resultsY", pattern=".tif",full.names=FALSE)[1:length(link2b)]
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

###lc
#futlgx2<-list.files("X:\\DRtemp\\temp\\input_masked\\",recursive=TRUE,pattern="lcXXX",full.names=TRUE)
#futlgx1<-futlgx2
#futlgx1<-gsub("XXX",";",futlgx1)
#futlgx1<-gsub(".tif","",futlgx1)
#futlgx1<-read.table(text=futlgx1,sep=";",stringsAsFactors=FALSE)
#futlgx1$filen<-futlgx2
#rm(futlgx2)
#futlg<-read.csv(file="X:\\DRtemp\\temp\\futlg.csv",stringsAsFactors=FALSE)

###
futlg<-read.csv(file="X:\\DRtemp\\temp\\futlg.csv",stringsAsFactors=FALSE)
futlg$filen<-gsub("C:","X:/DRtemp/",futlg$filen,fixed=TRUE)
#lc<-raster::stack(futlg$filen)

###fnf
#fnfx2<-list.files("X:\\DRtemp\\temp\\input_masked\\",recursive=TRUE,pattern="fnfXXX",full.names=TRUE)
#fnfx1<-fnfx2
#fnfx1<-gsub("XXX",";",fnfx1)
#fnfx1<-gsub(".tif","",fnfx1)
#fnfx1<-read.table(text=fnfx1,sep=";",stringsAsFactors=FALSE)
#fnfx1$filen<-fnfx2
#rm(fnfx2)

## forest not forest land
fnf<-raster("X:\\DRtemp\\temp\\fnf_map.txt")
#fnf2<-crop(fnf,template)
#fnf3<-disaggregate(fnf2,12,method="")
#rm(fnf2)
#fpop<-extend(fnf3,template)
#namx="fnfXXX"

###water/ice fraction in cell ## doesnt change with age unless do 1- all others.
##icewater<-raster("C:\\temp\\gicew.1700.txt")

##only keep need land-use & years
## gfsh1 wood harvested from secondary mature forest
## gfsh1 wood harvested from secondary young forest
## gfsh3 wood harvested from secondary non-forest
#futlg<-futlf[futlf$variable %in% c("gcrop","gsecd","gpast","gurbn","gothr","gfsh1","gfsh2","gfsh3")& futlf$year %in% years,]

###read in habitat preferences
###convert to CMIP5 land classes

##ecological data about species
spec_hab<-read.csv("X:\\DRtemp\\habitat_both_fin.csv",stringsAsFactors=FALSE)
spec_hab$Snow.and.ice[spec_hab$Snow.and.ice==0]<-NA
spec_hab$f<-(spec_hab$Evergreen.Needleleaf.forest+spec_hab$Evergreen.Broadleaf.forest+spec_hab$Deciduous.Needleleaf.forest+spec_hab$Deciduous.Broadleaf.forest+spec_hab$Mixed.forest)/5
spec_hab$nf<-(spec_hab$Closed.shrublands+spec_hab$Open.shrublands+spec_hab$Woody.savannas+spec_hab$Savannas+spec_hab$Grasslands)/5
spec_hab$water<-rowMeans(spec_hab[,c("Water","Permanent.wetlands","Snow.and.ice")],na.rm=TRUE)

spec_dist<-read.csv("X:\\DRtemp\\host_density_d_distance.csv",stringsAsFactors=FALSE)
spec_desp<-read.csv("X:\\DRtemp\\dispersal_distance_all.csv",stringsAsFactors=FALSE)
spec_desp$binom<-gsub("_"," ",spec_desp$binom)

###get current bioclim data
#spdata<-stack("X:\\DRtemp\\predictorsX.tif")
#names(spdata)<-read.csv("X:\\DRtemp\\names1.csv",stringsAsFactors=FALSE)$x
#spdata<-subset(spdata,c(8,10:length(names(spdata))))
#Altitude<-subset(spdata,1)
#names(Altitude)<-"Altitude"

#res1[res1$tag=="aedes communis_ochlerotatus communis sp_1_subs_21_29",]
#jamestown canyon x=22
# yellow x=15
#cache valley 
x=18

for (x in 1:length(dis1)){

	##choose one disease
	dis2<-diseases[diseases$disease==dis1[x],]
	if(length(dis2$tag[!is.na(dis2$tag)])==0){next}
	#print(dis2)
	#x=x+1
	
	##endemic region
	dis_trans<-d1[d1$name==dis1[x],]
	countr<-wrld_simpl[wrld_simpl$ISO2 %in% strsplit(dis_trans$countries,",")[[1]],]
	dis2$endemic=paste(sort(unique(countr$SUBREGION)),collapse="_")
	dis2$endemic2=paste(sort(unique(countr$REGION)),collapse="_")

	### find overlap of scenarios!!
	linkX<-link3[link3$species %in% dis2$name1,]
	linkX$codeX<-paste(linkX$year,linkX$RCP,sep="_")
	setDT(linkX)
	setDT(dis2)
	setkey(linkX,species)
	setkey(dis2,name1)
	linkX<-linkX[dis2]
	linkX<-linkX[!is.na(linkX$codeX),]
	lv<-unique(linkX$codeX[linkX$type=="vectors"])
	lh<-unique(linkX$codeX[linkX$type=="hosts"])
	if(length(lv[!is.na(lv)])==0 | length(lh[!is.na(lh)])==0){
	  bothmod<-unique(linkX$codeX)
	  }else{
	bothmod<-intersect(lh,lv)
	 	  }
	
	if(length(bothmod)==0){print("ERROR NO SCENARIO OVERLAP");next}
	
	res4<-NULL
	res5<-NULL	
    #system.time(
		##choose one host/vector
    st2<-Sys.time()
    
		for(z in 1:nrow(dis2)){

			dis3<-dis2[z,] #### SORT OUT MULTIPLIE FILES PROBLEM  ### neeed more masks
				
			### remove GLW dat			
			if(dis3$name1 %in% c("ducks","chickens","cattle","ducks","goats","human","pigs","sheep")){next}			

				##known locations
				load(file=as.vector(points3[points3$species==dis3$name1,"filen"])) ##called data1$data1
				dt1<-data.frame(coordinates(data1$data1),data1$data1@data)
				coordinates(dt1)<-~lon+lat
				projection(dt1)<-projection(countr)
				dt1$names2<-tolower(dt1$name)
				#dt1$names2<-gsub(" ","_",dt1$names2)
				dt1$names2<-gsub("_"," ",dt1$names2)
				species1<-unique(dt1$names2[!is.na(dt1$names2)])##species in gbif
				####function to count number of spaces in a string
				no_spaces<-function(x,what=" ") lengths(regmatches(x, gregexpr(what, x)))
				species1<-species1[no_spaces(species1)<2]
				num_spec<-length(species1)
				species1<-species1[sapply(species1,FUN=function (x) is.na(as.numeric(x)==((as.numeric(x))/1)))] ###species in gbif
				if(length(species1)>0){species2<-unique(read.table(text=species1,sep=" ",stringsAsFactors=FALSE)$V1)}else{species2<-c()} ##genus in gbif	

				###habitats etc
				ttt<-gsub("_1_subs_",";1;",dis3$name1,fixed=TRUE)
				ttt<-gsub("_2_subs_",";2;",ttt,fixed=TRUE)
				ttt<-gsub("_3_subs_",";3;",ttt,fixed=TRUE)
				ttt<-gsub("_4_subs_",";4;",ttt,fixed=TRUE)
				species1a<-strsplit(ttt,";")[[1]][1]
				species1a<-strsplit(species1a,"_")[[1]][1]
				if(length(strsplit(species1a," ")[[1]])>1){species2a<-strsplit(species1a," ")[[1]][1]}else{species2a<-species1a}
				species2<-unique(c(species2,species2a))##genus
				#species1a<-gsub(" ","_",species1a)
				species1a<-gsub("_"," ",species1a)
				species1<-unique(c(species1,species1a))###species
				#species1;species2
        
				##endemic region
				ttt2<-as.numeric(strsplit(strsplit(ttt,";")[[1]][3],"_")[[1]])
				#if(paste(paste(ttt2,collapse="_"),".tif",sep="") %in% list.files("X:\\DRtemp\\temp\\masks\\", pattern=".tif",full.names=FALSE)){next}
				ws2<-wrld_simpl[wrld_simpl$SUBREGION %in% ttt2,]
				ws3<-suppressWarnings(crop(ws2,remove_small_islands(ws2,min_val=10)))
				template2<-suppressWarnings(crop(template,ws3))

				### load a mask
				mask1<-raster(paste("X:\\DRtemp\\temp\\masks\\",paste(ttt2,collapse="_"),".tif",sep=""))
			
				#lc_mask<-stack(futlgx1[futlgx1$V2==(strsplit(ttt,";")[[1]][3]),"filen"])
				#fnf2<-raster(fnfx1[fnfx1$V2==(strsplit(ttt,";")[[1]][3]),"filen"])
				template2<-crop(template,mask1)
        
				#}else{species1<-dis3$name1;species2<-dis3$name1}

			spec_dist$lower2<-tolower(paste(spec_dist$Taxa,spec_dist$lower,sep=" "))	
			## get the species data
			spec_hab2<-spec_hab[spec_hab$species %in% species1,]
			spec_dist2<-spec_dist[spec_dist$lower2 %in% species1,]
			spec_desp2<-spec_desp[spec_desp$binom %in% species1,]
			if(nrow(spec_hab2)==0){spec_hab2<-spec_hab[spec_hab$species %in% species2,]}
			if(nrow(spec_dist2)==0){spec_dist2<-spec_dist[tolower(spec_dist$Taxa) %in% species2,]}
			if(nrow(spec_desp2)==0){spec_desp2<-spec_desp[spec_desp$binom %in% species2,]}
			
			if(nrow(spec_desp2)==0|nrow(spec_dist2)==0){
			
			ff<-gnr_resolve(species1)
			ff<-ff[ff$data_source_title=="NCBI",]
			ff<-ff[!duplicated(ff$user_supplied_name),"matched_name"]
			ff2<-get_uid_(ff)
			ff3<-do.call("rbind",ff2)
      
			c1<-classification(ff3$uid,db="ncbi", accepted = TRUE)
			NewVar<-function (x,y,z) if(!is.na(x)) {x[x$rank==y,z]}
			species1f<-base::suppressWarnings(tolower(unique(unlist(lapply(names(c1), function(x) NewVar(c1[[x]], "family" ,"name"))))))
			species1o<-base::suppressWarnings(tolower(unique(unlist(lapply(names(c1), function(x) NewVar(c1[[x]], "order" ,"name"))))))
			species1c<-base::suppressWarnings(tolower(unique(unlist(lapply(names(c1), function(x) NewVar(c1[[x]], "class" ,"name"))))))
			
			ffx<-gnr_resolve(species2)
			ffx<-ffx[ffx$data_source_title=="NCBI",]
			ffx<-ffx[!duplicated(ffx$user_supplied_name),"matched_name"]
			ffx2<-get_uid_(ffx)
			ffx3<-do.call("rbind",ffx2)
			
			c2<-classification(ffx3$uid,db="ncbi", accepted = TRUE)
			species2f<-base::suppressWarnings(tolower(unique(unlist(lapply(names(c1), function(x) NewVar(c1[[x]], "family" ,"name"))))))
			species2o<-base::suppressWarnings(tolower(unique(unlist(lapply(names(c1), function(x) NewVar(c1[[x]], "order" ,"name"))))))
			species2c<-base::suppressWarnings(tolower(unique(unlist(lapply(names(c1), function(x) NewVar(c1[[x]], "class" ,"name"))))))
			
			species3f<-unique(c(species1f,species2f))
			species3o<-unique(c(species1o,species2o))
			species3c<-unique(c(species1c,species2c))
			
			if(nrow(spec_hab2)==0){spec_hab2<-spec_hab[spec_hab$species %in% species3f,]}
			if(nrow(spec_dist2)==0){spec_dist2<-spec_dist[tolower(spec_dist$Taxa) %in% species3f,]}
			if(nrow(spec_desp2)==0){spec_desp2<-spec_desp[spec_desp$binom %in% species3f,]}
			if(nrow(spec_hab2)==0){spec_hab2<-spec_hab[spec_hab$species %in% species3o,]}
			if(nrow(spec_dist2)==0){spec_dist2<-spec_dist[tolower(spec_dist$Taxa) %in% species3o,]}
			if(nrow(spec_desp2)==0){spec_desp2<-spec_desp[spec_desp$binom %in% species3o,]}
			if(nrow(spec_hab2)==0){spec_hab2<-spec_hab[spec_hab$species %in% species3c,]}
			if(nrow(spec_dist2)==0){spec_dist2<-spec_dist[tolower(spec_dist$Taxa) %in% species3c,]}
			if(nrow(spec_desp2)==0){spec_desp2<-spec_desp[spec_desp$binom %in% species3c,]}
			rm(c1,c2,ff,ff2,ff3,ffx,ffx2,ffx3,species1o,species1c,species1f,species2o,species2c,species2f,species3o,species3c,species3f)
			}
			
			##put in mean if still none <- ADD MORE TO DATAFRAME
			##spec_hab ok
			if(nrow(spec_desp2)==0){print("no species data");break}#spec_desp2<-data.frame(binom=NA,min=mean(spec_desp$min,trim=0.1,na.rm=TRUE),max=mean(spec_desp$max,trim=0.1,na.rm=TRUE),mean=mean(spec_desp$mean,trim=0.1,na.rm=TRUE),min_mean=mean(spec_desp$min_mean,trim=0.1,na.rm=TRUE),max_mean=mean(spec_desp$max_mean,trim=0.1,na.rm=TRUE),stringsAsFactors = FALSE)}
			if(nrow(spec_dist2)==0){print("no species data");break}#spec_dist2<-data.frame(Taxa=NA,lower=NA,Average.of.Mass=mean(spec_dist$Average.of.Mass,trim=0.1,na.rm=TRUE),Average.of.walking.speed=mean(spec_dist$Average.of.walking.speed,trim=0.1,na.rm=TRUE),Average.of.d=mean(spec_dist$Average.of.d,trim=0.1,na.rm=TRUE),Average.of.density=mean(spec_dist$Average.of.density,trim=0.1,na.rm=TRUE),stringsAsFactors = FALSE)}

					### spec-dist fail gracefully
			### quantufy dispersal ability 
			disp=log(mean(spec_desp2$mean,na.rm=TRUE)+1)#/10   # ned a more complex function here
			### set of lots of if but - or sort out beforehand
			fw1<-focalWeight(mask1, d=c(disp, disp*2.8), "Gauss")### could make the shape a function of maximum to mean
			#dim(fw1)*5.6
			#sum(fw1)
			fw1<-fw1*1.1/max(fw1)#sum(dim(fw1))#(1/min(fw1))#spec_dist2$Average.of.density
			#sum(fw1)
			#plot((fw1))
	
			#### limit to each species
			link4<-link3[link3$species==dis3$name1,]
			link4$code<-paste(link4$model,link4$year,link4$RCP,sep="_")
			link4$codeX<-paste(link4$year,link4$RCP,sep="_")
			
			###remove those with no overlap for hosts and vectors
			###but what what are best models <-choose both here
			link4<-link4[link4$codeX %in% bothmod,]
			link4<-link4[order(link4$year,link4$RCP),]
			link4<-link4[c((1:nrow(link4))[link4$RCP==999],(1:nrow(link4))[!link4$RCP==999]),]
			link4$RCP[link4$RCP==999]<-2.6
			
			#summarise gas parameters
			link4$speed<-mean(spec_dist2$Average.of.walking.speed,na.rm=TRUE)
			link4$d<-mean(spec_dist2$Average.of.d,na.rm=TRUE)
			link4$density<-mean(spec_dist2$Average.of.density*num_spec,na.rm=TRUE)
			link4$density[link4$density>3000000]<-3000000

			##PRESENT DAY FIRST <--------------------------------------------------------START
			link5<-link4[1, ]		
			if(link5$year!=2010){print("no present day");next}
					
			###check has this been done before?
			if( paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,".tif",sep="") %in% list.files("X:/DRtemp/resultsZ/",pattern=".tif",full.names=TRUE)){
					  
					  IP<-raster(paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,".tif",sep=""))
					  clims2<-raster(paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"clims2.tif",sep=""))
					  lcsuit<-raster(paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"lcsuit.tif",sep=""))
					  k<-raster(paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"k.tif",sep=""))
					  load(file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"meds4.r",sep=""))						  
					       
					  
			}else{    
										
					###future climate data
					clims<-raster(link5$filen)
					clims2<-crop(clims,template2)
					rm(clims)
					clims2<-mask(clims2,mask1)

					#clims2[clims2>=0.5]<-1 ## good
					#clims2[clims2>0.15 & clims2<0.9]<-0.5 ## marginal
					clims2[clims2<=0.05]<-0 ## no
					
					###future land-use data
					futlh<-futlg[futlg$year==link5$year & futlg$RCP==link5$RCP,]
					#link3[link3$species==dis3$name1 & link3$year==years[y],"RCP"]
					#lc2<-subset(lc_mask,as.numeric(rownames(futlh)))
					lc2<-stack(futlh$filen)
					names(lc2)<-futlh$variable
					#if(years[y]==2005){lc<-calc(lc,mean,na.rm=TRUE)}##only for 2015 mean across all types doens work
					#lc2<-crop(lc,template2)
					#lc2<-resample(lc2,template2,method="ngb")
					#lc2<-mask(lc2,mask1)

					####chose all primary and make not forest
					f1<-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"])
					nf1<-f1
					f1[values(fnf)!=1]<-0
					nf1[values(fnf)==1]<-0
					
					####chose all secondary and make forest/not forest
					s1<-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"])
					#wood1<-max(subset(lc2,(1:nlayers(lc2))[futlh$variable %in% c("gfsh1","gfsh2","gfsh3")]),na.rm=TRUE)
					
					###pasture as grasslands
					### snow ice water
					water<-1-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"])
					#water2<-focal(water)
					
					## create habitat data for those that at missing it
					if(nrow(spec_hab2)==0){
					     ##only create first time around on current land-use
					    #if(y==1){
						    #tx<-raster(lc2)
						    #res(tx)<-2
						    #values(tx)<-1:ncell(tx)
						    #same1<-extract(tx,dt1)
						    #data2<-dt1[!duplicated(same1),]
						    # dtR<-sampleRandom(lc2, 5000, ext=extent(dt1), xy = TRUE, sp=TRUE, na.rm = TRUE)
						    #habs<-extract(subset(lc2,c(1,5:8)),dtR, method='bilinear') ## whats the null expectation? ## prortion of habitats in bounding box
					  	  habsX<-extract(subset(lc2,c(1,5:8)),dt1, method='bilinear')
					  	  #meds<-apply(habs,2,function(x) mean(x,trim=0,na.rm=TRUE))
						    meds2<-apply(habsX,2,function(x) mean(x,trim=0,na.rm=TRUE))
						    #meds3=(meds2-meds)/meds
						    #meds4=((meds3-min(meds3))/max(meds3-min(meds3))/2)+0.5 ### push towards urban due to bias
						    meds4=(meds2/2)+0.5 ### push towards urban due to bias
					     #}
					  
					  lcsuit<-(meds4[2]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"]))+(meds4[4]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"]))+(meds4[3]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(meds4[1]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(meds4[5]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))

						}else{
						###sum up suitability
						  #lcsuit<-(wood1*mean(spec_hab2$f,na.rm=TRUE))+(mean(spec_hab2$Grasslands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(mean(spec_hab2$Croplands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(mean(spec_hab2$Urban.and.built.up,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
						  lcsuit<-(f1*mean(spec_hab2$f,na.rm=TRUE))+(nf1*mean(spec_hab2$nf,na.rm=TRUE))+(water*mean(spec_hab2$water,na.rm=TRUE)) + (mean(spec_hab2$Grasslands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(mean(spec_hab2$Croplands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(mean(spec_hab2$Urban.and.built.up,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
						  meds4<-NA
						  				}
					
					lcsuit<-resample(lcsuit,clims2)
					###lcsuit<-focal(lcsuit,w=matrix(1/25,nrow=5,ncol=5))
					lcsuit<-mask(lcsuit,mask1)
					### create climate layer

					### mosuitos 1.21 -3 million km2 silver et al Mosquito ecology 2007
					### all mosquitos species ## several studies 40-60 species per site ### no references as yet
					### 20,000-75000 per species including the rare species
					### 150,000-300,000 of single common species in 4km2 [Constanini et al. Medical and Veterinary Entomology (1996) 10,203-219] 

					### carrying capacity layer (summed accross habitat preference)
					k<-lcsuit*link5$density ## make these number relavent to habitat type
					k[k>link5$density]<-link5$density
					k[k<0]<-0
					#k<-mask(k,mask1)
					#k[k==0]<-NA
					#k[is.na(k)]<-0 ## quicker but causes edge effects near the sea
					rm(f1,nf1,water)
					gc()
	
					### Initial populations
					IP<-raster(clims2)
					values(IP)<-0
					IP<-mask(IP,mask1)
					##set all cells current at maximum carrying capacity where samples are from
					##subset by known countries infected
					dt2<-dt1[is.na(over(dt1,countr)$NAME),]
					if(nrow(dt2)==0){dt2<-spsample(countr,100,type="regular")}		
					IP[cellFromXY(IP,dt1)]<-link5$density
					##sample from poor climate lowered
					IP<-IP*clims2
					IP[IP>k]<-k[IP>k]
					
					### test whether all countries have high values
					#countr$test1<-extract(IP,countr,mean,na.rm=TRUE)
					#if(length(countr$test1[countr$test1<1])>0){ 
					#  ############################################################################################<-<-<-----
					#  others2<-unique(link3[link3$species2==link5$species2,"region"])
					#  current<-link5$region
					#  others2<-others2[others2!=current]
					#  listed<-strsplit(others2,"_")
					#  listed2<-lapply(listed,function (x) match(x,current))
					 #					  }
					
					### start diffusion loop
					resx<-c()
						for (i in 1:150){

								resx[i]<-length(values(IP)[values(IP)>0& !is.na(values(IP))])
								print(resx[i])
								if(i>1){if(resx[i-1]/resx[i]==1){break}}
								##spread about using diffusion
								IP<-focal(IP,w=fw1,na.rm=TRUE)
								### stop going into sea or use mask below to remove step
								IP<-mask(IP,mask1)
								##make carry capacity k limit in different environments
								IP[IP>k]<-k[IP>k]
								##probability of occuring dependent on climate
								IP<-IP*clims2
								gc()
							}##end of I loop
						  writeRaster(IP,format="GTiff",file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,".tif",sep=""))				
						  writeRaster(k,format="GTiff",file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"k.tif",sep=""))				
						  writeRaster(lcsuit,format="GTiff",file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"lcsuit.tif",sep=""))				
						  writeRaster(clims2,format="GTiff",file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"clims2.tif",sep=""))				
						    
						  save(meds4,file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"meds4.r",sep="")	)					  
		  	
    			} ## if allready done for another disease

          
					##mask everything - I forget why
					IP<-mask(IP,mask1)
					lcsuit2<-mask(lcsuit,mask1)
					specX<-rep(link5$species,ncell(mask1))
					specX[is.na(values(mask1))]<-NA
					
					res1<-data.table(disease=dis3$disease,type=dis3$type,species=specX,speed=link5$speed,d=link5$d,density=link5$density,code=link5$code,cell.id=values(template2), clim=values(clims2),lc_suit=values(lcsuit2),abundance=values(k),realised=values(IP),i=1)
					res1<-res1[!is.na(species)]	
					names(res1)[9:12]<-paste(names(res1)[9:12],"present",sep="_")
					
					
				########################################### DO FUTURE	#################################					
				#res3<-NULL				
		
				## choose from link4 present day forward
				#yyy<-2:nrow(link4)#}
				### do relative change in this loop? -from present
				gc()
				rbl<-function(...) rbindlist(list(...))
				cl <- makeCluster(16)
				registerDoParallel(cl)
				
				st1<-Sys.time()
				res3 = foreach(ii=2:nrow(link4), .combine=rbl,.packages=c("raster","data.table")) %dopar%  {
				  #for(y in yyy){
				  
				  #ii=ii+1
				  ##choose one climate model
				  link5<-link4[ii, ]					
				  
				  
				  ###check has this been done before?
				  if( paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,".tif",sep="") %in% list.files("X:/DRtemp/resultsZ/",pattern=".tif",full.names=TRUE)){
				    
				    IP<-raster(paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,".tif",sep=""))
				    clims2<-raster(paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"clims2.tif",sep=""))
				    lcsuit<-raster(paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"lcsuit.tif",sep=""))
				    k<-raster(paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"k.tif",sep=""))
				    
				  }else{
				  
				    ###future climate data
				    clims<-raster(link5$filen)
				    clims2<-crop(clims,template2)
				    rm(clims)
				    clims2<-mask(clims2,mask1)
				  
				    #clims2[clims2>=0.5]<-1 ## good
				    #clims2[clims2>0.15 & clims2<0.9]<-0.5 ## marginal
				    clims2[clims2<=0.05]<-0 ## no
				  
				    #if(y==1){old_RCP=0 ;old_year=0}
				    #if(link5$RCP==old_RCP & link5$year==old_year){}else{
				    
				    ###future land-use data
				    futlh<-futlg[futlg$year==link5$year & futlg$RCP==link5$RCP,]
				    #link3[link3$species==dis3$name1 & link3$year==years[y],"RCP"]
				    #lc2<-subset(lc_mask,as.numeric(rownames(futlh)))
				    lc2<-stack(futlh$filen)
				    names(lc2)<-futlh$variable
				    #if(years[y]==2005){lc<-calc(lc,mean,na.rm=TRUE)}##only for 2015 mean across all types doens work
				    #lc2<-crop(lc,template2)
				    #lc2<-resample(lc2,template2,method="ngb")
				    #lc2<-mask(lc2,mask1)
				    
				    ####chose all primary and make not forest
				    f1<-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"])
				    nf1<-f1
				    f1[values(fnf)!=1]<-0
				    nf1[values(fnf)==1]<-0
				    
				    ####chose all secondary and make forest/not forest
				    s1<-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"])
				    #wood1<-max(subset(lc2,(1:nlayers(lc2))[futlh$variable %in% c("gfsh1","gfsh2","gfsh3")]),na.rm=TRUE)
				    
				    ###pasture as grasslands
				    ### snow ice water
				    water<-1-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"])
				    #water2<-focal(water)
				    
				    ## create habitat data for those that at missing it
				    if(nrow(spec_hab2)==0){
				      ##only create first time around on current land-use
				    	      
				      lcsuit<-(meds4[2]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"]))+(meds4[4]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"]))+(meds4[3]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(meds4[1]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(meds4[5]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
				      
				    }else{
				      ###sum up suitability
				      #lcsuit<-(wood1*mean(spec_hab2$f,na.rm=TRUE))+(mean(spec_hab2$Grasslands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(mean(spec_hab2$Croplands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(mean(spec_hab2$Urban.and.built.up,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
				      lcsuit<-(f1*mean(spec_hab2$f,na.rm=TRUE))+(nf1*mean(spec_hab2$nf,na.rm=TRUE))+(water*mean(spec_hab2$water,na.rm=TRUE)) + (mean(spec_hab2$Grasslands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(mean(spec_hab2$Croplands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(mean(spec_hab2$Urban.and.built.up,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
				      
				    }
				    
				    lcsuit<-resample(lcsuit,clims2)
				    ###lcsuit<-focal(lcsuit,w=matrix(1/25,nrow=5,ncol=5))
				    lcsuit<-mask(lcsuit,mask1)
				    ### create climate layer
				    
				    ### mosuitos 1.21 -3 million km2 silver et al Mosquito ecology 2007
				    ### all mosquitos species ## several studies 40-60 species per site ### no references as yet
				    ### 20,000-75000 per species including the rare species
				    ### 150,000-300,000 of single common species in 4km2 [Constanini et al. Medical and Veterinary Entomology (1996) 10,203-219] 
				    
				    ### carrying capacity layer (summed accross habitat preference)
				    k<-lcsuit*link5$density ## make these number relavent to habitat type
				    k[k>link5$density]<-link5$density
				    k[k<0]<-0
				    #k<-mask(k,mask1)
				    #k[k==0]<-NA
				    #k[is.na(k)]<-0 ## quicker but causes edge effects near the sea
				    rm(f1,nf1,water)
				    gc()
				      #} ## end on if the same year as last one no need to redo - just climate changes
				  
				    ### Initial populations
				    IP<-raster(clims2)
				    values(IP)<-0
				    IP<-mask(IP,mask1)
				    ##set all cells current at maximum carrying capacity where samples are from
				    ##subset by known countries infected
				    dt2<-dt1[!is.na(over(dt1,countr)$NAME),]
				    if(nrow(dt2)==0){dt2<-spsample(countr,100,type="regular")}		
				    IP[cellFromXY(IP,dt1)]<-link5$density
				    ##sample from poor climate lowered
				    IP<-IP*clims2
				    IP[IP>k]<-k[IP>k]
				  
				    ### test whether all countries have high values
				    #countr$test1<-extract(IP,countr,mean,na.rm=TRUE)
				    #if(length(countr$test1[countr$test1<1])>0){ 
				    #  ############################################################################################<-<-<-----
				    #  others2<-unique(link3[link3$species2==link5$species2,"region"])
				    #  current<-link5$region
				    #  others2<-others2[others2!=current]
				    #  listed<-strsplit(others2,"_")
				    #  listed2<-lapply(listed,function (x) match(x,current))
				    #					  }
				  
				    ### start diffusion loop
				    resx<-c()
			
				    for (i in 1:150){
				      
				      resx[i]<-length(values(IP)[values(IP)>0& !is.na(values(IP))])
				      print(resx[i])
				      if(i>1){if(resx[i-1]/resx[i]==1){break}}
				      ##spread about using diffusion
				      IP<-focal(IP,w=fw1,na.rm=TRUE)
				      ### stop going into sea or use mask below to remove step
				      IP<-mask(IP,mask1)
				      ##make carry capacity k limit in different environments
				      IP[IP>k]<-k[IP>k]
				      ##probability of occuring dependent on climate
				      IP<-IP*clims2
				      gc()
				      }##end of I loop
				    
				    ###write IP output so no need to keep doing
				    writeRaster(IP,format="GTiff",file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,".tif",sep=""))				
				    writeRaster(k,format="GTiff",file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"k.tif",sep=""))				
				    writeRaster(lcsuit,format="GTiff",file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"lcsuit.tif",sep=""))				
				    writeRaster(clims2,format="GTiff",file=paste("X:/DRtemp/resultsZ/",link5$species,"_",link5$code,"clims2.tif",sep=""))				
				    
				    				    
				  }  ##end of if this code and species has been done before
              ##if so just start here with IP
				  
				  IP<-mask(IP,mask1)
				  lcsuit2<-mask(lcsuit,mask1)
				  specX<-rep(link5$species,ncell(mask1))
				  specX[is.na(values(mask1))]<-NA
				  
				  #res1<-data.table(disease=dis3$disease,type=dis3$type,species=specX,speed=link5$speed,d=link5$d,density=link5$density,code=link5$code,cell.id=values(template2), clim=values(clims2),lc_suit=values(lcsuit2),abundance=values(k),realised=values(IP))
				  res_foreach<-data.table(species=specX,code=link5$code,cell.id=values(template2), clim_future=values(clims2),lc_suit_future=values(lcsuit2),abundance_future=values(k),realised_future=values(IP),i=ii)
				  return(res_foreach[!is.na(species)])
				  
				}## END OF FOREACH
				print("per species run")
				print(Sys.time()-st1)
				
				stopCluster(cl)
	  
			if(is.null(res4)){res4<-res1}else{res4<-rbindlist(list(res4,res1))}
			if(is.null(res5)){res5<-res3}else{res5<-rbindlist(list(res5,res3))}
			#rm(res1,res3)
				
					
				gc()
		}##end of z loop ## z is species/vectors/hosts
    
	  print("per disease time")
	 print(Sys.time()-st2)
	
  	### do we need to seperate by present and future?
    ###create merge columns
    res4[, merge1:= paste(species,cell.id,sep="-")]
    res5[, merge1:= paste(species,cell.id,sep="-")]
  
    ###make year by cell.id
    res5$code<-gsub("_2",";2",res5$code)
    res5$code<-gsub("0_","0;",res5$code)
    res5[,c("model","year","RCP"):=tstrsplit(code,";")]
    res5[,cell_by_year_by_RCP:=paste(cell.id,year,RCP,sep="_")]
    
    ##remove effects of model but add SD
    ## this is across models so need means not sums
    res5b<-res5[ ,.(mean(clim_future,na.rm=TRUE), mean(lc_suit_future,na.rm=TRUE), mean(abundance_future,na.rm=TRUE), mean(realised_future,na.rm=TRUE),sd(clim_future,na.rm=TRUE), sd(lc_suit_future,na.rm=TRUE), sd(abundance_future,na.rm=TRUE), sd(realised_future,na.rm=TRUE)), by = .(species,cell_by_year_by_RCP)]
    res5b[is.na(res5b)]<-0
    names(res5b)[3:ncol(res5b)]<-c(paste(names(res5)[4:7],"mean",sep="_"),paste(names(res5)[4:7],"sd",sep="_"))
    
    ###ADD IN TYPE SO YOU CAN THEM SUMMARISE BY TYPE
    add1<-res4[,c("species","type"),with=FALSE]
    add1<-add1[!duplicated(add1$species),]
    setkey(add1,species)
    setkey(res5b,species)
    res5b<-res5b[add1]
    
    ## what about mapping underlying species level patterns rather than disease level <- save raster?
    ## REMOVE EFFECT OF SPECIES
    ##present day and gas models params summarised for present day
    res4b<-res4[ ,.(mean(clim_present,na.rm=TRUE), mean(lc_suit_present,na.rm=TRUE), sum(abundance_present,na.rm=TRUE), sum(realised_present,na.rm=TRUE),mean(speed,na.rm=TRUE), mean(d,na.rm=TRUE), mean(density,na.rm=TRUE)), by = .(type,cell.id)]
    names(res4b)[3:ncol(res4b)]<-c(paste(names(res4)[c(9:12,4:6)],"mean",sep="_"))
    res4b[,merge1:=paste(cell.id,type,sep="_")]
    
    ##future by type cellby RCP 
    ## what about variation across species + mutliplying probabitiles?
    res5c<-res5b[ ,.(mean(clim_future_mean + (clim_future_sd*2),na.rm=TRUE), mean(clim_future_mean - (clim_future_sd*2),na.rm=TRUE), sum(realised_future_mean + (realised_future_sd*2),na.rm=TRUE), sum(realised_future_mean - (realised_future_sd*2),na.rm=TRUE)), by = .(type,cell_by_year_by_RCP)]
    names(res5c)[3:ncol(res5c)]<-c("clim_upper95","clim_lower95","realised_upper95","realised_lower95")
    
    ###make type by cell.id called merge1
    res5c$cell_by_year_by_RCP<-gsub("_2",";2",res5c$cell_by_year_by_RCP)
    res5c$cell_by_year_by_RCP<-gsub("0_","0;",res5c$cell_by_year_by_RCP)
    res5c[,c("cell.id","year","RCP"):=tstrsplit(cell_by_year_by_RCP,";")]
    res5c[,merge1:=paste(cell.id,type,sep="_")]
      
    ##merge present and future
    setkey(res4b,merge1)
    setkey(res5c,merge1)
    ###rem res5c is bigger as it has RCPs
    res6x<-res5c[res4b]
    #rm(res4,res5,res4b,res5b,res5c);gc()
    
    ### go to wide
    ### might have different RCPs so need NAs in that case <-----------PROBLEM
    ### what if no hosts or no vectors<-----
    res6h<-res6x[res6x$type=="hosts"]  
    res6v<-res6x[res6x$type=="vectors",c("cell_by_year_by_RCP","clim_present_mean","lc_suit_present_mean","abundance_present_mean","realised_present_mean","clim_upper95","clim_lower95","realised_upper95","realised_lower95","speed_mean","d_mean","density_mean"),with=FALSE]        
    #rm(res6x);gc()
    
    #### if human vec human
    if(nrow(res6h)==0 & dis_trans$Type == "HUMAN->VECTOR->HUMAN"){res6h<-res6v;rm(res6v)}
    
    #### if true absence - generate dataframe with same RCP year combinations but all zeros
    #### if false absence - fill in etc.
    if(nrow(res6v)==0 & dis_trans$Type %in% c("HOST->HUMAN","HOST->HUMAN->HUMAN") ){
          res6v<-data.table(cell_by_year_by_RCP=unique(res6h$cell_by_year_by_RCP),clim_present_mean=1,lc_suit_present_mean=1,abundance_present_mean=1,realised_present_mean=1,clim_upper95=1,clim_lower95=1,realised_upper95=1,realised_lower95=1,speed_mean=1,d_mean=1,density_mean=1)
          }
    
    #### if vector or host missing then break??
    #rm(res6h,res6v);gc()
    #VectHost<-c("HOST->VECTOR->HUMAN","HOST->VECTOR->HUMAN->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN")
    #noVect<-c("HOST->HUMAN","HOST->HUMAN->HUMAN")
    #noHost<-c("HUMAN->VECTOR->HUMAN")
    
    #### merge column wise for vectors or dummy vectors
    names(res6v)[2:ncol(res6v)]<-c(paste(names(res6v)[2:ncol(res6v)],"vector",sep="_"))
    setkey(res6h,cell_by_year_by_RCP)
    setkey(res6v,cell_by_year_by_RCP)
    res6<-res6h[res6v,nomatch=0]    
    
    ###make year by cell.id by year key
    res6[,cell_by_year:=paste(cell.id,year,sep="_")]
   
    ##merge with human pop
    setkey(res6,cell_by_year)
    res6<-res6[ft2,nomatch=0]
   
    ##merge with human pop and GDP
    res6<-res6[GP2,nomatch=0]
    
    ##merge with livestock
    lt3<-lt2[,c("cell_by_year",names(lt2)[tolower(names(lt2)) %in% dis2$name1]),with=FALSE]
    
    ##name a single column called secondary
    ###what about secondary??
    #### ADD SECONDARY TO HOSTS ONLY WHEN HOSTS>0 
    #### VECTOR=0 then hosts=0?
    
    if(ncol(lt3)>1){
      lt3[,secondary:=round(rowSums(lt3[,2:ncol(lt3)]),0)]
      lt3<-lt3[,c("cell_by_year","secondary")]
      res6<-res6[lt3,nomatch=0]
        } else{ res6[,secondary:=0]} ## should be zero because adding to hosts where host + vectors>1
    rm(lt3)#;gc()
    
    ##endemic_countries <- too validate <-----------------problem
    #pt1<-as.data.frame(xyFromCell(template,res6$cell.id),stringsAsFactors=FALSE)
    #coordinates(pt1)<-~x+y
    #projection(countr)<-projection(pt1)
    #res6$ec<-over(pt1,countr)$UN
    
    #### get rid of identical columns
    res6<-res6[,names(res6)[!names(res6) %in% c("i.id","i.ID","i.cell.id.2","i.time","i.time","i.type","i.cell.id","i.cell.id.1") ], with = FALSE ]
    ##how to bring hosts and vectors together? switch to wide format? <---- SORT OUT!!!
    gc()
    
    ## right estimate contact rate etc.
    
    ##calculate d and speed for humans
    ##proper gas model? <- two stage?
    
    res6[,present_realised_upper:= c((realised_upper95*speed_mean*d_mean) * (realised_upper95_vector*speed_mean_vector*d_mean_vector) *(secondary) * (humans * 0.0001 * 5)),]
    res6[,future_realised_upper:= c((realised_upper95*speed_mean*d_mean) * (realised_upper95_vector*speed_mean_vector*d_mean_vector) *(secondary) * (humans * 0.0001 * 5)),]
    res6[,future_realised_change_upper:= future_realised_upper-present_realised_change,]
    res6[,future_realised_change_upper_div:= present_realised/future_realised_change_upper,]
              
    res6[,present_realised_lower:= c((realised_lower95*speed_mean*d_mean) * (realised_lower95_vector*speed_mean_vector*d_mean_vector) *(secondary) * (humans * 0.0001 * 5)),]
    res6[,future_realised_lower:= c((realised_lower95*speed_mean*d_mean) * (realised_lower95_vector*speed_mean_vector*d_mean_vector) *(secondary) * (humans * 0.0001 * 5)),]
    res6[,future_realised_change_lower:= future_realised_lower-present_realised_change,]
    res6[,future_realised_change_lower_div:= present_realised/future_realised_change_lower,]
    
    res6[,future_realised_clim_upper:= clim_upper95-clim_present_mean,]
    res6[,future_realised_clim_upper_div:= clim_present_mean/clim_upper95,]

    res6[,future_realised_clim_lower:= clim_lower95-clim_present_mean,]
    res6[,future_realised_clim_lower_div:= clim_present_mean/clim_lower95,]    
    
    ### cut out population by GDP - add in future gini by region
    
    ### endemic country validation
    
    ## save
    
    ## run sensivity
    
    ###sumarise accross area
    
    ##save again for ggplots!
    
		
}##end of x loop ## per disease loop

	

####################

testmap1<-as.data.frame(res4[res4$species=="trachypithecus_presbytis_2_subs_35",])
pt2<-xyFromCell(template,unique(testmap1$cell.id),spatial=FALSE)
testmap1<-cbind(testmap1,pt2)
sp::coordinates(testmap1) <- ~x+y
pres1<-rasterize(testmap1,template,field="realised_present",sum,na.rm=TRUE)
pres1<-crop(pres1,pt1)
pres1[pres1==0]<-NA
values(pres1)[1]<-141
plot(pres1)
plot(wrld_simpl,col="light grey",border="white",add=TRUE)
plot(pres1,add=TRUE)
plot(countr,border="red",add=TRUE)

uni1<-unique(res5$species)
testmap1f<-res5[res5$species=="trachypithecus_presbytis_2_subs_35" & res5$code=="esm2m_2070_6",]
#testmap1f<-as.data.frame(res5[1:207003,])
pt2f<-xyFromCell(template,unique(testmap1f$cell.id),spatial=FALSE)
testmap1f<-cbind(testmap1f,pt2f)
sp::coordinates(testmap1f) <- ~x+y
pres1f<-rasterize(testmap1f,template,field="realised_future",sum,na.rm=TRUE)
pres1f<-crop(pres1f,pt1)
pres1f[pres1f==0]<-NA
plot(pres1f)
plot(wrld_simpl,col="light grey",border="white",add=TRUE)
plot(pres1f,add=TRUE)
plot(countr,border="red",add=TRUE)
