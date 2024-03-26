
require(raster)
require(rgdal)
require(dismo)
require(rJava)

#library(ncdf4)
#library(lubridate)
#library(grr)
#library(doParallel)
library(maptools)
#library(dismo)

setwd('/scratch/scratch/ucbtdw0/sdms')

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#source("./functions6.r")
#source("E:\\Dropbox\\R_scripts\\functions6.r")

data(wrld_simpl)

###look up all diseases
link1<-list.files("C:/temp/disease_analyses/", pattern=".csv",full.names=TRUE)

for (i in 1:length(link1)){
	tt<-read.csv(link1[i],stringsAsFactors=FALSE)
	l2<-gsub("C:/temp/disease_analyses/","",link1[i],fixed=TRUE)
	l2<-gsub(".csv","",l2,fixed=TRUE)
	l2<-gsub("_XXX","",l2,fixed=TRUE)
	l3<-strsplit(l2,"-")[[1]]
	if(length(l3)==2){tt$id=l3[2];tt$tag=l3[1]} else {tt$id=l3[1];tt$tag=NA}
	tt<-tt[,c("disease","type","name1","id","tag")]
	if(i==1){res1<-tt}else{res1<-rbind(res1,tt)}
	print(i)
	}
dis1<-unique(res1$disease)
res1$unique<-paste(res1$disease,res1$type,res1$name1,sep="_")
res1<-res1[!duplicated(res1$unique),]

##get all points data for each
points1a<-list.files("C:\\temp\\disease_analyses\\", pattern="all_points",full.names=TRUE)
points1<-list.files("C:\\temp\\disease_analyses\\", pattern="_points",full.names=TRUE)
points2b<-points1[!points1 %in% points1a]
points1a<-list.files("C:\\temp\\disease_analyses\\", pattern="all_points",full.names=FALSE)
points1<-list.files("C:\\temp\\disease_analyses\\", pattern="_points",full.names=FALSE)
points2<-points1[!points1 %in% points1a]
points2<-gsub("_points.r","",points2, fixed=TRUE)
points3<-data.frame(species=points2,filen=points2b)

####Get all future climate niches
link2b<-list.files("C:/temp/resultsY", pattern=".tif",full.names=TRUE)
link2<-list.files("C:/temp/resultsY", pattern=".tif",full.names=FALSE)
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
link2<-gsub("_present_",";present;2005;",link2, fixed=TRUE) ##check when in
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

link3<-read.table(text=link2,sep=";",stringsAsFactors=FALSE)
link3$filen<-link2b
names(link3)[1:4]<-c("species","model","year","RCP")

years=unique(link3$year)

####land covers this is for CMIP6  not CMIP5
futlc<-list.files("C:/LUH1/",recursive=TRUE,pattern=".txt",full.names=TRUE)
futld<-read.table(text=futlc,sep="/",stringsAsFactors=FALSE)
futle<-read.table(text=futld$V6,sep=".",stringsAsFactors=FALSE)
futlf<-cbind(futld,futle)
futlf$RCP<-2.6
futlf$RCP[futlf$V1=="LUHa_u2t1.v1_aim.v1.1"]<-6.0
futlf$RCP[futlf$V1=="LUHa_u2t1.v1_message.v1"]<-8.5
futlf$RCP[futlf$V1=="LUHa_u2t1.v1_minicam.v1"]<-4.5
futlf$filen<-futlc
#futlf$filen2<-futlc
#futlf$filen<-gsub(".txt",".asc",futlf$filen)
names(futlf)<-c("scrap6","scrap5","scrap4","model","scrap","scrap2","variable","year","scrap3","RCP","filen")
#unique(futlf$variable)
## forest not forest land
fnf<-raster("C:\\LUH1\\fnf_map.txt")

###water/ice fraction in cell ## doesnt change with age unless do 1- all others.
icewater<-raster("C:\\LUH1\\gicew.1700.txt")

##only keep need land-use 
## gfsh1 wood harvested from secondary mature forest
## gfsh1 wood harvested from secondary young forest
## gfsh3 wood harvested from secondary non-forest
futlg<-futlf[futlf$variable %in% c("gcrop","gsecd","gpast","gurbn","gothr","gfsh1","gfsh2","gfsh3"),]

###read in habitat preferences
###convert to CMIP5 land classes

spec_hab<-read.csv("E:\\Dropbox\\data\\habitat_both_fin.csv",stringsAsFactors=FALSE)
spec_hab$f<-(spec_hab$Evergreen.Needleleaf.forest+spec_hab$Evergreen.Broadleaf.forest+spec_hab$Deciduous.Needleleaf.forest+spec_hab$Deciduous.Broadleaf.forest+spec_hab$Mixed.forest)/5
spec_hab$nf<-(spec_hab$Closed.shrublands+spec_hab$Open.shrublands+spec_hab$Woody.savannas+spec_hab$Savannas+spec_hab$Grasslands)/5

spec_dist<-read.csv("E:\\Dropbox\\data\\host_density_d_distance.csv",stringsAsFactors=FALSE)

spdata<-stack("./predictorsX.tif")
names(spdata)<-read.csv("./names1.csv",stringsAsFactors=FALSE)$x
spdata<-subset(spdata,c(8,10:length(names(spdata))))
Altitude<-subset(spdata,1)
#names(Altitude)<-"Altitude"

#res1[res1$tag=="aedes communis_ochlerotatus communis sp_1_subs_21_29",]
#jamestown canyon x=22
# yellow x=15

for (x in 1:length(dis1)){

	##choose one disease
	dis2<-res1[res1$disease==dis1[x],]

		##choose one host/vector
		for(z in 1:length(dis2)){

			dis3<-dis2[z,] #### SORT OUT MULTIPLIE FILES PROBLEM

			###habitats etc
			ttt<-gsub("_1_subs_",";1;",dis3$name1,fixed=TRUE)
			ttt<-gsub("_2_subs_",";2;",ttt,fixed=TRUE)
			ttt<-gsub("_3_subs_",";3;",ttt,fixed=TRUE)
			ttt<-gsub("_4_subs_",";4;",ttt,fixed=TRUE)
			species2<-strsplit(ttt,";")[[1]][1]
			species1<-strsplit(species2,"_")[[1]][1]
			spec_hab2<-spec_hab[spec_hab$species==species1,]
			spec_dist2<-spec_dist[spec_dist$lower==species1,]

			##endemic region
			ttt2<-as.numeric(strsplit(strsplit(ttt,";")[[1]][3],"_")[[1]])
			if(paste(paste(ttt2,collapse="_"),".tif",sep="") %in% list.files("C:\\temp\\masks\\", pattern=".tif",full.names=FALSE)){next}
			ws2<-wrld_simpl[wrld_simpl$SUBREGION %in% ttt2,]
			ws3<-crop(ws2,remove_small_islands(ws2,min_val=10))
			template2<-suppressWarnings(crop(template,ws3))
		
			### make a mask
			mask1<-rasterize(ws3, template2)##make these beforehand
			names(mask1)<-"mask1"
			writeRaster(mask1,file=paste("C:\\temp\\masks\\",paste(ttt2,collapse="_"),".tif",sep=""))
			#mask1<-raster(paste("C:\\temp\\masks\\",paste(ttt2,collapse="_"),".tif",sep=""))
		}#end of z loop
	}#end of x loop

			##known locations
			load(file=as.vector(points3[points3$species==dis3$name1,"filen"])) ##called data1$data1
			dt1<-data.frame(coordinates(data1$data1),data1$data1@data)
			coordinates(dt1)<-~lon+lat
			
			
				for(y in 1:length(years)){

					clims<-stack(link3[link3$species==dis3$name1 & link3$year==years[y],"filen"])
					clims2<-crop(clims,ws3)
					clims2<-mask(clims2,mask1)

					clims2[clims2>=0.5]<-1 ## good
					clims2[clims2>0.15 & clims2<0.9]<-0.5 ## marginal
					clims2[clims2<=0.15]<-NA ## no
					plot(clims2)
					points(dt1)
					
					futlh<-futlg[futlg$year==years[y],]
					#link3[link3$species==dis3$name1 & link3$year==years[y],"RCP"]
					lc<-stack(futlh[,"filen"])
					lc2<-crop(lc,ws3)
					f1<-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"])
					nf1<-f1
					f1[fnf!=1]<-0
					nf1[fnf==1]<-0
					
					if(nrow(spec_hab2)==0){
						tx<-raster(lc2)
						res(tx)<-2
						values(tx)<-1:ncell(tx)
						same1<-extract(tx,dt1)
						data2<-dt1[!duplicated(same1),]
						habs<-extract(subset(lc,c(1,5:8)),dt1, method='bilinear') ## whats the null expectation? ## prortion of habitats in bounding box
						habsX<-extract(subset(lc,c(1,5:8)),data2)
						boxplot(log(habs+1))
						##sumarise by category
						habs2<-aggregate(habs,by=list(habs$
						

					### create climate layer
					### carrying capacity layer (summed accross habitat preference)
					k<-raster(clims2)
					values(k)<-spec_dist2$Average.of.density ## make these number relavent to habitat type
					k<-mask(k,mask1)
					#k[is.na(k)]<-0 ## quicker but causes edge effects near the sea
					
					### create a habitat layer and multiply by preferences
					(spec_hab$nf*nf1)+(spec_hab$f*f1)+(lc$gurbn*spec_hab2$Urban.and.built.up)+(lc$gcrop*spec_hab2$Croplands)

					### Initial populations
					IP<-raster(clims2)
					values(IP)<-0
					IP<-mask(IP,mask1)
					IP[cellFromXY(IP,dt1)]<-spec_dist2$Average.of.density
					IP[IP>k]<-k[IP>k]
					### 
						for (i in 1:1000){

							print(i)
							print(length(values(IP)[values(IP)>50& !is.na(values(IP))]))
							##spread about

							##spec_dist2
     							#Taxa lower Average.of.Mass Average.of.walking.speed Average.of.d Average.of.density
							#Culex culex             0.1                        1     1.74e-05             150000


							## make matrix = to dispersal ability
							IP<-focal(IP,w=,mean,na.rm=TRUE)
							### stop going into sea or use mask below to remove step
							IP<-mask(IP,mask1)
							##make carry capacity k limit in different environments
							IP[IP>k]<-k[IP>k]
							##probability of occuring dependent on climate
							IP<-IP*clims2
						
						
						}##end of I loop

						plot(IP,colNA="blue")
					

				}##end of y loopn

		}##end of z loop

	#https://htmlpreview.github.io/?https://github.com/MCMaurer/D_RUG_Winter_18_Talk/blob/master/Talk_Slides.html#1
	#write.feather(disease,etc.

}##end of x loop


