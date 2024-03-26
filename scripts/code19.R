
require(raster)
require(rgdal)
require(dismo)
#require(rJava)
library(data.table)
library(sp)

#library(ncdf4)
#library(lubridate)
#library(grr)
#library(doParallel)
library(maptools)
#library(dismo)

#setwd('/scratch/scratch/ucbtdw0/sdms')

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
values(template)<-1:ncell(template)

#source("./functions6.r")
#source("E:\\Dropbox\\R_scripts\\functions6.r")
source("C:\\Users\\xxxx\\Dropbox\\R_scripts\\functions6.r")

##data.table function
cbind2 <- function(...) (setattr(do.call(c,list(...)),"class",c("data.table","data.frame")))


years=c(2030,2050,2070,2080,2010)
RCPS<-c(2.6,8.5,6.0,4.5)
lcs<-c("gcrop","gsecd","gpast","gurbn","gothr","gfsh1","gfsh2","gfsh3")

data(wrld_simpl)

###read disease data
d1<-read.csv("C:\\Users\\xxxx\\Dropbox\\disease_niche_analysis\\disease_data\\DRs_data\\disease_table28.csv",stringsAsFactors=FALSE)

###look up all diseases
for (i in 1:length(link1)){
	tt<-read.csv(link1[i],stringsAsFactors=FALSE)
	l2<-gsub("C:/temp/disease_analyses/","",link1[i],fixed=TRUE)
	l2<-gsub(".csv","",l2,fixed=TRUE)
	l2<-gsub("_XXX","",l2,fixed=TRUE)
	l3<-strsplit(l2,"-")[[1]]
	if(length(l3)==2){tt$id=l3[2];tt$tag=l3[1]} else {tt$id=l3[1];tt$tag=NA}
	tt<-tt[,c("disease","type","name1","id","tag")]
	if(i==1){diseases<-tt}else{diseases<-rbind(diseases,tt)}
	print(i)
}
diseases$name2<-paste(diseases$name1,";",sep="")
diseases$name2<-gsub("_subs_",";",diseases$name2)
diseases$name2<-gsub(";;",";",diseases$name2)
dd<-read.table(text=diseases$name2,sep=";",stringsAsFactors = F)
dd$V1<-gsub("_1","",dd$V1)
dd$V1<-gsub("_2","",dd$V1)
dd$V1<-gsub("_3","",dd$V1)
dd$V1<-gsub("_4","",dd$V1)

dis1<-unique(diseases$disease)
diseases$unique<-paste(diseases$disease,diseases$type,diseases$name1,sep="_")
diseases<-diseases[!duplicated(diseases$unique),]

##get all points data for each
points1a<-list.files("C:\\temp\\disease_analyses\\", pattern="all_points",full.names=TRUE)
points1<-list.files("C:\\temp\\disease_analyses\\", pattern="_points",full.names=TRUE)
points2b<-points1[!points1 %in% points1a]
points1a<-list.files("C:\\temp\\disease_analyses\\", pattern="all_points",full.names=FALSE)
points1<-list.files("C:\\temp\\disease_analyses\\", pattern="_points",full.names=FALSE)
points2<-points1[!points1 %in% points1a]
points2<-gsub("_points.r","",points2, fixed=TRUE)
points3<-data.frame(species=points2,filen=points2b)

#load gridded livestock of the world
#load("C:\\temp\\livestock_future_2030_2050_2070_2080.r")
#fwrite(lt2,file="C:\\temp\\livestock_future_2030_2050_2070_2080.csv")
lt2<-fread("C:\\temp\\livestock_future_2030_2050_2070_2080.csv")


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

###convert to masked

##poppop
#gdppop<-stack("C:\\POP_GDP_SSP13\\gdp_pop_ssp1_3.tif")
gpt<-read.csv(file="C:\\POP_GDP_SSP13\\gdp_pop_ssp1_3.csv",stringsAsFactors=FALSE)$x
gpt2<-gpt
gpt2<-gsub("op2","op;2",gpt2)
gpt2<-gsub("dp2","dp;2",gpt2)
gpt2<-gsub("ssp",";ssp",gpt2)
gptx<-read.table(text=gpt2,sep=";",stringsAsFactors=FALSE)
gptx$gpt<-gpt

gdppop2<-list.files("C:\\temp\\input_masked\\",recursive=TRUE,pattern="gdppopXXX",full.names=TRUE)
gdppopx<-gdppop2
gdppopx<-gsub("XXX",";",gdppopx)
gdppopx<-gsub(".tif","",gdppopx)
gdppopx1<-read.table(text=gdppopx,sep=";",stringsAsFactors=FALSE)
gdppopx1$filen<-gdppop2

###futpop
futpopx2<-list.files("C:\\temp\\input_masked\\",recursive=TRUE,pattern="fpopXXX",full.names=TRUE)
futpopx1<-futpopx2
futpopx1<-gsub("XXX",";",futpopx1)
futpopx1<-gsub(".tif","",futpopx1)
futpopx1<-read.table(text=futpopx1,sep=";",stringsAsFactors=FALSE)
futpopx1$filen<-futpopx2

futpop<-read.csv(file="C:\\temp\\futpop.csv",stringsAsFactors=FALSE)

###lc
futlgx2<-list.files("C:\\temp\\input_masked\\",recursive=TRUE,pattern="lcXXX",full.names=TRUE)
futlgx1<-futlgx2
futlgx1<-gsub("XXX",";",futlgx1)
futlgx1<-gsub(".tif","",futlgx1)
futlgx1<-read.table(text=futlgx1,sep=";",stringsAsFactors=FALSE)
futlgx1$filen<-futlgx2

futlg<-read.csv(file="C:\\temp\\futlg.csv",stringsAsFactors=FALSE)

###fnf
fnfx2<-list.files("C:\\temp\\input_masked\\",recursive=TRUE,pattern="fnfXXX",full.names=TRUE)
fnfx1<-fnfx2
fnfx1<-gsub("XXX",";",fnfx1)
fnfx1<-gsub(".tif","",fnfx1)
fnfx1<-read.table(text=fnfx1,sep=";",stringsAsFactors=FALSE)
fnfx1$filen<-fnfx2

## forest not forest land
fnf<-raster("C:\\temp\\fnf_map.txt")

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
spec_hab<-read.csv("E:\\Dropbox\\data\\habitat_both_fin.csv",stringsAsFactors=FALSE)
spec_hab$Snow.and.ice[spec_hab$Snow.and.ice==0]<-NA
spec_hab$f<-(spec_hab$Evergreen.Needleleaf.forest+spec_hab$Evergreen.Broadleaf.forest+spec_hab$Deciduous.Needleleaf.forest+spec_hab$Deciduous.Broadleaf.forest+spec_hab$Mixed.forest)/5
spec_hab$nf<-(spec_hab$Closed.shrublands+spec_hab$Open.shrublands+spec_hab$Woody.savannas+spec_hab$Savannas+spec_hab$Grasslands)/5
spec_hab$water<-rowMeans(spec_hab[,c("Water","Permanent.wetlands","Snow.and.ice")],na.rm=TRUE)

spec_dist<-read.csv("E:\\Dropbox\\data\\host_density_d_distance.csv",stringsAsFactors=FALSE)
spec_desp<-read.csv("E:\\Dropbox\\data\\dispersal_distance_all.csv",stringsAsFactors=FALSE)
spec_desp$binom<-gsub("_"," ",spec_desp$binom)

###get current bioclim data
spdata<-stack("C:\\temp\\predictorsX.tif")
names(spdata)<-read.csv("E:\\Dropbox\\names1.csv",stringsAsFactors=FALSE)$x
spdata<-subset(spdata,c(8,10:length(names(spdata))))
Altitude<-subset(spdata,1)
#names(Altitude)<-"Altitude"

#res1[res1$tag=="aedes communis_ochlerotatus communis sp_1_subs_21_29",]
#jamestown canyon x=22
# yellow x=15
#cache valley x=132

for (x in 1:length(dis1)){

	##choose one disease
	dis2<-diseases[diseases$disease==dis1[x],]
	if(length(dis2$tag[!is.na(dis2$tag)])==0){next}
	
	##endemic region
	dis_trans<-d1[d1$name==dis1[x],]
	countr<-wrld_simpl[wrld_simpl$ISO2 %in% strsplit(dis_trans$countries,",")[[1]],]
		

		##choose one host/vector
		for(z in 1:nrow(dis2)){

			dis3<-dis2[z,] #### SORT OUT MULTIPLIE FILES PROBLEM
			
			
			### ARGHHHH LOAD GLW			
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
				#if(paste(paste(ttt2,collapse="_"),".tif",sep="") %in% list.files("C:\\temp\\masks\\", pattern=".tif",full.names=FALSE)){next}
				ws2<-wrld_simpl[wrld_simpl$SUBREGION %in% ttt2,]
				ws3<-crop(ws2,remove_small_islands(ws2,min_val=10))
				template2<-suppressWarnings(crop(template,ws3))
		
				### load a mask
				mask1<-raster(paste("C:\\temp\\masks\\",paste(ttt2,collapse="_"),".tif",sep=""))
			
				lc_mask<-stack(futlgx1[futlgx1$V2==(strsplit(ttt,";")[[1]][3]),"filen"])
				gdp_mask<-stack(gdppopx1[gdppopx1$V2==(strsplit(ttt,";")[[1]][3]),"filen"])
				futpop_mask<-stack(futpopx1[futpopx1$V2==(strsplit(ttt,";")[[1]][3]),"filen"])
				fnf2<-raster(fnfx1[fnfx1$V2==(strsplit(ttt,";")[[1]][3]),"filen"])
				template2<-crop(template,mask1)

				#}else{species1<-dis3$name1;species2<-dis3$name1}

			## get the species data
			spec_hab2<-spec_hab[spec_hab$species %in% species1,]
			spec_dist2<-spec_dist[spec_dist$lower %in% species1,]
			spec_desp2<-spec_desp[spec_desp$binom %in% species1,]
			if(nrow(spec_hab2)==0){spec_hab2<-spec_hab[spec_hab$species %in% species2,]}
			if(nrow(spec_dist2)==0){spec_dist2<-spec_dist[spec_dist$lower %in% species2,]}
			if(nrow(spec_desp2)==0){spec_desp2<-spec_desp[spec_desp$binom %in% species2,]}

			##put in mean
			#if(nrow(spec_desp2)==0){
			group1<-c("host_type","secondary_type","vector_type")[dis3$type==c("hosts","secondary","vectors")] 
			spec_desp2<-spec_desp[spec_desp$binom %in% dis_trans[,group1],]

			### spec_desp might not have any values - all NAs for instance  - sort out
			### spec-dist fail gracefully

			### quantufy dispersal ability 
			disp=mean(spec_desp2$mean,na.rm=TRUE)/10   # ned a more complex function here
			### set of lots of if but - or sort out beforehand
			fw1<-focalWeight(mask1, d=c(disp, disp*2.8), "Gauss")### could make the shape a function of maximum to mean
			#dim(fw1)*5.6
			#sum(fw1)
			fw1<-fw1*1.1/max(fw1)#sum(dim(fw1))#(1/min(fw1))#spec_dist2$Average.of.density
			#sum(fw1)
			#plot(raster(fw1))
	
			#### limit to each species
			link4<-link3[link3$species==dis3$name1,]
			link4$code<-paste(link4$model,link4$year,link4$RCP,sep="_")
			link4<-link4[order(link4$year,link4$RCP),]
			link4<-link4[c((1:nrow(link4))[link4$RCP==999],(1:nrow(link4))[!link4$RCP==999]),]
			link4$RCP[link4$RCP==999]<-2.6
			#summarise gas parameters
			link4$speed<-mean(spec_dist2$Average.of.walking.speed,na.rm=TRUE)
			link4$d<-mean(spec_dist2$Average.of.d,na.rm=TRUE)
			link4$density<-mean(spec_dist2$Average.of.density*num_spec,na.rm=TRUE)
			link4$density[link4$density>3000000]<-3000000

				res3<-NULL				
				NUM1<-5 ### number of models sampled
				if(nrow(link4)>NUM1){yyy<-c(1,sample(2:nrow(link4),NUM1-1,replace=FALSE))}else{yyy<-1:nrow(link4)}
				### do relative change in this loop? -from present
				for(y in yyy){
					
					##choose one climate model
					link5<-link4[y, ]					
					
					###future climate data
					clims<-raster(link5$filen)
					clims2<-crop(clims,template2)
					clims2<-mask(clims2,mask1)

					#clims2[clims2>=0.5]<-1 ## good
					#clims2[clims2>0.15 & clims2<0.9]<-0.5 ## marginal
					clims2[clims2<=0.15]<-0 ## no
					
					if(y==1){old_RCP=0 ;old_year=0}
					if(link5$RCP==old_RCP & link5$year==old_year){}else{
					
					###future land-use data
					futlh<-futlg[futlg$year==link5$year & futlg$RCP==link5$RCP,]
					#link3[link3$species==dis3$name1 & link3$year==years[y],"RCP"]
					lc2<-subset(lc_mask,as.numeric(rownames(futlh)))
					names(lc2)<-futlh$variable
					#if(years[y]==2005){lc<-calc(lc,mean,na.rm=TRUE)}##only for 2015 mean across all types doens work
					#lc2<-crop(lc,template2)
					#lc2<-resample(lc2,template2,method="ngb")
					#lc2<-mask(lc2,mask1)

					####chose all primary and make not forest
					f1<-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"])
					nf1<-f1
					f1[values(fnf2)!=1]<-0
					nf1[values(fnf2)==1]<-0
					
					####chose all secondary and make forest/not forest
					s1<-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"])
					wood1<-max(subset(lc2,(1:nlayers(lc2))[futlh$variable %in% c("gfsh1","gfsh2","gfsh3")]),na.rm=TRUE)
					
					###pasture as grasslands
					### snow ice water
					water<-1-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"])-subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"])
					
					## create habitat data for those that at missing it
					if(nrow(spec_hab2)==0){
					   ##only create first time around on current land-use
					   if(y==1){
						#tx<-raster(lc2)
						#res(tx)<-2
						#values(tx)<-1:ncell(tx)
						#same1<-extract(tx,dt1)
						#data2<-dt1[!duplicated(same1),]
						#dtR<-sampleRandom(lc2, 5000, ext=extent(dt1), xy = TRUE, sp=TRUE, na.rm = TRUE)
						#habs<-extract(subset(lc2,c(1,5:8)),dtR, method='bilinear') ## whats the null expectation? ## prortion of habitats in bounding box
						habsX<-extract(subset(lc2,c(1,5:8)),dt1, method='bilinear')
						#meds<-apply(habs,2,function(x) mean(x,trim=0,na.rm=TRUE))
						meds2<-apply(habsX,2,function(x) mean(x,trim=0,na.rm=TRUE))
						#meds3=(meds2-meds)/meds
						#meds4=((meds3-min(meds3))/max(meds3-min(meds3))/2)+0.5 ### push towards urban due to bias
						meds4=(meds2/2)+0.5 ### push towards urban due to bias
							}
						lcsuit<-(meds4[2]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"]))+(meds4[4]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"]))+(meds4[3]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(meds4[1]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(meds4[5]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))

						}else{
						###sum up suitability
						lcsuit<-(f1*spec_hab2$f)+(nf1*spec_hab2$nf)+(water*spec_hab2$water)+((s1-wood1)*spec_hab2$nf)+((wood1)*spec_hab2$f)+(spec_hab2$Grasslands*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(spec_hab2$Croplands*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(spec_hab2$Urban.and.built.up*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
						}
					
					### lcsuit<-mask(lcsuit,mask1)
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
					rm(f1,nf1,water,wood1)
					gc()
						} ## end on if the same plot

					old_RCP=link5$RCP ;old_year=link5$year

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
					
					### start diffusion loop
					resx<-c()
					###check has this been done before?
					if( paste("X:/resultsZ/",link5$species,"_",link5$code,".tif",sep="") %in% list.files("X:/resultsZ/",pattern=".tif",full.names=TRUE)){

						IP<-raster(paste("X:/resultsZ/",link5$species,"_",link5$code,".tif",sep=""))
						}else{
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
						writeRaster(IP,format="GTiff",file=paste("X:/resultsZ/",link5$species,"_",link5$code,".tif",sep=""))				

						}

					IP<-mask(IP,mask1)
					lcsuit2<-mask(lcsuit,mask1)
					specX<-rep(link5$species,ncell(mask1))
					specX[is.na(values(mask1))]<-NA
					
				gc()		
				##create data frame use data.table or feather
				## read up
				## add the year coloumns onto side
				if(y==1){res1<-data.table(disease=dis3$disease,type=dis3$type,species=specX,speed=link5$speed,d=link5$d,density=link5$density,speed=link5$speed,cell.id=values(template2), clim_present=values(clims2),lc_suit_present=values(lcsuit2),abundance_present=values(k),realised_present=values(IP));res1<-res1[!is.na(species)];next}
				gc()
				if(is.null(res3)){res3<-data.table(species=specX,code=link5$code, clim_future=values(clims2),lc_suit_future=values(lcsuit2),abundance_future=values(k),realised_future=values(IP));res3<-res3[!is.na(species)];next}
	   			gc()
				res3<-rbindlist(list(res3,data.table(species=specX,code=link5$code, clim_future=values(clims2),lc_suit_future=values(lcsuit2),abundance_future=values(k),realised_future=values(IP))))
				res3<-res3[!is.na(species)]	   			
				gc()
			print("#########")
			print(y)
			print("#########")
			}##end of y loopn ### some point sumarise by disease?? need to account for vectors? set a updating field for host presence? or vector presence or just create a bifg dataframe?

			###summarise by species here? or just add together and deal with later??
			### some of this has to be done several times e.g. vector spread by aedes over NA
			### save it and redo

			if(z==1){res4<-res1}else{res4<-rbindlist(res4,res1)}
			if(z==1){res5<-res3}else{res5<-rbindlist(res5,res3)}
			rm(res1,res3)
			gc()
		}##end of z loop ## z is species/vectors/hosts

		##aggregate vectors by model code
		##aggregate hosts
		## add secondary
		## add human
		
		dis2a<-dis2[is.na(dis2$tag)]
		cols1<-c((1:6)[dis2a$name %in% names(lt2)[1:6]],(1:6)[dis2a$name %in% names(lt2)[1:6]]+7,(1:6)[dis2a$name %in% names(lt2)[1:6]]+13,(1:6)[dis2a$name %in% names(lt2)[1:6]]+19,(1:6)[dis2a$name %in% names(lt2)[1:6]]+25)
		lt3<-lt2[res1$cell.id,cols1]

		## use dis3 to collate human and secondary and livestock
		##minus present day
		##divide by present day
		##have column for country reporting the disease

		#dataset2$contactsx=ceiling(1*(((dataset2$susceptible)*0.67*0.8)*(dataset2$host_density)*5.6*transx2$d*sqrt(transx2$host_distance^2+mean(dataset2$daily_walking_distance)^2)))

		#res1<-res1[!is.na(res1$maxent)  , ] ## and is na any of the pops ## need to be the same do to relative change
					
		#save per disease data frame
		##address contact rates here using an average of d and distance as otherwise there is just too many columns ( should be approximate)

		need t match year to each

		##Future populations						
		futpopX<-futpop[futpop$year==link5$year & !futpop$location=="total",]
		futpop2<-subset(futpop_mask,as.numeric(rownames(futpopX)))
		names(futpop2)<-paste(futpopX$location,futpopX$SSP,futpopX$year,sep="_") #### NEED TO CHECK THESE!!!!!

		##Future gdp						
		gptx2<-gptx[gptx$V2==link5$year,]
		gdppop2<-subset(gdp_mask,as.numeric(rownames(gptx2)))
		names(gdppop2)<-gptx2$gpt

		
	#https://htmlpreview.github.io/?https://github.com/MCMaurer/D_RUG_Winter_18_Talk/blob/master/Talk_Slides.html#1
	#write.feather(disease,etc.

}##end of x loop ## per disease loop

	##save summarize by disease database


