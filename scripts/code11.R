
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

#source("./functions6.r")
source("E:\\Dropbox\\R_scripts\\functions6.r")

years=c(2030,2050,2070,2080,2010)
RCPS<-c(2.6,8.5,6.0,4.5)
lcs<-c("gcrop","gsecd","gpast","gurbn","gothr","gfsh1","gfsh2","gfsh3")

data(wrld_simpl)

###read disease data
d1<-read.csv("E:\\Dropbox\\disease_niche_analysis\\disease_data\\DRs_data\\disease_table25.csv")

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
link2b<-list.files("X:/resultsY", pattern=".tif",full.names=TRUE)
link2<-list.files("X:/resultsY", pattern=".tif",full.names=FALSE)
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
#cache valley x=24

for (x in 1:length(dis1)){

	##choose one disease
	dis2<-res1[res1$disease==dis1[x],]

		##choose one host/vector
		for(z in 1:nrow(dis2)){

			dis3<-dis2[z,] #### SORT OUT MULTIPLIE FILES PROBLEM

			### ARGHHHH LOAD GLW			
			if(dis3$name1 %in% c("ducks","chickens","cattle","ducks","goats","human","pigs","sheep")){next} ## add in sensbile way


			##known locations
			load(file=as.vector(points3[points3$species==dis3$name1,"filen"])) ##called data1$data1
			dt1<-data.frame(coordinates(data1$data1),data1$data1@data)
			coordinates(dt1)<-~lon+lat
			dt1$names2<-tolower(dt1$name)
			#dt1$names2<-gsub(" ","_",dt1$names2)
			dt1$names2<-gsub("_"," ",dt1$names2)
			species1<-unique(dt1$names2[!is.na(dt1$names2)])##species in gbif
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

			## get the species data
			spec_hab2<-spec_hab[spec_hab$species %in% species1,]
			spec_dist2<-spec_dist[spec_dist$lower %in% species1,]
			spec_desp2<-spec_desp[spec_desp$binom %in% species1,]
			if(nrow(spec_hab2)==0){spec_hab2<-spec_hab[spec_hab$species %in% species2,]}
			if(nrow(spec_dist2)==0){spec_dist2<-spec_dist[spec_dist$lower %in% species2,]}
			if(nrow(spec_desp2)==0){spec_desp2<-spec_desp[spec_desp$binom %in% species2,]}

			### spec_desp might not have any values - all NAs for instance  - sort out
			### spec-dist fail gracefully

			### aggregate multirow tables to species rather than add another loop??

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


			### loop through years to create dataframe ## create years/RCP combinations ## get rid of combinations that don't exist
			### can parallelise then?

			#### limit to each species
			link4<-link3[link3$species==dis3$name1,]
			link4$code<-paste(link4$model,link4$year,link4$RCP,sep="_")
			link4<-link4[c((1:nrow(link4))[link4$RCP==999],(1:nrow(link4))[!link4$RCP==999]),]

				for(y in 1:nrow(link4)){
					
					##choose one climate model
					link5<-link4[y, ]					
					
					RCP1<-link5$RCP
					RCP2<-link5$RCP

					if(link5$year==2010){RCP1=999;RCP2=2.6}

					###future climate data
					clims<-raster(link5$filen)
					clims2<-crop(clims,template2)
					clims2<-mask(clims2,mask1)

					#clims2[clims2>=0.5]<-1 ## good
					#clims2[clims2>0.15 & clims2<0.9]<-0.5 ## marginal
					clims2[clims2<=0.15]<-0 ## no
					#plot(clims2)
					#points(dt1)
					
					###future land-use data
					futlh<-futlg[futlg$year==link5$year & futlg$RCP==RCP1,]
					#link3[link3$species==dis3$name1 & link3$year==years[y],"RCP"]
					lc2<-subset(lc_mask,as.numeric(rownames(futlh)))
					names(lc2)<-futlh$variable
					#if(years[y]==2005){lc<-calc(lc,mean,na.rm=TRUE)}##only for 2015 mean across all types doens work
					#lc2<-crop(lc,template2)
					#lc2<-resample(lc2,template2,method="ngb")
					#lc2<-mask(lc2,mask1)

					###SORT OUT
					#fnf2<-crop(fnf,lc2)
					#fnf2<-resample(fnf2,template2,method="ngb")
					#fnf2<-subset(fnf_mask,as.numeric(rownames(fnf_mask)))

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
						tx<-raster(lc2)
						res(tx)<-2
						values(tx)<-1:ncell(tx)
						same1<-extract(tx,dt1)
						data2<-dt1[!duplicated(same1),]
						dtR<-sampleRandom(lc2, 5000, ext=extent(dt1), xy = TRUE, sp=TRUE, na.rm = TRUE)
						habs<-extract(subset(lc2,c(1,5:8)),dtR, method='bilinear') ## whats the null expectation? ## prortion of habitats in bounding box
						habsX<-extract(subset(lc2,c(1,5:8)),data2, method='bilinear')
						meds<-apply(habs,2,function(x) mean(x,trim=0,na.rm=TRUE))
						meds2<-apply(habsX,2,function(x) mean(x,trim=0,na.rm=TRUE))
						meds3=(meds2-meds)/meds
						meds4=((meds3-min(meds3))/max(meds3-min(meds3))/2)+0.5 ### push towards urban due to bias
						
						lcsuit<-(meds4[2]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gothr"]))+(meds4[4]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gsecd"]))+(meds4[3]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(meds4[1]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(meds4[5]*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))

						}else{
						###sum up suitability
						lcsuit<-(f1*spec_hab2$f)+(nf1*spec_hab2$nf)+(water*spec_hab2$water)+((s1-wood1)*spec_hab2$nf)+((wood1)*spec_hab2$f)+(spec_hab2$Grasslands*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(spec_hab2$Croplands*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(spec_hab2$Urban.and.built.up*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
						}

					### create climate layer
					### carrying capacity layer (summed accross habitat preference)
					k<-lcsuit*spec_dist2$Average.of.density ## make these number relavent to habitat type
					k[k>spec_dist2$Average.of.density]<-spec_dist2$Average.of.density
					k[k<0]<-0
					#k<-mask(k,mask1)
					#k[k==0]<-NA
					#k[is.na(k)]<-0 ## quicker but causes edge effects near the sea
					
					### Initial populations
					IP<-raster(clims2)
					values(IP)<-0
					IP<-mask(IP,mask1)
					IP[cellFromXY(IP,dt1)]<-spec_dist2$Average.of.density
					IP[IP>k]<-k[IP>k]
					
					### quantufy dispersal ability 
					disp=mean(spec_desp2$mean,na.rm=TRUE)#/10   # ned a more complex function here
					### set of lots of if but - or sort out beforehand
					fw1<-focalWeight(IP, d=c(disp, disp*2.8), "Gauss")### could make the shape a function of maximum to mean
					dim(fw1)*5.6
					#sum(fw1)
					fw1<-fw1*1.1/max(fw1)#sum(dim(fw1))#(1/min(fw1))#spec_dist2$Average.of.density
					#sum(fw1)
					plot(raster(fw1))	
					
					resx<-c()
					system.time(
					for (i in 1:10){

							resx[i]<-length(values(IP)[values(IP)>0& !is.na(values(IP))])
							print(resx[i])

							if(i>1){if(resx[i-1]/resx[i]>0.999){break}}

							##spread about using diffusion
							IP<-focal(IP,w=fw1,na.rm=TRUE)
							### stop going into sea or use mask below to remove step
							IP<-mask(IP,mask1)
							##make carry capacity k limit in different environments
							IP[IP>k]<-k[IP>k]
							##probability of occuring dependent on climate
							IP<-IP*clims2
											
						}##end of I loop
						)
						#plot(IP,colNA="blue")

					##Future populations						
					futpopX<-futpop[futpop$year==link5$year & !futpop$location=="total",]
					futpop2<-subset(futpop_mask,as.numeric(rownames(futpopX)))
					names(futpop2)<-paste(futpopX$location,futpopX$SSP,futpopX$year,sep="_") #### NEED TO CHECK THESE!!!!!

					##Future gdp						
					gptx2<-gptx[gptx$V2==link5$year,]
					gdppop2<-subset(gdp_mask,as.numeric(rownames(gptx2)))
					names(gdppop2)<-gptx2$gpt

					
				##create data frame use data.table or feather
				## read up
				## add the year coloumns onto side
				res1<-data.table(disease=dis3$disease,type=dis3$type,species=link5$species,year=link5$year,RCP=link5$RCP,model=link5$model, maxent=values(clims2),forest=values(f1),nonforest=values(nf1),water=values(water),lc_suit=values(lcsuit),abundance=values(k),realised=values(IP),values(futpop2),values(gdppop2))
				res1<-res1[!is.na(res1$maxent)  , ] ## and is na any of the pops

				}##end of v loop

			}##end of y loopn ### some point sumarise by disease?? need to account for vectors? set a updating field for host presence? or vector presence or just create a bifg dataframe?

		}##end of z loop

	#https://htmlpreview.github.io/?https://github.com/MCMaurer/D_RUG_Winter_18_Talk/blob/master/Talk_Slides.html#1
	#write.feather(disease,etc.

}##end of x loop


