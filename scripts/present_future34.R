

###run disease validation
library(biomod2)
library(dismo)
library(maptools)
library(data.table)
library(sp)

data(wrld_simpl)
wrld_simpl[wrld_simpl$SUBREGION==5,"REGION"]<-20

library(raster)
library(fasterize)
library(dismo)
library(sf)
library(rgdal)

library(viridis)
library(sp)

library(ggplot2)
library(velox)

load(file="X:\\DRtemp\\babesia_problem.r")
source("C:\\Users\\xxxx\\Dropbox\\R_scripts\\TSS.r")


###missing GINI present
gin<-stack("V:\\GINI_STACK.tif")
names(gin)<-read.csv("V:\\GINI_STACK_NAMES.csv",stringsAsFactors = FALSE)$x
gin<-velox(gin)
##load pop dgp ssp
#popgdp<-stack("V:\\POP_GDP_SSP13\\gdp_pop_ssp1_3.tif")
#names(popgdp)<-read.csv("V:\\POP_GDP_SSP13\\gdp_pop_ssp1_3.csv",stringsAsFactors = FALSE)$x
#popgdp<-velox(popgdp)

##make dummy <- nned better  soluation
zz<-data.frame(x=-180,y=90)
coordinates(zz)<-~x+y
projection(zz)<-projection(raster())

##read in admins to move away from GRID
ad1<-readOGR("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\admin1\\ne_10m_admin_1_states_provinces.shp","ne_10m_admin_1_states_provinces")

###load GP4 - pop gdp
#load(file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\GP4.r")
#GP5=gin$extract(ad1,fun=function(x) mean(x,na.rm=TRUE))
#GP5<-as.data.frame(GP5,stringsAsFactors=FALSE)
#names(GP5)<-read.csv("V:\\GINI_STACK_NAMES.csv",stringsAsFactors = FALSE)$x
#setDT(GP5)
#GP5[,diss_me:=ad1$diss_me]
#GPpres<-GP5[,c("X2010_GINI_PRESENT_C","diss_me")]
#GPpres[,year_diss_me:=paste(2010,diss_me,sep="_")]
#setkey(GPpres,diss_me)
load(file="V:\\GPpres.r")
setkey(GPpres,diss_me)
#GP5<-GP5[,mget(names(GP5)[2:ncol(GP5)])]
#GP6<-melt(GP5,id.vars="diss_me", measure.vars=names(GP5)[1:12],variable.factor=FALSE,value.factor=FALSE)
#GP6$variable<-gsub("_C","",GP6$variable)
#GP6$variable<-gsub("X","",GP6$variable)
#GP6$variable<-gsub("_GINI","",GP6$variable)
#GP6[,c("year","SSP"):=tstrsplit(variable,"_")]
#setkey(GP6,diss_me)
#names(GP6)[2:7]<-gsub("-","_",names(GP4)[2:7],fixed=TRUE)
#GP6[,year_diss_me:=paste(year,diss_me,sep="_")]
#GP6[,type_ssp:=paste("gini",SSP,sep="_")]
#GP7<-dcast.data.table(GP6,year_diss_me~type_ssp,value.var = "value")
#save(GP7,file="V:\\GP7.r")
load(file="V:\\GP7.r")
setkey(GP7,year_diss_me)

##merge with ad1#
#ad1@data<-cbind(ad1@data,as.data.frame(GP4))

##make popgdp
#popgdp5=popgdp$extract(ad1,fun=function(x) sum(x,na.rm=TRUE))
#popgdp5<-as.data.frame(popgdp5,stringsAsFactors=FALSE)
#names(popgdp5)<-read.csv("V:\\POP_GDP_SSP13\\gdp_pop_ssp1_3.csv",stringsAsFactors = FALSE)$x
#setDT(popgdp5)
#popgdp5[,diss_me:=ad1$diss_me]
#popgdp6<-melt(popgdp5,id.vars="diss_me", measure.vars=names(popgdp5)[1:30],variable.factor=FALSE,value.factor=FALSE)
#popgdp6$variable<-gsub("p20","p_20",popgdp6$variable)
#popgdp6$variable<-gsub("ssp","_ssp",popgdp6$variable)
#popgdp6$variable<-gsub("_GINI","",popgdp6$variable)
#popgdp6[,c("type","year","SSP"):=tstrsplit(variable,"_")]
#popgdp6[,year_diss_me:=paste(year,diss_me,sep="_")]
#popgdp6[,type_ssp:=paste(type,SSP,sep="_")]
#popgdp7<-dcast.data.table(popgdp6,year_diss_me~type_ssp,value.var = "value")
#save(popgdp7,file="V:\\popgdp7.r")
load(file="V:\\popgdp7.r")
setkey(popgdp7,year_diss_me)
#popgdp2010<-popgdp6[year==2010,.(mean(value)) ,by=.(type,diss_me)]
#popgdp2010b<-dcast.data.table(popgdp2010,diss_me~type,value.var = "V1")
#names(popgdp2010b)[2:3]<-c("gdp2010","humans2010")
#save(popgdp2010b,file="V:\\popgdp2010b.r")
load(file="V:\\popgdp2010b.r")
setkey(popgdp2010b,diss_me)

##set up breaks for risk calculation
breaks1<-seq(from=0,to=max(log(popgdp7$pop_ssp3+1),na.rm=TRUE),by=(max(log(popgdp7$pop_ssp3+1),na.rm=TRUE)/10))
breaks2<-seq(from=min(log(1/(popgdp7$gdp_ssp1/popgdp7$pop_ssp1)),na.rm=TRUE),to=max(log(1/(popgdp7$gdp_ssp1/popgdp7$pop_ssp1)),na.rm=TRUE),by=((max(log(1/(popgdp7$gdp_ssp1/popgdp7$pop_ssp1)),na.rm=TRUE)-(min(log(1/(popgdp7$gdp_ssp1/popgdp7$pop_ssp1)),na.rm=TRUE)))/10))

##base raster template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

##create rasterized admin
pres_ad1<-fasterize(st_as_sf(ad1),template,field="diss_me")

##read in two country points
lesstwo<-read.csv("C:\\Users\\xxxx\\Dropbox\\data\\disease_in_two_or_fewer_countries4.csv",stringsAsFactors=FALSE)
lesstwo$name2<-paste(" ",lesstwo$name,sep="")
lesstwo$dummy<-1
coordinates(lesstwo)<-~Longitude+Latitude


##read in empress-i point data
point_data<-read.csv("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\all_empress_i_data4.csv",stringsAsFactors = FALSE)
point_data<-point_data[!is.na(point_data$SumCases),]
point_data$LU<-paste(" ",point_data$name_LU,sep="")

###read disease data
d1<-read.csv("C:\\Users\\xxxx\\Dropbox\\disease_table32b.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4


###find all present day
pres1<-list.files("V:\\per_disease2\\",pattern=" ALL1",full.names=TRUE)
pres2<-list.files("V:\\per_disease2\\",pattern=" ALL1",full.names=FALSE)
pres2<-gsub("  ALL1.r","",pres2)
pres2<-gsub(" ALL1.r","",pres2)

###find all present day
pres1b<-list.files("V:\\per_disease2\\",pattern=" ALL2",full.names=TRUE)
pres2b<-list.files("V:\\per_disease2\\",pattern=" ALL2",full.names=FALSE)
pres2b<-gsub("  ALL2.r","",pres2b)
pres2b<-gsub(" ALL2.r","",pres2b)

##get just the correct ones
pres1<-pres1[!pres2 %in% pres2b]
pres2<-pres2[!pres2 %in% pres2b]
pres1<-c(pres1,pres1b)
pres2<-c(pres2,pres2b)

##change projecion??
wrld_simpl2<-spTransform(wrld_simpl,CRS=projection(template))
wrld_simpl2[wrld_simpl2$NAME=="Western Sahara","AREA"]<-999
#wrld_simpl2[wrld_simpl2$NAME=="Luxembourg","AREA"]
wrld_simpl2<-wrld_simpl2[wrld_simpl2$NAME!="Antarctica" & wrld_simpl2$AREA>10,]

###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl2),template)

###create new raster
t2<-ws2
t2[!is.na(t2)]<-0

#load optimum weigtings
#save(res100,file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\res100b.r")
load("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\res100b.r")
res100$uni<-paste(res100$AUC_best.stat,res100$name,sep="_")

##fin best AUC
best1<-aggregate(res100$AUC_best.stat,by=list(res100$name),max)
best1$uni<-paste(best1$x,best1$Group.1,sep="_")
res100b<-res100[res100$uni %in% best1$uni,]
#hist(res100b$AUC_best.stat)

##loop 104 problems
res11<-NULL
i=0
dev.off()
#i=46 #dengue
for (i in (i+1):length(pres1)){
  
    ##load a disease
  load(pres1[i])
  
  ##optimums
  res100c<-res100b[res100b$name2==pres2[i],]

  #remove NAs  
  res6<-res6[!is.na(clim_present_mean),]
  if(nrow(res6)==0){print("no data");break}
  
    ### just need present day at the moment? Do all?
  #res6<-res6[,c(1:25,40:43)]
  res6<-cbind(res6,xyFromCell(template,res6$cell.id))

  ###get popGDP
  ##turnfrom wide to long
  ## test res6 is ok  
  if(sd(res6$clim_present_mean,na.rm=TRUE)==0){print("no data");break}
  if(nrow(res6)==0){print(pres1[i]);print("no raw data");next}
  if(!"cell.id" %in% names(res6)){next}
  
  ##check point data
  pd<-point_data[point_data$LU %in% pres2[i] & point_data$Status!="Denied" & point_data$SumCases>0,]
  if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd);projection(pd)<-projection(raster())} else {rw1<-1000}
  
  ##bv problem
  if(pres2[i]==paste(" ",bv,sep="")){pres2[i]<-" babesia venatorum"}
  
  ##create endemic region
  dis_trans<-d1[d1$name2==pres2[i],]
  CFR=mean(dis_trans$CFR.low,dis_trans$CFR.high,na.rm=TRUE)
  ccc<-strsplit(dis_trans$countries,",")[[1]]
  if(is.na(ccc[1])){print("no recognised countries");print(pres2[i]);next}
  countr<-wrld_simpl2[wrld_simpl2$ISO2 %in% ccc ,]
  
  
  ###test for county level data
  if(nrow(countr)!=length(ccc)){print("issue with countries");print(ccc)}
  if(nrow(countr)==0){print("no recognised countries");print(pres2[i]);next}
  regn<-wrld_simpl2[wrld_simpl2$SUBREGION %in% countr$SUBREGION,] ##psuedoabsence region
  regn2<-wrld_simpl2[wrld_simpl2$REGION %in% countr$REGION,] ##psuedoabsence region
  
  ###if coveerage of reporting countries is too high switch to larger region  <---- artribrary
  if(length(unique(countr$REGION))>1){regn<-regn2}
  regnPA<-regn ##region for psuego absences
  if(nrow(countr)/nrow(regn)>0.5){regnPA<-regn2}
  if(nrow(countr)/nrow(regn)>0.5){regnPA<-wrld_simpl2}
  
  ## but if not countries in subregion 
  
  ###if less than two countries
  if(nrow(countr)<3){
    l2<-lesstwo[lesstwo$name2 == pres2[i],]
    projection(l2)<-projection(ad1)
    ##overlay
    overpoints<-over(ad1,l2)
    ad2<-ad1[!is.na(overpoints$dummy), ] ### used to set present points samples
    ad2$dummy<-1:nrow(ad2)
    
    #regn<-countr
    
    #plot(ad2,border="red",add=TRUE)
    }else{
    ##if not less than two ad2 the "TRUE" sample locations are just the countries
    ad2<-countr    ### used to set present points samples
    ad2$dummy<-1:nrow(ad2)
  }
  
  ###determine if in regn and regn2
	pres_ad2<-fasterize(st_as_sf(ad2),template,field="dummy")
	regn11<-fasterize(st_as_sf(regn),template,field = "UN") ##either  subregion or region depending
	regn22<-fasterize(st_as_sf(regnPA),template,field = "UN")###region for PA needed?
	#diss_me2<-fasterize(st_as_sf(ad1),template)
	
	res6[,countr:=extract(pres_ad2,res6$cell.id)]
	res6[,sub_region:=extract(regn11,res6$cell.id)] ##whats used to determine projection area
	res6[,region:=extract(regn22,res6$cell.id)]  ##not used
	res6[,diss_me:=extract(pres_ad1,res6$cell.id)]

	####load rasters of popgdp
	#popgdp<-extract(popdgp,res6$cell.id)
	#popgdp<-melt(popgdp)
	
	### rid of any lines no inside	
	### regn
	### countries
	res6<-res6[!is.na(res6$region),]
	#res6b<-res6[!is.na(res6$region),]
	#res6<-res6[!is.na(res6$countr),] ##run again but for countries that have reported

	
	##create unique
	res6[,year_RCP:=paste(RCP,year,sep="_")]
	#res6[,cell_by_year:=paste(cell.id,year,sep="_")]
	#setkey(res6,cell_by_year)
	#year_RCP<-unique(res6$year_RCP)
	
	##spill over rate cases per year over population at risk 2005
	#spillover_rate<-as.numeric(dis_trans$cases_per_year)/(nrow(res6)/length(unique(res6$year_RCP)))
	spillover_rate<-as.numeric(dis_trans$cases_per_year)/sum(countr$POP2005,na.rm=TRUE)
	

  #### CREATE HAZARD DIFFERENCES
  
  ### Hazard options
  ### with and without land use - will affect area by threshold
  ### w/wo spill over rate - per person probability - not affect area just units
  ### incidence rate - how to calculate - with vulnerability? - doesn't make sense if then calculated differently in the end
  ### what defines the proportion at risk - all incidence is rough - gross population
  ### does it matter as just needs to relative?? Only if working out change per disease rather than comparing sums
  
  ### comes back to contact model
  
  ### w/wo CFR/burden (CFR risk of death - rather than risk of infection)
  ### what to display? Raw climate difference, LU climate difference, spillover weighted, CFR weighted
  
  ### Exposure - just numbers of people
  
  ### Vulnerability - proportion in poverty from gini * GDP -how does this then intersect with spillover rate?
  
    ## DON'T NEED TO RUN CELLS WITH 0 CLIM
    ## (MAYBE) DON'T NEED TO CREATE GAS MOEL WITH THESHOLDS AS ZERO OR 1 ONLY VARIATION COMES FROM lAND-USE
    ## USE UPPER AND LOWER AS MEASURES OF CERTAINANTY
    ## THOSE WITHIN TWO SD OF ZERO - EFFECTIVELY 0??
    
    ## What about cuts? no hosts, medium, high (0,0.5,0.75) (use cut (midpoint)=TRUE)

    ##abandon dispersal?
    ## create new value from clim * mean density
    ## include or not
    
	### if H->V then:
	### make new "host" from vector
	### create concept of 'contact species'
	##probability of meeting a host or vector
    ##but what if vector == 1
    if(sd(res6$clim_present_mean_vector,na.rm=TRUE)>0){
    res6[,clim_present_mean2:= c(1-((1-clim_present_mean) *(1-clim_present_mean_vector))),]  
    res6[,lc_suit_present_mean2:= lc_suit_present_mean_vector,]  
    res6[,density_mean:= density_mean_vector,]  
    res6[,clim_future_mean2:= c(1-((1-clim_future_mean) *(1-clim_future_mean_vector))),]  
    res6[,lc_suit_future_mean2:= lcsuit_future_mean_vector,]  
    #res6[,density_mean:= density_mean_vector,] 

   }else{
    res6[,clim_present_mean2:= clim_present_mean,]### change to fit in below
     res6[,clim_future_mean2:= clim_future_mean,]### change to fit in below
     res6[,lc_suit_present_mean2:= lc_suit_present_mean,]  
     res6[,lc_suit_future_mean2:= lcsuit_future_mean,]  
     
       
    }

	### BRING IN OPTIMISED LIMITS FROM AUC11 RUN  SAVED DATA
	#res6$clim_present_mean2[res6$clim_present_mean2>=AUC_THRESH]<-1
	#res6$clim_present_mean2[res6$clim_present_mean2<AUC_THRESH]<-0
  ##run gas model - any point in variation not going to influence spatial contact patterns
  #res6b[,present_realisedT:= c(((realised_present_meanH*ttmean)*(speed_mean*ttspeed)*(d_mean*ttd)) * (realised_present_meanV*abs(speed_mean_vector*ttspeedv))*(d_mean_vector*ttdv) *(secondary_present) * ((humans2010/5.6) * (0.0001*ttdh) * (5*ttspeedh))),]

  ###create Hazard
  ### intersection of host|vector|pathogen - spatial only no intensity (assume homogenous)
  ### pathogeen by TSS of pathogen given host/vector
  ### ROC curve plot(e,'ROC')
  ### clim should have been timesed by LULC when it was per species - re run!
  ### continous raster function using country centroids to find contiguous areas
  
  ### work out threshold for there or not
  ### clim etc. used to say is disease in this grid cell?
  ### then risk is infection or chance of death at a given incidence rate
  ### population at risk * yearly cases (* CFR)
  
  ### not looking at intensity - only way to examine that is through vlunerability (or habitat preference?)
  ### can measure intensity for only 10 or so diseases  

	###create realised present mean - abudance from habitat????
	res6[,Hazard := clim_present_mean2 * lc_suit_present_mean2 ]
	
	###create realised present mean - abudance from habitat????
	res6[,Hazard_future:= clim_future_mean2 * lc_suit_future_mean2 ]

	###create realised present mean - abudance from habitat????
	#res6[,Hazard_delta:= Hazard_future-Hazard  ]

	###Hazard outside reported subregion is 0 <-big assumption
	#res6[,Hazard2:=Hazard]
	#res6$Hazard2[is.na(res6$sub_region)]
	## final hazard
	## CFR needs to be greater than 0 - to approximate dalys
	## can we estimate DALYS (proportion of life lost * 1-CFR etc.)
	#res6[,Hazard:=HazardUW * spillover_rate2] ## or without CFR and spill over and divide up into cost groups
	## final hazard
	#res6[,Hazard_future:=Hazard_future_UW * spillover_rate2] ## or without CFR and spill over and divide up into cost groups

	#res8<-res8[,lapply(.SD,function (x) mean(x,trim=0.1,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=c("Hazard","HazardUW","Hazard_future","Hazard_future_UW","RiskP","RiskSSP1","RiskSSP2","RiskSSP3","delta_riskSSP1","delta_riskSSP2","delta_riskSSP3")]
	res7<-res6[!is.na(diss_me),lapply(.SD,function (x) mean(x,trim=0.1,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=c("Hazard","Hazard_future","clim_present_mean2","clim_future_mean2")]
  res7[,c("RCP","year"):=tstrsplit(year_RCP,"_")]
	res7[,year_diss_me:=paste(year,diss_me,sep="_")]
	setkey(res7,year_diss_me)
  rm(res6);gc()
	
	#### END OF HAZARD ####
	### SUMMARIZE AND COMBINE INVARIANT DATA
	### combine pop gdp gini and hazard
	res8<-res7[popgdp7,nomatch=0]
	res8<-res8[GP7,nomatch=0]
	setkey(res8,diss_me)
	res8<-res8[GPpres]
	setkey(res8,diss_me)
	res8<-res8[popgdp2010b]
	#rm(res6);gc()
	
	#res8<-res8[,-c("i.year_diss_me","i.gdp2010", "i.humans2010" )]
	## final exposure
	## could include cattle??
	## could cut by at risk groups
	#res8[,pop_ssp1b:=pop_ssp1]
	#res8[,pop_ssp2b:=pop_ssp2]
	#res8[,pop_ssp3b:=pop_ssp3]
	#res8$pop_ssp1b[res8$pop_ssp1b>quantile(res8$pop_ssp1,0.9,na.rm=TRUE)]<-quantile(res8$pop_ssp1,0.90,na.rm=TRUE)
	#res8$pop_ssp2b[res8$pop_ssp2b>quantile(res8$pop_ssp2,0.9,na.rm=TRUE)]<-quantile(res8$pop_ssp2,0.90,na.rm=TRUE)
	#res8$pop_ssp3b[res8$pop_ssp3b>quantile(res8$pop_ssp3,0.9,na.rm=TRUE)]<-quantile(res8$pop_ssp3,0.90,na.rm=TRUE)
	
	#res8[,ExposureSSP1:=pop_ssp1]
	#res8[,ExposureSSP2:=pop_ssp2]
	#res8[,ExposureSSP3:=pop_ssp3]

	##should this be made relative - so examine change multipled by Hazard as we cannot test intensity
	### new aproach - hazard then validate and THRESHOLD - first part of paper.
	### then we examine the fates of population in the predicted endemic area (and just reporting countries)
	### set hazard to 1 then (exposure* vulnerability to see new value)
	### centre exposure + vulnerability to 1 - say waht is the fate of people in endemic areas?
	### then weight by CFR and spill over rate to see different patterns (for 2050?)
	### 1) Hazard 2) raw RISK 3) spill-over weigthed RISK 4) spill +CFR Burden RISK (have effective diseases per cell to show how general)
	### use scenarios to set varaiblity/agreement
	### then do 3) across RCP?
		
	###threshold HAZARD not risk
	### sort out region etc
	### more then one continent then whole continent not subregion right?
	
	quanttrim<-function(x,thresh=0.99){
		  x[x>quantile(x,thresh,na.rm=TRUE)]<-quantile(x,thresh,na.rm=TRUE)
	  x[x<quantile(x,1-thresh,na.rm=TRUE)]<-quantile(x,1-thresh,na.rm=TRUE)
	  return(x)
	}
	
	res8[,Exposure:=log(humans2010+1)]
	res8[,ExposureSSP1:=log(pop_ssp1+1)]
	res8[,ExposureSSP2:=log(pop_ssp2+1)]
	res8[,ExposureSSP3:=log(pop_ssp3+1)]
  ##what does this mean?
	#res8[,ExposureSSP1:=quanttrim((log(pop_ssp1+1)-Exposure)/Exposure)]
	#res8[,ExposureSSP2:=quanttrim((log(pop_ssp2+1)-Exposure)/Exposure)]
	#res8[,ExposureSSP3:=quanttrim((log(pop_ssp3+1)-Exposure)/Exposure)]

	res8[,Vulnerability:=log((gdp2010/humans2010)+1)]
	#res8$Vulnerability<-max(res8$Vulnerability,na.rm=TRUE)-res8$Vulnerability
	res8[,VulnerabilitySSP1:=log((gdp_ssp1/pop_ssp1)+1)]
	res8[,VulnerabilitySSP2:=log((gdp_ssp2/pop_ssp2)+1)]
	res8[,VulnerabilitySSP3:=log((gdp_ssp3/pop_ssp3)+1)]
	
	#res8[,VulnerabilitySSP1:=log((gdp_ssp1/pop_ssp1)+1)]
	#res8[,VulnerabilitySSP2:=log((gdp_ssp2/pop_ssp2)+1)]
	#res8[,VulnerabilitySSP3:=log((gdp_ssp3/pop_ssp3)+1)]
	#res8$VulnerabilitySSP1<-max(res8$VulnerabilitySSP1,na.rm=TRUE)-res8$VulnerabilitySSP1
	#res8$VulnerabilitySSP2<-max(res8$VulnerabilitySSP2,na.rm=TRUE)-res8$VulnerabilitySSP2
	#res8$VulnerabilitySSP3<-max(res8$VulnerabilitySSP3,na.rm=TRUE)-res8$VulnerabilitySSP3
	#res8$VulnerabilitySSP1<-quanttrim(res8$VulnerabilitySSP1/res8$Vulnerability,thresh=0.99)
	#res8$VulnerabilitySSP2<-quanttrim(res8$VulnerabilitySSP2/res8$Vulnerability,thresh=0.99)
	#res8$VulnerabilitySSP3<-quanttrim(res8$VulnerabilitySSP3/res8$Vulnerability,thresh=0.99)
	
	vw=res100c$vw
	ew=res100c$ew
	
	res8[,RiskP:= (Hazard) + (vw*(Vulnerability)) + (ew*(Exposure))]
	res8[,RiskSSP1:= (Hazard) + (vw*(VulnerabilitySSP1))* (ew*(ExposureSSP1))]
	res8[,RiskSSP2:= (Hazard) + (vw*(VulnerabilitySSP2))* (ew*(ExposureSSP2))]
	res8[,RiskSSP3:= (Hazard) + (vw*(VulnerabilitySSP3))* (ew*(ExposureSSP3))]
	
	##combine with dis_trans
	res9<-cbind(res8,res100c[,32:ncol(res100c)])
	
	if(is.null(res11)){res11<-res9}else{res11<-rbind(res11,res9)}
	
  rm(res8,res9,res100c)#,adc_data,ad1c)
	gc()
  print(i)
}

##save results
fwrite(res11,file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\final_risk4.csv")

res11<-fread(file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\final_risk4.csv")

res11[,yearRCPdissmename:=paste(RCP,i.year_diss_me,name,sep="_")]


res11[,RiskP:= (Hazard) + (vw*(Vulnerability)) + (ew*(Exposure))]
res11[,RiskSSP1:= (Hazard) + (vw*(VulnerabilitySSP1))* (ew*(ExposureSSP1))]
res11[,RiskSSP2:= (Hazard) + (vw*(VulnerabilitySSP2))* (ew*(ExposureSSP2))]
res11[,RiskSSP3:= (Hazard) + (vw*(VulnerabilitySSP3))* (ew*(ExposureSSP3))]


res12<-res11[!is.na(year_RCP),names(res11)[c(100,3:6,23:34,81:97)],with=F]
res12[,c("year","RCP","dissme","name"):=tstrsplit(yearRCPdissmename,"_")] ###order wrong
res12[,year_RCP_name:=paste(year,RCP,name,sep="_")]
res12[,year_RCP_dissme:=paste(year,RCP,dissme,sep="_")]


##for maps
res13<-aggregate(res12[,2:33],by=list(res12$year_RCP_dissme),mean,na.rm=TRUE)
res13[,c("year","RCP","dissme"):=tstrsplit(Group.1,"_")]

##for graphs
res14<-aggregate(res12[,2:33],by=list(res12$year_RCP_name),mean,na.rm=TRUE)
setDT(res14)
res14[,c("year","RCP","name"):=tstrsplit(Group.1,"_")]

## figure 1 ggplot 

ggplot(data=res14,aes(x=year,y=RiskSSP2,group=name))+
  geom_point()

  




###after optimum 
ad1c$RiskPT<-ad1c$Risk1
ad1c$RiskPT[ad1c$RiskPT>thr1[[2]]]<-1
ad1c$RiskPT[ad1c$RiskPT!=1]<-0


##when best breaks are found 
jpeg(file=paste("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\risk_maps\\",pres2[i],"_",sample(1:1000,1),"Hazard.jpg",sep=""),width=2500,height=1000)
print(
  ggplot(data = ad1c) +
    geom_sf(aes(fill = Hazard),colour=NA) +
    scale_fill_viridis_c(option = "plasma")+
    geom_sf(data=countr2,aes(),colour="white", alpha=0,size=1)+
    geom_sf(data=pd2,aes(),cex=0.1)+
    labs(title=paste(pres2[i]," - AUC ",round(e@auc,2),"  AUC points ",round(epd_auc,2)," ",names(thresh)," ",thresh[[1]]," TSS ",TSSv, sep=""))+
    theme(plot.title=element_text(size=32),axis.text = element_text(size=28))
  #+	facet_wrap(.~year_RCP)
)
dev.off()

jpeg(file=paste("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\risk_maps\\",pres2[i],"_",sample(1:1000,1),"Threshold_Hazard.jpg",sep=""),width=2500,height=1000)
print(
  ggplot(data = ad1c) +
    geom_sf(aes(fill = RiskPT),colour=NA) +
    scale_fill_viridis_c(option = "plasma")+
    geom_sf(data=countr2,aes(),colour="grey", alpha=0,size=1)+
    geom_sf(data=pd2,aes(),cex=0.1)+
    labs(title=paste(pres2[i]," - AUC ",round(e@auc,2),"  AUC points ",round(epd_auc,2)," ",names(thresh)," ",thresh[[1]]," TSS ",TSSv, sep=""))+
    theme(plot.title=element_text(size=32),axis.text = element_text(size=28))
  #+	facet_wrap(.~year_RCP)
)
dev.off()




##create final line plots (Hazard and Risk)

## create global summed risk map by aggregating by adl3 (Hazard and Risk)

##divide hazard (non weighted) into

###
#Emerging
#0-100
#Established
#100-10,000
#Endemic
#10000-2,500,000
#Hyperendemic
#2,500,000+



