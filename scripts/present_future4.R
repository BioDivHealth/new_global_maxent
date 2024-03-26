

###run disease validation
library(dismo)
library(maptools)
library(data.table)
library(sp)

data(wrld_simpl)

library(raster)
library(fasterize)
library(dismo)
library(sf)
library(rgdal)

library(viridis)
library(sp)

load(file="X:\\DRtemp\\babesia_problem.r")

#ft2<-fread("X:/DRtemp/temp\\future_pop_maskedb.csv")
#setkey(ft2,cell_by_year) ##create cell.id by code??

##read in admins to move away from GRID
ad1<-readOGR("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\admin1\\ne_10m_admin_1_states_provinces.shp","ne_10m_admin_1_states_provinces")

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

GP2<-fread("X:/DRtemp/temp\\future_gdppop_maskedb.csv")
setkey(GP2,cell_by_year) ##create cell.id by code??

#GP2[,c("x","y"):=as.data.frame(xyFromCell(template,GP2$cell.id))]
#coordinates(GP2)<-~x+y
#projection(GP2)<-projection(ad1)


##get gini
gini<-stack("V:\\GINI_STACK.tif")
names(gini)<-read.csv(file="V:\\GINI_STACK_NAMES.csv",stringsAsFactors = FALSE)$x
gini2<-extract(gini,xyFromCell(template,GP2$cell.id),df=TRUE)
gini2$cell.id<-GP2$cell.id
gini2<-gini2[,-15]
gini2<-gini2[,-1]
gini3<-melt(gini2)
gini3$cell.id<-GP2$cell.id
gini3$year<-read.table(text=as.character(gini3$variable),sep="_",stringsAsFactors = FALSE)[,1]
gini3$year<-as.numeric(gsub("X","",gini3$year))
setDT(gini3)
gini3[,cell_by_year:=paste(cell.id,year,sep="_")]
setkey(gini3,cell_by_year)
#setkey(GP2,time)
GP2<-GP2[gini3,nomatch=0]
rm(gini3);gc()

##addin missing values done need to do that as is all the same just need na.rm=TRUE
##GP2$value[is.na(GP2$value)]<-median(GP2$value,na.rm=TRUE)

##create rasterized admin
pres_ad1<-fasterize(st_as_sf(ad1),template,field="diss_me")
fff<-unique(GP2$cell.id)
diss_me<-extract(pres_ad1,fff)
fff2<-data.frame(cell.id=fff,diss_me=diss_me)
setDT(fff2)
setkey(fff2,cell.id)

##merge together
setkey(GP2,cell.id)
GP2<-GP2[fff2]
#GP2<-GP2[,-i.cell.id]

#GP3<-aggregate(GP2[,c(1:6,8,9)],by=list(GP2[,diss_me]),FUN=sum)
###data.table fast aggregate
GP3<-GP2[!is.na(diss_me),lapply(.SD,sum,na.rm=TRUE),by=diss_me,.SDcols=c("gdp-ssp1","gdp-ssp2","gdp-ssp3","pop-ssp1","pop-ssp2","pop-ssp3","gdp2010","humans2010")]
GP3b<-GP2[!is.na(diss_me),lapply(.SD,mean,na.rm=TRUE),by=diss_me,.SDcols=c("value")]

setkey(GP3,diss_me)
setkey(GP3b,diss_me)
GP4<-GP3[GP3b]

###match to
GP4<-GP4[match(ad1$diss_me,GP4$diss_me),]

##lauren's equation - check
## work out vulnerable per admin
GP4[,pov_frac_SSP1:= 1 - (1-(693.50/(gdp-ssp1*10000)))^((1+value)/(1-value))]
GP4[,pov_frac_SSP2:= 1 - (1-(693.50/(gdp-ssp2*10000)))^((1+value)/(1-value))]
GP4[,pov_frac_SSP3:= 1 - (1-(693.50/(gdp-ssp3*10000)))^((1+value)/(1-value))]

##merge with ad1
ad1@data<-cbind(ad1@data,as.data.frame(GP4))

##write new file
writeOGR("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\admin3_POPGDP.shp")



#TEST ITS ALL OK
#t2<-pres_ad1
#t2[!is.na(t2)]<-0
#t2[GP2[year==2030,]$cell.id]<-GP2[year==2030,]$'pop-ssp1'
#t2[GP2[year==2030,]$cell.id]<-GP2[year==2030,]$value
#plot(t2)



##aggregate to


extract(GP2)



#sum(GP2[GP2$time==2030,"pop-ssp3"],na.rm=TRUE)

#t2<-template
#values(t2)<-NA
#t2[GP2$cell.id]<-GP2$X2010_GINI_PRESENT_C
#plot(t2)
#data(wrld_simpl)
#plot(wrld_simpl,add=TRUE)
#gc()

#test1<-over(GP2,ad1)
#head(test1)

##CREATE VULNERABILITY
	
     #October 2015 world bank $1.9 dollars a day
	 # $693.50 per year
	 
	 ##what about GDP??
	 
	 ## instead of fraction  could have function of number with poverty as driver of vulnerability
	##where 1 perfect equatility and 0 all vulnerable 

	 ##ACTUALLY COULD SUMMERISE ALL THE POPULATION AND GDP STUFF INTO ADMINS AND SAVE!!



#### need to load lt2 etc. and create series of ad1 columns that capture mean(?) values for
#### each polygon
##humans
#th<-t2
#th[res6b$cell.id]<-res6b$humans2010

##can use fasterize then zonal
#ad1$dummy<-1:nrow(ad1)
#ad2<-fasterize(st_as_sf(ad1),template,field="dummy")
#ad1$meanH<-raster::zonal(th,ad2,fun=mean,na.rm=TRUE)
#ad1$meanH<-raster::extract(th,ad1,fun=mean,na.rm=TRUE)
#ad1$meanH<-log(ad1$meanH+1)

#mycolours <- brewer.pal(9, "YlOrRd")
#spplot(ad1,"meanH",  cuts = 8, col.regions = mycolours) 

##read in two country points
lesstwo<-read.csv("C:\\Users\\xxxx\\Dropbox\\data\\disease_in_two_or_fewer_countries4.csv",stringsAsFactors=FALSE)
lesstwo$name2<-paste(" ",lesstwo$name,sep="")
lesstwo$dummy<-1
coordinates(lesstwo)<-~Longitude+Latitude


##read in empress-i point data
point_data<-read.csv("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\all_empress_i_data2.csv",stringsAsFactors = FALSE)
point_data$LU<-paste(" ",point_data$name_LU,sep="")

###read disease data
d1<-read.csv("X:\\DRtemp\\disease_table30.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)

###find all present day
pres1<-list.files("V:\\per_disease2\\",pattern=" ALL1",full.names=TRUE)
pres2<-list.files("V:\\per_disease2\\",pattern=" ALL1",full.names=FALSE)
pres2<-gsub("  ALL1.r","",pres2)
pres2<-gsub(" ALL1.r","",pres2)
#pres3<-as.Date(file.mtime(pres1))
#pres1<-pres1[pres3>=as.Date("2019-04-10")]
#pres2<-pres2[pres3>=as.Date("2019-04-10")]


template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
values(template)<-1:ncell(template)

##change projecion??
wrld_simpl2<-spTransform(wrld_simpl,CRS=projection(template))

###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl2),template)
#plot(ws2)

###create new raster
t2<-ws2
t2[!is.na(t2)]<-0

##read in population and gdp and gin
#gdp<-stack(list.files("V:\\global_population_and_gdp\\data",full.names=TRUE)[c(4,6,8,10,11,15,17,19,20,24,26,28,29)])
#pop<-stack(list.files("V:\\global_population_and_gdp\\data",full.names=TRUE)[c(4,6,8,10,11,15,17,19,20,24,26,28,29)+31])
#gdp<-gdp/pop
#gini<-stack("V:\\GINI_STACK.tif")
#ginames<-read.csv("V:\\GINI_STACK_NAMES.csv",stringsAsFactors = FALSE)$x
#ginames<-gsub("X","gini_",ginames)
#ginames<-gsub("_C","",ginames)
#ginames<-gsub("_GINI_SSP","_",ginames)
#ginames<-gsub("_GINI_PRESENT","_",ginames)
#names(gini)<-ginames

##load up optimised thresholds
#AUCs<-read.csv("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\AUC3.csv",stringsAsFactors = FALSE)

res3<-NULL

##loop
i=24
for (i in i:length(pres1)){
  
    ##load a disease
  load(pres1[i])
  
  #remove NAs  
  res6<-res6[!is.na(clim_present_mean),]
  if(nrow(res6)==0){print("no data");break}
  
    ### just need present day at the moment? Do all?
  #res6<-res6[,c(1:25,40:43)]
  res6<-cbind(res6,xyFromCell(template,res6$cell.id))
  
  if(sd(res6$clim_present_mean,na.rm=TRUE)==0){print("no data");break}

  if(nrow(res6)==0){print(pres1[i]);print("no raw data");next}
  if(!"cell.id" %in% names(res6)){next}
  
  pd<-point_data[point_data$LU %in% pres2[i] & point_data$Status!="Denied",]
  if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd)} else {rw1<-1000}
  
  ##create endemic region
  if(pres2[i]==paste(" ",bv,sep="")){pres2[i]<-" babesia venatorum"}
  
  dis_trans<-d1[d1$name2==pres2[i],]
  ccc<-strsplit(dis_trans$countries,",")[[1]]
  if(is.na(ccc[1])){print("no recognised countries");print(pres2[i]);next}
  countr<-wrld_simpl[wrld_simpl$ISO2 %in% ccc ,]
  if(nrow(countr)!=length(ccc)){print("issue with countries");print(ccc);next}
  if(nrow(countr)==0){print("no recognised countries");print(pres2[i]);next}
  regn<-wrld_simpl[wrld_simpl$SUBREGION %in% countr$SUBREGION,] ##psuedoabsence region
  regn2<-wrld_simpl[wrld_simpl$REGION %in% countr$REGION,] ##psuedoabsence region
  
  ###if coveerage of reporting countries is too high switch to larger region  <---- artribrary
  if(nrow(countr)/nrow(regn)>0.75){regn<-regn2}
  if(nrow(countr)/nrow(regn)>0.75){regn<-wrld_simpl}
  
  ###if less than two countries
  if(nrow(countr)<3){
    l2<-lesstwo[lesstwo$name2 == pres2[i],]
    projection(l2)<-projection(ad1)
    ##overlay
    overpoints<-over(ad1,l2)
    ad2<-ad1[!is.na(overpoints$dummy), ]
    regn<-countr
    #plot(ad2,border="red",add=TRUE)
    
  }else{
    ##if not less than two ad2 the "TRUE" sample locations are just the countries
    ad2<-countr  
  }
  
  ###determine if in regn and regn2
  
	pres_ad2<-fasterize(st_as_sf(ad2),template)
	res6$countr<-extract(pres_ad2,res6$cell.id)
	res6$regn<-extract(regn,res6$cell.id)
	res6$regn2<-extract(regn2,res6$cell.id)
	res6$adl3<-extract(ad1,res6$cell.id)

	### rid of any lines no inside	
	### regn
	### countries
	res6<-res6[!is.na(res6$regn),]
	#res6<-res6[!is.na(res6$countr),] ##run again but for countries that have reported

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
    res6b[,clim_present_mean2:= c(1-((1-clim_present_mean) *(1-clim_present_mean_vector))),]  
    res6b[,lc_suit_present_mean:= lc_suit_present_mean_vector,]  
    res6b[,density_mean:= density_mean_vector,]  

   }else{
    res6b[,clim_present_mean2:= clim_present_mean,]### change to fit in below
	
    }
    	
	### WHAT ABOUT FUTURE??????????????????
		
	### BRING IN OPTIMISED LIMITS FROM AUC11 RUN  SAVED DATA
	#res6$clim_present_mean2[res6$clim_present_mean2>=AUC_THRESH]<-1
	#res6$clim_present_mean2[res6$clim_present_mean2<AUC_THRESH]<-0
		
	###create realised present mean - abudance from habitat????
	res6[,realised_present_mean2:= clim_present_mean2 * lc_suit_present_mean * density_mean]
       
    ##rid of zeros
    ## make a  raster of arena
    res6b<-res6[!is.na(realised_present_mean2),]### WHAT ABOUT NEW AREAS IN THE FUTURE
    res6b<-res6b[!(realised_present_mean2==0),]
    
    ##run gas model - any point in variation not going to influence spatial contact patterns
    #res6b[,present_realisedT:= c(((realised_present_meanH*ttmean)*(speed_mean*ttspeed)*(d_mean*ttd)) * (realised_present_meanV*abs(speed_mean_vector*ttspeedv))*(d_mean_vector*ttdv) *(secondary_present) * ((humans2010/5.6) * (0.0001*ttdh) * (5*ttspeedh))),]
    
	#### END OF HAZARD ####
	
	### SUMMARIZE AND COMBINE INVARIANT DATA

	### summarise to admin level 3 units
	### is mean correct - any other measure of change?
	res8<- aggregate(res7,by=list(res7$adml3,res$RCP,res7$year),mean)
	
	### combine with per admin level 3 EXPOSURE (population size) AND VULNERABILITY (POVERTY FRACTION)
	setkey(res8,adl3)
	ad2d<-ad2$data
	setkey(ad2d,xxx)
	
	res8<-res8[ad2d]
	
	### need CFR per disease
	### need spillover rate per disease
	### get disease_table30 cbinded on
	dis1<-diseases[,]
	res8<-cbind(res8,dis1)
	
	###res 8 is in long format with years and rcp differences PER ADMIN AREAS (SO ONLY A COUPLE OF HUNDRED LINES)
	
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

	## final hazard
	res8[,res8$Hazard:=realised_future_mean2 * CFR * spillover_rate] ## or without CFR and spill over and divide up into cost groups
	
	## final exposure
	## could include cattle??
	## could cut by at risk groups
	res8[,res8$Exposure:=sum_human_population_future]

	## final Vulnerability
	res8[,res8$Vulnerability:=pov_frac_future]
	
	## final risk (units are deaths with CFR - cases without)
	res8[,res8$Risk:=Hazard * Vulnerability * Exposure]

	###do again for present day
	res8[,res8$RiskP:=realised_present_mean2 * CFR * spillover_rate * sum_human_population * pov_frac]
	
	## delta risk
	res8[,delta_risk:= RiskP-Risk]
	
	##create unique
	res8[,year_RCP=paste(RCP,year,sep="_")]
	year_RCP<-unique(res8$year_RCP

					
	  ##when best breaks are found 
    pdf(file=paste("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\risk_maps\\risk_map_",pres2[i],"_",sample(1:1000,1),".pdf",sep=""),width=15,height=10)
    
	##create multi panel facet plot of maps - cut to endemic area??
	ggplot(ad1b,aes(
			ggpolygon()
	
	
	dev.off()
	
	### or create panel plot
	
	### also create database of all diseases
	
	if(i==1) {res10<-res9} else {res10<-rbind(res10,res9)}

	 
  #gc()
  print(i)
}

##save results
write.csv(res10,file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\final_risk1.csv")

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



