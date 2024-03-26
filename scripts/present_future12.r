

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

library(ggplot2)
library(velox)

load(file="X:\\DRtemp\\babesia_problem.r")

##read in popGDP
popGDP<-stack(list.files("C:\\Users\\xxxx\\Downloads\\global_population_and_gdp\\data", pattern=".tif",full.names=TRUE))
yyy<-read.table(text=names(popGDP),sep="_",stringsAsFactors = FALSE)
yyy$layers<-1:nrow(yyy)
yyy$names<-names(popGDP)
popGDP<-subset(popGDP,yyy[yyy$V2 %in% c(2010,2030,2050,2070,2080),"layers"])
popGDP<-velox(popGDP)
  

###get popGDP
df1<-data.frame(cell.id=1:ncell(pres_ad1),sort=values(pres_ad1))
df1<-df1[!is.na(df1$sort),]
sp1<-SpatialPoints(xyFromCell(template,df1$cell.id))
popGDP2<-popGDP$extract_points(sp1)
popGDP2<-as.data.frame(popGDP2)
names(popGDP2)<-yyy[yyy$V2 %in% c(2010,2030,2050,2070,2080),"names"]
setDT(popGDP2)
popGDP3<-melt(popGDP2)
popGDP3[,c("SSP","year"):=tstrsplit(variable,split="_",fixed=TRUE)]
popGDP3$cell.id<-df1$cell.id

popGDPfuture<-dcast(popGDP3,SSP~year,value.var = "value")

  
##read in admins to move away from GRID
ad1<-readOGR("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\admin1\\ne_10m_admin_1_states_provinces.shp","ne_10m_admin_1_states_provinces")

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

##create rasterized admin
pres_ad1<-fasterize(st_as_sf(ad1),template,field="diss_me")

###load GP4 - pop gdp
load(file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\GP4.r")
##merge with ad1
ad1@data<-cbind(ad1@data,as.data.frame(GP4))


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

###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl2),template)

###create new raster
t2<-ws2
t2[!is.na(t2)]<-0

res3<-NULL

##loop
i=1
for (i in i:length(pres1)){
  
    ##load a disease
  load(pres1[i])
  
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
  pd<-point_data[point_data$LU %in% pres2[i] & point_data$Status!="Denied",]
  if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd)} else {rw1<-1000}
  
  ##bv problem
  if(pres2[i]==paste(" ",bv,sep="")){pres2[i]<-" babesia venatorum"}
  
  ##create endemic region
  dis_trans<-d1[d1$name2==pres2[i],]
  CFR=mean(dis_trans$CFR.low,dis_trans$CFR.high,na.rm=TRUE)
  ccc<-strsplit(dis_trans$countries,",")[[1]]
  if(is.na(ccc[1])){print("no recognised countries");print(pres2[i]);next}
  countr<-wrld_simpl[wrld_simpl$ISO2 %in% ccc ,]
  
  ##spill over rate cases per year over population at risk 2005
  spillover_rate<-as.numeric(dis_trans$cases_per_year)/sum(countr$POP2005,na.rm=TRUE)
    
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
    }#else{
    ##if not less than two ad2 the "TRUE" sample locations are just the countries
    #ad2<-countr  
  #}
  
  ###determine if in regn and regn2
  
	pres_ad2<-fasterize(st_as_sf(countr),template,field="UN")
	regn11<-fasterize(st_as_sf(regn),template,field = "UN")
	regn22<-fasterize(st_as_sf(regn2),template,field = "UN")
	#diss_me2<-fasterize(st_as_sf(ad1),template)
	
	res6[,countr:=extract(pres_ad2,res6$cell.id)]
	res6[,sub_region:=extract(regn11,res6$cell.id)]
	res6[,region:=extract(regn22,res6$cell.id)]
	res6[,diss_me:=extract(pres_ad1,res6$cell.id)]

	####load rasters of popgdp
	#popgdp<-extract(popdgp,res6$cell.id)
	#popgdp<-melt(popgdp)
	
	### rid of any lines no inside	
	### regn
	### countries
	res6<-res6[!is.na(res6$sub_region),]
	#res6b<-res6[!is.na(res6$region),]
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
    res6[,clim_present_mean2:= c(1-((1-clim_present_mean) *(1-clim_present_mean_vector))),]  
    res6[,lc_suit_present_mean:= lc_suit_present_mean_vector,]  
    res6[,density_mean:= density_mean_vector,]  
    res6[,clim_future_mean2:= c(1-((1-clim_future_mean) *(1-clim_future_mean_vector))),]  
    res6[,lc_suit_future_mean:= lcsuit_future_mean_vector,]  
    #res6[,density_mean:= density_mean_vector,] 

   }else{
    res6[,clim_present_mean2:= clim_present_mean,]### change to fit in below
     res6[,clim_future_mean2:= clim_future_mean,]### change to fit in below
     
    }
    	
	
	### WHAT ABOUT FUTURE??????????????????
		
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
	res6[,HazardUW := clim_present_mean2 * lc_suit_present_mean * density_mean]
	
	###create realised present mean - abudance from habitat????
	res6[,Hazard_future_UW:= clim_future_mean2 * lcsuit_future_mean * density_mean]
	
	## final hazard
	## CFR needs to be greater than 0 - to approximate dalys
	## can we estimate DALYS (proportion of life lost * 1-CFR etc.)
	res6[,Hazard:=HazardUW * (CFR+0.001) * spillover_rate] ## or without CFR and spill over and divide up into cost groups
	
	## final hazard
	res6[,Hazard_future:=Hazard_future_UW * (CFR+0.001) * spillover_rate] ## or without CFR and spill over and divide up into cost groups
	
	
	##create unique
	res6[,year_RCP:=paste(RCP,year,sep="_")]
	#year_RCP<-unique(res6$year_RCP)
	
	#Can do the admin3 at the same time as RCP
	#### AGGREGATE TO ADMIN3
	#### DIFFERENT VERSIONS??
	
	### this doesn't work as admin might mean to zero
	### turn 0 to NA?
	### if some up and some down then no change??
	### unless we use threshold...???
	
	### does it make sense to do hazard at grid cell and exposure at admin??
	
	## is na country for just reported country trends
	##do this for after grid not now and take trimmed mean

	res7<-res6[!is.na(sub_region),lapply(.SD,function (x) mean(x,trim=0.1,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=c("Hazard","HazardUW","Hazard_future","Hazard_future_UW")]
	setkey(res7,diss_me)

	###match to
	res7<-res7[match(ad1$diss_me,res7$diss_me),]
	
	#### END OF HAZARD ####
	### SUMMARIZE AND COMBINE INVARIANT DATA
	res8<-cbind(res7,GP4)

	
	## final exposure
	## could include cattle??
	## could cut by at risk groups
	res8[,Exposure:=sum_human_population_future]

	## final Vulnerability
	res8[,Vulnerability:=pov_frac_future]
	
	## final risk (units are deaths with CFR - cases without)
	res8[,Risk:=Hazard * Vulnerability * Exposure]

	###do again for present day
	res8[,RiskP:=HazardUW * CFR * spillover_rate * sum_human_population * pov_frac]
	
	## delta risk
	res8[,delta_risk:= RiskP-Risk]
	
	##combine with dis_trans
	res8<-cbind(res8,dis_trans)
	
  ## plot result
	ad1b<-ad1
  ad1b@data<-cbind(ad1@data,as.data.frame(res8))
  
  ##prepare for ggplot
  #ad1b@data$id<-rownames(ad1b@data)
  #ad1bpoints<-fortify(ad1b,region="id")
  #ad1bDF<-merge(ad1bpoints,ad1b@data,by="id")
  #ad1bDF<-ad1bDF[,-(107:116)]
  
  ###just convert to sf and then can used natively
  
  ad1c<-st_as_sf(ad1b)
  
					
	##when best breaks are found 
  pdf(file=paste("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\risk_maps\\risk_map_",pres2[i],"_",sample(1:1000,1),".pdf",sep=""),width=15,height=10)
    
	  ggplot(data = ad1c) +
    geom_sf(aes(fill = delta_risk)) +
    scale_fill_viridis_c(option = "plasma", trans = "sqrt")+
	facet_wrap(.~year_RCP)
  	
	dev.off()
	
	### or create panel plot
	
	### also create database of all diseases
	
	if(i==1) {res9<-res8} else {res9<-rbind(res9,res8)}

	 
  #gc()
  print(i)
}

##save results
write.csv(res9,file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\final_risk1.csv")

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



