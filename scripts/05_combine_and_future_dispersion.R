
###run disease validation

# Load libraries -----
library(biomod2)
library(dismo)
library(maptools)
library(data.table)
library(sp)
library(rgeos)
library(raster)
library(fasterize)
library(dismo)
library(sf)
library(rgdal)
library(viridis)
library(sp)
library(ggplot2)
library(velox)
library(smoothr)

# Import coastlines and livestock data ----
download.file("https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/ne_10m_coastline.zip", 
              destfile = 'coastlines.zip')

# unzip the file
unzip(zipfile = "coastlines.zip", 
      exdir = 'ne-coastlines-10m')

# load the data - rgdal version
#coastlines <- readOGR("data/ne_10m_coastline/ne_10m_coastline.shp")
# load the data - sf version
coastlines <- st_read("data/ne_10m_coastline/ne_10m_coastline.shp")
plot(coastlines)

# load function files
source("scripts/functions/do_auc.R")
source("scripts/functions/rast_to_range.r")

# read livestock data

# List all the TIFF files in the directory
tif_files <- list.files("data/global_IUCN", pattern="_Da.tif", full.names = TRUE)

# Check if there are any files
if (length(tif_files) == 0) {
  stop("No TIFF files found matching the pattern.")
}

lvst <- stack(list.files("data/global_IUCN",pattern="_Da.tif",full.names=TRUE))
names(lvst) <- c("buffalo","chickens","cattle","ducks","goats","horses","pigs","sheep")

# read in admins to move away from GRID
ad1 <- st_read("data/admin1/ne_10m_admin_1_states_provinces.shp",
             "ne_10m_admin_1_states_provinces")

# Read in PREDICTS ----
# base raster template
template <- raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),
                   crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# defaunation surface
predicts <- raster("data/global_IUCN/predicts_resampled_template.tif")

# create rasterized admin
pres_ad1 <- fasterize(st_as_sf(ad1),template,field="diss_me")

# read in two country points
lesstwo <- read.csv("data/disease_in_two_or_fewer_countries4.csv",stringsAsFactors=FALSE)
lesstwo$name2 <- paste(" ",lesstwo$name,sep="")
lesstwo$dummy <- 1
coordinates(lesstwo) <- ~Longitude+Latitude


# Read in empress-i point data ----
point_data <- read.csv("data/all_empress_i_data4.csv",stringsAsFactors = FALSE)
point_data <- point_data[!is.na(point_data$SumCases),]
#point_data$LU<-paste(" ",point_data$name_LU,sep="")

# Read disease data ----
d1<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\disease_table32c.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4

###do or not
d2<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\diseases4a.csv",stringsAsFactors=FALSE)
d2a<-d2[!is.na(d2$disease),]
d2a$level[is.na(d2a$level)]<-999
d2b<-d2a[d2a$level>=1 & d2a$level!=999,]
nums<-unique(d2b$combine)

###find all present day
pres1<-list.files("D:\\Users\\Public\\Documents\\per_disease3\\",pattern=" ALL4",full.names=TRUE)
pres2<-list.files("D:\\Users\\Public\\Documents\\per_disease3\\",pattern=" ALL4",full.names=FALSE)
pres2<-gsub("  ALL4.r","",pres2)
pres2<-gsub(" ALL4.r","",pres2)


###find all present day
iucn1<-list.files("D:\\Users\\xxxx\\Documents\\disease_host_ranges\\",pattern="3.tif",full.names=TRUE)
iucn2<-list.files("D:\\Users\\xxxx\\Documents\\disease_host_ranges\\",pattern="3.tif",full.names=FALSE)
iucn2<-gsub("3.tif","",iucn2)

##load world simple
wrld_simpl3<-readOGR("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\wrld_simpl3.shp","wrld_simpl3")
#wrld_simpl3[wrld_simpl3$NAME %in% c("Russia","France"),]
wrld_simpl3@data[wrld_simpl3$NAME=="Russia","SUBREGION"]<-wrld_simpl3@data[wrld_simpl3$NAME=="Ukraine","SUBREGION"]
#wrld_simpl3[wrld_simpl3$NAME %in% c("Siberia","Mongolia"),]
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","SUBREGION"]<-30
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","REGION"]<-142
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","ISO2"]<-"XX"
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","FIPS"]<-"XX"
wrld_simpl3@data[wrld_simpl3$NAME=="Chile","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Paraguay","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Uruguay","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Argentina","SUBREGION"]<-4

###find neighbours
neigh1a<-read.csv(file="C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\country_borders.csv",stringsAsFactors = FALSE)
neigh1a<-neigh1a[!is.na(neigh1a$country_border_code),]

###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl3),template)

###create new raster
t2<-ws2
t2[!is.na(t2)]<-0

###start loop
for (i in sample(nums)){

  #get data
  d2c<-d2b[d2b$combine==i,]
  
  ##alter name
  if(d2c$disease[1]=="angiostrongylus costaricensis "){d2c$disease[1]<-"angiostrongylus costaricensis"}
  
  ###check in dataframes
  d2c<-d2c[paste0(d2c$disease) %in% pres2,]
  if(nrow(d2c)==0){next}
  
    ##nake name from muliple
  nam_dis<-paste(d2c$disease,collapse = "_")


  ##check what has been processed
   been_done<-list.files("D:\\Users\\xxxx\\Documents\\gainlossresults\\",pattern="_2.csv")
  been_done<-gsub("_X_",";",been_done)
  been_done<-read.table(text=been_done,sep=";")[,1]
   
   if(nam_dis %in% been_done){next}#else {break}}
  

  ## make livstock
  lvsth<-subset(lvst,(1:nlayers(lvst))[names(lvst) %in% (strsplit(d2c$all_host,",")[[1]])])
  if(nlayers(lvsth)>1){lvsth<-calc(lvsth,sum,na.rm=TRUE)}
  if(nlayers(lvsth)>0){lvsth2<-resample(lvsth,template,method='ngb');lvsth[lvsth<10]<-NA;lvsth[lvsth>10]<-1}
  
  ##get IUCN ranges
  sp1<-iucn1[ iucn2 %in% d2c$disease]
  
  if(length(sp1)==0){
    
    if(nlayers(lvsth)>0){iucn4<-lvsth2}else{print("problem with IUCN");print(d2c$disease);iucn4<-ws2}
    
  }else{
    
    if(length(sp1)==1){iucn4<-raster(sp1)
    }else{
      iucn3<-stack(sp1)
      iucn4<-calc(iucn3,sum,na.rm=TRUE)
    }
    
    if(nlayers(lvsth)>0){iucn4<-calc(stack(iucn4,lvsth2),sum,na.rm=TRUE)}
    
  }
  
    #load a disease
  if(nrow(d2c)==1){load(pres1[pres2==paste0("",d2c$disease)[1]]);res6[,name:=d2c$disease[1]]} else {
    
    for (rr in 1:nrow(d2c)){
      
      if(d2c$disease[rr]=="babesia venatorum"){d2c$disease[rr]<-"babesia venatorum";load(pres1[pres2==d2c$disease[rr]])}
      
      load(pres1[pres2==paste0("",d2c$disease[rr])[1]])
      
      
      res6[,name:=d2c$disease[rr]]
      
      if(rr==1) {res6z<-res6} else {res6z=rbindlist(list(res6z,res6))}
      rm(res6);gc()
    }
    
    ##agragte by name to combine
    res6<-res6z[,                       
                lapply(.SD, function (x) mean(x,na.rm=TRUE)),  ## compute the mean
                by = .(year,RCP,cell.id),           ## for every 'year,RCP,cell.id'
                .SDcols = names(res6z)[!names(res6z) %in% c("year","RCP","cell.id","name")]] ## for just those specified in .SDcols
    rm(res6z);gc()
  }
  
  
  
  #remove NAs  
  res6<-res6[!is.na(clim_present_mean),]
  ## test res6 is ok  
  if(sd(res6$clim_present_mean,na.rm=TRUE)==0){print("no data");break}
  if(nrow(res6)==0){print("no data set present");break}
  if(!"cell.id" %in% names(res6)){print("no cell.id present");next}

  ### just need present day
  res6<-cbind(res6,xyFromCell(template,res6$cell.id))
  

  ##check point data # point_data$LU  %in% d2c$disease
  pd<-point_data[(point_data$name_LU  %in% d2c$disease) & point_data$Status!="Denied" & point_data$SumCases>0,]
  if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd);projection(pd)<-projection(raster())} else {rw1<-1000}
  
  ##create endemic region
  dis_trans<-d1[d1$name2 %in% paste0(" ",d2c$disease),]
  CFR=mean(c(dis_trans$CFR.low,dis_trans$CFR.high),na.rm=TRUE)
  ccc<-strsplit(paste(dis_trans$countries,collapse = ","),",")[[1]] ##deal with mulitple ROWS ##unique
  if(is.na(ccc[1])){print("no recognised countries");print(pres2[i]);next}
  countr<-wrld_simpl3[wrld_simpl3$ISO2 %in% ccc ,]
  
  ##distance of points
  centrd<-colMeans(coordinates(countr))
  cent2<-template
  values(cent2)<-NA
  cent2[cellFromXY(template,centrd)]<-1
  cent3<-distance(cent2)
  
  ###test for county level data
  if(nrow(countr)!=length(ccc)){print("issue with countries");print(ccc)}
  if(nrow(countr)==0){print("no recognised countries");print(pres2[i]);next}
  
  neigh2<-neigh1a[neigh1a$country_code %in% as.character(countr$ISO2) , ]

  regn<-wrld_simpl3[wrld_simpl3$SUBREGION %in% countr$SUBREGION,] ##psuedoabsence region
  regn2<-wrld_simpl3[wrld_simpl3$REGION %in% countr$REGION,] ##psuedoabsence region
  regn3<-wrld_simpl3[wrld_simpl3$ISO2 %in% c(neigh2$country_border_code,as.character(countr$ISO2)),] ##region of interest
  
  ###if coveerage of reporting countries is too high switch to larger region  <---- artribrary
  regnPA<-regn ##region for psuego absences
  if(nrow(countr)/nrow(regn)>=(2/3)){regnPA<-regn2}

  ## but if not countries in subregion 
  ###if less than two countries
  
  l2<-lesstwo[lesstwo$name2  %in% paste0(" ",d2c$disease),]
  if(nrow(l2)>0){
    projection(l2)<-projection(ad1)
    ##overlay
    overpoints<-over(ad1,l2)
    ad2<-ad1[!is.na(overpoints$dummy), ] ### used to set present points samples
    ad2$dummy<-1:nrow(ad2)
    
  }else{
    ##if not less than two ad2 the "TRUE" sample locations are just the countries
    ad2<-countr    ### used to set present points samples
    ad2$dummy<-1:nrow(ad2)
  }
  
  ##intersection between hosts and wildlife
  ad2r<-fasterize(st_as_sf(ad2),template)
  host_and_present<-mask(iucn4,ad2r)
  
  ###determine if in regn and regn2
  pres_ad2<-fasterize(st_as_sf(ad2),template,field="dummy")
  presc<-fasterize(st_as_sf(countr),template,field = "UN") ##either  subregion or region depending
  
  regn11<-fasterize(st_as_sf(regn),template,field = "UN") ##either  subregion or region depending
  regn22<-fasterize(st_as_sf(regnPA),template,field = "UN")###region for PA needed?
  regn33<-fasterize(st_as_sf(regn3),template,field = "UN")###region for PA needed?
  
  ### get country and region info
  res6[,subcountr:=extract(pres_ad2,res6$cell.id)]
  res6[,countr:=extract(presc,res6$cell.id)]
  
  res6[,sub_region:=extract(regn11,res6$cell.id)] ##largest
  res6[,regionPA:=extract(regn22,res6$cell.id)]  ## PA area 
  res6[,surrounding:=extract(regn33,res6$cell.id)]  ## - surrounding countries
  #res6[,diss_me:=extract(pres_ad1,res6$cell.id)]
  
  ###cut do to coarses area
  res6<-res6[!is.na(res6$sub_region),]
  
  ###test host range
  res6[,host_range:=extract(iucn4,res6$cell.id)]

  ##distance to centre of known country
  res6[,dist_cent:=extract(cent3,res6$cell.id)]
  res6[,dist_cent:=1-(dist_cent/max(dist_cent))]


  ##create unique
  res6[,year_RCP:=paste(RCP,year,sep="_")]

  ##spill over rate cases per year over population at risk 2005
  #spillover_rate<-as.numeric(dis_trans$cases_per_year)/(nrow(res6)/length(unique(res6$year_RCP)))
  spillover_rate<-as.numeric(sum(dis_trans$cases_per_year),na.rm=TRUE)/sum(countr$POP2005,na.rm=TRUE)
  
  #### CREATE HAZARD DIFFERENCES
  
  if(sd(res6$clim_present_mean_vector,na.rm=TRUE)>0){
    res6$clim_present_mean_vector[res6$clim_present_mean<0.001]<-1 ##do we want this?
    res6[,clim_present_mean2:= c(1-((1-clim_present_mean) *(1-clim_present_mean_vector))),]  
    
    res6[,clim_future_mean2:= c(1-((1-clim_future_mean) *(1-clim_future_mean_vector))),]  
  
      ##there is no vector land-use data so use host land-se
    res6[,lc_suit_present_mean2:= lc_suit_present_mean_vector * lc_suit_present_mean ,]
    res6[,lc_suit_future_mean2:= lcsuit_future_mean_vector * lcsuit_future_mean,] 
    
  }else{
    res6[,clim_present_mean2:= clim_present_mean,]### change to fit in below
    res6[,clim_future_mean2:= clim_future_mean,]### change to fit in below
    res6[,lc_suit_present_mean2:= lc_suit_present_mean ,]  
    res6[,lc_suit_future_mean2:= lcsuit_future_mean,]
    }
  
   #set real range as intersection of country data and host ranges
  res6[,countr2:=ifelse(!is.na(subcountr) & !is.na(host_range),1,NA)]
  
  ###set dummy
  res6[,dummy:=1]
  
  ##reduce to one instance
  res6a<-res6[res6$year_RCP==unique(res6$year_RCP)[1] & !is.na(res6$regionPA),]
 
  ##threshold again
 
      for(ww in 1:50){
      
        ##set randoms
        timesP=sample(c(100,500,1500,5000,10000),1)
        timesA=sample(c(500,1500,5000,10000),1)
    
        weighting=sample(seq(from=0.1, to=0.9, by=0.1),1)
        
        ##make Hazard
        res6a[,Hazard := (weighting*clim_present_mean2) +((1-weighting) * lc_suit_present_mean2) ]
        
        thr1<-do_auc(data1=res6a,test_column="Hazard",group_column="countr2",weights="dist_cent",timesP=timesP,timesA=timesA,summarise=FALSE)
    
        rt1<-data.frame(thr1,timesA=timesA,timesP=timesP,weighting=weighting,column="centre_weighted")
        #rt1b<-data.frame(thr2,timesA=timesA,timesP=timesP,weighting=weighting,column="random")
    
       if(ww==1){rt2<-rt1} else {rt2<-rbind(rt2,rt1)}
       #print(ww)
      }
   
    ##get best values
    dt1<-data.frame(weighting=aggregate(rt2$AUC_best.stat,by=list(rt2$weighting),mean)$Group.1,mean=aggregate(rt2$AUC_best.stat,by=list(rt2$weighting),mean)[,"x",drop=TRUE],sd=aggregate(rt2$AUC_best.stat,by=list(rt2$weighting),sd)[,"x",drop=TRUE])
    dt1=cbind(dt1,dt1[dt1$mean==max(dt1$mean),c("mean","sd")])
    names(dt1)[4:5]<-c("maxmean","maxsd")
    dt1$twosdless<-dt1$maxmean-(2*dt1$maxsd)
    dt2<-dt1[dt1$mean>dt1$twosdless,]
    
    ##set final weighting
    dis_trans$final_weighting=mean(dt2$weighting)
    
    ##make Hazard
    res6[,Hazard := (dis_trans$final_weighting*clim_present_mean2) +((1-dis_trans$final_weighting) * lc_suit_present_mean2) ]
    
    ##remove very bad cut-offs
    rt3<-rt2[rt2$AUC_best.stat>quantile(rt2$AUC_best.stat,0.50) & rt2$weighting %in% dt2$weighting,]
    
    ##put in main data frame
    dis_trans$AUC_cutoff<-mean(rt3$AUC_cutoff)
    dis_trans$AUC<-mean(rt3$AUC_best.stat)
 
    ##put in main data frame
    dis_trans$AUC_cutoffmax<-max(rt3$AUC_cutoff)
    dis_trans$AUCmax<-max(rt3$AUC_best.stat)
  
    ##put in main data frame
    dis_trans$AUC_cutoffmin<-min(rt3$AUC_cutoff)
    dis_trans$AUCmin<-min(rt3$AUC_best.stat)
  

    #real points were here
    dis_trans$AUC_cutoff_real<-NA;dis_trans$AUC_real<-NA # end of if pb
  
  #without future land-use change
  res6[,Hazard_future := (dis_trans$final_weighting*clim_future_mean2) +((1-dis_trans$final_weighting) * lc_suit_present_mean2) ]

  ##ensure the cutoffs are not too strict
  if(dis_trans$AUC_cutoffmax>quantile(res6$Hazard,0.975,na.rm=TRUE)){print("high bounded");dis_trans$AUC_cutoffmax<-quantile(res6$Hazard,0.95,na.rm=TRUE)}
  if(dis_trans$AUC_cutoffmin<quantile(res6$Hazard,0.025,na.rm=TRUE)){print("low bounded");dis_trans$AUC_cutoffmin<-quantile(res6$Hazard,0.05,na.rm=TRUE)}
  
  ### percell thresholds
  res6[,HazardT:=ifelse(Hazard>=dis_trans$AUC_cutoff,1,0)]
  res6[,Hazard_futureT:=ifelse(Hazard_future>=dis_trans$AUC_cutoff,1,0)]
  res6[,HazardTmin:=ifelse(Hazard>=dis_trans$AUC_cutoffmin,1,0)]
  res6[,Hazard_futureTmin:=ifelse(Hazard_future>=dis_trans$AUC_cutoffmin,1,0)]
  res6[,HazardTmax:=ifelse(Hazard>=dis_trans$AUC_cutoffmax,1,0)]
  res6[,Hazard_futureTmax:=ifelse(Hazard_future>=dis_trans$AUC_cutoffmax,1,0)]
  
  ### compare to combined scenario
  unis<-unique(res6$year_RCP)

  ##what columns?
  cols1<-c("climT","clim_futureT","HazardT","Hazard_futureT")
  titls<-c(paste(cols1[1:2],round(mean(dis_trans$AUC),2),sep=" "),paste(cols1[3:4],round(mean(dis_trans$AUC),2),sep=" "))
    
    yy=1
    zz=3
    res9<-res6[res6$year_RCP==unis[yy] ,]
    
    res9<-as.data.frame(res9)

    pl1<-ws2
    pl1[pl1==1]<-0
    
    ##remove cells outside intersection of countries and host ranges
    res9b<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding),]
    
    ##populate with HazardT
    pl1[res9b$cell.id]<-res9b[ ,names(res9b)[names(res9b)=="HazardT"],drop=TRUE]
    
    #if(zz==9){xxx<-(log(res9$Hazard+1))-(log(res9$Hazard_future+1));xxx[xxx==0]=NA;pl1[res9$cell.id]<-xxx}
    #crop to general area
    pl1<-crop(pl1,extent( min(res9b$x)-5, max(res9b$x)+5,min(res9b$y)-5, max(res9b$y)+5))
    pl1[pl1==0]<-NA
    
    ##make range
    poly1<-rast_to_range(x=pl1,present=host_and_present,crumb_size=1e9)
    pl2<-pl1
    rm(pl1)
    
    ##plot to see where the data are MINNNNNNN
    pl1<-ws2
    pl1[pl1==1]<-0
    
    ##remove cells outside intersection of countries and host ranges
    res9b<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding),]
    
    ##populate with HazardT
    pl1[res9b$cell.id]<-res9b[ ,names(res9b)[names(res9b)=="HazardTmin"],drop=TRUE]
    
    #crop to general area
    pl1<-crop(pl1,extent( min(res9b$x)-5, max(res9b$x)+5,min(res9b$y)-5, max(res9b$y)+5))
    pl1[pl1==0]<-NA
    
    ##make range
    ee<-simpleError("no polygon overlap")
    
    ##make range
    poly1min<-tryCatch(rast_to_range(x=pl1,present= host_and_present,crumb_size=1e9),error=function(ee) ee)
    
    if(class(poly1min)[1]=="simpleError"|is.na(poly1min)){poly1min<-poly1}
    rm(pl1)
    
    ##plot to see where the data are MAXXXX
    pl1<-ws2
    pl1[pl1==1]<-0
    
    ##remove cells outside intersection of countries and host ranges
    res9b<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding),]
    
    ##populate with HazardT
    pl1[res9b$cell.id]<-res9b[ ,names(res9b)[names(res9b)=="HazardTmax"],drop=TRUE]
    
    #crop to general area
    pl1<-crop(pl1,extent( min(res9b$x)-5, max(res9b$x)+5,min(res9b$y)-5, max(res9b$y)+5))
    pl1[pl1==0]<-NA
    
    ee<-simpleError("no polygon overlap")
    
    ##make range
    
    poly1max<-tryCatch(rast_to_range(x=pl1,present= host_and_present,crumb_size=1e9),error=function(ee) ee)
    
    if(class(poly1max)[1]=="simpleError"|is.na(poly1max)){poly1max<-poly1}
    rm(pl1)

    ##plot to see where the data are BLANK
    #plot(wrld_simpl)
    pl1<-ws2
    pl1[pl1==1]<-0
    
    #crop to general area
    pl1<-crop(pl1,extent( min(res9b$x)-5, max(res9b$x)+5,min(res9b$y)-5, max(res9b$y)+5))
    pl1[pl1==0]<-NA
    
    ## add buffer
    poly2030<-gBuffer(poly1,width=1/7)
    poly2050<-gBuffer(poly1,width=4/7)
    poly2070<-gBuffer(poly1,width=6/7)
    poly2080<-gBuffer(poly1,width=1)
    poly2030i<-gBuffer(poly1,width=-1/7)
    if(length(poly2030i)==0){poly2030i<-poly2030}
    poly2050i<-gBuffer(poly1,width=-4/7)
    if(length(poly2050i)==0){poly2050i<-poly2030i}
    poly2070i<-gBuffer(poly1,width=-6/7)
    if(length(poly2070i)==0){poly2070i<-poly2050i}
    poly2080i<-gBuffer(poly1,width=-1)
    if(length(poly2080i)==0){poly2080i<-poly2070i}
    
    
    ## add buffer but how far???
    poly2030min<-gBuffer(poly1min,width=1/7)
    poly2050min<-gBuffer(poly1min,width=4/7)
    poly2070min<-gBuffer(poly1min,width=6/7)
    poly2080min<-gBuffer(poly1min,width=1)
    poly2030i<-gBuffer(poly1min,width=-1/7)
    poly2030mini<-gBuffer(poly1min,width=-1/7)
    if(length(poly2030mini)==0){poly2030mini<-poly2030min}
    poly2050mini<-gBuffer(poly1min,width=-4/7)
    if(length(poly2050mini)==0){poly2050mini<-poly2030mini}
    poly2070mini<-gBuffer(poly1min,width=-6/7)
    if(length(poly2070mini)==0){poly2070mini<-poly2050mini}
    poly2080mini<-gBuffer(poly1min,width=-1)
    if(length(poly2080mini)==0){poly2080mini<-poly2070mini}
    
    ## add buffer but how far???
    poly2030max<-gBuffer(poly1max,width=1/7)
    poly2050max<-gBuffer(poly1max,width=4/7)
    poly2070max<-gBuffer(poly1max,width=6/7)
    poly2080max<-gBuffer(poly1max,width=1)
    poly2030maxi<-gBuffer(poly1max,width=-1/7)
    if(length(poly2030maxi)==0){poly2030maxi<-poly2030max}
    poly2050maxi<-gBuffer(poly1max,width=-4/7)
    if(length(poly2050maxi)==0){poly2050maxi<-poly2030maxi}
    poly2070maxi<-gBuffer(poly1max,width=-6/7)
    if(length(poly2070maxi)==0){poly2070maxi<-poly2050maxi}
    poly2080maxi<-gBuffer(poly1max,width=-1)
    if(length(poly2080maxi)==0){poly2080maxi<-poly2070maxi}
    
########### <---------    
  ##based on thomas et al (science 2011) its mean 16.9km or four cells per decade so 7 decades times four is 28 but that is for 2080 #should prevent from going too quick - what about each host at time?

    ##do for all decades
    resx<-data.frame(x=res6$x,y=res6$y,dummy=1,stringsAsFactors = FALSE)
    coordinates(resx)<-~x+y
    projection(resx)<-projection(poly1)

    ##extract for each point
    res6[,present2010:=over(resx,poly1)[,1]]
    res6[,present2010min:=over(resx,poly1min)[,1]]
    res6[,present2010max:=over(resx,poly1max)[,1]]

    for (ee in c(2030,2050,2070,2080)){
      
      res6t<-res6[res6$year==ee,]
      
      res6t[,presentfuture:=over(resx[res6$year==ee,],get(paste0("poly",ee)))]
      res6t[,presentfuturemin:=over(resx[res6$year==ee,],get(paste0("poly",ee,"min")))]
      res6t[,presentfuturemax:=over(resx[res6$year==ee,],get(paste0("poly",ee,"max")))]
      res6t[,presentfuturei:=over(resx[res6$year==ee,],get(paste0("poly",ee,"i")))]
      res6t[,presentfuturemini:=over(resx[res6$year==ee,],get(paste0("poly",ee,"min","i")))]
      res6t[,presentfuturemaxi:=over(resx[res6$year==ee,],get(paste0("poly",ee,"max","i")))]
      
            
      if(ee==2030){resF=res6t}else{resF<-rbindlist(list(resF,res6t))}
      gc()
    }  
    
    ##write main data
    fwrite(resF,file=paste("D:\\Users\\xxxx\\Documents\\datasets1\\",nam_dis,"_",sample(1:1000,1),"_all_data3.csv",sep=""))
    
    ##write gain loss
    fwrite(dis_trans,file=paste("D:\\Users\\xxxx\\Documents\\gainlossresults\\",nam_dis,"_",sample(1:1000,1),"_2.csv",sep=""))
    
   print("DONE")
} ##end of nums loop