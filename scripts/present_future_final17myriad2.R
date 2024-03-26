
###run disease validation
#library(biomod2)
library(dismo)
library(maptools)
library(data.table)
library(sp)
library(rgeos)

library(raster)
library(fasterize)
#library(dismo)
library(sf)
library(rgdal)

#library(viridis)
#library(sp)

library(ggplot2)
#library(velox)

library(smoothr)

#download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", 
#              destfile = 'coastlines.zip')

# unzip the file
#unzip(zipfile = "coastlines.zip", 
#      exdir = 'ne-coastlines-10m')

# load the data 
coastlines <- readOGR("/lustre/scratch/scratch/ucbtdw0/ne-coastlines-10m/ne_10m_coastline.shp")

source("/home/ucbtdw0/Scratch/New_Global_MAXENT/do_auc.R")
source("/home/ucbtdw0/Scratch/New_Global_MAXENT/rast_to_range.R")


##read livestock data
#CHECK
lvst<-stack(list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/data/",pattern="_Da.tif",full.names=TRUE))
names(lvst)<-c("buffalo","chickens","cattle","ducks","goats","horses","pigs","sheep")

##read in admins to move away from GRID
ad1<-readOGR("/home/ucbtdw0/Scratch/New_Global_MAXENT/admin1/ne_10m_admin_1_states_provinces.shp","ne_10m_admin_1_states_provinces")


#load(file="X:/DRtemp/babesia_problem.r")
#source("/home/ucbtdw0/Scratch/R_scripts/TSS.r")

##load popgdp
#load(file="V:/GPpres.r")
#setkey(GPpres,diss_me)
#load(file="V:/GP7.r")
#setkey(GP7,year_diss_me)
#load(file="C:/Users/Public/Documents/globalIUCN/popgdp7.r")

###make 2010
#popgdp7[,c("year","dissme"):=tstrsplit(year_diss_me,"_")]
#popgdp2010b<-popgdp7[year==2010,c(1,2,5)]
#names(popgdp2010b)[2:3]<-c("gdp2010","humans2010")
#max1<-quantile(popgdp2010b$gdp2010/popgdp2010b$humans2010,0.95,na.rm=TRUE)
#popgdp2010b[,gdp_per_person:=(gdp2010/humans2010)/max1]
#popgdp2010b$gdp_per_person[popgdp2010b$gdp_per_person>1]<-1
#setkey(popgdp2010b,year_diss_me)

##future
#setkey(popgdp7,year_diss_me)
#popgdp7[,gdp_per_personssp1:=(gdp_ssp1/pop_ssp1)/max1]
#popgdp7$gdp_per_personssp1[popgdp7$gdp_per_personssp1>1]<-1
#popgdp7[,gdp_per_personssp2:=(gdp_ssp2/pop_ssp2)/max1]
#popgdp7$gdp_per_personssp2[popgdp7$gdp_per_personssp2>1]<-1
#popgdp7[,gdp_per_personssp3:=(gdp_ssp3/pop_ssp3)/max1]
#popgdp7$gdp_per_personssp3[popgdp7$gdp_per_personssp3>1]<-1


##base raster template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

##defaunation surface
#predicts3<-raster("/home/ucbtdw0/Scratch/global_IUCN/data/final-abund-bii-isl-main.tif")
#predicts2<-aggregate(predicts3,fact=5,min,na.rm=TRUE)
#predicts2b<-extend(predicts2,template)
#predicts<-resample(predicts2,template,method="bilinear")
#writeRaster(predicts2b,format="GTiff",file="/home/ucbtdw0/Scratch/global_IUCN/data/predicts_resampled_template.tif",overwrite=TRUE);gc()
predicts<-raster("/home/ucbtdw0/Scratch/New_Global_MAXENT/data/predicts_resampled_template.tif")

##create rasterized admin
pres_ad1<-fasterize(st_as_sf(ad1),template,field="diss_me")

#gini<-stack("C:/Users/Public/Documents/globalIUCN/GINI_STACK.tif")
#names(gini)<-read.csv("C:/Users/Public/Documents/globalIUCN/GINI_STACK_NAMES.csv")$x
#gini2<-velox(gini)
#gini3<-as.data.frame(gini2$extract(ad1,fun=function (x) mean(x,na.rm=TRUE)))
#gini3<-extract(gini2,ad1,mean,na.rm=TRUE)
#names(gini3)<-read.csv("C:/Users/Public/Documents/globalIUCN/GINI_STACK_NAMES.csv")$x
#gini3$diss_me<-ad1$diss_me
#gini4<-melt(gini3,id.vars="diss_me")
#gini4$value[is.na(gini4$value)]<-median(gini4$value,na.rm=TRUE)
#setDT(gini4)
#gini4[,c("year","scrap","SSP","scrap3"):= tstrsplit(variable,"_")]
#gini4[,year:=as.numeric(gsub("X","",year))]
#gini4[,c("variable", "scrap","scrap3"):=NULL]
#gini4[,year_dissme:=paste(year,diss_me,sep="_")]
#gini5<-dcast(gini4,year_dissme+year+diss_me~SSP,value.var = "value")
#setkey(gini5,year_dissme)
#save(gini5,file="C:/Users/Public/Documents/globalIUCN/GINI5.r")
#load(file="C:/Users/Public/Documents/globalIUCN/GINI5.r")

##read in two country points
lesstwo<-read.csv("/home/ucbtdw0/Scratch/New_Global_MAXENT/disease_in_two_or_fewer_countries4.csv",stringsAsFactors=FALSE)
lesstwo$name2<-paste(" ",lesstwo$name,sep="")
lesstwo$dummy<-1
coordinates(lesstwo)<-~Longitude+Latitude

##read in empress-i point data
point_data<-read.csv("/home/ucbtdw0/Scratch/New_Global_MAXENT/all_empress_i_data4.csv",stringsAsFactors = FALSE)
point_data<-point_data[!is.na(point_data$SumCases),]
#point_data$LU<-paste(" ",point_data$name_LU,sep="")

###read disease data
d1<-read.csv("/home/ucbtdw0/Scratch/New_Global_MAXENT/disease_table32c.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4

###do or not
d2<-read.csv("/home/ucbtdw0/Scratch/New_Global_MAXENT/diseases4a.csv",stringsAsFactors=FALSE)
d2a<-d2[!is.na(d2$disease),]
d2a$level[is.na(d2a$level)]<-999
d2b<-d2a[ d2a$level!=999,]
nums<-unique(d2b$combine)

###find all present day
#pres1<-list.files("V:/per_disease2/",pattern=" ALL1",full.names=TRUE)
#pres2<-list.files("V:/per_disease2/",pattern=" ALL1",full.names=FALSE)
#pres2<-gsub("  ALL1.r","",pres2)
#pres2<-gsub(" ALL1.r","",pres2)

###find all present day
#pres1b<-list.files("V:/per_disease2/",pattern=" ALL2",full.names=TRUE)
#pres2b<-list.files("V:/per_disease2/",pattern=" ALL2",full.names=FALSE)
#pres2b<-gsub("  ALL2.r","",pres2b)
#pres2b<-gsub(" ALL2.r","",pres2b)

###find all present day
pres1<-list.files("/home/ucbtdw0/Scratch/per_disease3/",pattern=" ALL4",full.names=TRUE)
pres2<-list.files("/home/ucbtdw0/Scratch/per_disease3/",pattern=" ALL4",full.names=FALSE)
pres2<-gsub("  ALL4.r","",pres2)
pres2<-gsub(" ALL4.r","",pres2)


###find all present day
iucn1<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/disease_host_ranges/",pattern="3.tif",full.names=TRUE)
iucn2<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/disease_host_ranges/",pattern="3.tif",full.names=FALSE)
iucn2<-gsub("3.tif","",iucn2)


##get just the correct ones
#pres1<-pres1[!pres2 %in% pres2b]
#pres2<-pres2[!pres2 %in% pres2b]
#pres1<-c(pres1,pres1b)
#pres2<-c(pres2,pres2b)

## same with latest
#pres1<-pres1[!pres2 %in% pres2c]
#pres2<-pres2[!pres2 %in% pres2c]
#pres1<-c(pres1,pres1c)
#pres2<-c(pres2,pres2c)


##change projecion??
##change projecion??
#load(file="/home/ucbtdw0/Scratch/New_Global_MAXENT/wrld_simpl3.r")
#wrld_simpl3<-spTransform(wrld_simpl3,CRS=projection(template))
#writeOGR(wrld_simpl3,driver="ESRI Shapefile",dsn="/home/ucbtdw0/Scratch/New_Global_MAXENT/wrld_simpl3.shp","wrld_simpl3")
wrld_simpl3<-readOGR("/home/ucbtdw0/Scratch/New_Global_MAXENT/wrld_simpl3.shp","wrld_simpl3")
#wrld_simpl3[wrld_simpl3$NAME %in% c("Russia","France"),]
wrld_simpl3@data[wrld_simpl3$NAME=="Russia","SUBREGION"]<-wrld_simpl3@data[wrld_simpl3$NAME=="Ukraine","SUBREGION"]
#wrld_simpl3[wrld_simpl3$NAME %in% c("Siberia","Mongolia"),]
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","SUBREGION"]<-30
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","REGION"]<-142
wrld_simpl3@data$ISO2<-as.character(wrld_simpl3@data$ISO2)
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","ISO2"]<-"XX"
wrld_simpl3@data$FIPS<-as.character(wrld_simpl3@data$FIPS)
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","FIPS"]<-"XX"
wrld_simpl3@data[wrld_simpl3$NAME=="Chile","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Paraguay","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Uruguay","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Argentina","SUBREGION"]<-4
###find neighbours

neigh1a<-read.csv(file="/home/ucbtdw0/Scratch/New_Global_MAXENT/country_borders.csv",stringsAsFactors = FALSE)
neigh1a<-neigh1a[!is.na(neigh1a$country_border_code),]

#wrld_simpl3[wrld_simpl3$NAME=="Western Sahara","AREA"]<-999
#wrld_simpl2[wrld_simpl2$NAME=="Luxembourg","AREA"]
#wrld_simpl3<-wrld_simpl3[wrld_simpl2$NAME!="Antarctica" & wrld_simpl2$AREA>10,]

###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl3),template)

###create new raster
t2<-ws2
t2[!is.na(t2)]<-0

#load optimum weigtings
#save(res100,file="/home/ucbtdw0/Scratch/New_Global_MAXENT/res100b.r")
#load("/home/ucbtdw0/Scratch/New_Global_MAXENT/res100d.r")
#res100$uni<-paste(res100$AUC_best.stat,res100$Group.1,sep="_")
#res100$Group.1<-as.character(res100$Group.1)

##fin best AUC
#res100b<-aggregate(res100[,2:ncol(res100)],by=list(res100$Group.1),mean)
#res100b$Group.1[res100b$Group.1=="angiostrongylus costaricensis "]<-"angiostrongylus costaricensis"
#best1$uni<-paste(best1$x,best1$Group.1,sep="_")
#res100b<-res100[res100$uni %in% best1$uni,]
#hist(res100b$AUC_best.stat)
#plot(jitter(res100b$vw),jitter(res100b$ew))
#res100b<-aggregate(res100b[,2:ncol(res100b)],by=list(res100b$Group.1),mean)
#plot(best1$AUC_cutoff,res100b$AUC_cutoff[res100b$Group.1 %in% best1$Group.1])

#plot(best1$ew,res100b$ew[res100b$Group.1 %in% best1$Group.1])
#plot(best1$vw,res100b$vw[res100b$Group.1 %in% best1$Group.1])

#i=11
##start loop
##problem with calabazo

### run some test for highest cut off

for (i in sample(nums)){
  #gc()
  #if(i<90){next}
  
  #get data
  d2c<-d2b[d2b$combine==i,]
  
  ##alter name
  if(d2c$disease[1]=="angiostrongylus costaricensis "){d2c$disease[1]<-"angiostrongylus costaricensis"}
  
  ###check in dataframes
  d2c<-d2c[paste0(d2c$disease) %in% pres2,]
  if(nrow(d2c)==0){next}
  
    ##nake name from muliple
  nam_dis<-paste(d2c$disease,collapse = "_")

   print(nam_dis)
  
    ##first writing and test
   #  if( paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/being_done//",nam_dis,"_1.csv",sep="") %in% list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/being_done/",full.names=TRUE) ){next}
  #write.csv(nam_dis,file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/being_done//",nam_dis,"_1.csv",sep=""))

   been_done<-list.files("D:\\Users\\xxxx\\Documents\\gainlossresults\\",pattern="_2.csv")
  been_done<-gsub("_X_",";",been_done)
  been_done<-read.table(text=been_done,sep=";")[,1]
   
   if(nam_dis %in% been_done){next}
  
#}
  ## make livstock
  lvsth<-subset(lvst,(1:nlayers(lvst))[names(lvst) %in% (strsplit(d2c$all_host,",")[[1]])])
  if(nlayers(lvsth)>1){lvsth<-calc(lvsth,sum,na.rm=TRUE)}
  if(nlayers(lvsth)>0){lvsth2<-resample(lvsth,template,method='ngb');lvsth[lvsth<10]<-NA;lvsth[lvsth>10]<-1}
  
  
  
  ##get IUCN ranges
  sp1<-iucn1[ iucn2 %in% d2c$disease]
  ##replace with wrld_simpl
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
  
  
  ##check multiple if so do best if one 
  #if(nrow(d2c)>1){d2c<-d2c[d2c$level==max(d2c$level,na.rm=TRUE),]}
  
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
    #res6<-aggregate(res6z,by=list(res6z$year,res6z$RCP,res6z$cell.id),mean,na.rm=TRUE)
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

  ### just need present day at the moment? Do all?
  #res6<-res6[,c(1:25,40:43)]
  res6<-cbind(res6,xyFromCell(template,res6$cell.id))
  
  ##plot to see where the data are 
  #plot(wrld_simpl)
  #pl1<-ws2
  #pl1[pl1==1]<-0
  #pl1[res6$cell.id[res6$year==2030]]<-res6$clim_present_mean[res6$year==2030]
  #plot(pl1,add=TRUE)
  #points(xyFromCell(template,unique(res6$cell.id)),col="red",pch=20)
  
  ##check point data # point_data$LU  %in% d2c$disease
  pd<-point_data[(point_data$name_LU  %in% d2c$disease) & point_data$Status!="Denied" & point_data$SumCases>0,]
  if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd);projection(pd)<-projection(raster())} else {rw1<-1000}
  
  ##bv problem
  #if(pres2[i]==" babesia venatorum"){pres2[i]<-" babesia venatorum"}
  
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
  #neigh3<-neigh2[apply(neigh2,1,function(x) any(x)==TRUE),]###not working
  
  
  regn<-wrld_simpl3[wrld_simpl3$SUBREGION %in% countr$SUBREGION,] ##psuedoabsence region
  regn2<-wrld_simpl3[wrld_simpl3$REGION %in% countr$REGION,] ##psuedoabsence region
  regn3<-wrld_simpl3[wrld_simpl3$ISO2 %in% c(neigh2$country_border_code,as.character(countr$ISO2)),] ##region of interest
  
  ###if coveerage of reporting countries is too high switch to larger region  <---- artribrary
  #if(length(unique(countr$REGION))>2){regn<-regn2}
  regnPA<-regn ##region for psuego absences
  if(nrow(countr)/nrow(regn)>=(2/3)){regnPA<-regn2}
  #if(nrow(countr)/nrow(regnPA)>=(3/4)){regnPA<-wrld_simpl2}
  
  ## but if not countries in subregion 
  ###if less than two countries
  
  l2<-lesstwo[lesstwo$name2  %in% paste0(" ",d2c$disease),]
  if(nrow(l2)>0){
    projection(l2)<-projection(ad1)
    ##overlay
    overpoints<-over(ad1,l2)
    ad2<-ad1[!is.na(overpoints$dummy), ] ### used to set present points samples
    ad2$dummy<-1:nrow(ad2)
    
    #regnPA<-regn3
    #regn3<-countr
    
    #plot(ad2,border="red",add=TRUE)
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
  
  #diss_me2<-fasterize(st_as_sf(ad1),template)
  
  ##COULD THIS BE WHAT IS CAUSING THE PATCHYNESS?
  res6[,subcountr:=extract(pres_ad2,res6$cell.id)]
  res6[,countr:=extract(presc,res6$cell.id)]
  
  res6[,sub_region:=extract(regn11,res6$cell.id)] ##largest
  res6[,regionPA:=extract(regn22,res6$cell.id)]  ## PA area 
  res6[,surrounding:=extract(regn33,res6$cell.id)]  ## - surrounding countries
  #res6[,diss_me:=extract(pres_ad1,res6$cell.id)]
  
  ###cut do to coarses area
  res6<-res6[!is.na(res6$sub_region),]
  
  ##intersection of IUCN4 and known area limit
  #zzz<-spsample(ad2,1000,type="regular")
  #zzz2<-as(raster(extent(zzz)), "SpatialPolygons")
  #zzz3<-rasterize(zzz2,template)
  #iucn5<-mask(iucn4,zzz3)
  
  ###test host range
  res6[,host_range:=extract(iucn4,res6$cell.id)]
  #res6[,host_range2:=extract(iucn5,res6$cell.id)]
  
  ##defunation
  #res6[,defaun:=extract(predicts,res6$cell.id)]

  ##distance to centre of known country
  res6[,dist_cent:=extract(cent3,res6$cell.id)]
  res6[,dist_cent:=1-(dist_cent/max(dist_cent))]
  
  ####load rasters of popgdp
  #popgdp<-extract(popdgp,res6$cell.id)
  #popgdp<-melt(popgdp)
  
  ### rid of any lines no inside	
  ### regn
  ### countries
  #res6<-res6[!is.na(res6$regionPA),]
  #res6<-res6[!is.na(res6$sub_region),]
  #res6<-res6[!is.na(res6$countr),] ##run again but for countries that have reported
  
  
  ##create unique
  res6[,year_RCP:=paste(RCP,year,sep="_")]
  #res6[,cell_by_year:=paste(cell.id,year,sep="_")]
  #setkey(res6,cell_by_year)
  #year_RCP<-unique(res6$year_RCP)
  
  ##spill over rate cases per year over population at risk 2005
  #spillover_rate<-as.numeric(dis_trans$cases_per_year)/(nrow(res6)/length(unique(res6$year_RCP)))
  spillover_rate<-as.numeric(sum(dis_trans$cases_per_year),na.rm=TRUE)/sum(countr$POP2005,na.rm=TRUE)
  
  
  #### CREATE HAZARD DIFFERENCES
  
  ##second waymarker
  #write.csv(nam_dis,file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/being_done//",nam_dis,"_2.csv",sep=""))
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
    res6$clim_present_mean_vector[res6$clim_present_mean<0.001]<-1 ##do we want this?
    res6[,clim_present_mean2:= c(1-((1-clim_present_mean) *(1-clim_present_mean_vector))),]  
    
    #res6[,density_mean:= density_mean_vector,]  
    res6[,clim_future_mean2:= c(1-((1-clim_future_mean) *(1-clim_future_mean_vector))),]  
    ##there is no vector land-use data so use host land-se
    res6[,lc_suit_present_mean2:= lc_suit_present_mean_vector * lc_suit_present_mean ,]
    #res6[,lc_suit_present_mean3:= lc_suit_present_mean2 *density_mean ,]
    res6[,lc_suit_future_mean2:= lcsuit_future_mean_vector * lcsuit_future_mean,] 
    #res6[,lc_suit_future_mean3:= lc_suit_future_mean2 *density_mean ,]
    #res6[,density_mean:= density_mean_vector,] 
    
  }else{
    res6[,clim_present_mean2:= clim_present_mean,]### change to fit in below
    res6[,clim_future_mean2:= clim_future_mean,]### change to fit in below
    res6[,lc_suit_present_mean2:= lc_suit_present_mean ,]  
    res6[,lc_suit_future_mean2:= lcsuit_future_mean,]
    #res6[,lc_suit_present_mean3:= lc_suit_present_mean *defaun *density_mean ,]  
    #res6[,lc_suit_future_mean3:= lcsuit_future_mean *defaun *density_mean ,]
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
  
  ### SHOULD IT BE LOG CLIMATE MEAN TOP END MEANS LITTLE????
  #res6[,Hazard := clim_present_mean2 * lc_suit_present_mean2 ]
  #res6[,HazardD := clim_present_mean2 * lc_suit_present_mean3 ]
  
  ###set real range as intersection of country data and host ranges
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
        
        #weighting
        #mean(res6$lc_suit_present_mean2)
        #zero_one<-sample(c(1,0),1)
  
        #if(zero_one==1){
        thr1<-do_auc(data1=res6a,test_column="Hazard",group_column="countr2",weights="dist_cent",timesP=timesP,timesA=timesA,summarise=FALSE)
        #thr2<-do_auc(data1=res6a,test_column="clim_present_mean2",group_column="countr2",weights="dist_cent",timesP=timesP,timesA=timesA,summarise=FALSE)
       #        }else{
       #thr2<-do_auc(data1=res6a,test_column="Hazard",group_column="countr2",weights="dummy",timesP=timesP,timesA=timesA,summarise=FALSE)
       #thr2<-do_auc(data1=res6a,test_column="clim_present_mean2",group_column="countr2",weights="dummy",timesP=timesP,timesA=timesA,summarise=FALSE)
         #  }
    
        rt1<-data.frame(thr1,timesA=timesA,timesP=timesP,weighting=weighting,column="centre_weighted")
        #rt1b<-data.frame(thr2,timesA=timesA,timesP=timesP,weighting=weighting,column="random")
    
       if(ww==1){rt2<-rt1} else {rt2<-rbind(rt2,rt1)}
       #print(ww)
      }
    
    ##seperate out
    #rt2<-rt2A[rt2A$column=="Hazard",]
    #rt2C<-rt2A[rt2A$column!="Hazard",]
    
    #ggplot(rt2,aes(x=AUC_cutoff,y=AUC_best.stat,col=column))+
    #  geom_point()+
    #  theme_classic()
    
    #ggplot(rt2,aes(x=as.factor(weighting),y=AUC_best.stat,col=column))+
    #  geom_boxplot()+
    #  theme_classic()
    
    ##plot if any relationship
    #plot(rt2$AUC_cutoff,rt2$AUC_best.stat,main=d2c$disease[1])
    #lm1<-lm(rt2$AUC_best.stat~rt2$AUC_cutoff)
    #abline(a=lm1$coefficients[1],b=lm1$coefficients[2])
    
    
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
  


    
    dis_trans$AUC_cutoff_real<-NA;dis_trans$AUC_real<-NA # end of if pb
  
  ###create realised present mean - abudance from habitat????
  #res6[,Hazard_future:= clim_future_mean2 * lc_suit_future_mean2 ]
  #res6[,Hazard_futureD:= clim_future_mean2 * lc_suit_future_mean3 ]

  #without future land-use change
  res6[,Hazard_future := (dis_trans$final_weighting*clim_future_mean2) +((1-dis_trans$final_weighting) * lc_suit_present_mean2) ]

  #with land-use change
  #res6[,Hazard_future := (dis_trans$final_weighting*clim_future_mean2) +((1-dis_trans$final_weighting) * lc_suit_future_mean2) ]

  
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
  
  
  ##third waymarker
  write.csv(nam_dis,file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/being_done//",nam_dis,"_3.csv",sep=""))
  
  #res6[,climT:=ifelse(clim_present_mean2>=dis_trans$AUC_cutoffclim,1,0)]
  #res6[,clim_futureT:=ifelse(clim_future_mean2>=dis_trans$AUC_cutoffclim,1,0)]
  
  ##make estimated density
  #res6[,HazardD:=HazardT * lc_suit_present_mean2 *defaun *density_mean]
  #res6[,Hazard_futureD:=Hazard_futureT * lc_suit_future_mean2 *defaun *density_mean]
  
  ##make uni
  #res6[,year_diss_me:=paste(year,diss_me,sep="_")]
  
  #res6<-res6[,c("cell.id","year","RCP","x","y","countr","sub_region","regionPA","surrounding","host_range","host_range2","defaun","diss_me","year_RCP","clim_present_mean2","clim_future_mean2","Hazard","Hazard_future","HazardT","Hazard_futureT","year_diss_me")]
  
  #fwrite(res6,file=paste("C:\\Users\\xxxx\\Documents\\datasets1\\",nam_dis,"_",sample(1:1000,1),"_all_data2.csv",sep=""))
  
  
  ##do we want to focal it a bit? To spread out? Could do it with NAs to get rid of small distance parts?
  
  ### large ggplot for each disease with difference as once columns and scenario
  ### compare to combined scenario
  unis<-unique(res6$year_RCP)

  ##what columns?
  cols1<-c("climT","clim_futureT","HazardT","Hazard_futureT")
  titls<-c(paste(cols1[1:2],round(mean(dis_trans$AUC),2),sep=" "),paste(cols1[3:4],round(mean(dis_trans$AUC),2),sep=" "))
    
      #par(mfrow=c(3,3))
    #for(yy in 1:length(unis)){
    yy=1
    zz=3
    res9<-res6[res6$year_RCP==unis[yy] ,]
    
    res9<-as.data.frame(res9)
  #titls<-paste(c("HazardT","HazardL","HazardFLSSP3","Hazard Change"),round(mean(thr2[,1]),2),sep=" ")
  #cols1<-c("clim_present_mean2","Hazard","HazardT","HazardTD","Hazard_future","Hazard_futureTD","host_range","host_range2","RiskL")
 
    ##plot to see where the data are MEANNNNN
    #plot(wrld_simpl)
    pl1<-ws2
    pl1[pl1==1]<-0
    
    ##remove cells outside intersection of countries and host ranges
    #res9b<-res9[ !is.na(res9$countr2),]
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
    #plot(wrld_simpl)
    pl1<-ws2
    pl1[pl1==1]<-0
    
    ##remove cells outside intersection of countries and host ranges
    #res9b<-res9[ !is.na(res9$countr2),]
    res9b<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding),]
    
    ##populate with HazardT
    pl1[res9b$cell.id]<-res9b[ ,names(res9b)[names(res9b)=="HazardTmin"],drop=TRUE]
    
    #if(zz==9){xxx<-(log(res9$Hazard+1))-(log(res9$Hazard_future+1));xxx[xxx==0]=NA;pl1[res9$cell.id]<-xxx}
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
    #plot(wrld_simpl)
    pl1<-ws2
    pl1[pl1==1]<-0
    
    ##remove cells outside intersection of countries and host ranges
    #res9b<-res9[ !is.na(res9$countr2),]
    res9b<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding),]
    
    ##populate with HazardT
    pl1[res9b$cell.id]<-res9b[ ,names(res9b)[names(res9b)=="HazardTmax"],drop=TRUE]
    
    #if(zz==9){xxx<-(log(res9$Hazard+1))-(log(res9$Hazard_future+1));xxx[xxx==0]=NA;pl1[res9$cell.id]<-xxx}
    #crop to general area
    pl1<-crop(pl1,extent( min(res9b$x)-5, max(res9b$x)+5,min(res9b$y)-5, max(res9b$y)+5))
    pl1[pl1==0]<-NA
    
    ee<-simpleError("no polygon overlap")
    
    ##make range
    
    > poly1max
    <simpleError in st_area.sfc(x[1]): package lwgeom required, please install it first>
    poly1max<-tryCatch(rast_to_range(x=pl1,present= host_and_present,crumb_size=1e9),error=function(ee) ee)
    
    if(class(poly1max)[1]=="simpleError"|is.na(poly1max)){poly1max<-poly1}
    rm(pl1)

    ##plot to see where the data are BLANK
    #plot(wrld_simpl)
    pl1<-ws2
    pl1[pl1==1]<-0
    
    ##remove cells outside intersection of countries and host ranges
    #res9b<-res9[ !is.na(res9$countr2),]
    #res9b<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding),]
    
    ##populate with HazardT
    #pl1[res9b$cell.id]<-res9b[ ,names(res9b)[names(res9b)=="lc_suit_present_mean2"],drop=TRUE]
    
    #if(zz==9){xxx<-(log(res9$Hazard+1))-(log(res9$Hazard_future+1));xxx[xxx==0]=NA;pl1[res9$cell.id]<-xxx}
    #crop to general area
    pl1<-crop(pl1,extent( min(res9b$x)-5, max(res9b$x)+5,min(res9b$y)-5, max(res9b$y)+5))
    pl1[pl1==0]<-NA
    
    ##fourth waymarker
    write.csv(nam_dis,file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/being_done//",nam_dis,"_4.csv",sep=""))
    
    #pdf(file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/risk_maps2/",nam_dis,"_",sample(1:1000,1),"all_plots.pdf",sep=""),width=8,height=8)
    #  plot(pl1,main=paste(nam_dis,round(dis_trans$AUC,2),round(dis_trans$AUCmax,2),round(dis_trans$AUCmin,2),sep=" "),legend=FALSE)
    #  plot(coastlines,add=TRUE,cex=2)
    #  #if(nrow(pd)>0){points(pd$Longitude,pd$Latitude,pch=20,cex=1)}
    #  plot(poly1,add=TRUE,border="grey")
    #  plot(poly1max,add=TRUE,border="green")
    #  plot(poly1min,add=TRUE,border="blue")
    #  plot(ad2,add=TRUE,cex=2,border="red")
    #  #plot(coastlines,add=TRUE,cex=2,border="green")
    
     # legend("bottomright", legend=c("Endemic","Max", "Min","Recorded"),
     #        col=c("grey","green", "blue","red"), lty=1, cex=1)
      
    #dev.off()
    
    
    ## add buffer but how far???
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
    
    
    ### can we just do this far all the dataframe at once - all future cells in or out?
    ### no need to do future polygon
    ###also define which are lost (not hazard_futureT==1) 
    ###do a xtabs in or out 2010 vs in or out 20X0 (but only for cells in correct buffer)
    ###get rid of failed models
    ###due to climate and due to land-use
    
    ##based on thomas et al (science 2011) its mean 16.9km or four cells per decade so 7 decades times four is 28 but that is for 2080 #should prevent from going too quick - what about each host at time?
    #pl2<-suppressWarnings(buffr(pl1,distance=28,units="cell",mask=FALSE))
    #pl2030<-suppressWarnings(buffr(pl1,distance=4*((2030-2010)/10),units="cell",mask=FALSE))
    #pl2050<-suppressWarnings(buffr(pl1,distance=4*((2050-2010)/10),units="cell",mask=FALSE))
    #pl2070<-suppressWarnings(buffr(pl1,distance=4*((2070-2010)/10),units="cell",mask=FALSE))
    #pl2080<-suppressWarnings(buffr(pl1,distance=4*((2080-2010)/10),units="cell",mask=FALSE))
    
    ##fifth waymarker
    write.csv(nam_dis,file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/being_done//",nam_dis,"_5.csv",sep=""))
    
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
    
   ##sixth waymarker
    write.csv(nam_dis,file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/being_done//",nam_dis,"_6.csv",sep=""))
    
   

    ##write main data
    fwrite(resF,file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/datasets/",nam_dis,"_X_",sample(1:1000,1),"_all_data4.csv",sep=""))
    
    ##write gain loss
    fwrite(dis_trans,file=paste("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",nam_dis,"_X_",sample(1:1000,1),"_2.csv",sep=""))
    
   print("DONE")
} ##end of nums loop