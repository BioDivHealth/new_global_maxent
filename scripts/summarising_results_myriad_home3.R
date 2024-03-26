library(data.table)
library(raster)
#library(velox)
#library(ggplot2)
#library(ggpubr)
#library(cowplot)
#library(ggthemes)
library(sp)
library(rgdal)
library(sf)
library(fasterize)

##base raster template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

source("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\df_to_rast.R")
#source("/home/ucbtdw0/Scratch/New_Global_MAXENT/df_to_rast.R")

##change projecion??
##change projecion??
#load(file="/home/ucbtdw0/Scratch/New_Global_MAXENT/wrld_simpl3.r")
#wrld_simpl3<-spTransform(wrld_simpl3,CRS=projection(template))
#writeOGR(wrld_simpl3,driver="ESRI Shapefile",dsn="/home/ucbtdw0/Scratch/New_Global_MAXENT/wrld_simpl3.shp","wrld_simpl3")
#wrld_simpl3<-readOGR("New_Global_MAXENT/wrld_simpl3.shp","wrld_simpl3")
#wrld_simpl3<-readOGR("/home/ucbtdw0/Scratch/New_Global_MAXENT/wrld_simpl3.shp","wrld_simpl3")
wrld_simpl3<-readOGR("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\wrld_simpl3.shp","wrld_simpl3")
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

#neigh1a<-read.csv(file="/home/ucbtdw0/Scratch/New_Global_MAXENT/country_borders.csv",stringsAsFactors = FALSE)
#neigh1a<-neigh1a[!is.na(neigh1a$country_border_code),]

#wrld_simpl3[wrld_simpl3$NAME=="Western Sahara","AREA"]<-999
#wrld_simpl2[wrld_simpl2$NAME=="Luxembourg","AREA"]
#wrld_simpl3<-wrld_simpl3[wrld_simpl2$NAME!="Antarctica" & wrld_simpl2$AREA>10,]

###rasterize world simple
###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl3),template)

##make blank
ws3<-ws2
ws2[ws2==1]<-0

### plot and summarising resutls
### areas that have increased due to climate versus those that increased due to
### land use. Size of dots - burden estimate?
### panels for each year - facet across year and RCP?
### if over/under threshold due climate versus 

#gini5<-fread(file="C:\\Users\\Public\\Documents\\GINI_by_cell1.csv")
#setkey(gini5,cell_by_year)
#gini5<-dcast(gini4,year_dissme+year+diss_me~SSP,value.var = "value")
#setkey(gini5,year_dissme)


#ad1<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/cellids/",pattern="cellids",full.names=TRUE)
#ad2<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/cellids/",pattern="cellids",full.names=FALSE)
ad1<-list.files("C:/Users/david.redding/Documents/cellids/",pattern="cellids",full.names=TRUE)
ad2<-list.files("C:/Users/david.redding/Documents/cellids/",pattern="cellids",full.names=FALSE)
ad3<-gsub("_1_cellids.csv","",ad2)


#gl1<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",pattern="_1.csv",full.names=TRUE)
#gl2<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",pattern="_1.csv",full.names=FALSE)
gl1<-list.files("C:/Users/david.redding/Documents/gainlossresults/",pattern="_1.csv",full.names=TRUE)
gl2<-list.files("C:/Users/david.redding/Documents/gainlossresults/",pattern="_1.csv",full.names=FALSE)
gl2<-gsub("_1.csv","",gl2)
gl2<-gsub("_X_",";",gl2)
gl3<-read.table(text=gl2,sep=";",stringsAsFactors = F)$V1

###get gain loss
for (i in 1:length(gl1)){
  
  ##gainloss
  glF<-tryCatch(fread(gl1[i]), error=function(err2) err2)
  if(class(glF)[1]=="simpleError"){ next }
  
  if(i==1) {glF2<-glF} else {glF2<-rbind(glF2,glF)}
  
}

##
#hist(glF2$AUC,breaks=25)
#length(glF2$AUC[glF2$AUC>0.6])

###read disease data
#d1<-read.csv("/home/ucbtdw0/Scratch/New_Global_MAXENT/disease_table32c.csv",stringsAsFactors=FALSE)

d1<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT/disease_table32c.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4


glF3<-merge(glF2,d1,by.x="name",by.y="name",all.x=TRUE,all.y=FALSE)

##popdgp
#popgdp<-fread("C:\\Users\\Public\\Documents\\future_gdppop_maskedb.csv")
#setkey(popgdp,cell_by_year)

##popdgp
#pop<-fread("/home/ucbtdw0/Scratch/New_Global_MAXENT/data/future_pop_maskedb.csv")
pop<-fread("C:\\Users\\Public\\Documents\\future_pop_maskedb.csv")
names(pop)[9]<-"humans_present"
setkey(pop,cell_by_year)

##gini model (called lm1)
#load(file="C:\\Users\\Public\\Documents\\fitted_gini_10th.r")

err = simpleError("Error in read1")
err2 = simpleError("Error in read2")

#res5<-NULL
for (j in 68:length(ad1)){
 
  ##load up each disease
  resF2<-tryCatch(fread(ad1[j]), error=function(err) err)
  if(class(resF2)[1]=="simpleError"){ next }
  
  if(length(glF3$name[glF3$name==resF2$name[1]])==0){next}
  
  ###disease details
  disdet<-glF3[glF3$name==resF2$name[1],]
  
  ## remove if AUC is less than 0.6
  if(disdet$AUCmax<0.6){next}
  
  ## do for some subsections
  #"HOST->VECTOR->HUMAN"                       "HOST->HUMAN"                               "HOST->HUMAN->HUMAN"                       
  #"HUMAN->VECTOR->HUMAN"                      "HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN" "HOST->VECTOR->HUMAN->VECTOR->HUMAN"       
  #"HUMAN->VECTOR->HUMAN->HUMAN" 
  typex="ALL"
  #if(!disdet$Type.x %in% c("HOST->VECTOR->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN")){typex="HOST->VECTOR->HUMAN";next}
  #if(!disdet$Type.x %in% c("HOST->HUMAN","HOST->HUMAN->HUMAN")){typex="HOST->HUMAN";next}
  #if(!disdet$Type.x %in% c("HUMAN->VECTOR->HUMAN","HUMAN->VECTOR->HUMAN->HUMAN")){typex="HUMAN->VECTOR->HUMAN";next}
  
  
  ##make cell by year
  resF2[,c("RCP","year"):=tstrsplit(year_RCP,"_")]
  resF2[,cell_by_year:=paste(cell.id,year,sep="_")]
  setkey(resF2,cell_by_year)
  
  ###merge
  resF<-pop[resF2]
  setkey(resF,cell_by_year)
  #resF<-popgdp[resF]
  #setkey(resF,cell_by_year)
  
  rm(resF2);gc()
  #merge gini
  #resF<-gini5[resF]
  
  #remove duplicate columns
  #resF[,c("i.cell.id","i.cell.id.2","i.time","i.year","RCP","i.cell.id.1"):=NULL]
  resF<-resF[,!duplicated(as.list(resF)),with=FALSE]
  #resF[,c("i.year","V1"):=NULL]
  
  ##xy coorindates
  resF[,c("lon","lat"):=as.data.frame(xyFromCell(template,resF$cell.id))]
  

  ## calculate Exposure
  #"gdp-ssp1" "gdp-ssp2" "gdp-ssp3"
  #"pop-ssp1" "pop-ssp2"  "pop-ssp3"
  #"gdp2010" 
  #"humans2010" "ssp1"  "ssp2" "ssp3" "ssp4" "ssp5"
  
  ##make proportion of people in poverty GINI
  ## according to gini assuming a lognormal this is the proporation of people earning less than the per capita GDP
  #resF$present_value[is.na(resF$present_value)]<-0
  #resF$future_value[is.na(resF$future_value)]<-0
  #resF[,prop_in_pov:=stats::predict.lm(newdata=data.frame(gini=PRESENT/100),lm1)]
  #resF[,prop_in_pov_futureSSP1:=stats::predict.lm(newdata=data.frame(gini=SSP1/100),lm1)]
  #resF[,prop_in_pov_futureSSP2:=stats::predict.lm(newdata=data.frame(gini=SSP2/100),lm1)]
  #resF[,prop_in_pov_futureSSP3:=stats::predict.lm(newdata=data.frame(gini=SSP3/100),lm1)]
  #resF[,prop_in_pov_futureSSP4:=stats::predict.lm(newdata=data.frame(gini=SSP4/100),lm1)]
  #resF[,prop_in_pov_futureSSP5:=stats::predict.lm(newdata=data.frame(gini=SSP5/100),lm1)]
  
  ## calculate Exposure
  #"humans2010" "ssp1"  "ssp2" "ssp3" "ssp4" "ssp5"
  #resF[,people_in_poverty:=prop_in_pov*humans_present]
  #resF[,people_in_povertySSP1:=(prop_in_pov_futureSSP1*ssp1)-people_in_poverty]
  #resF[,people_in_povertySSP2:=(prop_in_pov_futureSSP2*ssp2)-people_in_poverty]
  #resF[,people_in_povertySSP3:=(prop_in_pov_futureSSP3*ssp3)-people_in_poverty]
  #resF[,people_in_povertySSP4:=(prop_in_pov_futureSSP4*ssp4)-people_in_poverty]
  #resF[,people_in_povertySSP5:=(prop_in_pov_futureSSP5*ssp5)-people_in_poverty]
  
  #resF[,peopleSSP1:=ssp1-humans_present]
  #resF[,peopleSSP2:=ssp2-humans_present]
  #resF[,peopleSSP3:=ssp3-humans_present]
  #resF[,peopleSSP4:=ssp4-humans_present]
  #resF[,peopleSSP5:=ssp5-humans_present]
  
  ###
  resF[,c("cell_by_year"):=NULL]
  
  ##split
  pres<-resF[resF$time1=="present",]
  pres[,c("cell_by_year","time1","time"):=NULL]
  pres<-pres[!duplicated(cell.id),]
  pres[,c("ssp1","ssp2","ssp3","ssp4","ssp5","year_RCP","year","RCP"):=NULL]
  
  ##cases per year divide by population at risk
  spill_over_rate<-disdet$cases_per_year.x/sum(pres$humans_present,na.rm=TRUE)
  
  ###future 
  fut<-resF[resF$time1!="present",]
  fut[,c("cell_by_year","time1","time"):=NULL]
  
  rm(resF);gc()
  
  ##which are not in one but in other
  #in present not any future
  pres2<-pres[!cell.id %in% fut$cell.id,]
  
  ##in future not in present
  fut2<-fut[!cell.id %in% pres$cell.id,]
  
  rm(pres,fut);gc()
  
  # raster of present cases
  pres4<-df_to_rast(ws2,pres2,"humans_present")
  
  ##times spill over rate?
  pres_spill<-pres4*spill_over_rate
  #assign(combs$uni[qq],pres_spill)
  if(j==1){assign(paste(typex,"present",sep="_"),pres_spill)} else {assign(paste(typex,"present",sep="_"),get(paste(typex,"present",sep="_"))+pres_spill)}
  
  #print(j)}
  
  
  ##all future combs
  combs<-expand.grid(names(fut2)[2:6],unique(fut2$year_RCP))
  combs$uni<-paste(combs$Var1,combs$Var2,typex,sep="X")
  if(j==1){combs$done<-0} else {combs$done<-0}
  
  ### loop through
  #for (qq in 1:nrow(combs)){
  
  for (qq in 18:22){
      
    
    #create blank
    if(combs$done[qq]==0){assign(combs$uni[qq],ws3);combs$done[qq]<-1}
    
    #create rast
    fut3<-df_to_rast(ws2,fut2[fut2$year_RCP==combs[qq,2],],combs[qq,1])
    
    assign(combs$uni[qq],get(combs$uni[qq])+(fut3*spill_over_rate))
    #assign(paste(typex,"present",sep="_"),get(paste(typex,"present",sep="_"))+pres_spill)
    
    #print(qq)
  }
  
  
  #writeRaster(xxxx,"/home/ucbtdw0/Scratch/New_Global_MAXENT/data/future_pop_maskedb.csv")
  
  
  print(j)

}# end of j

