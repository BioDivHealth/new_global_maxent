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
library(rgeos)
library(maptools)
library(smoothr)
library(viridis)

##base raster template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

source("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\df_to_rast.R")
#source("/home/ucbtdw0/Scratch/New_Global_MAXENT/df_to_rast.R")


no_disease<-c("asama","ashriver","azagny","babesia meri","banna","bear canyon","neembucu","bloodland lake","blue river","bokeloh bat lyssavirus","borna","calabazo","campripley","cercopithecine herpes","el moro canyon","eyach","fort morgan","gumbo limbo","imjin","ippy","isla vista","jemezsprings","lee","ljungan","luxi","magboi","mobala","mopeia","muju","muleshoe","tioman virus","pergamino","pirital","reston","rickettsia helvetica","rio mamore","rio mearim","rio segundo","sangassou","tacaribe","tamiami","tanganya","tensaw","thailand","thottapalayam","topografov","una","hojo","prospect hill","khujand","histoplasma capsulatum")


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
ad1<-list.files("C:/Users/david.redding/Documents/cellids/",pattern="_2_cellids",full.names=TRUE)
ad2<-list.files("C:/Users/david.redding/Documents/cellids/",pattern="_2_cellids",full.names=FALSE)
ad3<-gsub("_2_cellids.csv","",ad2)


#gl1<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",pattern="_1.csv",full.names=TRUE)
#gl2<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",pattern="_1.csv",full.names=FALSE)
gl1<-list.files("C:/Users/david.redding/Dropbox/New_Global_MAXENT/gainlossresults/",pattern="_2.csv",full.names=TRUE)
gl2<-list.files("C:/Users/david.redding/Dropbox/New_Global_MAXENT/gainlossresults/",pattern="_2.csv",full.names=FALSE)
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

#d1<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT/disease_table32c.csv",stringsAsFactors=FALSE)
#d1$name2<-paste(" ",d1$name,sep="")
#d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
#d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4
#glF3<-merge(glF2,d1,by.x="name",by.y="name",all.x=TRUE,all.y=FALSE)

###make spill over rate
glF3<-glF2

###create meta type
glF3$Type2<-glF3$Type
glF3$Type2[glF3$Type %in% c("HOST->VECTOR->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN")]<-"HOST_VECTOR_HUMAN"
glF3$Type2[glF3$Type %in% c("HOST->HUMAN","HOST->HUMAN->HUMAN")]<-"HOST_HUMAN"
glF3$Type2[glF3$Type %in% c("HUMAN->VECTOR->HUMAN","HUMAN->VECTOR->HUMAN->HUMAN")]<-"HUMAN_VECTOR_HUMAN"

###create spillover categories
glF3$spillover_rate2<-cut(glF3$cases_per_year,breaks=c(0,4.99,99,9999,1e9),labels=FALSE)
#glF3$spillover_rate2<-cut(log(glF3$cases_per_year),breaks=getJenksBreaks(log(glF3$cases_per_year), 5),labels=FALSE)
table(glF3$spillover_rate2)
glF3$dummy=1
#meancases<-aggregate(glF3$cases_per_year,by=list(glF3$spillover_rate2),function (x) mean(x,trim=0.1,na.rm=TRUE))
#round(meancases,0)
mediancases<-aggregate(glF3$cases_per_year,by=list(glF3$spillover_rate2),median)
round(mediancases,0)
names(mediancases)<-c("spillover_rate2","cases_per_year2")

glF3<-merge(glF3,mediancases,by="spillover_rate2")
#rowSums(xtabs(data=glF3,dummy~spillover_rate2+Type2))
##popdgp
#popgdp<-fread("C:\\Users\\Public\\Documents\\future_gdppop_maskedb.csv")
#setkey(popgdp,cell_by_year)

##popdgp
#pop<-fread("/home/ucbtdw0/Scratch/New_Global_MAXENT/data/future_pop_maskedb.csv")
pop<-fread("C:/Users/david.redding/Dropbox/New_Global_MAXENT/future_pop_maskedb.csv")
names(pop)[9]<-"humans_present"
setkey(pop,cell_by_year)

##gini model (called lm1)
#load(file="C:\\Users\\Public\\Documents\\fitted_gini_10th.r")

err = simpleError("Error in read1")
err2 = simpleError("Error in read2")

###########################################################################
###########################################################################
###########################################################################
bef1<-ls()
j=1
#res5<-NULL
for (j in j:length(ad1)){

    ##load up each disease
  resF2<-tryCatch(fread(ad1[j]), error=function(err) err)
  if(class(resF2)[1]=="simpleError"){ next }
  
  #if(length(glF3$name[glF3$name==resF2$name[1]])==0){next}
  
  ###disease details
  disdet<-glF3[glF3$name==resF2$name[1],]
  ccc<-strsplit(paste(disdet$countries[1],collapse = ","),",")[[1]] ##deal with mulitple ROWS ##unique
  #if(is.na(ccc[1])){print("no recognised countries");print(pres2[i]);next}
  countr2<-wrld_simpl3[wrld_simpl3$ISO2 %in% ccc ,]
  countr2$dummy=1
  countr<-unionSpatialPolygons(countr2,ID=countr2$dummy)
  countr<-drop_crumbs(countr,threshold=1e9)
  #plot(countr)
  nnn<-resF2$name[1]
  
  ## remove if AUC is less than 0.6
  if(disdet$AUCmax<0.6){next}
  
  ##no non diseases
  #if(nnn %in% no_disease){next}
   
     ##make cell by year
     resF2[,c("RCP","year"):=tstrsplit(year_RCP,"_")]
      resF2[,cell_by_year:=paste(cell.id,year,sep="_")]
     setkey(resF2,cell_by_year)
  
     ###merge
     resF<-pop[resF2]
     #setkey(resF,cell_by_year)
     #resF<-popgdp[resF]
     #setkey(resF,cell_by_year)
     
     rm(resF2);gc()
     #merge gini
     #resF<-gini5[resF]
     
     #remove duplicate columns
     #resF[,c("i.cell.id","i.cell.id.2","i.time","i.year","RCP","i.cell.id.1"):=NULL]
     #resF<-resF[,!duplicated(as.list(resF)),with=FALSE]
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
     #resF[,c("cell_by_year"):=NULL]
     ##create dummy
     resF[,dummy:=1]
     
     #k=1
     for (k in 1:4){
       
       #skip=0
       #if(k==1) {typex="ALL"}
       #if(k==2) {typex="HOST_VECTOR_HUMAN";if(!disdet$Type %in% c("HOST->VECTOR->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN")){skip=1}}
       #if(k==3) {typex="HOST_HUMAN";if(!disdet$Type %in% c("HOST->HUMAN","HOST->HUMAN->HUMAN")){skip=1}}
       #if(k==4) {typex="HUMAN_VECTOR_HUMAN";if(!disdet$Type %in% c("HUMAN->VECTOR->HUMAN","HUMAN->VECTOR->HUMAN->HUMAN")){skip=1}}
       
       #if(skip==1){next}
       typex=NULL
       skip=0
       if(k==1) {typex="Emerging";if(disdet$spillover_rate2 %in% c(2,3,4)){skip=1}}
       if(k==2) {typex="Established";if(disdet$spillover_rate2 %in% c(1,3,4)){skip=1}}
       if(k==3) {typex="Endemic";if(disdet$spillover_rate2 %in% c(1,2,4)){skip=1}}
       if(k==4) {typex="Hyperendemic";if(disdet$spillover_rate2 %in% c(1,2,3)){skip=1}}
       
       if(skip==1){next}
       
     
     ##split
     pres<-resF[resF$time1=="present",]
     pres[,c("cell_by_year","time1","time"):=NULL]
     pres<-pres[!duplicated(cell.id),]
     pres[,c("ssp1","ssp2","ssp3","ssp4","ssp5","year_RCP","year","RCP"):=NULL]
     
     ##cases per year divide by population at risk
     spill_over_rate<-disdet$cases_per_year2/sum(pres$humans_present,na.rm=TRUE)
     
     ###future 
     fut<-resF[resF$time1!="present",]
     fut[,c("cell_by_year","time1","time"):=NULL]
     
     # raster of present cases
     pres4<-df_to_rast(template=ws2,data=pres,column="dummy")
     
     #assign(combs$uni[qq],pres_spill)
     if(class(tryCatch(get(paste(typex,"endemic",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(typex,"endemic",sep="_"),pres4)} else {assign(paste(typex,"endemic",sep="_"),get(paste(typex,"endemic",sep="_"))+pres4)}
     
     pres4b<-df_to_rast(ws2,pres,"humans_present")
     
     #assign(combs$uni[qq],pres_spill)
     if(class(tryCatch(get(paste(typex,"humans_interactions",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(typex,"humans_interactions",sep="_"),pres4b)} else {assign(paste(typex,"humans_interactions",sep="_"),get(paste(typex,"humans_interactions",sep="_"))+pres4b)}
     
     pres_spill<-log((pres4b*spill_over_rate)+1)
     
     #assign(combs$uni[qq],pres_spill)
     if(class(tryCatch(get(paste(typex,"human_cases",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(typex,"human_cases",sep="_"),pres_spill)} else {assign(paste(typex,"human_cases",sep="_"),get(paste(typex,"human_cases",sep="_"))+pres_spill)}
     
     pres_spill2<-pres_spill*((disdet$CFR.high+disdet$CFR.low)/2)
     
     #assign(combs$uni[qq],pres_spill)
     if(class(tryCatch(get(paste(typex,"deaths",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(typex,"deaths",sep="_"),pres_spill2)} else {assign(paste(typex,"deaths",sep="_"),get(paste(typex,"deaths",sep="_"))+pres_spill2)}
     
     #pdf(file=paste("C:\\Users\\david.redding\\Documents\\risk_maps_pres\\",nnn,"_",sample(1:1000,1),"_endemic_area.pdf",sep=""),width=12,height=6)
     #plot(pres4,main=paste(nnn,round(disdet$AUCmax,2),sep=" "),legend=FALSE)
     #plot(countr,add=TRUE,border="red",lwd=1)
     #dev.off()
     
        
     ##times spill over rate?
     #pres_spill<-pres4*spill_over_rate*((disdet$CFR.high+disdet$CFR.low)/2)
     #print(j)}
     
     ##all future combs
     combs2<-expand.grid(names(fut)[2:6],unique(fut$year_RCP))
     combs<-combs2[combs2$Var1=="ssp3" & combs2$Var2 %in% c("2.6_2080","4.5_2080","6_2080","8.5_2080"),] ##reduce to 1
     if(nrow(combs)==0){print("skipped no 2080");next}
      # combs<-combs2[combs2$Var1=="ssp3" & combs2$Var2 %in% c("2.6_2070","4.5_2070","6_2070","8.5_2070"),]     
       # }
     combs$uni<-paste(combs$Var1,combs$Var2,typex,sep="X")
     #combs$uni2<-paste(combs$Var1,combs$Var2,typex,"endemic",sep="X")
     
     #if(j==1){combs$done<-0} else {combs$done<-0}
     
     ### loop through
     for (qq in 1:nrow(combs)){
     #  for (qq in 1){
         
       fut4<-df_to_rast(ws2,fut[fut$year_RCP==combs[qq,2],],"dummy")
       
       
       #assign(combs$uni[qq],fut_spill)
       if(class(tryCatch(get(paste(combs$uni[qq],"future_endemic",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_endemic",sep="_"),fut4)} else {assign(paste(combs$uni[qq],"future_endemic",sep="_"),get(paste(combs$uni[qq],"future_endemic",sep="_"))+fut4)}
       
              
       change1<-fut4-pres4
        change2<-change1
        change2[change2>0]<-0
        change1[change1<0]<-0
       
       
        if(class(tryCatch(get(paste(combs$uni[qq],"future_endemic_increase",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_endemic_increase",sep="_"),change1)} else {assign(paste(combs$uni[qq],"future_endemic_increase",sep="_"),get(paste(combs$uni[qq],"future_endemic_increase",sep="_"))+change1)}
        
        if(class(tryCatch(get(paste(combs$uni[qq],"future_endemic_decrease",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_endemic_decrease",sep="_"),change2)} else {assign(paste(combs$uni[qq],"future_endemic_decrease",sep="_"),get(paste(combs$uni[qq],"future_endemic_decrease",sep="_"))+change2)}
        rm(change1,change2) 
         
       fut4b<-df_to_rast(ws2,fut[fut$year_RCP==combs[qq,2],],combs[qq,1])
       
       #assign(combs$uni[qq],fut_spill)
       if(class(tryCatch(get(paste(combs$uni[qq],"future_humans_interactions",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_humans_interactions",sep="_"),fut4b)} else {assign(paste(combs$uni[qq],"future_humans_interactions",sep="_"),get(paste(combs$uni[qq],"future_humans_interactions",sep="_"))+fut4b)}
       
       change1<-fut4b-pres4b
       change2<-change1
       change2[change2>0]<-0
       change1[change1<0]<-0
       
       
       if(class(tryCatch(get(paste(combs$uni[qq],"future_humans_interactions_increase",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_humans_interactions_increase",sep="_"),change1)} else {assign(paste(combs$uni[qq],"future_humans_interactions_increase",sep="_"),get(paste(combs$uni[qq],"future_humans_interactions_increase",sep="_"))+change1)}
       
       if(class(tryCatch(get(paste(combs$uni[qq],"future_humans_interactions_decrease",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_humans_interactions_decrease",sep="_"),change2)} else {assign(paste(combs$uni[qq],"future_humans_interactions_decrease",sep="_"),get(paste(combs$uni[qq],"future_humans_interactions_decrease",sep="_"))+change2)}
       rm(change1,change2) 
       
       fut_spill<-log((fut4b*spill_over_rate)+1)
       
       #assign(combs$uni[qq],fut_spill)
       if(class(tryCatch(get(paste(combs$uni[qq],"future_human_cases",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_human_cases",sep="_"),fut_spill)} else {assign(paste(combs$uni[qq],"future_human_cases",sep="_"),get(paste(combs$uni[qq],"future_human_cases",sep="_"))+fut_spill)}

       
       change1<-fut_spill-pres_spill
       change2<-change1
       change2[change2>0]<-0
       change1[change1<0]<-0
       
       
       if(class(tryCatch(get(paste(combs$uni[qq],"future_human_cases_increase",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_human_cases_increase",sep="_"),change1)} else {assign(paste(combs$uni[qq],"future_human_cases_increase",sep="_"),get(paste(combs$uni[qq],"future_human_cases_increase",sep="_"))+change1)}
       
       if(class(tryCatch(get(paste(combs$uni[qq],"future_human_cases_decrease",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_human_cases_decrease",sep="_"),change2)} else {assign(paste(combs$uni[qq],"future_human_cases_decrease",sep="_"),get(paste(combs$uni[qq],"future_human_cases_decrease",sep="_"))+change2)}
       rm(change1,change2) 
              
       fut_spill2<-fut_spill*((disdet$CFR.high+disdet$CFR.low)/2)
       
       #assign(combs$uni[qq],fut_spill)
       if(class(tryCatch(get(paste(combs$uni[qq],"future_deaths",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_deaths",sep="_"),fut_spill2)} else {assign(paste(combs$uni[qq],"future_deaths",sep="_"),get(paste(combs$uni[qq],"future_deaths",sep="_"))+fut_spill2)}
       
       change1<-fut_spill2-pres_spill2
       change2<-change1
       change2[change2>0]<-0
       change1[change1<0]<-0
       
       
       if(class(tryCatch(get(paste(combs$uni[qq],"future_deaths_increase",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_deaths_increase",sep="_"),change1)} else {assign(paste(combs$uni[qq],"future_deaths_increase",sep="_"),get(paste(combs$uni[qq],"future_deaths_increase",sep="_"))+change1)}
       
       if(class(tryCatch(get(paste(combs$uni[qq],"future_deaths_decrease",sep="_")), error=function(err2) err2))[1]=="simpleError"){assign(paste(combs$uni[qq],"future_deaths_decrease",sep="_"),change2)} else {assign(paste(combs$uni[qq],"future_deaths_decrease",sep="_"),get(paste(combs$uni[qq],"future_deaths_decrease",sep="_"))+change2)}
       rm(change1,change2) 
      
      if(qq %in% seq(from=10,to=80,by=10)) {gc()}
      
       #print(qq)
     }
     
     rm(pres,fut)
  
  } #end k loop
  #writeRaster(david.redding,"/home/ucbtdw0/Scratch/New_Global_MAXENT/data/future_pop_maskedb.csv")
     rm(resF);gc()
     
  #dim((fut[fut$year_RCP==combs[qq,2],"cell.id"]))

  print(j)

}# end of j

plot(Emerging_endemic,main="Emerging",col=viridis(20,option="A"))
plot(Established_endemic,main="Established",col=viridis(20,option="A"))
plot(Endemic_endemic,main="Endemic",col=viridis(20,option="A"))
plot(Hyperendemic_endemic,main="Hyperendemic",col=viridis(20,option="A"))

plot(ssp3X8.5_2080XEmerging_future_endemic,main="Emerging",col=viridis(20,option="A"))
plot(ssp3X8.5_2080XEmerging_future_endemic_increase/ssp3X8.5_2080XEmerging_future_endemic,main="Emerging",col=viridis(20,option="A"))
plot((-1*ssp3X8.5_2080XEmerging_future_endemic_decrease)/ssp3X8.5_2080XEmerging_future_endemic,main="Emerging",col=viridis(20,option="A"))

plot(ssp3X8.5_2080XEstablished_future_endemic,main="Established_future",col=viridis(20,option="A"))
plot(ssp3X8.5_2080XEstablished_future_endemic_increase/ssp3X8.5_2080XEstablished_future_endemic,main="Established_future",col=rev(viridis(20,option="A")))
plot((-1*ssp3X8.5_2080XEstablished_future_endemic_decrease)/ssp3X8.5_2080XEstablished_future_endemic,main="Established_future",col=rev(viridis(20,option="A")))

plot(ssp3X8.5_2080XEndemic_future_endemic,main="Endemic_future",col=viridis(20,option="A"))
plot(ssp3X8.5_2080XEndemic_future_endemic_increase/ssp3X8.5_2080XEndemic_future_endemic,main="Endemic_future",col=rev(viridis(20,option="A")))
plot((-1*ssp3X8.5_2080XEndemic_future_endemic_decrease)/ssp3X8.5_2080XEndemic_future_endemic,main="Endemic_future",col=rev(viridis(20,option="A")))


plot(ssp3X8.5_2080XHyperendemic_future_endemic)
plot(ssp3X8.5_2080XHyperendemic_future_endemic_increase/ssp3X8.5_2080XHyperendemic_future_endemic)
plot(-1*(ssp3X8.5_2080XHyperendemic_future_endemic_decrease/ssp3X8.5_2080XHyperendemic_future_endemic))



names1<-ls()[c(7:10,30:41,55:246)]

for (w in 1:length(names1)){
  
  rast1<-get(names1[w])
  names(rast1)<-names1[w]
  
  writeRaster(rast1,format="GTiff",file=paste("C:\\Users\\david.redding\\Documents\\more_results\\",gsub(">","",names1[w]),"_",sample(1:1e8,1),"_no_nondisease.tif",sep=""))
  
  print(w)
  
}
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")

test1<-aggregate(HOST_HUMAN_endemic,100,method='bilinear')

test2<-aggregate(ssp3X8.5_2080XHOST_HUMAN_future_endemic_increase,100,method='bilinear')
test2[test2==0]<-NA

plot(test1,col=rev(viridis(20,option="magma")))

plot(test2,col=rev(viridis(20)),add=TRUE,alpha=0.5)
plot(coastlines,add=TRUE)

names1[1:10]
