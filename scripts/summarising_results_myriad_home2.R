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



##change projecion??
##change projecion??
#load(file="/home/ucbtdw0/Scratch/New_Global_MAXENT/wrld_simpl3.r")
#wrld_simpl3<-spTransform(wrld_simpl3,CRS=projection(template))
#writeOGR(wrld_simpl3,driver="ESRI Shapefile",dsn="/home/ucbtdw0/Scratch/New_Global_MAXENT/wrld_simpl3.shp","wrld_simpl3")
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


### plot and summarising resutls
### areas that have increased due to climate versus those that increased due to
### land use. Size of dots - burden estimate?
### panels for each year - facet across year and RCP?
### if over/under threshold due climate versus 

gini5<-fread(file="C:\\Users\\Public\\Documents\\GINI_by_cell1.csv")
setkey(gini5,cell_by_year)
#gini5<-dcast(gini4,year_dissme+year+diss_me~SSP,value.var = "value")
#setkey(gini5,year_dissme)


ad1<-list.files("C:/Users/david.redding/Documents/cellids/",pattern="cellids",full.names=TRUE)
ad2<-list.files("C:/Users/david.redding/Documents/cellids/",pattern="cellids",full.names=FALSE)
ad3<-gsub("_1_cellids.csv","",ad2)


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
hist(glF2$AUC)

###read disease data
d1<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT/disease_table32c.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4


glF3<-merge(glF2,d1,by.x="name",by.y="name",all.x=TRUE,all.y=FALSE)

##popdgp
#popgdp<-fread("C:\\Users\\Public\\Documents\\future_gdppop_maskedb.csv")
#setkey(popgdp,cell_by_year)

##popdgp
pop<-fread("C:\\Users\\Public\\Documents\\future_pop_maskedb.csv")
names(pop)[9]<-"humans_present"
setkey(pop,cell_by_year)

##gini model (called lm1)
load(file="C:\\Users\\Public\\Documents\\fitted_gini_10th.r")

err = simpleError("Error in read1")
err2 = simpleError("Error in read2")

res5<-NULL
for (j in sample(1:length(ad1))){
 
  ##load up each disease
  resF2<-tryCatch(fread(ad1[j]), error=function(err) err)
  if(class(resF2)[1]=="simpleError"){ next }
  
  if(length(glF3$name[glF3$name==resF2$name[1]])==0){next}
  
  ###disease details
  disdet<-glF3[glF3$name==resF2$name[1],]
  
  ##make cell by year
  resF2[,c("RCP","year"):=tstrsplit(year_RCP,"_")]
  resF2[,cell_by_year:=paste(cell.id,year,sep="_")]
  setkey(resF2,cell_by_year)
  
  ###merge
  resF<-pop[resF2]
  setkey(resF,cell_by_year)
  #resF<-popgdp[resF]
  #setkey(resF,cell_by_year)
  
  #rm(resF2):gc()
  #merge gini
  resF<-gini5[resF]
  
  #remove duplicate columns
  #resF[,c("i.cell.id","i.cell.id.2","i.time","i.year","RCP","i.cell.id.1"):=NULL]
  resF<-resF[,!duplicated(as.list(resF)),with=FALSE]
  resF[,c("i.year","V1"):=NULL]
  
  ##make proportion of people in poverty GINI
  ## according to gini assuming a lognormal this is the proporation of people earning less than the per capita GDP
  resF$present_value[is.na(resF$present_value)]<-0
  resF$future_value[is.na(resF$future_value)]<-0
  resF[,prop_in_pov:=stats::predict.lm(newdata=data.frame(gini=present_value/100),lm1)]
  resF[,prop_in_pov_future:=stats::predict.lm(newdata=data.frame(gini=future_value/100),lm1)]
  
  ## calculate Exposure
  #"gdp-ssp1" "gdp-ssp2" "gdp-ssp3"
  #"pop-ssp1" "pop-ssp2"  "pop-ssp3"
  #"gdp2010" 
  #"humans2010" "ssp1"  "ssp2" "ssp3" "ssp4" "ssp5"
  
  resF[,Risk1:=1]
  
  
  ##cells in present day ## add in AUC score
  mean_pres_size<-nrow(resF[HazardT==1 & year_RCP==unique(resF$year_RCP)[1],])
  min_pres_size<-nrow(resF[HazardTmin==1 & year_RCP==unique(resF$year_RCP)[1],])
  max_pres_size<-nrow(resF[HazardTmax==1 & year_RCP==unique(resF$year_RCP)[1],])
  pres_size<-data.frame(term=c("mean","min","max"),pres_size=c(mean_pres_size,min_pres_size,max_pres_size))
      
  ##work out cells with change
  resF[,Clim_round := round(clim_present_mean2,1) ]
  resF[,LU_round := round(lc_suit_present_mean2,1) ]
  resF[,Clim_future_round := round(clim_future_mean2,1) ]
  resF[,LU_future_round := round(lc_suit_future_mean2,1) ]
  
  ###record if increase or decrease in the future
  resF[,clim_up:=ifelse(Clim_future_round>Clim_round,1,0)]
  resF[,clim_down:=ifelse(Clim_future_round<Clim_round,1,0)]
  resF[,lu_up:=ifelse(LU_future_round>LU_round,1,0)]
  resF[,lu_down:=ifelse(LU_future_round<LU_round,1,0)]
  
  
  
  
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
  #plot(poly2080)
  #resF3<-resF2[resF2$year==2080,]
  #points(resF3$x,resF3$y)
  #resF2[,dummy:=1]
  
  ## what about clim up clim down in just areas that are experiencing an increase/decrease - why the change rather than how is changing over whole area
  
  ## also present should just be in present poly - rest should be NA
  
  ## bring in humans here for each and humans in poverty
  
  termx<-c("dummy~HazardT+Hazard_futureT+","dummy~HazardTmin+Hazard_futureTmin+","dummy~HazardTmax+Hazard_futureTmax+","clim_up~HazardT+Hazard_futureT+","clim_up~HazardTmin+Hazard_futureTmin+","clim_up~HazardTmax+Hazard_futureTmax+","clim_down~HazardT+Hazard_futureT+","clim_down~HazardTmin+Hazard_futureTmin+","clim_down~HazardTmax+Hazard_futureTmax+","lu_up~HazardT+Hazard_futureT+","lu_up~HazardTmin+Hazard_futureTmin+","lu_up~HazardTmax+Hazard_futureTmax+","lu_down~HazardT+Hazard_futureT+","lu_down~HazardTmin+Hazard_futureTmin+","lu_down~HazardTmax+Hazard_futureTmax+")
 
  minmax<-rep(c("mean","min","max"),length(termx)/3)
  
  for(ww in 1:length(termx)){
    
    ###summarise gains and losses
    results1<-as.data.frame(xtabs(data=resF,paste0(termx[ww],"year+RCP")))
    names(results1)[1:2]<-c("HazardT","Hazard_futureT")
    termx2<-gsub("HazardT","",termx[ww],fixed=TRUE)
    termx2<-gsub("Hazard_futureT","",termx2,fixed=TRUE)
    termx2<-gsub("dummy","hazard_NA",termx2,fixed=TRUE)
    termx2<-gsub("~++","_mean",termx2,fixed=TRUE)
    
    #termx2<-gsub("~","_",termx2,fixed=TRUE)
    termx2<-gsub("+","_",termx2,fixed=TRUE)
    termx2<-gsub("~min_min_","_min",termx2,fixed=TRUE)
    #termx2<-gsub("~mean_","",termx2,fixed=TRUE)
    termx2<-gsub("~max_max_","_max",termx2,fixed=TRUE)
    
    results1$type=termx2#paste0(minmax[ww],"_",termx2)

    results1$year_RCP<-paste(results1$year,results1$RCP,sep="_")
    
    if(ww==1){res2<-results1}  else {res2<-rbind(res2,results1)}
      ###combine
      
    }
  res2$year_RCP_group<-paste(res2$type,res2$year,res2$RCP,sep="_")
  
  ###combine
  
  ##remove very wrong ones - failed models
  res2x<-res2[res2$type %in% c("hazard_NA_mean","hazard_NA_max","hazard_NA_min")  & res2$HazardT==1 ,]
  results_test<-aggregate(res2x$Freq,by=list(res2x$year_RCP_group,res2x$year_RCP),sum)
  results_test$prop<-results_test$x/endsize
  results_test<-results_test[results_test$prop>0.15 , ]
  
  if(nrow(results_test)==0){print("all too small");print(ad2[j]);next}  
  
  ##remove same for present and future
  results2<-res2[res2$HazardT!=res2$Hazard_futureT & res2$year_RCP %in% results_test$Group.2 ,]
  results2$Freq2<-results2$Freq
  results2$Freq2[results2$Hazard_futureT==0]<-results2$Freq2[results2$Hazard_futureT==0]*-1
  results3a<-cbind(results2,read.table(text=results2$type,sep="_",stringsAsFactors = FALSE,col.names=c("haz","updown","minmax")))
  
  ##sum together negative (losses) and gains (positive)
  results3<-aggregate(results3a$Freq2,by=list(results3a$haz,results3a$minmax,results3a$RCP,results3a$year),sum)
  names(results3)<-c("Type","MinMax","RCP","Year","Gain")
  
  ##cut oout clim and lu
  rs3b<-results3[results3$Type=="hazard",]
  rs3b$year_RCP_minmax<-paste(rs3b$Year,rs3b$RCP,rs3b$MinMax,sep="_")
  
  library(reshape2)
  rs3c<-results3[results3$Type!="hazard",]
  rs3d<-dcast(rs3c,MinMax+RCP+Year~Type,value.var="Gain")
  rs3d$year_RCP_minmax<-paste(rs3d$Year,rs3d$RCP,rs3d$MinMax,sep="_")
  
  
  
  ##combine
  results4a<-merge(rs3b,rs3d[,c("year_RCP_minmax","clim","lu")],by.x="year_RCP_minmax",by.y="year_RCP_minmax")
  results4<-merge(results4a,pres_size,by.x="MinMax",by.y="term")
  results4$change<-(results4$Gain+results4$pres_size)/results4$pres_size
  results4$changeClim<-(results4$clim+results4$pres_size)/results4$pres_size
  results4$changeLU<-(results4$lu+results4$pres_size)/results4$pres_size
  results4$disease<-resF$name[1]
    
  fwrite(results4,file=paste0("/home/ucbtdw0/Scratch/New_Global_MAXENT/summaries/data_summary1_",ad2[j],"_",sample(1:1e6,1),".csv"))
  
  #if(is.null(res5)){res5<-results4} else {res5<-rbind(res5,results4)}
  
  
} #### end of j lop

#fwrite(res5,file=paste0("/home/ucbtdw0/Scratch/New_Global_MAXENT/data_summary1_",sample(1:1e6,1),".csv"))

##add disease info
#d1<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\disease_table32c.csv",stringsAsFactors=FALSE)
#d1$name2<-paste(" ",d1$name,sep="")
#d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
#d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4

##read in other data
#d2<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\diseases4a.csv",stringsAsFactors=FALSE)

  ###remove failed models - how? remove and predict using MICE - see how different - chisq
  #lm1<-lm(data=results3,Gain~Year)
  #results3$resid<-scale(lm1$residuals)
  #results3$dn<-dnorm(results3$resid,mean=mean(results3$resid),sd=sd(results3$resid))
  
  ##long to wide
  #res6<-dcast(res5[res5$Type=="hazard" & res5$MinMax =="max",],disease+RCP+Year~Type,value.var="change")
  
  ##merge with all data
  #res6<-merge(res5,d1,by.x="disease",by.y="name",all.x=TRUE,all.y=FALSE)
  #res7<-merge(res6,d2,by.x="disease",by.y="disease",all.x=TRUE,all.y=FALSE)
  
  ##make plots
 # ggplot(res7,aes(x=changeClim,y=changeLU,col=host_type.1))+
    #geom_point(size=log(log(res7$cases_per_year+1)+1)+1)+
 #   geom_point(size=log(log(res7$cases_per_year+1)+1))+
#    theme_cowplot()+
 #   facet_grid(RCP~Year,scales="free")+
 #   xlim(0,2)+
 #   ylim(0,2)+
 #   geom_hline(yintercept=1,lty=2)+
 #   geom_vline(xintercept=1,lty=2)+
 #   theme(panel.grid.major = element_blank(),
 #         panel.grid.minor = element_blank(),
 #         panel.border = element_rect(colour = "black"))+
 #   xlab("Climate Change Effect")+
 #   ylab("Land-use Change Effect")
  
  

