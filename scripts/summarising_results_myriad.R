library(data.table)
library(raster)
library(reshape2)
#library(velox)
#library(ggplot2)
#library(ggpubr)
#library(cowplot)
#library(ggthemes)

##base raster template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


### plot and summarising resutls
### areas that have increased due to climate versus those that increased due to
### land use. Size of dots - burden estimate?
### panels for each year - facet across year and RCP?
### if over/under threshold due climate versus 

#gini<-stack("C:\\Users\\Public\\Documents\\GINI_STACK.tif")
#names(gini)<-read.csv("C:\\Users\\Public\\Documents\\GINI_STACK_NAMES.csv")$x

#gini3<-extract(gini,xyFromCell(template,1:ncell(template)))
#fwrite(gini3,file="C:\\Users\\Public\\Documents\\GINI_STACK.csv")
#names(gini3)<-read.csv("C:\\Users\\Public\\Documents\\GINI_STACK_NAMES.csv")$x
#gini4<-fread(file="C:\\Users\\Public\\Documents\\GINI_STACK.csv")
#gini4[,cell.id:=1:ncell(template)]

#gini5<-melt(gini4,id.vars="cell.id")
#setkey(gini5,value)
#gc()
#gini5<-gini5[!is.na(value),]
#gc()
#gini5[,c("year","scrap","SSP","scrap3"):= tstrsplit(variable,"_")]
#gc()
#gini5[,c("variable", "scrap","scrap3"):=NULL]
#gc()
#gini5[,year:=as.numeric(gsub("X","",year))]
#gc()
#gini5[,cell_by_year:=paste(cell.id,year,sep="_")]


#gini6<-gini5[year==2010,]
#gini6<-gini6[!duplicated(cell.id),]
#gini6[,c("year","SSP","cell_by_year"):=NULL]
#names(gini6)[2]<-"present_value"
#setkey(gini6,cell.id)
#gini5<-gini5[year!=2010,]
#gini5[, future_value2:=value]

#setkey(gini5,cell.id)
#gc()
#gini7<-gini6[gini5]
#gc()
#gini7[,future_value2:=NULL]
#fwrite(gini7,file="C:\\Users\\Public\\Documents\\GINI_by_cell.csv")

gini5<-fread(file="/home/ucbtdw0/Scratch/New_Global_MAXENT/data/GINI_by_cell1.csv")
setkey(gini5,cell_by_year)
#gini5<-dcast(gini4,year_dissme+year+diss_me~SSP,value.var = "value")
#setkey(gini5,year_dissme)


ad1<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/datasets1/",pattern="all_data3",full.names=TRUE)
ad2<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/datasets1/",pattern="all_data3",full.names=FALSE)
ad2<-gsub("_all_data3.csv","",ad2)
ad2<-gsub("_X_",";",ad2)
ad3<-read.table(text=ad2,sep=";",stringsAsFactors = F)$V1

gl1<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",pattern="_1.csv",full.names=TRUE)
gl2<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",pattern="_1.csv",full.names=FALSE)
gl2<-gsub("_1.csv","",gl2)
gl2<-gsub("_X_",";",gl2)
gl3<-read.table(text=gl2,sep=";",stringsAsFactors = F)$V1

##popdgp
#popgdp<-fread("C:\\Users\\Public\\Documents\\future_gdppop_maskedb.csv")
#setkey(popgdp,cell_by_year)

##popdgp
pop<-fread("/home/ucbtdw0/Scratch/New_Global_MAXENT/data/future_pop_maskedb.csv")
names(pop)[9]<-"humans_present"
setkey(pop,cell_by_year)

##gini model (called lm1)
load(file="/home/ucbtdw0/Scratch/New_Global_MAXENT/data/fitted_gini_10th.r")

err = simpleError("Error in read1")
err2 = simpleError("Error in read2")

#res5<-NULL

for (j in sample(1:length(ad1))){
 
  if(paste0("data_summary1_",ad2[j],"_",1,".csv") %in% list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/summaries/",full.names=FALSE)){next}
  
  ##load up each disease
  resF2<-tryCatch(fread(ad1[j]), error=function(err) err)
  if(class(resF2)[1]=="simpleError"){ next }
  
  if(length(gl1[gl3==ad3[j]])==0){next}
  
  ##make cell by year
  resF2[,cell_by_year:=paste(cell.id,year,sep="_")]
  setkey(resF2,cell_by_year)
  
  ###merge
  resF<-pop[resF2]
  setkey(resF,cell_by_year)
  #resF<-popgdp[resF]
  #setkey(resF,cell_by_year)
  
  rm(resF2);gc()
  #merge gini
  resF<-gini5[resF]
  
  ##gainloss
  glF<-tryCatch(fread(gl1[gl3==ad3[j]][1]), error=function(err2) err2)
  if(class(glF)[1]=="simpleError"){ next }
  
  #table(resF$Clim_round)
  #table(resF$LU_round)
  
  ##mean value
  mean_clim<-median(resF$clim_present_mean2[!is.na(resF$present2010) & resF$year_RCP==unique(resF$year_RCP)[1]],na.rm=TRUE)
  mean_LU<-median(resF$lc_suit_present_mean2[!is.na(resF$present2010) & resF$year_RCP==unique(resF$year_RCP)[1]],na.rm=TRUE)
  
  ##set up constants
  thresh1<-glF$AUC_cutoff[1]
  thresh2<-glF$AUC_cutoffmin[1]
  thresh3<-glF$AUC_cutoffmax[1]
  #names(resF)[ !names(resF) %in% c("name","density_mean","density_mean_vector","year","RCP","year_RCP","secondary_future","secondary_present","x","y","subcountr","countr","sub_region","regionPA","surrounding","host_range","defaun","countr2","dummy","Hazard","Hazard_future","HazardT","Hazard_futureT","HazardTmin","Hazard_futureTmin","HazardTmax","Hazard_futureTmax","present2010","present2010min","presentfuturemax","HazardTC","Hazard_futureTC","HazardTminC","Hazard_futureTminC","HazardTmaxC","Hazard_futureTmaxC")]
  

  ##cells in endemic area
  endsize<-nrow(resF[!is.na(countr) & year_RCP==unique(resF$year_RCP)[1],])

# "present2010"                 "present2010min"             
# "present2010max"              "presentfuture"               "presentfuturemin"            "presentfuturemax"            "presentfuturei"             
# "presentfuturemini"           "presentfuturemaxi"
  
  
  ##make RCP
  resF[,c("RCP","scrap"):=tstrsplit(year_RCP,"_")]
  
  #remove unneeded columns
  resF[,c("i.cell.id","i.year","i.cell.id.1","time","scrap"):=NULL]
  
  ###remove all not in any polygons - WORK OUT WITH NEW INNERS AND OUTERS
  resF$HazardT[is.na(resF$present2010)]<-NA
  resF$Hazard_futureT[is.na(resF$presentfuture) & !is.na(resF$presentfuturei)]<-NA
  
  resF$HazardTmin[is.na(resF$present2010min)]<-NA
  resF$Hazard_futureTmin[is.na(resF$presentfuturemin) & !is.na(resF$presentfuturemini)]<-NA
  
  resF$HazardTmax[is.na(resF$present2010max)]<-NA
  resF$Hazard_futureTmax[is.na(resF$presentfuturemax) & !is.na(resF$presentfuturemaxi)]<-NA
  
  resF<-resF[!is.na(resF$HazardT)|!is.na(resF$HazardTmin)|!is.na(resF$HazardTmax)|!is.na(resF$Hazard_futureT)|!is.na(resF$Hazard_futureTmin)|!is.na(resF$Hazard_futureTmax)|!is.na(resF$HazardT), ]
  
  ##remember: res6[,Hazard := clim_present_mean2 * lc_suit_present_mean2 ]
  
  ##make proportion of people in poverty GINI
  ## according to gini assuming a lognormal this is the proporation of people earning less than the per capita GDP
  #resF$present_value[is.na(resF$present_value)]<-0
  #resF$future_value[is.na(resF$future_value)]<-0
  resF[,prop_in_pov:=stats::predict.lm(newdata=data.frame(gini=PRESENT/100),lm1)]
  resF[,prop_in_pov_futureSSP1:=stats::predict.lm(newdata=data.frame(gini=SSP1/100),lm1)]
  resF[,prop_in_pov_futureSSP2:=stats::predict.lm(newdata=data.frame(gini=SSP2/100),lm1)]
  resF[,prop_in_pov_futureSSP3:=stats::predict.lm(newdata=data.frame(gini=SSP3/100),lm1)]
  resF[,prop_in_pov_futureSSP4:=stats::predict.lm(newdata=data.frame(gini=SSP4/100),lm1)]
  resF[,prop_in_pov_futureSSP5:=stats::predict.lm(newdata=data.frame(gini=SSP5/100),lm1)]
  
  ## calculate Exposure
  #"humans2010" "ssp1"  "ssp2" "ssp3" "ssp4" "ssp5"
  resF[,people_in_poverty:=prop_in_pov*humans_present]
  resF[,people_in_povertySSP1:=(prop_in_pov_futureSSP1*ssp1)-people_in_poverty]
  resF[,people_in_povertySSP2:=(prop_in_pov_futureSSP2*ssp2)-people_in_poverty]
  resF[,people_in_povertySSP3:=(prop_in_pov_futureSSP3*ssp3)-people_in_poverty]
  resF[,people_in_povertySSP4:=(prop_in_pov_futureSSP4*ssp4)-people_in_poverty]
  resF[,people_in_povertySSP5:=(prop_in_pov_futureSSP5*ssp5)-people_in_poverty]
  
  resF[,peopleSSP1:=ssp1-humans_present]
  resF[,peopleSSP2:=ssp2-humans_present]
  resF[,peopleSSP3:=ssp3-humans_present]
  resF[,peopleSSP4:=ssp4-humans_present]
  resF[,peopleSSP5:=ssp5-humans_present]
  
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
  
  #plot(poly2080)
  #resF3<-resF2[resF2$year==2080,]
  #points(resF3$x,resF3$y)
  #resF2[,dummy:=1]
  
  ## what about clim up clim down in just areas that are experiencing an increase/decrease - why the change rather than how is changing over whole area
  
  ## also present should just be in present poly - rest should be NA
  
  ## bring in humans here for each and humans in poverty
  
  termx<-c("people_in_povertySSP5~HazardT+Hazard_futureT+","people_in_povertySSP5~HazardTmin+Hazard_futureTmin+","people_in_povertySSP5~HazardTmax+Hazard_futureTmax+","people_in_povertySSP4~HazardT+Hazard_futureT+","people_in_povertySSP4~HazardTmin+Hazard_futureTmin+","people_in_povertySSP4~HazardTmax+Hazard_futureTmax+","people_in_povertySSP3~HazardT+Hazard_futureT+","people_in_povertySSP3~HazardTmin+Hazard_futureTmin+","people_in_povertySSP3~HazardTmax+Hazard_futureTmax+","people_in_povertySSP2~HazardT+Hazard_futureT+","people_in_povertySSP2~HazardTmin+Hazard_futureTmin+","people_in_povertySSP2~HazardTmax+Hazard_futureTmax+","people_in_povertySSP1~HazardT+Hazard_futureT+","people_in_povertySSP1~HazardTmin+Hazard_futureTmin+","people_in_povertySSP1~HazardTmax+Hazard_futureTmax+","peopleSSP5~HazardT+Hazard_futureT+","peopleSSP5~HazardTmin+Hazard_futureTmin+","peopleSSP5~HazardTmax+Hazard_futureTmax+","peopleSSP4~HazardT+Hazard_futureT+","peopleSSP4~HazardTmin+Hazard_futureTmin+","peopleSSP4~HazardTmax+Hazard_futureTmax+","peopleSSP3~HazardT+Hazard_futureT+","peopleSSP3~HazardTmin+Hazard_futureTmin+","peopleSSP3~HazardTmax+Hazard_futureTmax+","peopleSSP2~HazardT+Hazard_futureT+","peopleSSP2~HazardTmin+Hazard_futureTmin+","peopleSSP2~HazardTmax+Hazard_futureTmax+","peopleSSP1~HazardT+Hazard_futureT+","peopleSSP1~HazardTmin+Hazard_futureTmin+","peopleSSP1~HazardTmax+Hazard_futureTmax+","dummy~HazardT+Hazard_futureT+","dummy~HazardTmin+Hazard_futureTmin+","dummy~HazardTmax+Hazard_futureTmax+","clim_up~HazardT+Hazard_futureT+","clim_up~HazardTmin+Hazard_futureTmin+","clim_up~HazardTmax+Hazard_futureTmax+","clim_down~HazardT+Hazard_futureT+","clim_down~HazardTmin+Hazard_futureTmin+","clim_down~HazardTmax+Hazard_futureTmax+","lu_up~HazardT+Hazard_futureT+","lu_up~HazardTmin+Hazard_futureTmin+","lu_up~HazardTmax+Hazard_futureTmax+","lu_down~HazardT+Hazard_futureT+","lu_down~HazardTmin+Hazard_futureTmin+","lu_down~HazardTmax+Hazard_futureTmax+")
  
  #minmax<-rep(c("mean","min","max"),length(termx)/3)
  
  for(ww in 1:length(termx)){
    
    ###summarise gains and losses
    results1<-as.data.frame(xtabs(data=resF,paste0(termx[ww],"year+RCP")))
    names(results1)[1:2]<-c("HazardT","Hazard_futureT")
    type1<-strsplit(termx[ww],"~")[[1]][1]
    
    termx2<-strsplit(termx[ww],"~")[[1]][2]
    termx2<-gsub("HazardT","",termx2,fixed=TRUE)
    termx2<-gsub("Hazard_futureT","",termx2,fixed=TRUE)

    #termx2<-gsub("~","_",termx2,fixed=TRUE)
    termx2<-gsub("++","mean",termx2,fixed=TRUE)
    termx2<-gsub("min+min+","min",termx2,fixed=TRUE)
    #termx2<-gsub("~mean_","",termx2,fixed=TRUE)
    termx2<-gsub("max+max+","max",termx2,fixed=TRUE)
    
    results1$type1=type1
    results1$type2=termx2
    
    results1$year_RCP<-paste(results1$year,results1$RCP,sep="_")
    
    if(ww==1){res2<-results1}  else {res2<-rbind(res2,results1)}
      ###combine
      
  }
  
  res2$year_RCP_group<-paste(res2$type1,res2$type2,res2$year,res2$RCP,sep=";")
  
  ###combine
  
  ##remove very wrong ones - failed models
  res2x<-res2[ res2$type1=="dummy" ,]
  results_test<-aggregate(res2x$Freq,by=list(res2x$year_RCP),sum)
  #results_test$prop<-results_test$x/endsize
  results_test<-results_test[results_test$x!=0 , ]
  
  #if(nrow(results_test)==0){print("all too small");print(ad2[j]);next}  
  
  ##remove same for present and future
  results2<-res2[res2x$year_RCP %in% results_test$Group.1  ,]
  results2$Freq2<-results2$Freq
  results2$Freq2[results2$Hazard_futureT==0]<-results2$Freq2[results2$Hazard_futureT==0]*-1
  
  pres1<-aggregate(res2$Freq[res2$HazardT==1],by=list(res2$year_RCP_group[res2$HazardT==1]),sum)
  fut1<-aggregate(res2$Freq[res2$Hazard_futureT==1],by=list(res2$year_RCP_group[res2$Hazard_futureT==1]),sum)
  
  ##make difference
  fut1$change1=fut1$x-pres1$x
  results3<-cbind(fut1,read.table(text=fut1$Group.1,sep=";",stringsAsFactors = FALSE))
  
  #results3a<-cbind(results2,read.table(text=results2$type,sep="_",stringsAsFactors = FALSE,col.names=c("haz","updown","minmax")))
  
  ##sum together negative (losses) and gains (positive)
  #results3<-aggregate(results2$Freq2,by=list(results2$type1,results2$type2,results2$RCP,results2$year),sum)
  names(results3)<-c("year_RCP_group","future_value","Gain","Type","MinMax","Year","RCP")
  
  ##merge
  #results4<-merge(results3,pres_size,by.x="MinMax",by.y="term")
  #results4$change<-(results4$Gain+results4$pres_size)/results4$pres_size
  results3$disease<-resF$name[1]
  
  ##remove duds
  results4<-results3[results3$future_value!=0,]
  
  #fwrite(results4,file=paste0("/home/ucbtdw0/Scratch/New_Global_MAXENT/summaries/data_summary1_",ad2[j],"_",1,".csv"))
  
  
  #fwrite(results4,file=paste0("C:\\Users\\Public\\Documents\\summaries/data_summary1_",ad2[j],"_",sample(1:1e6,1),".csv"))
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
  
  

