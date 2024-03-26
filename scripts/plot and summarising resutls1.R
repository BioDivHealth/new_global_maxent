library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggthemes)

### plot and summarising resutls
### areas that have increased due to climate versus those that increased due to
### land use. Size of dots - burden estimate?
### panels for each year - facet across year and RCP?
### if over/under threshold due climate versus 

ad1<-list.files("C:\\Users\\david.redding\\Documents\\myriad_data",pattern="all_data3",full.names=TRUE)
ad2<-list.files("C:\\Users\\david.redding\\Documents\\myriad_data",pattern="all_data3",full.names=FALSE)
ad2<-gsub("_all_data3.csv","",ad2)

gl1<-list.files("C:\\Users\\david.redding\\Documents\\myriad_data",pattern="_1.csv",full.names=TRUE)
gl2<-list.files("C:\\Users\\david.redding\\Documents\\myriad_data",pattern="_1.csv",full.names=FALSE)
gl2<-gsub("_1.csv","",gl2)

res5<-NULL
for (j in 12:length(ad1)){
  
  ##load up each disease
  resF<-fread(ad1[j])
  
  if(length(gl1[gl2==resF$name[1]])==0){next}
  
  ##gainloss
  glF<-fread(gl1[gl2==resF$name[1]])
  
  ##remember: res6[,Hazard := clim_present_mean2 * lc_suit_present_mean2 ]
  
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

  ##cells in present day ## add in AUC score
  mean_pres_size<-nrow(resF[HazardT==1 & year_RCP==unique(resF$year_RCP)[1],])
  min_pres_size<-nrow(resF[HazardTmin==1 & year_RCP==unique(resF$year_RCP)[1],])
  max_pres_size<-nrow(resF[HazardTmax==1 & year_RCP==unique(resF$year_RCP)[1],])
  pres_size<-data.frame(term=c("mean","min","max"),pres_size=c(mean_pres_size,min_pres_size,max_pres_size))
      
  ###remove all not in any polygons
  resF2<-resF[!is.na(present2010) | !is.na(presentfuture), ]
  resF2min<-resF[!is.na(present2010min) | !is.na(presentfuturemin), ]
  resF2max<-resF[!is.na(present2010max) | !is.na(presentfuturemax), ]
  
  rm(resF);gc()
  #plot(poly2080)
  #resF3<-resF2[resF2$year==2080,]
  #points(resF3$x,resF3$y)
  #resF2[,dummy:=1]
  
  termx<-c("dummy~HazardT+Hazard_futureT+","clim_up~","clim_down~","lu_up~","lu_down~")
  
  for(ww in 1:length(termx)){
  
    ###summarise gains and losses
  results1<-as.data.frame(xtabs(data=resF2,paste0(termx[ww],"year+RCP")))
  results1$type=paste0("mean_",termx[ww])
  results1min<-as.data.frame(xtabs(data=resF2min,paste0(termx[ww],"year+RCP")))
  results1min$type=paste0("min_",termx[ww])
  results1max<-as.data.frame(xtabs(data=resF2max,paste0(termx[ww],"year+RCP")))
  results1max$type=paste0("max_",termx[ww])
  
  results1b<-rbind(results1,results1min,results1max)

  results1b$year_RCP<-paste(results1b$year,results1b$RCP,sep="_")
  
     if(ww==1){res2<-results1b} else {
    
      if(ww==2){res3<-results1b} else {res3<-rbind(res3,results1b)}
  ###combine
  
    }
  
  }
  
  ## unique(res3$type)
  
  setDT(res3)
  ups<-c("mean_clim_up~","min_clim_up~","max_clim_up~","mean_lu_up~","min_lu_up~","max_lu_up~") 
  downs<-c("mean_lu_down~","min_lu_down~","max_lu_down~","mean_clim_down~","min_clim_down~","max_clim_down~")
  res3[,times1:=ifelse(type %in% ups,1,-1)]
  res3[,c("type2","clim_lu","up_down"):=tstrsplit(type,"_")]
  res3[,Freq2:=Freq*times1]
  
  ###aggregate up and downs due to climate and lu
  res4<-res3[, lapply(.SD, sum, na.rm=TRUE), by=list(year,RCP,type2,clim_lu), .SDcols=c("Freq2") ] 
  names(res4)<-c("Year","RCP","Type","Type2","Gain")
  res4<-as.data.frame(res4)
  
  ##remove very wrong ones - failed models
  results_test<-aggregate(res2$Freq,by=list(res2$type,res2$year,res2$RCP),sum)
  results_test$prop<-results_test$x/endsize
  results_test$year_RCP<-paste(results_test$Group.2,results_test$Group.3,sep="_")
  results_test<-results_test[results_test$prop>0.1 & results_test$prop<1.9, ]
  
  if(nrow(results_test)==0){next}  
  
  ##remove same for present and future
  results2<-res2[res2$HazardT!=res2$Hazard_futureT & res2$year_RCP %in% results_test$year_RCP ,]
  results2<-results2[results2$Freq!=0,]###check not throwing out good data
  results2$Freq2<-results2$Freq
  results2$Freq2[results2$Hazard_futureT==0]<-results2$Freq2[results2$Hazard_futureT==0]*-1
  results3<-aggregate(results2$Freq2,by=list(results2$type,results2$RCP,results2$year),sum)
  names(results3)<-c("Type","RCP","Year","Gain")
  results3$Type<-gsub("_dummy~HazardT+Hazard_futureT+","",results3$Type,fixed=TRUE)
  results3$Type2="hazard"

  
  ##combine
  results4<-rbind(results3,res4)
  resutls4<-merge(results4,pres_size,by.x="Type",by.y="term")
  results4$change=(resutls4$Gain+resutls4$pres_size)/resutls4$pres_size
  results4$disease<-resF2$name[1]
    
  
  if(is.null(res5)){res5<-results4} else {res5<-rbind(res5,results4)}
  
  
} #### end of j lop

fwrite(res5,file=paste0("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\data_summary1_",sample(1:1e6,1),".csv"))

  ###remove failed models - how? remove and predict using MICE - see how different - chisq
  #lm1<-lm(data=results3,Gain~Year)
  #results3$resid<-scale(lm1$residuals)
  #results3$dn<-dnorm(results3$resid,mean=mean(results3$resid),sd=sd(results3$resid))
  
  res6<-dcast(res5[res5$Type2!="hazard" & res5$Type=="max",],disease+RCP+Year~Type2,value.var="change")
  
  ##make plots
  ggplot(res6,aes(x=clim,y=lu,col=disease))+
    geom_point()+
   theme_cowplot()+
    facet_grid(Year~RCP,scales="free")+
    xlim(-2,2)+
    ylim(-2,2)+
    geom_hline(yintercept=0,lty=2)+
    geom_vline(xintercept=0,lty=2)
  
  

