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
for (j in 1:length(ad1)){
  
  ##load up each disease
  resF<-fread(ad1[j])
  
  if(length(gl1[gl2==resF$name[1]])==0){next}
  
  ##gainloss
  glF<-fread(gl1[gl2==resF$name[1]])
  
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
  resF$HazardT[is.na(resF$present2010)]<-NA
  resF$Hazard_futureT[is.na(resF$presentfuture)]<-NA
  resF$HazardTmin[is.na(resF$present2010min)]<-NA
  resF$Hazard_futureTmin[is.na(resF$presentfuturemin)]<-NA
  resF$HazardTmax[is.na(resF$present2010max)]<-NA
  resF$Hazard_futureTmax[is.na(resF$presentfuturemax)]<-NA
  resF<-resF[!is.na(resF$HazardT)|!is.na(resF$HazardTmin)|!is.na(resF$HazardTmax)|!is.na(resF$Hazard_futureT)|!is.na(resF$Hazard_futureTmin)|!is.na(resF$Hazard_futureTmax)|!is.na(resF$HazardT), ]
  
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
  
  #plot(poly2080)
  #resF3<-resF2[resF2$year==2080,]
  #points(resF3$x,resF3$y)
  #resF2[,dummy:=1]
  
  ## what about clim up clim down in just areas that are experiencing an increase/decrease - why the change rather than how is changing over whole area
  
  ## also present should just be in present poly - rest should be NA
  
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
  results_test<-aggregate(res2$Freq,by=list(res2$type,res2$year,res2$RCP),sum)
  results_test$prop<-results_test$x/endsize
  results_test$year_RCP_group<-paste(results_test$Group.1,results_test$Group.2,results_test$Group.3,sep="_")
  results_test<-results_test[results_test$prop>0.1 , ]
  
  if(nrow(results_test)==0){next}  
  
  ##remove same for present and future
  results2<-res2[res2$HazardT!=res2$Hazard_futureT & res2$year_RCP_group %in% results_test$year_RCP_group ,]
  results2$Freq2<-results2$Freq
  results2$Freq2[results2$Hazard_futureT==0]<-results2$Freq2[results2$Hazard_futureT==0]*-1
  results3a<-cbind(results2,read.table(text=results2$type,sep="_",stringsAsFactors = FALSE,col.names=c("haz","updown","minmax")))
  
  results3<-aggregate(results3a$Freq2,by=list(results3a$haz,results3a$minmax,results3a$RCP,results3a$year),sum)
  names(results3)<-c("Type","MinMax","RCP","Year","Gain")

  
  ##combine
  #results4<-rbind(results3,res4)
  results4<-merge(results3,pres_size,by.x="MinMax",by.y="term")
  results4$change<-(results4$Gain+results4$pres_size)/results4$pres_size
  results4$disease<-resF2$name[1]
    
  
  if(is.null(res5)){res5<-results4} else {res5<-rbind(res5,results4)}
  
  
} #### end of j lop

fwrite(res5,file=paste0("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\data_summary1_",sample(1:1e6,1),".csv"))

##add disease info
d1<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\disease_table32c.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4


  ###remove failed models - how? remove and predict using MICE - see how different - chisq
  #lm1<-lm(data=results3,Gain~Year)
  #results3$resid<-scale(lm1$residuals)
  #results3$dn<-dnorm(results3$resid,mean=mean(results3$resid),sd=sd(results3$resid))
  
  ##long to wide
  res6<-dcast(res5[res5$Type2!="hazard" & res5$Type=="max",],disease+RCP+Year~Type2,value.var="change")
  
  ##merge with all data
  res7<-merge(res6,d1,by.x="disease",by.y="name",all.x=TRUE,all.y=FALSE)
  
  ##make plots
  ggplot(res7,aes(x=clim,y=lu,col=Type))+
    #geom_point(size=log(log(res7$cases_per_year+1)+1)+1)+
    geom_point(size=2)+
    theme_cowplot()+
    facet_grid(.~Year,scales="free")+
    xlim(0,2)+
    ylim(0,2)+
    geom_hline(yintercept=1,lty=2)+
    geom_vline(xintercept=1,lty=2)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"))+
    xlab("Climate Change Effect")+
    ylab("Land-use Change Effect")
  
  

