library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(sinaplot)

### plot and summarising results
### areas that have increased due to climate versus those that increased due to
### land use. Size of dots - burden estimate?
### panels for each year - facet across year and RCP?
### if over/under threshold due climate versus 


res1<-list.files("C:\\Users\\david.redding\\Documents\\more_results2\\",pattern="ALL.csv",full.names=TRUE)
res2<-list.files("C:\\Users\\david.redding\\Documents\\more_results2\\",pattern="ALL.csv",full.names=FALSE)
res2<-gsub("per_disease_results_","",res2)
res2<-gsub(".csv","",res2)
res2b<-read.table(text=res2,sep="_",stringsAsFactors = FALSE)

res3<-NULL
for (i in 1:length(res1)){
  res2x<-read.csv(res1[i],stringsAsFactors = F)
  if(length(names(res2x)[names(res2x) %in% c("disease")])==0){res2x$disease=res2b$V1[i]}
  if(is.null(res3)){res3<-res2x} else {res3<-rbind(res3,res2x)}
}

res3a<-res3[res3$model=="present",]
names(res3a)<-paste0(names(res3a),"ALL")
  
res3b<-res3[res3$model!="present",]
res3b2<-read.table(text=res3b$model,sep="X",stringsAsFactors = FALSE)
res3b3<-read.table(text=res3b2$V2,sep="_",stringsAsFactors = FALSE)
names(res3b2)[c(1,3)]<-c("SSP","DISEASE_GROUP")
names(res3b3)<-c("RCP","Year")


res3c<-cbind(res3b[,2:3],(res3b[,4:ncol(res3b)])-res3a[,4:ncol(res3a)],res3b3,res3b2[c(1,3)])

hist(res3c$richness[res3c$disease=="dengue"])

hist(res3c$richness[res3c$disease=="rocio"])

hist(res3c$richness[res3c$disease=="leishmania infantum chagasi"])

hist(res3c$richness[res3c$disease=="fasciola gigantica"])




# get rid of 10 extreme values either side

res5<-res3c[res3c$richness>(-4.65e+05) & res3c$richness<(4.65e+05),]#dcast(res3,MinMax+RCP+Year+disease+pres_size~Type,value.var=c("Gain","change"))#rbind(res4,res3)


##add disease info
d1<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\disease_table32c.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4
d1$Vectored<-1
d1$Outbreak<-0
d1$Vectored[d1$Type %in% c("HOST->HUMAN->HUMAN","HOST->HUMAN")]<-0
d1$Outbreak[d1$Type %in% c("HOST->VECTOR->HUMAN->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN","HOST->HUMAN->HUMAN","HUMAN->VECTOR->HUMAN->HUMAN")]<-1


##read in other data
#d2<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\diseases4a.csv",stringsAsFactors=FALSE)

  ###remove failed models - how? remove and predict using MICE - see how different - chisq
  #lm1<-lm(data=results3,Gain~Year)
  #results3$resid<-scale(lm1$residuals)
  #results3$dn<-dnorm(results3$resid,mean=mean(results3$resid),sd=sd(results3$resid))
  
  ##long to wide
  #res6<-dcast(res5[res5$Type=="hazard" & res5$MinMax =="max",],disease+RCP+Year~Type,value.var="change")
  
  ##merge with all data
  res6<-merge(res5,d1,by.x="disease",by.y="name",all.x=TRUE,all.y=FALSE)
 # res7<-merge(res6,d2,by.x="disease",by.y="disease",all.x=TRUE,all.y=FALSE)
  

  ##make SSP numeric
  res6$SSP<-as.numeric(gsub("ssp","",res6$SSP))
  
  res6$dummy=1
  
  #xtabs(dummy~page+model,res6)
  
  xtabs(dummy~page+disease,res6)
  
  table(res6$page)
  
  unique(paste(res6$page,res6$disease,sep="_"))
  
  
  ##aggregate by page
  res7<-aggregate(res6[,c("richness","disease_richness","humans","cases","deaths","x","y","max_x","max_y","min_x","min_y","RCP","Year","SSP","CFR.low","CFR.high","cases_per_year","spillover_rate2","Outbreak","Vectored")],by=list(res6$page,res6$model),FUN=function(x) mean(x,trim=0.1))
  
  
  
  ggplot(res7[ res7$richness!=0 ,],aes(x=Year,y=richness)) +
    #stat_summary(fun =  mean, geom = "point",size=2,position = position_dodge(3)) + 
    stat_summary(fun =  median, geom = "line",size=2,position = position_dodge(3),col="red") + 
    stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    
    #stat_summary(fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_smooth(method="lm")+
   #geom_point(aes(fill=as.factor(RCP)),size=3,alpha=0.5) +
    #geom_violin()+
    #scale_x_log10()+
    theme_cowplot()+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("Year")+
    ylab("Change in endemic area")+
    geom_abline(intercept=0,slope=0,lty=2)+
    #xlim(2,9)+
    facet_grid(Vectored~RCP,scales="free")
    #facet_grid(Vectored~SSP,scales="free")
  
    
  ggplot(res7[ res7$disease_richness!=0 ,],aes(x=Year,y=disease_richness)) +
    #stat_summary(fun =  mean, geom = "point",size=2,position = position_dodge(3)) + 
    stat_summary(fun =  median, geom = "line",size=2,position = position_dodge(3),col="red") + 
    stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_smooth(method="lm")+
    #geom_point(aes(fill=as.factor(RCP)),size=3,alpha=0.5) +
    #geom_violin()+
    #scale_x_log10()+
    theme_cowplot()+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("Year")+
    ylab("Change in endemic area")+
    geom_abline(intercept=0,slope=0,lty=2)+
    #xlim(2,9)+
    facet_grid(Vectored~RCP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")

  

  ggplot(res7[ res7$cases!=0 ,],aes(x=Year,y=cases)) +
    #stat_summary(fun = mean, geom = "line",size=2,position = position_dodge(3)) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #stat_summary(fun =  median, geom = "line",size=2,position = position_dodge(3),col="red") + 
    #stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    
    geom_smooth()+
    #geom_point(aes(fill=as.factor(RCP)),size=3,alpha=0.5) +
    #geom_violin()+
    #scale_x_log10()+
    theme_cowplot()+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("Year")+
    ylab("Change in Expected Cases")+
    #geom_abline(intercept=0,slope=1,lty=2)
    #xlim(2,9)+
    facet_grid(RCP~SSP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")
  
  
  
  
  ggplot(res7[ res7$deaths!=0 ,],aes(x=Year,y=deaths)) +
    stat_summary(fun = function(x) mean(x,trim=0.1), geom = "line",size=2,position = position_dodge(3)) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_smooth(method="lm")+
    #geom_point(aes(fill=as.factor(RCP)),size=3,alpha=0.5) +
    #geom_violin()+
    #scale_x_log10()+
    theme_cowplot()+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("RCP")+
    ylab("Change in Expected Deaths")+
    #geom_abline(intercept=0,slope=1,lty=2)
    #xlim(2,9)+
    facet_grid(Vectored~RCP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")
  
  
  
  
  