# Load libraries ----
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggthemes)

### plot and summarising results
### areas that have increased due to climate versus those that increased due to
### land use. Size of dots - burden estimate?
### panels for each year - facet across year and RCP?
### if over/under threshold due climate versus 

# Read in data ----
res1<-list.files("data/summaries/",pattern="_1.csv",full.names=TRUE)
res2<-list.files("data/summaries/",pattern="_1.csv",full.names=FALSE)
res2<-gsub("data_summary1_","",res2)
res2<-gsub(".csv","",res2)
res2b<-read.table(text=res2,sep=";",stringsAsFactors = FALSE)

res3<-NULL
for (i in 1:length(res1)){
  res2x<-read.csv(res1[i],stringsAsFactors = F)
  if(length(names(res2x)[names(res2x) %in% c("disease")])==0){res2x$disease=res2b$V1[i]}
  if(is.null(res3)){res3<-res2x} else {res3<-rbind(res3,res2x)}
}

#only the first time - temp remove
res5<-res3#dcast(res3,MinMax+RCP+Year+disease+pres_size~Type,value.var=c("Gain","change"))#rbind(res4,res3)


##add disease info
d1<-read.csv("data/disease_table32c.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4

##read in other data
d2<-read.csv("data/diseases4a.csv",stringsAsFactors=FALSE)

  ###remove failed models - how? remove and predict using MICE - see how different - chisq
  #lm1<-lm(data=results3,Gain~Year)
  #results3$resid<-scale(lm1$residuals)
  #results3$dn<-dnorm(results3$resid,mean=mean(results3$resid),sd=sd(results3$resid))
  
  ##long to wide
  #res6<-dcast(res5[res5$Type=="hazard" & res5$MinMax =="max",],disease+RCP+Year~Type,value.var="change")
  
  ##merge with all data
  res6<-merge(res5,d1,by.x="disease",by.y="name",all.x=TRUE,all.y=FALSE)
  res7<-merge(res6,d2,by.x="disease",by.y="disease",all.x=TRUE,all.y=FALSE)
  
  ##add present back in
  res7$present_value=(res7$Gain*-1)+res7$future_value
  
  #"peopleSSP5"            "peopleSSP4"            "clim_down"             "lu_down"               "dummy"                 "lu_up"                
  #"people_in_povertySSP1" "people_in_povertySSP2" "peopleSSP2"            "peopleSSP3"            "people_in_povertySSP3" "clim_up"              
  #"people_in_povertySSP4" "people_in_povertySSP5" "peopleSSP1"           
  
  ##cases per year per unit
  res7$caseperunit=res7$cases_per_year/res7$present_value
  
  
  ##make plots
  ggplot(res7[ res7$Type.x %in% c("people_in_povertySSP1","people_in_povertySSP2","people_in_povertySSP3","people_in_povertySSP4","people_in_povertySSP5")  & !is.na(res7$host_type.1)   ,],aes(x=log(present_value+1),y=log(future_value+1),col=host_type.1))+
    #geom_point(size=log(log(res7$cases_per_year+1)+1)+1)+
    stat_summary(fun = mean, geom = "point",size=1.5,position = position_dodge(3)) + 
    stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_point(size=(log(res7$cases_per_year[ res7$MinMax=="max" & res7$Type.x=="dummy"  & !is.na(res7$host_type.1)]+1)/3))+
    theme_cowplot()+
    facet_grid(RCP~Year,scales="free")+
    xlim(2.5,22.5)+
    ylim(2.5,22.5)+
    geom_abline(intercept=0,slope=1,lty=2)+
    #geom_hline(yintercept=1,lty=2)+
    #geom_vline(xintercept=1,lty=2)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"))+
    ylab("Future")+
    xlab("Present")
  
  
  res7$extra_cases<-(res7$cases_per_year*res7$changeLU)-res7$cases_per_year
  
  quantile(res7$extra_cases,0.01)
  quantile(res7$extra_cases,0.99)
  
  #res7<-res7[res7$extra_cases>(-500) , ]
  #res7<-res7[res7$extra_cases<(1000) , ]
  
    
  ##make plots
  ggplot(res7[ res7$MinMax=="max" & res7$Type.x=="dummy"  & !is.na(res7$host_type.1) ,],aes(x=Year,y=Gain*cases_per_year,col=as.factor(disease)))+
    stat_summary(fun = mean, geom = "point",size=1.5,position = position_dodge(3)) + 
    stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    theme_cowplot()+
    facet_wrap(RCP~host_type.1,scales="free",ncol=3)+
    #geom_hline(yintercept=1,lty=2)+
    #geom_vline(xintercept=1,lty=2)+
    theme(legend.position = "none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"))+
    xlab("Year")+
    ylab("change in cases per year")

