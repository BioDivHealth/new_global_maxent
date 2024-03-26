library(maptools)
library(data.table)
library(sp)

data(wrld_simpl)
wrld_simpl[wrld_simpl$SUBREGION==5,"REGION"]<-20

library(raster)
library(fasterize)
library(dismo)
library(sf)
library(rgdal)

library(viridis)
library(sp)

library(ggplot2)
library(velox)

library(stringr)

#read in table of data about diseases
d1<-read.csv("C:\\Users\\David\\Dropbox\\disease_table32b.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")
d1$name2<-gsub(" angiostrongylus costaricensis "," angiostrongylus costaricensis",d1$name2)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.001,33,99999,999999999999),labels=FALSE)/4
setDT(d1)
d1[,burden:=ifelse(spillover_rate2==0.25,"potential/n0 cases",ifelse(spillover_rate2==0.5,"rare/n1-10 cases",ifelse(spillover_rate2==0.75,"endemic/n10-10000 cases","hyperendemic/n10000-100000000 cases")))]
setkey(d1,name)

###read in results file
res11<-fread(file="C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\final_risk4.csv")

res11[,yearRCPdissmename:=paste(RCP,year,diss_me,name,sep="_")]


res11[,RiskP:= (Hazard) + (vw*(Vulnerability)) + (ew*(Exposure))]
res11[,RiskSSP1:= (Hazard) + (vw*(VulnerabilitySSP1))* (ew*(ExposureSSP1))]
res11[,RiskSSP2:= (Hazard) + (vw*(VulnerabilitySSP2))* (ew*(ExposureSSP2))]
res11[,RiskSSP3:= (Hazard) + (vw*(VulnerabilitySSP3))* (ew*(ExposureSSP3))]


res12<-res11[!is.na(year),names(res11)[c(100,3:6,23:34,81:97)],with=F]
res12[,c("RCP","year","dissme","name"):=tstrsplit(yearRCPdissmename,"_")]
res12[,year_RCP_name:=paste(year,RCP,name,sep="_")]
res12[,year_RCP_dissme:=paste(year,RCP,dissme,sep="_")]


### change in number of admins??
res12[,area:=1]
#aggregate(res15$area,by=list(res15$uni2),sum)



##for maps
#res13<-aggregate(res12[,2:34],by=list(res12$year_RCP_dissme),mean,na.rm=TRUE)
#setDT(res13)
#res13[,c("year","RCP","dissme"):=tstrsplit(Group.1,"_")]


##for graphs
res14<-aggregate(res12[,c(2:34,41)],by=list(res12$year_RCP_name),sum,na.rm=TRUE)
setDT(res14)
res14[,c("year","RCP","name"):=tstrsplit(Group.1,"_")]
setkey(res14,name)
res14$year<-as.numeric(res14$year)


##combine with disease data table
res14<-d1[res14]

##make broad grouping
res14[,group3:=group,]
res14[str_detect(res14[["group"]],"vir"),"group3"]<-"virus"
res14[,Type2:="HOST",]
res14[str_detect(res14[["Type"]],"VECTOR"),"Type2"]<-"VECTOR"

##create single SSP column
res14a<-res14
res14a$RiskF<-res14a$RiskSSP1
res14a$SSP<-"SSP1"
res14b<-res14
res14b$RiskF<-res14b$RiskSSP2
res14b$SSP<-"SSP2"
res14c<-res14
res14c$RiskF<-res14c$RiskSSP3
res14c$SSP<-"SSP3"

res15<-rbind(res14a,res14b,res14c)
#res15<-res15[,-c(54:60,62:64),with=F]

###create 2010
resP<-res15[year==2030 ,]
resP$year=2010
resP$RiskF<-resP$RiskP
resP$Hazard_future<-resP$Hazard
res15<-rbind(res15,resP)


#threshold?? or at res12 stage???

res15[,burd:= cases_per_year * (RiskF/RiskP)]


##need to redo name so powassan is the same and make decisions accross all of them
coefs_area<-res15[AUC_best.stat>0.65,list(intercept=coef(lm(log(area+1)~year))[1], coef=coef(lm(log(area+1)~year))[2],prob=(summary(lm(area~year)))$coefficients[2,4]),by=name]

hist(coefs_area$coef[coefs_area$prob<0.05])

coefs_hazard<-res15[AUC_best.stat>0.65,list(intercept=coef(lm(Hazard_future~year))[1], coef=coef(lm(Hazard_future~year))[2],prob=(summary(lm(Hazard_future~year)))$coefficients[2,4]),by=name]

setorder(coefs_hazard,coef)

hist(coefs_hazard$coef[coefs_hazard$prob<0.05])

coefs_risk<-res15[AUC_best.stat>0.65,list(intercept=coef(lm(RiskF~year))[1], coef=coef(lm(RiskF~year))[2],prob=(summary(lm(RiskF~year)))$coefficients[2,4]),by=name]

hist(coefs_risk$coef[coefs_risk$prob<0.05])


coefs_risk<-res15[AUC_best.stat>0.65,list(intercept=coef(lm(RiskF~year))[1], coef=coef(lm(RiskF~year))[2],prob=(summary(lm(RiskF~year)))$coefficients[2,4]),by=name]

hist(coefs_risk$coef[coefs_risk$prob<0.05])



#set(res11, i=which(str_detect(res11[["group"]],"vir")), j=group3, value="virus")

### force through 1 by dividing hazardP by itself? create a new data.frame with RCP as grouping and 2010 as year?

## figure 1 ggplot 
ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi"& SSP=="SSP2" & RCP==8.5,],aes(x=as.numeric(year),y=Hazard_future/Hazard,group=name))+
  geom_smooth(method="loess",se=FALSE,span=1)+
  facet_grid(group3~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi" & SSP=="SSP2",],aes(x=as.numeric(year),y=(Hazard_future/Hazard),group=name))+
  geom_smooth(method="loess",se=FALSE,span=1)+
  facet_grid(RCP~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi" & SSP=="SSP2" & RCP==8.5, ],aes(x=as.numeric(year),y=Hazard_future/Hazard,group=name))+
  geom_smooth(method="loess",se=FALSE,span=1)+
  facet_grid(Type2~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+ 
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) , ],aes(x=as.numeric(year),y=Hazard_future/Hazard,col=RCP))+
  geom_smooth(method="loess")+
  facet_grid(group3~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) , ],aes(x=as.numeric(year),y=RiskF/RiskP,col=SSP))+
  geom_smooth(method="loess")+
  facet_grid(group3~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())


ggplot(data=res15[!is.na(spillover_rate2) , ],aes(x=as.numeric(year),y=burd,col=SSP))+
  geom_smooth(method="loess")+
  facet_grid(group3~.,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())



ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi" & RCP==8.5, ],aes(x=as.numeric(year),y=RiskF/RiskP,group=name))+
  geom_smooth(method="loess",se=FALSE)+
  facet_grid(Type2~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi", ],aes(x=as.numeric(year),y=RiskF/RiskP,group=name))+
  geom_smooth(method="loess",se=FALSE)+
  facet_grid(RCP~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi" & RCP==8.5,],aes(x=as.numeric(year),y=RiskF/RiskP,group=name))+
  geom_smooth(method="loess",se=FALSE)+
  facet_grid(SSP~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi" & RCP==8.5,],aes(x=as.numeric(year),y=area,group=name))+
  geom_smooth(method="loess",se=FALSE)+
  facet_grid(Type2~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())






###after optimum 
ad1c$RiskPT<-ad1c$Risk1
ad1c$RiskPT[ad1c$RiskPT>thr1[[2]]]<-1
ad1c$RiskPT[ad1c$RiskPT!=1]<-0


##when best breaks are found 
jpeg(file=paste("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\risk_maps\\",pres2[i],"_",sample(1:1000,1),"Hazard.jpg",sep=""),width=2500,height=1000)
print(
  ggplot(data = ad1c) +
    geom_sf(aes(fill = Hazard),colour=NA) +
    scale_fill_viridis_c(option = "plasma")+
    geom_sf(data=countr2,aes(),colour="white", alpha=0,size=1)+
    geom_sf(data=pd2,aes(),cex=0.1)+
    labs(title=paste(pres2[i]," - AUC ",round(e@auc,2),"  AUC points ",round(epd_auc,2)," ",names(thresh)," ",thresh[[1]]," TSS ",TSSv, sep=""))+
    theme(plot.title=element_text(size=32),axis.text = element_text(size=28))
  #+	facet_wrap(.~year_RCP)
)
dev.off()

jpeg(file=paste("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\risk_maps\\",pres2[i],"_",sample(1:1000,1),"Threshold_Hazard.jpg",sep=""),width=2500,height=1000)
print(
  ggplot(data = ad1c) +
    geom_sf(aes(fill = RiskPT),colour=NA) +
    scale_fill_viridis_c(option = "plasma")+
    geom_sf(data=countr2,aes(),colour="grey", alpha=0,size=1)+
    geom_sf(data=pd2,aes(),cex=0.1)+
    labs(title=paste(pres2[i]," - AUC ",round(e@auc,2),"  AUC points ",round(epd_auc,2)," ",names(thresh)," ",thresh[[1]]," TSS ",TSSv, sep=""))+
    theme(plot.title=element_text(size=32),axis.text = element_text(size=28))
  #+	facet_wrap(.~year_RCP)
)
dev.off()




##create final line plots (Hazard and Risk)

## create global summed risk map by aggregating by adl3 (Hazard and Risk)

##divide hazard (non weighted) into

###
#Emerging
#0-100
#Established
#100-10,000
#Endemic
#10000-2,500,000
#Hyperendemic
#2,500,000+



