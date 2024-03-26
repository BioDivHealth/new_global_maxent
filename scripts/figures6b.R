library(data.table)
library(sp)

library(gridExtra)
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
res11<-fread(file="C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\final_risk5.csv")
#res11<-res11[!is.na(res11$Hazard),]
res11[,yearRCPdissmename:=paste(year_RCP,diss_me,Group.1,sep="_")]

res11b<-fread(file="C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\final_risk5b.csv")
#res11<-res11[!is.na(res11$Hazard),]
res11b[,yearRCPdissmename:=paste(year_RCP,diss_me,Group.1,sep="_")]

#remove any overlap
res11b<-res11b[!res11b$Group.1 %in% res11$Group.1,]
res11<-rbindlist(list(res11,res11b))

###lose old columns
##set(res11,,10:34,NULL)

quanttrim<-function(x,thresh=0.99){
  x[x>quantile(x,thresh,na.rm=TRUE)]<-quantile(x,thresh,na.rm=TRUE)
  x[x<quantile(x,1-thresh,na.rm=TRUE)]<-quantile(x,1-thresh,na.rm=TRUE)
  return(x)
}

res11[,Exposure:=log(humans2010+1)]
res11[,ExposureSSP1:=log(pop_ssp1+1)]
res11[,ExposureSSP2:=log(pop_ssp2+1)]
res11[,ExposureSSP3:=log(pop_ssp3+1)]
##what does this mean?
#res11[,ExposureSSP1:=quanttrim((log(pop_ssp1+1)-Exposure)/Exposure)]
#res11[,ExposureSSP2:=quanttrim((log(pop_ssp2+1)-Exposure)/Exposure)]
#res11[,ExposureSSP3:=quanttrim((log(pop_ssp3+1)-Exposure)/Exposure)]

res11[,Vulnerability:=log((gdp2010/(humans2010+0.01))+1)]
#res8$Vulnerability<-max(res8$Vulnerability,na.rm=TRUE)-res8$Vulnerability
res11[,VulnerabilitySSP1:=log((gdp_ssp1/(pop_ssp1+0.01))+1)]
res11[,VulnerabilitySSP2:=log((gdp_ssp2/(pop_ssp2+0.01))+1)]
res11[,VulnerabilitySSP3:=log((gdp_ssp3/(pop_ssp3+0.01))+1)]


res12<-res11[!is.na(year),]
#res12[,c("RCP","year","dissme","name"):=tstrsplit(yearRCPdissmename,"_")]
res12[,year_RCP_name:=paste(year,RCP,Group.1,sep="_")]
res12[,year_RCP_dissme:=paste(year,RCP,dissme,sep="_")]
res12[,name_dissme:=paste(Group.1,dissme,sep="_")]

res12$year<-as.numeric(res12$year)

### change in number of admins??
res12[,area:=1]
#aggregate(res15$area,by=list(res15$uni2),sum)

names(res12)

##create risk present
#res12[,RiskPEW:= (8*Hazard) + (1*(gdp2010)) + (1*(humans2010))]
#res12[,RiskP:= HazardT * scale(Vulnerability) * scale(Exposure)]



##for maps
#res13<-aggregate(res12[,2:34],by=list(res12$year_RCP_dissme),mean,na.rm=TRUE)
#setDT(res13)
#res13[,c("year","RCP","dissme"):=tstrsplit(Group.1,"_")]

##make data. frame
res12b<-as.data.frame(res12)

###create 2010
resP<-res12b[!duplicated(res12b$name_dissme),]
resP$year<-2010
#resP$RiskF<-resP$RiskP
resP$Hazard_future<-resP$Hazard
resP$Hazard_futureT<-resP$HazardT
#resP$Vulnerability<-resP$gdp2010
#resP$Exposure<-resP$humans2010
resP$RCP<-rep("0",nrow(resP))
resP$SSP<-rep("present",nrow(resP))

##create single SSP column
res14a<-res12b
#res14a$RiskF<-res14a$RiskSSP1
res14a$Vulnerability<-res14a$VulnerabilitySSP1
res14a$Exposure<-res14a$ExposureSSP1
res14a$SSP<-rep("SSP1",nrow(res14a))

res14b<-res12b
#res14b$RiskF<-res14b$RiskSSP2
res14b$Vulnerability<-res14b$VulnerabilitySSP2
res14b$Exposure<-res14b$ExposureSSP2
res14b$SSP<-rep("SSP2",nrow(res14b))

res14c<-res12b
#res14c$RiskF<-res14c$RiskSSP3
res14c$Vulnerability<-res14c$VulnerabilitySSP3
res14c$Exposure<-res14c$ExposureSSP3
res14c$SSP<-rep("SSP3",nrow(res14c))


res14d<-rbind(resP,res14a,res14b,res14c)
setDT(res14d)
#res15<-res15[,-c(54:60,62:64),with=F]

#rmake this
res14d[,year_RCP_SSP_name:=paste(year,RCP,SSP,Group.1,sep="_")]


#remake this column
#res14[year_RCP_name:=paste(year,RCP,name,sep="_")]
res14d[,Vuln:=1-(Vulnerability/max(Vulnerability))]

res14d[,Risk:= Hazard_futureT * Vuln * Exposure]
#res14d[,RiskFEW:= (8*Hazard_future) + (1*(Vulnerability)) + (1*(Exposure))]

####create cutoff column
#res14d[,thresh:=ifelse(RiskF>AUC_cutoff,1,0)]

##threshold is key - compare threholded and non-thredholded or just current and see the different
##difdference in expdected contact - same per grid cell contact rate?
## theshold using weighted background points - less more more in vulnerable areas - more in high population areas - postive and negative weights


### agregate to find mean

### graphs of per disease

##for graphs
res14<-aggregate(res14d[,c("Hazard","Hazard_future","HazardT","Hazard_futureT","Exposure","Vulnerability","Risk","area")],by=list(res14d$year_RCP_SSP_name),sum,na.rm=TRUE)
setDT(res14)
res14[,c("year","RCP","SSP","name"):=tstrsplit(Group.1,"_")]
res14$year<-as.numeric(res14$year)

##combine with disease data table
setkey(res14,name)
res15<-d1[res14]

#threshold?? or at res12 stage???
res15[,burd:= cases_per_year * (Risk) * ((CFR.low+CFR.high)/2)]

##make broad grouping
res15[,group3:=group,]
res15[str_detect(res15[["group"]],"vir"),"group3"]<-"virus"
res15[,Type2:="HOST",]
res15[str_detect(res15[["Type"]],"VECTOR"),"Type2"]<-"VECTOR"

res15[,RCP_SSP:=paste(RCP,SSP,sep="_")]

resx<-res15
uuu<-unique(res14$name)

for (i in 1:length(uuu)){
  
  res15<-resx
  res15<-res15[res15$name==uuu[i],]
  
  
  haz<-ggplot(data=res15,aes(x=year,y=Hazard_future,color=RCP))+
    geom_smooth(method="loess",se=TRUE,span=1)+
    geom_point()+
    #facet_grid(group3~burden,scales="free_y")+
    #geom_hline(yintercept = 1,lty=2)+
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.title = element_blank())+
    labs(title="Hazard")
  
  ris1<-ggplot(data=res15,aes(x=year,y=Risk,col=SSP))+
    geom_smooth(method="loess",se=TRUE,span=1)+
    geom_point()+
    #facet_grid(group3~burden,scales="free_y")+
    geom_hline(yintercept = 1,lty=2)+
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.title = element_blank())+
    labs(title="Risk")
  
  area1<-ggplot(data=res15,aes(x=year,y=Vulnerability,col=SSP))+
    geom_smooth(method="loess",se=TRUE,span=1)+
    geom_point()+
    #facet_grid(group3~burden,scales="free_y")+
    #geom_hline(yintercept = 1,lty=2)+
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.title = element_blank())+
    labs(title="Vulnerability")
  
  burd1<-ggplot(data=res15,aes(x=year,y=Exposure,col=SSP))+
    geom_smooth(method="loess",se=TRUE,span=1)+
    geom_point()+
    #facet_grid(group3~burden,scales="free_y")+
    #geom_hline(yintercept = 1,lty=2)+
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.title = element_blank())+
    labs(title="Exposure")
  
png(paste("C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\test_plots\\plot4_",uuu[i],".png",sep=""))

grid.arrange(haz,ris1,area1,burd1,nrow=2,top=paste(uuu[i],sep=" ") ,bottom="Year")

dev.off()
  
  
  print(i)
  
  
}##end i loop


#set(res11, i=which(str_detect(res11[["group"]],"vir")), j=group3, value="virus")

### force through 1 by dividing hazardP by itself? create a new data.frame with RCP as grouping and 2010 as year?

res15<-res

## figure 1 ggplot 
ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi",],aes(x=as.numeric(year),y=Risk,group=name,col=SSP))+
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

ggplot(data=res15[!is.na(spillover_rate2)  & group3!="fungi" , ],aes(x=as.numeric(year),y=Hazard_future/Hazard,col=RCP))+
  geom_smooth(method="loess")+
  facet_grid(group3~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) , ],aes(x=as.numeric(year),y=Risk,col=SSP))+
  geom_smooth(method="loess")+
  facet_grid(group3~burden,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())


ggplot(data=res15[!is.na(spillover_rate2) , ],aes(x=as.numeric(year),y=burd,col=SSP))+
  geom_smooth(method="loess")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggplot(data=res15[!is.na(spillover_rate2) , ],aes(x=as.numeric(year),y=burd,col=group3))+
  geom_smooth(method="loess")+
  facet_grid(.~group3,scales="free_y")+
  geom_hline(yintercept = 1,lty=2)+
  theme_minimal()+
  theme(panel.grid = element_blank())


ggplot(data=res15[!is.na(spillover_rate2) & group3!="fungi" & RCP==8.5, ],aes(x=as.numeric(year),y=RiskF/RiskP,group=name,col=group3))+
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
jpeg(file=paste("C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\risk_maps\\",pres2[i],"_",sample(1:1000,1),"Hazard.jpg",sep=""),width=2500,height=1000)
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

jpeg(file=paste("C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\risk_maps\\",pres2[i],"_",sample(1:1000,1),"Threshold_Hazard.jpg",sep=""),width=2500,height=1000)
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



