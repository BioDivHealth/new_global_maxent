

###run disease validation

library(dismo)

library(maptools)

library(sp)

data(wrld_simpl)

library(raster)
library(fasterize)
library(sf)
library(data.table)


##read in empress-i point data
point_data<-read.csv("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\all_empress_i_data.csv",stringsAsFactors = FALSE)
point_data$LU<-paste(" ",point_data$name_LU,sep="")

###read disease data
d1<-read.csv("X:\\DRtemp\\disease_table28.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")

###find all present day
pres1<-list.files("X:\\DRtemp\\per_disease\\",pattern=" full",full.names=TRUE)
pres2<-list.files("X:\\DRtemp\\per_disease\\",pattern=" full",full.names=FALSE)
pres2<-gsub(" full.r","",pres2)

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
values(template)<-1:ncell(template)

##change projecion??
wrld_simpl2<-spTransform(wrld_simpl,CRS=projection(template))

###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl2),template)

###create new raster
t2<-ws2
t2[!is.na(t2)]<-0

res1<-NULL

##loop
for (i in 1:length(pres1)){
  
  #t3<-t2
  
  load(pres1[i])
  
  if(nrow(res7)==0){print(pres1[i]);print("no raw data");next}
  if(!"cell.id" %in% names(res7)){next}
  
  pd<-point_data[point_data$LU %in% pres2[i] & point_data$Status=="Confirmed",]
  if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd)} else {rw1<-1000}
  
  ##endemic region
  dis_trans<-d1[d1$name2==pres2[i],]
  countr<-wrld_simpl[wrld_simpl$ISO2 %in% strsplit(dis_trans$countries,",")[[1]],]
  
  ###rasterize endemic countries
  ws3<-fasterize(st_as_sf(countr),template,field="UN")
  res7$reported<-extract(ws3,res7$cell.id)
  
  
  
  
  
  ###sub
  #dis_trans$SUBREGION=paste(sort(unique(countr$SUBREGION)),collapse="_")
  #dis_trans$REGION=paste(sort(unique(countr$REGION)),collapse="_")
  
    #res6b<-res6[res6$year==2015,]
  #t3[res6$cell.id]<-res6$clim_present_mean
  #t3[res6$cell.id]<-res6$realised_present_mean
  
    
  #pdf(file=paste("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\Present_maps\\Present_maps",pres2[i],"_",sample(1:1000,1),".pdf",sep=""),width=15,height=10)
  #print(plot(t3))
  #print(plot(wrld_simpl[wrld_simpl$UN %in% unique(res6$ec),],border="red",add=TRUE))
  #print(plot(countr,border="blue",add=TRUE))
  #if(nrow(pd)>0){points(pd,pch=20,cex=0.5)}
  #dev.off()
  
  #summarise by year and RCP and ssp
  #currently summing over all cells - WHAT ABOUT AREA? means need median and sd ===== CHECK
  #res8<-res7[,.(sum(present_realised,na.rm=TRUE),sum(future_realised_upper,na.rm=TRUE),sum(future_realised_lower,na.rm=TRUE),mean(future_realised_clim_upper,na.rm=TRUE),mean(future_realised_clim_upper_div,na.rm=TRUE),mean(future_realised_clim_lower,na.rm=TRUE),mean(future_realised_clim_lower_div,na.rm=TRUE),sum(future_realised_change_upper,na.rm=TRUE),sum(future_realised_change_lower,na.rm=TRUE),sum(future_realised_change_lower_div,na.rm=TRUE)),.(year,RCP,ssp)]
  #names(res8)[4:ncol(res8)]<-c("present_realised","future_realised_upper","future_realised_lower","future_realised_clim_upper","future_realised_clim_upper_div","future_realised_clim_lower","future_realised_clim_lower_div","future_realised_change_upper","future_realised_change_lower","future_realised_change_lower_div")
  
  res7<-res7[!is.na(res7$reported),]
  
  if(nrow(res7)==0){print("no points in country");next}
  
  res8<-res7[,.(sum(future_realised_upper,na.rm=TRUE),sum(future_realised_lower,na.rm=TRUE),mean(future_realised_clim_upper,na.rm=TRUE),mean(future_realised_clim_lower,na.rm=TRUE),sum(future_realised_change_upper,na.rm=TRUE),sum(future_realised_change_lower,na.rm=TRUE)),.(year,RCP,ssp)]
  names(res8)[4:ncol(res8)]<-c("future_realised_upper","future_realised_lower","future_realised_clim_upper","future_realised_clim_lower","future_realised_change_upper","future_realised_change_lower")
  rm(res7)
  
  ##add in disease data
  res9<-cbind(res8,dis_trans)
  rm(res8)
  save(res9,file=paste("X:/DRtemp/per_disease/",pres2[i],"summary_incountry.r"))
  rm(res9)
  
  print(i)
}


#hist(res1$AUC,breaks=20,main="AUC scores of all diseases")
#abline(v=0.5,lty=2)

#write.csv(res1,file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\AUC1.csv")


