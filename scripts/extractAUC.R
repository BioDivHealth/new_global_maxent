

###run disease validation

library(dismo)

library(maptools)

data(wrld_simpl)

###find all present day
pres1<-list.files("W:\\DRtemp\\per_disease\\",pattern=" present",full.names=TRUE)


template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
values(template)<-1:ncell(template)

t2<-raster(template)

res1<-NULL

##loop
for (i in i:length(pres1)){
  
  
  load(pres1[i])
  
  res6b<-res6[res6$year==2015,]
  
  if(nrow(res6)==0){print(pres1[i]);print("no raw data");next}
    
  
  if(mean(res6$ea,na.rm=TRUE)==0){print(pres1[i]);print("all zero");next}
  
  if(mean(res6$ea,na.rm=TRUE)==1){print(pres1[i]);print("all one");next}
  
  
    if(nrow(res6)<500000){res7<-res6}else{
    
    res7<-res6[sample(1:nrow(res6),100000,replace=FALSE),]
  
  }
    
    res7<-res7[!is.na(res7$present_realised),]
  
    res7<-res7[!is.na(res7$ea),]
  
  
    values(t2)[row.names(res7)]<-res7$present_realised    
    
  
  
  if(nrow(res7)==0){print(pres1[i]);print("all NA");next}
    
    res7<-as.data.frame(res7)
    
  ff<-evaluate(p=res7[res7$ea==1,"present_realised",drop=TRUE],a=res7[res7$ea==0,"present_realised"])@auc
  
  if(is.null(res1)){res1<-data.frame(disease=gsub("X:\\DRtemp\\per_disease\\ ","",pres1[i]),AUC=ff)}else{res1<-rbind(res1,data.frame(disease=gsub("X:\\DRtemp\\per_disease\\ ","",pres1[i]),AUC=ff))}  
  
  print(i)
}


hist(res1$AUC,breaks=20,main="AUC scores of all diseases")
abline(v=0.5,lty=2)

write.csv(res1,file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\AUC1.csv")


