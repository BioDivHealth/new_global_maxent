

###run disease validation

library(dismo)

library(maptools)

library(sp)

data(wrld_simpl)

library(raster)
library(fasterize)
library(dismo)

##get gini data
## move cases in poor country
## worse reporting in poor countries
gini<-fread("X:/DRtemp/ssp_gini_countries2.csv",header=TRUE)
cdata<-merge(wrld_simpl@data,gini,by.x="ISO3",by.y="ISO3",all.x=TRUE)
cdata<-cdata[order(match(cdata$NAME,wrld_simpl$NAME)),]
cdata[,(15:ncol(cdata))] <- lapply(cdata[,(15:ncol(cdata))], function(x) ifelse(is.na(x), mean(x, trim=0.1,na.rm = TRUE), x))
identical(cdata$NAME,wrld_simpl$NAME)
wrld_simpl@data<-cdata

##read in empress-i point data
point_data<-read.csv("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\all_empress_i_data.csv",stringsAsFactors = FALSE)
point_data$LU<-paste(" ",point_data$name_LU,sep="")

###read disease data
d1<-read.csv("X:\\DRtemp\\disease_table28.csv",stringsAsFactors=FALSE)
d1$name2<-paste(" ",d1$name,sep="")

###find all present day
pres1<-list.files("X:\\DRtemp\\per_disease\\",pattern=" present",full.names=TRUE)
pres2<-list.files("X:\\DRtemp\\per_disease\\",pattern=" present",full.names=FALSE)
pres2<-gsub(" present.r","",pres2)

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
values(template)<-1:ncell(template)

##change projecion??
wrld_simpl2<-spTransform(wrld_simpl,CRS=projection(template))

###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl2),template)
#plot(ws2)

###create new raster
t2<-ws2
t2[!is.na(t2)]<-0


res1<-NULL

##loop
for (i in 1:length(pres1)){
  

  load(pres1[i])

  ##create baseline raster
  t3<-t2

  if(nrow(res6)==0){print(pres1[i]);print("no raw data");next}
  if(!"cell.id" %in% names(res6)){next}
  
  pd<-point_data[point_data$LU %in% pres2[i] & point_data$Status=="Confirmed",]
  if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd)} else {rw1<-1000}
  
  ##endemic region
  dis_trans<-d1[d1$name2==pres2[i],]
  countr<-wrld_simpl[wrld_simpl$ISO2 %in% strsplit(dis_trans$countries,",")[[1]],]
 
  ##get gini layer
  gini3<-fasterize(st_as_sf(countr),template,field="2010_GINI_PRESENT_C")
  gini3<-gini3/100
  res6$gini<-extract(gini3,res6$cell.id)
  res6$gini[is.na(res6$gini)]<-median(res6$gini,na.rm=TRUE)
  
  ###rasterize endemic countries
  ws3<-fasterize(st_as_sf(countr),template,field="UN")
  res6$reported<-extract(ws3,res6$cell.id)
  
  ##test for regions
  res6$reported<-extract(ws3,res6$cell.id)
  
  
  ###sub
  dis_trans$SUBREGION=paste(sort(unique(countr$SUBREGION)),collapse="_")
  dis_trans$REGION=paste(sort(unique(countr$REGION)),collapse="_")
  
  #reg1<-wrld_simpl[wrld_simpl$REGION %in% strsplit(dis_trans$REGION,"_")[[1]],]
  #reg2<-fasterize(st_as_sf(wrld_simpl),template,field="REGION")
  #res6$region<-extract(reg2,res6$cell.id)
  #res6<-res6[!is.na(res6$region),]
  #aggregate(res6$clim_present_mean,by=list(res6$region),sum)
  
  
    ##seven different thresholds
  ### what about sensitivity testing - find best values?
    ### try with just clim not realised
    ### turn into just arena of contact
    ### what about GINI?
  system.time(
  for (xx in 1:1){
    
    ## DON'T NEED TO RUN CELLS WITH 0 CLIM
    ## (MAYBE) DON'T NEED TO CREATE GAS MOEL WITH THESHOLDS AS ZERO OR 1 ONLY VARIATION COMES FROM lAND-USE
    ## USE UPPER AND LOWER AS MEASURES OF CERTAINANTY
    ## THOSE WITHIN TWO SD OF ZERO - EFFECTIVELY 0??
    
    ## What about cuts? no hosts, medium, high (0,0.5,0.75) (use cut (midpoint)=TRUE)
  
    breaks1=sample(seq(0,max(res6$clim_present_mean,na.rm=TRUE)-0.1,by=0.01),replace=FALSE,1)
    breaks2=sample(seq(0,max(res6$clim_present_mean_vector,na.rm=TRUE)-0.1,by=0.01),replace=FALSE,1)
    
    ##cutoff #choose threshold
    #tt<-sample(1:4,1)
    res6[,clim_present_mean2:= ifelse(clim_present_mean<breaks1,0,clim_present_mean)]
    res6[,clim_present_mean_vector2:= ifelse(clim_present_mean_vector<breaks2,0,clim_present_mean)]
   
    ##rid of zeros
    ## make a  raster of arena
    res6b<-res6[!(clim_present_mean2==0),]
    
    ##abandon dispersal?
    ## create new value from clim * mean density
    res6b[,realised_present_meanH:= clim_present_mean2 * lc_suit_present_mean * density_mean]
    res6b[,realised_present_meanV:= clim_present_mean_vector2 * lc_suit_present_mean_vector * density_mean_vector]
    
    ### present day humans
    ### gass model parameters
    ttmean<-rnorm(1,0,2.5)
    ttspeed<-rnorm(1,0,2.5)
    ttd<-rnorm(1,0,2.5)
    ttmeanv<-rnorm(1,0,2.5)
    ttspeedv<-rnorm(1,0,2.5)
    ttdv<-rnorm(1,0,2.5)
    ttspeedh<-rnorm(1,0,2.5)
    ttdh<-rnorm(1,0,2.5)
    
    ##run gas model - any point in variation not going to influence spatial contact patterns
    res6b[,present_realisedT:= c(((realised_present_mean2*ttmean)*abs(speed_mean*ttspeed)*abs(d_mean*ttd)) * (realised_present_mean_vector2*abs(speed_mean_vector*ttspeedv))*abs(d_mean_vector*ttdv) *(secondary_present) * ((humans2010/5.6) * abs(0.0001*ttdh) * abs(5*ttspeedh))),]
    
    #impact of gini
    res6b[,present_realisedT2:=present_realisedT*((gini^2)+(abs(rnorm(1,mean(gini^2),1))))]
    
    ##get rid of NAs
    res6b$present_realisedT2[is.na(res6$present_realisedT2)]<-0
    
    #res6b<-res6[res6$year==2015,]
    #t3[res6$cell.id]<-res6$clim_present_mean
    t3[res6b$cell.id]<-res6b$abundance_present_mean
  
    #pdf(file=paste("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\Present_maps\\Present_maps",pres2[i],"_",sample(1:1000,1),".pdf",sep=""),width=15,height=10)
    #print(plot(t3))
    #print(plot(wrld_simpl[wrld_simpl$UN %in% unique(res6$ec),],border="red",add=TRUE))
    #print(plot(countr,border="blue",add=TRUE))
    #if(nrow(pd)>0){points(pd,pch=20,cex=0.5)}
    #dev.off()

    ffR<-c();  ffR2<-c()
    for (z  in 1:5){
    
        rand1<-spsample(wrld_simpl[wrld_simpl$UN.x %in% unique(res6$ec),],rw1*10,type="random")
        
        
        rand2<-extract(t3,rand1)
        real2<-spsample(countr,rw1,type="random")
        real3<-extract(t3,real2)

      if(nrow(pd)>0){
    
      real1<-extract(t3,pd)
      res6c<-rbind(data.frame(ea=rep("0",length(rand1)),value=rand2),data.frame(ea=rep("1",rw1),value=real1))
      ffR2[z]<-evaluate(p=res6c[res6c$ea==1,"value",drop=TRUE],a=res6c[res6c$ea==0,"value"])@auc
  
      }else {ffR2[z]<-NA}
      
    res6b<-rbind(data.frame(ea=rep("0",length(rand1)),value=rand2),data.frame(ea=rep("1",rw1),value=real3))
    ffR[z]<-evaluate(p=res6b[res6b$ea==1,"value",drop=TRUE],a=res6b[res6b$ea==0,"value"])@auc
    #print(z) 
    }
  
  }##end of xx loop  
  )
  
  resX<-data.frame(disease=pres2[i],AUC=mean(ffR),AUC2=mean(ffR2,na.rm=TRUE),ttmean,ttspeed,ttd,ttmeanv,ttspeedv,ttdv,ttspeedh,ttdh,breaks1,breaks2)
  
  if(is.null(res1)){res1<-resX}else{res1<-rbind(res1,resX)}  
  
  print(i)
}


hist(res1$AUC,breaks=50,main="AUC scores of all diseases")
abline(v=0.5,lty=2)
abline(v=0.65,lty=1)


write.csv(res1,file="C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT\\AUC1.csv")


