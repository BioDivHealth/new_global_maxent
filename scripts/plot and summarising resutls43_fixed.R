library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(sinaplot)
library(raster)
library(maptools)
data(wrld_simpl)
library(fasterize)
library(sf)
library(waffle)
library(doParallel)
library(gstat)
#library(raster)
library(RColorBrewer)
library(viridis)

source("C:\\Users\\david.redding\\Dropbox\\R_scripts\\bivariate_plot_functions2.R")
source("C:\\Users\\david.redding\\Dropbox\\R_scripts\\HDPI.R")


##read in empress-i point data
point_data<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\all_empress_i_data4.csv",stringsAsFactors = FALSE)
point_data<-point_data[!is.na(point_data$SumCases),]
#point_data$LU<-paste(" ",point_data$name_LU,sep="")
point_data<-point_data[!point_data$name_LU %in% c("borrelia burgdorferi","babesia microti","borrelia miyamotoi"),]


##change projecion??
#load(file="C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\wrld_simpl3.r")
#wrld_simpl3<-spTransform(wrld_simpl3,CRS=projection(template))
#writeOGR(wrld_simpl3,driver="ESRI Shapefile",dsn="C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\wrld_simpl3.shp","wrld_simpl3")
wrld_simpl3<-readOGR("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\wrld_simpl3.shp","wrld_simpl3")
#wrld_simpl3[wrld_simpl3$NAME %in% c("Russia","France"),]
wrld_simpl3@data[wrld_simpl3$NAME=="Russia","SUBREGION"]<-wrld_simpl3@data[wrld_simpl3$NAME=="Ukraine","SUBREGION"]
#wrld_simpl3[wrld_simpl3$NAME %in% c("Siberia","Mongolia"),]
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","SUBREGION"]<-30
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","REGION"]<-142
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","ISO2"]<-"XX"
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","FIPS"]<-"XX"
wrld_simpl3@data[wrld_simpl3$NAME=="Chile","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Paraguay","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Uruguay","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Argentina","SUBREGION"]<-4

##base raster template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

template2<-template
values(template2)<-0

template1<-aggregate(template,25)

#make area per sqaure
ar1<-mean(values(area(template2)),na.rm=TRUE)

###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl3),template)

ws2b<-fasterize(st_as_sf(wrld_simpl3),template1)

regionsA<-readOGR("C:\\Users\\david.redding\\Dropbox\\ethnologue\\subregions_layer.shp","subregions_layer")
regionsA<-regionsA[regionsA$IMAGE24 !="Antarctica",]


###cell_ids
#ci1<-list.files("C:\\Users\\david.redding\\Documents\\cellids/",pattern="_2_cellids",full.names=TRUE)
#ci2<-list.files("C:\\Users\\david.redding\\Documents\\cellids/",pattern="_2_cellids",full.names=FALSE)
#ci2<-gsub("_2_cellids.csv","",ci2)

ci1<-list.files("D:\\cellids/",pattern="_2_cellids",full.names=TRUE)
ci2<-list.files("D:\\cellids/",pattern="_2_cellids",full.names=FALSE)
ci2<-gsub("_2_cellids.csv","",ci2)


##make wrld simpl mask
mask1<-fasterize(st_as_sf(wrld_simpl[wrld_simpl$NAME!="Antartica",]),template)

### plot and summarising results
### areas that have increased due to climate versus those that increased due to
### land use. Size of dots - burden estimate?
### panels for each year - facet across year and RCP?
### if over/under threshold due climate versus 


#gl1<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",pattern="_1.csv",full.names=TRUE)
#gl2<-list.files("/home/ucbtdw0/Scratch/New_Global_MAXENT/gainlossresults/",pattern="_1.csv",full.names=FALSE)

gl1<-list.files("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\newdata\\gainlossresults2/",pattern="_2.csv",full.names=TRUE)
gl2<-list.files("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\newdata\\gainlossresults2/",pattern="_2.csv",full.names=FALSE)

#gl1<-list.files("C:\\Users\\david.redding\\Documents\\gainlossresults2/",pattern="_2.csv",full.names=TRUE)
#gl2<-list.files("C:\\Users\\david.redding\\Documents\\gainlossresults2/",pattern="_2.csv",full.names=FALSE)
gl2<-gsub("_2.csv","",gl2)
gl2<-gsub("_X_",";",gl2)
gl3<-read.table(text=gl2,sep=";",stringsAsFactors = F)$V1

###get gain loss
for (i in 1:length(gl1)){
  
  ##gainloss
  glF<-tryCatch(fread(gl1[i]), error=function(err2) err2)
  
  if(nrow(glF)>1){
    
      glF2a<-glF[1,]
      glF2a$cases_per_year=sum(glF$cases_per_year,na.rm=TRUE)
      glF2a$CFR.low=weighted.mean(glF$CFR.low,w=glF$cases_per_year,na.rm=TRUE)
      glF2a$CFR.high=weighted.mean(glF$CFR.high,w=glF$cases_per_year,na.rm=TRUE)
      glF2a$countries=paste(glF$countries,collapse = ",")
      glF2a$spillover_rate2=weighted.mean(glF$spillover_rate2,w=glF$spillover_rate2,na.rm=TRUE)
      
      } else {glF2a<-glF}
  
  #if(class(glF)[1]=="simpleError"){  }
  
  if(i==1) {glF2<-glF2a} else {glF2<-rbind(glF2,glF2a)}
  
}

#"rickettsia africae"  "rickettsia honei"                  "rickettsia japonica"               "rickettsia rickettsii"            



#res1<-list.files("D:\\more_results2\\",pattern="PRESENT.csv",full.names=TRUE)
#res2<-list.files("D:\\more_results2\\",pattern="PRESENT.csv",full.names=FALSE)

res1<-list.files("D:\\Users\\xxxx\\Documents\\more_results5\\",pattern=".csv",full.names=TRUE)
res2<-list.files("D:\\Users\\xxxx\\Documents\\more_results5\\",pattern=".csv",full.names=FALSE)
#res2<-gsub("_PRESENT","",res2)
res2<-gsub(".csv","",res2)
res2b<-read.table(text=res2,sep="_",stringsAsFactors = FALSE)

res3<-NULL
for (i in 1:length(res1)){
  res2x<-read.csv(res1[i],stringsAsFactors = F)
  res2x$threshold<-res2b[i,"V5"]
  if(length(names(res2x)[names(res2x) %in% c("disease")])==0){res2x$disease=res2b$V1[i]}
  if(is.null(res3)){res3<-res2x} else {res3<-rbind(res3,res2x)}
}





res3a<-res3[res3$model=="present",]
names(res3a)<-paste0(names(res3a),"PRESENT")
res3a$model<-gsub("present","presentX0_2010XPRESENT",res3a$model)
res3a$mer1<-paste0(res3a$diseasePRESENT,res3a$thresholdPRESENT)
  
res3b<-res3[res3$model!="present",]
res3b2<-read.table(text=res3b$model,sep="X",stringsAsFactors = FALSE)
res3b3<-read.table(text=res3b2$V2,sep="_",stringsAsFactors = FALSE)
names(res3b2)[c(1,3)]<-c("SSP","DISEASE_GROUP")
names(res3b3)<-c("RCP","Year")

res3d<-cbind(res3b,res3b2[c(1,3)],res3b3)

res3d$mer1<-paste0(res3d$disease,res3d$threshold)

##merge with 2010
res3z<-merge(res3d,res3a[,c(4:14,32)],by.x="mer1",by.y="mer1")


##
#res3z<-cbind(res3b[,2:3],(res3c[,4:14])-res3c[,31:41],res3b3,res3b2[c(1,3)],res3c[,15:30])
#hist(res3z$richness[res3z$SSP=="ssp3"],breaks=50)

##merge with gainloss
res3c<-res3z#merge(res3z,glF2[,46:54],by.x="disease",by.y="name2")

#hist(res3c$AUC)


#hist(res3c$richness[res3c$disease=="dengue"])

#hist(res3c$richness[res3c$disease=="rocio"])

#hist(res3c$richness[res3c$disease=="leishmania infantum chagasi"])

#hist(res3c$richness[res3c$disease=="fasciola gigantica"])


##remove disease with very patchy known endemics areas
#res3c<-res3c[!res3c$disease %in% c("angiostrongylus cantonensis","fasciola gigantica","paragonimus kellicotti"),]# res3c$richness>quantile(res3c$richness,0.01) & res3c$richness<quantile(res3c$richness,0.99),]#dcast(res3,MinMax+RCP+Year+disease+pres_size~Type,value.var=c("Gain","change"))#rbind(res4,res3)


#hist(res3c$richness,breaks=50)

#remove extreme values
res5<-res3c#[res3c$richness>quantile(res3c$richness,0.05) & res3c$richness<quantile(res3c$richness,0.95),]#dcast(res3,MinMax+RCP+Year+disease+pres_size~Type,value.var=c("Gain","change"))#rbind(res4,res3)

#hist((res5$richness),breaks=50)


#ebs1<-res3c[res3c$disease=="dengue",]
#boxplot(ebs1$richness~ebs1$Year)

#ebs1<-res5[res5$disease=="dengue",]
#boxplot(ebs1$richness~ebs1$Year)


#ebs1<-res3c[res3c$disease=="ebola",]
#boxplot(ebs1$richness~ebs1$Year)

#ebs1<-res5[res5$disease=="ebola",]
#boxplot(ebs1$richness~ebs1$Year)



##add disease info
d1<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\disease_table33c.csv",stringsAsFactors=FALSE)
#d1$name<-paste(" ",d1$name,sep="")
d1$name<-gsub("angiostrongylus costaricensis ","angiostrongylus costaricensis",d1$name)
d1$spillover_rate2<-cut(d1$cases_per_year,breaks=c(0,0.99,99,99999,999999999999),labels=FALSE)/4
d1$Vectored<-1
d1$Vectored[d1$Type == c("HOST->HUMAN->HUMAN")]<-0
d1$Vectored[d1$Type == c("HOST->HUMAN")]<-0
d1$Vectored[d1$Type == c("HUMAN->VECTOR->HUMAN")]<-2
d1$Vectored[d1$Type == c("HUMAN->VECTOR->HUMAN->HUMAN")]<-2

##outbreak diseases
d1$Outbreak<-0
#d1$Outbreak[d1$Type %in% c("HOST->VECTOR->HUMAN->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN","HOST->HUMAN->HUMAN","HUMAN->VECTOR->HUMAN->HUMAN")]<-1
#d1$Outbreak[d1$Type %in% c("HOST->VECTOR->HUMAN")]<-0

d1$Outbreak[d1$Type == c("HOST->VECTOR->HUMAN->HUMAN")]<-1
d1$Outbreak[d1$Type == c("HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN")]<-1
d1$Outbreak[d1$Type == c("HOST->HUMAN->HUMAN")]<-1
d1$Outbreak[d1$Type == c("HUMAN->VECTOR->HUMAN->HUMAN")]<-1



##subset - change to absolute
d1<-d1[,c(1:20,46:48)]

d1[d1$name=="omsk","countries"]<-"RU"

for (w in 1:nrow(d1)){
  
  d2<-d1[w,]
  
  cc<-strsplit(d2$countries,",")[[1]]
  
  ws2b<-wrld_simpl3[wrld_simpl3$ISO2 %in% cc, ]
  
  r2<-data.table(disease=d2$name,popn=sum(ws2b$POP2005,na.rm=TRUE), ws2b@data,stringsAsFactors = F)
  
  if(w==1){r3<-r2} else {r3<-rbind(r3,r2)}
  
}

r3$dummy=1

regions<-as.data.frame.matrix(xtabs(data=r3,dummy~disease+REGION))
regions[regions>0]<-1
names(regions)<-c("Antarctica","Africa","Oceania","North America","South America","Asia","Europe")
regions$name<-row.names(regions)
regions$range<-paste0(regions$Antarctica,regions$Africa,regions$Oceania,regions$America,regions$Asia,regions$Europe)
#table(regions$range)

##pop at risk
r4<-r3[,c("disease","popn")]
r4<-r4[!duplicated(r3$disease),]

##put known pop at risk
regions<-merge(regions,r4,by.x="name",by.y="disease",all.x=T,all.y=F)

##merge with original data
d1<-merge(d1,regions,by.x="name",by.y="name")

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
  
  #res6$dummy=1
  #xtabs(dummy~page+model,res6)
  #xtabs(dummy~page+disease,res6)
  #table(res6$page)
  #unique(paste(res6$page,res6$disease,sep="_"))
  
  ##make new group for bacteria
  res6$group2[res6$group2=="bacteria"]<-"bartonella"
  res6$group2[res6$disease=="germiston"]<-"orthobunyavirus"
  res6$group2[res6$disease=="guaroa"]<-"orthobunyavirus"
  
   
  #table(res6[,"group2"])
  table(is.na(res6[,"group2"]))
  
  #res6$uni2<-paste(res6$disease,res6$SSP,res6$RCP,sep="_")
  #res6a<-res6[!duplicated(res6$uni2),]
  #res6a$Year<-2010
  #res6a[,3:13]<-0
  #res6<-rbind(res6,res6a)
  
  
  
  ###aggregate by page
  #res7<-aggregate(res6[res6$AUCmax>0.6 & res6$disease.y!="None reported",c("richness","disease_richness","humans","cases","deaths","x","y","max_x","max_y","min_x","min_y","RCP","Year","SSP","CFR.low","CFR.high","cases_per_year","Outbreak","Vectored")],by=list(res6$group2[res6$AUCmax>0.6 &res6$disease.y!="None reported"],res6$RCP[res6$AUCmax>0.6 &res6$disease.y!="None reported"],res6$Year[res6$AUCmax>0.6 &res6$disease.y!="None reported"]),FUN=median)
  #res7$uni<-paste(res7$Group.1,res7$Group.2,sep="_")

  #res6x<-res6[ res6$Vectored==1 & res6$AUCmax>0.6 & res6$group2!="angiostrongylus" & res6$disease.y!="None reported" ,]

  ##sort out missing data
  res6[res6$disease.y=="" & res6$disease=="germiston" ,"disease.y"]<-"Febrile illness with rash (arbocat)"
  res6[res6$disease.y=="" & res6$disease=="guaroa" ,"disease.y"]<-"Febrile illness (arbocat)"
  res6[res6$disease.y=="XXX"& res6$disease=="rocio" ,"disease.y"]<-"Encephalitis (arbocat)"
  res6[is.na(res6$disease.y) & res6$disease=="topografov"   ,"disease.y"] <-"None reported"
  
  
  ##make main kindgoms
  res6$Kingdom=res6$group
  res6$Kingdom[!res6$Kingdom %in% c("bacteria","parasite","fungi")]<-"virus"
  
  ##make main geography
  res6$Tropical="Temperate"
  res6$Tropical[res6$y>(-23.5) & res6$y<(23.5)]<-"Tropical"
  
  ##make main heimsphere
  res6$Hemisphere="North"
  res6$Hemisphere[res6$y<0]<-"South"
  
  ##create quadrant
  res6$Region<-paste(res6$Tropical,res6$Hemisphere,sep=" ")
  res6$Region<-factor(res6$Region,levels=c("Temperate North","Tropical North","Tropical South","Temperate South"))
  
    ##
  ##difference from 2010
  res6$change_x<-res6$x-res6$xPRESENT
  res6$change_y<-res6$y-res6$yPRESENT
  
  res6q2<-res6[,c("richness","disease_richness","humans","cases","poor_cases","x","y","max_x","max_y","min_x","min_y")]-res6[,paste0(c("richness","disease_richness","humans","cases","poor_cases","x","y","max_x","max_y","min_x","min_y"),"PRESENT")]
  names(res6q2)<-paste0(names(res6q2),"diff")
  
  ####add on differece
  res6<-cbind(res6,res6q2)
  
  ###rid na richness
  #res6<-res6[!is.na(res6$richness),]
  
  res6$TempTrop<-"Tropical"
  res6$TempTrop[res6$y<(-23.5) | res6$y>23.5]<-"Temperate"
  
  res6$Pole<-"North"
  res6$Pole[res6$y<0]<-"South"
  
  res6$Region<-paste(res6$TempTrop,res6$Pole,sep=" ")
  
  res6$Region<-factor(res6$Region,levels=c("Temperate North","Tropical North","Tropical South","Temperate South"))
  
  
  res6$Disease<-"Disease"
  
  res6$Disease[res6$disease %in% unique(res6[res6$disease.y=="None reported" ,"disease" ])]<-"No Disease"
  
  res6$disease2<-as.numeric(as.factor(res6$disease))
  
  #####MAPS???
  res6$model<-gsub("0X","0_",res6$model)
  
  #res6[1,]
  
  ##set response to be mapped
  res6$Response=res6$richnessdiff
  #res6$Response=res6$casesdiff
  #res6$Response=res6$poor_casesdiff
  
  #res6$Response[res6$Response<0]<-rnorm(n=nrow(res6[res6$Response<0,]),mean=(-1),sd=0.001)
  #res6$Response[res6$Response>0]<-rnorm(n=nrow(res6[res6$Response>0,]),mean=(1),sd=0.001)
  
  mods<-unique(res6$model)#[1:20]
  
  mods2<-read.table(text=mods,sep="_")
  mods2$mods<-mods
  mods3<-read.table(text=mods2$V1,sep="X")
  mods2<-cbind(mods2,mods3)
  names(mods2)[c(2:3,5:6)]<-c("Year","Threshold","SSP","RCP")
  mods2<-mods2[mods2$Threshold=="ALL",]
  
  
  for(pp in c(2030,2050,2070,2080)){
  
  
  for (jj in c(2.6,4.5,6.0,8.5)){
  
  
  ###choose subsection
  mods4<-mods2[mods2$Year %in% pp & mods2$SSP=="ssp3" & mods2$RCP==jj,]# & mods2$SSP=="ssp2" ,]
  
#  mods4<-mods2[mods2$RCP==6.0,]
  
  ff2=0
  sum1=NULL
  ##loop to create map
  #for (f in 1:250){
  
    #  res6q<-res6[ res6$AUC>0.5 & res6$Disease=="Disease" & res6$threshold=="ALL",]

    
    res6x<-res6[res6$model %in% sample(mods4$mods,1) &  res6$AUC>0.5 & res6$Disease=="Disease",]
  
    #res6x<-res6[res6$Year==2070 & res6$threshold=="ALL" & res6$RCP==4.5 & res6$SSP==1 , ]
    
    ##choose all options of that disease
    resN<-res6x#[res6x$disease %in% resNtemp$disease,]
    
    ##make spatial object
    aq<-resN
    coordinates(aq)<-~x+y
    aq$dummy=1
    
    aq$Response[aq$Response>0]<-1
    aq$Response[aq$Response<0]<-(-1)
    
    #aq$Response<-scale(aq$Response)
    
    r<-raster()
    res(r)<-12.5#sample(c(10,12.5,15,17.5,20),1)

    
    lunq<-function(...) length(unique(...))
    
    ### mean response
    ras2a<-rasterize(aq,r,field="Response",mean)
    #plot(ras2a)
  
    ##confidence
    ras2b<-rasterize(aq,r,field="group",lunq)
    #ras2b<-rasterize(aq,r,field="dummy",sum)
    #plot(ras2b)

    ### divide by confidence
    #ras2<-ras2a/max(values(ras2b),na.rm=TRUE)
    ras2<-ras2a*ras2b
    
    
    ##make NA zero
    ras2[is.na(ras2)]<-0
    #plot(ras2)
    
    #r2<-raster()
    #res(r2)<-15#sample(c(5,10,15,20),1)
    
    
    #ras2<-raster::resample(ras2,r2,method="bilinear")
    
    #plot(ras2)
    
    
    if(is.null(sum1)) {sum1=ras2} else{sum1=sum1+ras2}
    
    ff2=ff2+1
    print(ff2)
    #}
    
  sum2<-sum1/ff2
  #sum2[1,1]<-0.883
  
  sum3<-raster::resample(sum2,template,method="ngb")
  
  sum4<-mask(sum3,ws2)
  
  #plot(sum4)
  
  #sum4a<-focal(sum4,w=matrix(1/25,5,5),na.rm=TRUE)
  #sum4a<-focal(sum4a,w=matrix(1/25,5,5),na.rm=TRUE)
  #sum4a<-focal(sum4a,w=matrix(1/25,5,5),na.rm=TRUE)
  
  #sum4<-focal(sum4 ,w=matrix(1/9,3,3))
  #sum4<-focal(sum4 ,w=matrix(1/9,3,3))
  
  #pr1<-"+proj=moll"
  #sum4<-sum4a#projectRaster(from=sum4a,crs = pr1)
  #sum4[1,1]<-max(values(sum4),na.rm=T)*-1
  #sum4[1,2]<-min(values(sum4),na.rm=T)*-1
  
  breakpoints <- c(seq(from=-1*(max(abs(values(sum4)),na.rm=TRUE)),to=1*(max(abs(values(sum4)),na.rm=TRUE)),by=c(1*(max(abs(values(sum4)),na.rm=TRUE))-(-1*(max(abs(values(sum4)),na.rm=TRUE))))/16))
  #breakpoints <- c(seq(from=-1,to=1,by=0.125))
  rr<-8#(length(breakpoints)-1)/2
  colors <- c(rev(RColorBrewer::brewer.pal(rr, "Blues")),rep("White",1),(RColorBrewer::brewer.pal(rr, "Oranges")))
  
  
  png(paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\results_maps\\20maps_",paste(c(paste(unique(mods4$RCP,collapse="_")),paste(unique(mods4$Year,collapse="_")),paste(unique(mods4$SSP,collapse="_"))),collapse="_"),".png",sep=""),width=1900,height=900)
  plot(sum4, breaks = breakpoints, col = colors,colNA="light grey",main=paste(c(paste(unique(mods4$RCP,collapse=" ")),paste(unique(mods4$Year,collapse=" ")),paste(unique(mods4$SSP,collapse=" "))),collapse=" "))
  #sum2<-focal(sum2, w=matrix(1/9, nc=3, nr=3))
  dev.off()
  #sum3<-as.data.frame(sum2,xy=TRUE)
  #sum3<-sum3[!is.na(sum3$layer),]
  
  #sum3$layer[sum3$layer<(-0.05)]<-(-0.05)
  #sum3$layer[sum3$layer>(0.05)]<-0.05
  
  } ## end of jj
    
  } ## end of pp

  #plot(sum2)#,col=viridis(20,option = "E"))  
  
  ##auc RESUILTS
  nrow(res6[res6$AUC>=0.5,])/ nrow(res6)
  nrow(res6[res6$AUC>=0.6,])/ nrow(res6)
  nrow(res6[res6$AUC>=0.7,])/ nrow(res6)
  
  res6[1,]
  
  
  ####choose only good
  res6q<-res6[ res6$AUC>0.5 & res6$Disease=="Disease" & res6$threshold=="ALL",]
  
  ##
  
  ##overal mean difference
  zz<-aggregate(log(res6q[,c(5:15,36:46)]),by=list(res6q$Year),median)#function (x) mean(x,trim=0.1,na.rm=TRUE))
  ((exp(zz$richness[zz$Group.1==2080 ]) -exp(zz$richness[zz$Group.1==2030 ]))*ar1)/50
  
  zz<-aggregate(log(res6q[,c(5:15,36:46)]),by=list(res6q$Year),HPDI,prob=0.99)#function (x) mean(x,trim=0.1,na.rm=TRUE))
  ((exp(zz$richness [zz$Group.1==2080 ]) -exp(zz$richness[zz$Group.1==2030 ]))*ar1)/50
  
  
  
  zz<-aggregate(log(res6q[,c(5:15,36:46)]),by=list(res6q$Year,res6q$TempTrop),median)
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="Temperate"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="Temperate"]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="Tropical"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="Tropical"]))*ar1)/50
  
  
  zz<-aggregate(log(res6q[,c(5:15,36:46)]),by=list(res6q$Year,res6q$TempTrop),HPDI,prob=0.95)
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="Temperate"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="Temperate"]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="Tropical"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="Tropical"]))*ar1)/50
  
  
  # 0 host only, 1 vector and host, 2 vector only
  zz<-aggregate(log(res6q[,c(5:15,36:46)]),by=list(res6q$Year,res6q$Vectored),median)
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2==0]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2==0]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2==1]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2==1]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2==2]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2==2]))*ar1)/50
  

  zz<-aggregate(log(res6q[,c(5:15,36:46)]),by=list(res6q$Year,res6q$Vectored),HPDI,prob=0.95)
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2==0]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2==0]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2==1]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2==1]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2==2]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2==2]))*ar1)/50
  
  
  # 0 host only, 1 vector and host, 2 vector only
  zz<-aggregate(log(res6q[,c(5:15,36:46)]),by=list(res6q$Year,res6q$Kingdom),median)
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="bacteria"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="bacteria"]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="parasite"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="parasite"]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="virus"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="virus"]))*ar1)/50
  
  
  zz<-aggregate(log(res6q[,c(5:15,36:46)]),by=list(res6q$Year,res6q$Kingdom),HPDI,prob=0.95)
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="bacteria"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="bacteria"]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="parasite"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="parasite"]))*ar1)/50
  ((exp(zz$richness[zz$Group.1==2080 & zz$Group.2=="virus"]) -exp(zz$richness[zz$Group.1==2030 & zz$Group.2=="virus"]))*ar1)/50
  
  
  ggplot(res6q,aes(x=log(richness),col=as.factor(Year)))+
    geom_density(adjust=2)+
    facet_wrap(.~Vectored)+
    theme_cowplot()
  
  
  ggplot(res6q[res6q$Kingdom!="fungi", ],aes(x=(richness/richnessPRESENT),col=as.factor(Year)))+
    geom_density(adjust=2)+
    facet_wrap(.~Kingdom,scales="free")+
    theme_cowplot()+
    geom_vline(xintercept = 1,lty=2)
  
  
  ggplot(res6q[res6q$Kingdom!="fungi", ],aes(x=(richness/richnessPRESENT),fill=as.factor(Year)))+
    geom_density(adjust=2,alpha=0.5)+
    facet_grid(Vectored~Kingdom,scales="free")+
    theme_cowplot()+
    geom_vline(xintercept = 1,lty=2)
  
  
  xt<-xtabs(data=res6q[!duplicated(res6q$disease) & res6q$Kingdom!="fungi",],dummy~Kingdom+Vectored)
  chi1<-chisq.test(xt)
  chi1$observed
  xt/chi1$expected
  chi1$stdres
  
  aggregate(res6q[,c(5:15,36:46,82:94,)],by=list(res6q$Year,res6q$RCP,res6q$SSP),median,na.rm=TRUE)

  ##make subset
  res7<-res6[ res6$AUC>0.5 & res6$RCP %in% c(6.0) & res6$Year %in% c(2070) & res6$threshold=="ALL" & res6$SSP==3,]
  
  ##make mean by panel
  res8<-aggregate(res7[,c("change_x","change_y")],by=list(res7$Region),mean,na.rm=TRUE)
  res8sd<-aggregate(res7[,c("change_x","change_y")],by=list(res7$Region),sd,na.rm=TRUE)
  res8n<-aggregate(res7[,c("change_x","change_y")],by=list(res7$Region),function(x) sum( !is.na(x) ) )
  res8sem<-res8sd[,c("change_x","change_y")]/sqrt(res8n[,c("change_x","change_y")])
  res8upper<-res8[,c("change_x","change_y")]+res8sem
  res8lower<-res8[,c("change_x","change_y")]-res8sem
  res9<-cbind(res8,res8upper,res8lower)
  names(res9)<-c("Region","Mean_x","Mean_y","Max_x","Max_y","Min_x","Min_y")
  
  ##make plot
  
  #positive plot
  res7$posy<-res7$dummy
  res7$posy[round(res7$change_y,2)<=0]=0
  
  
  xt<-xtabs(data=res7,dummy~Region+posy)
  chi1<-chisq.test(xt)
  chi1$observed
  xt/chi1$expected
  chi1$stdres
  
  #Region              0         1
  #Temperate North 0.4477124 1.1823085
  #Tropical North  1.6117647 0.7980583
  #Tropical South  1.2820856 0.9068844
  #Temperate South 2.1696833 0.6138910

  #X-squared = 17.405, df = 3, p-value = 0.0005834

  
  se <- function(x) sqrt(var(x) / length(x))
  
  # 0 host only, 1 vector and host, 2 vector only
  zx<-aggregate(res7$change_y,by=list(res7$Region),mean)
  
  
  zx2<-aggregate(res7$change_y,by=list(res7$Region),se)
  
  zx$Higher=zx$x+zx2$x
  zx$Lower=zx$x-zx2$x
  
  zx*111
  
  #       Group.1            x     Higher      Lower
  #1 Temperate North  0.870944137  1.0034632  0.7384251
  #2  Tropical North  0.145473105  0.4409287 -0.1499825
  #3  Tropical South -0.004403239  0.1491136 -0.1579201
  #4 Temperate South -0.650990616 -0.2851792 -1.0168020
  
  
  
  
  ##Figure 2
  gg <- ggplot(data=res7)
  gg <- gg + geom_point( aes(x=change_x, y=change_y, group=disease,shape=as.factor(Vectored),fill=as.factor(Kingdom)),size=2)
  #gg <- gg + geom_point( aes(x=(Longitude), y=change_max_y),fill="seagreen",size=2,pch=21)
  #gg <- gg + geom_point( aes(x=(Longitude), y=change_min_y),fill="orange",size=2,pch=21)
  #gg <- gg + geom_line( aes(x=(change_x), y=change_y, group=disease),lty=2)
  gg <- gg+theme_cowplot()
  gg <- gg + geom_hline(yintercept = 0)
  gg <- gg + geom_vline(xintercept = 0)
  gg <- gg + xlim(-5,5)
  gg <- gg + ylim(-5.5,5.5)
  gg <- gg + xlab("Longitude")
  gg <- gg + ylab("Change in latitude")
  #gg <- gg + scale_fill_distiller(palette = "YlOrRd")
  gg <- gg + facet_wrap(~Region,ncol=2,scales="free")
  gg <- gg + theme(legend.position="right")
  gg <- gg + scale_shape_manual(values=c(21,24,22))
  
  gg
  #set_palette(gg,"npg")
  
  ##make map
  ###PRESENT DATA
  
  
  res3$model<-gsub("present","presentX0_2010XPRESENT",res3$model)
  res3b2<-read.table(text=res3$model,sep="X",stringsAsFactors = FALSE)
  res3b3<-read.table(text=res3b2$V2,sep="_",stringsAsFactors = FALSE)
  names(res3b2)[c(1,3)]<-c("SSP","DISEASE_GROUP")
  names(res3b3)<-c("RCP","Year")
  
  res3e<-cbind(res3,res3b2[c(1,3)],res3b3)
  
  ###make another unique
  res3e$uni9<-paste0(res3e$disease,res3e$RCP,res3e$Year)
  
  ##filter only different coordinates
  dat<-res3e[!duplicated(res3e$uni9),]
  
  
  library(ggplot2)  # FYI you need v2.0
  library(dplyr)    # yes, i could have not done this and just used 'subset' instead of 'filter'
  library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
  library(ggthemes) # theme_map and tableau colors
  
  world <- map_data("world")
  world <- world[world$region != "Antarctica",] # intercourse antarctica
  
  #dat <- res3# having factors here by default isn't a bad thing
  
  gg <- ggplot()
  gg <- gg + geom_map(data=world, map=world,
                      aes(x=long, y=lat, map_id=region),
                      color="white", fill="#7f7f7f", size=0.05, alpha=1/4)
  gg <- gg + geom_line(data=dat[dat$AUC>0.5 & dat$RCP %in% c(6.0,0.0) & dat$Year %in% c(2010,2070),], 
                       aes(x=x, y=y, group=disease), 
                       size=0.1, alpha=1)
  gg <- gg + geom_point(data=dat[dat$AUC>0.5 & dat$RCP %in% c(6.0,0.0) & dat$Year %in% c(2010,2070),], 
                        aes(x=x, y=y, group=disease,col=as.factor(Year)), 
                        size=0.1, alpha=1)
  gg <- gg + scale_color_tableau()
  #gg <- gg + coord_proj("+proj=wintri")
  #gg <- gg + facet_wrap(~RCP)
  gg <- gg + theme_map()
  gg <- gg + theme(strip.background=element_blank())
  gg <- gg + theme(legend.position="bottom")
  gg <- gg + geom_hline(yintercept = 0,lwd=0.1)
  gg <- gg + geom_hline(yintercept = 23.5,linetype = "dashed",lwd=0.1)
  gg <- gg + geom_hline(yintercept = -23.5,linetype = "dashed",lwd=0.1)
  gg <- gg +coord_sf(xlim = c(-165, 175), ylim = c(-55, 85), expand = FALSE)
  gg
  
  ###
  write.csv(dat,file="C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\north_south.csv")
  
  dat2<-dat[dat$AUC>0.5 & dat$RCP %in% c(4.5,0.0) & dat$Year %in% c(2010,2070),]
  
  dat2$TempTrop<-"Tropical"
  dat2$TempTrop[dat2$y<(-23.5) | dat2$y>23.5]<-"Temperate"
  
  dat2$Pole<-"North"
  dat2$Pole[dat2$y<0]<-"South"
  
  dat2$Region<-paste(dat2$TempTrop,dat2$Pole,sep=" ")
  
  dat2$Region<-factor(dat2$Region,levels=c("Temperate North","Tropical North","Tropical South","Temperate South"))
  
  datP<-dat2[dat2$Year==2010,]
  datF<-dat2[dat2$Year!=2010,]
  
  datPF<-merge(datF,datP,by="disease",all.x=T,all.y=F)
  
  ####make a decions about poleward
  datPF$change_x<-datPF$x.x-datPF$x.y
  datPF$change_y<-(datPF$y.x+90)-(datPF$y.y+90)
  datPF$change_max_y<-(datPF$max_y.x+90)-(datPF$max_y.y+90)
  datPF$change_min_y<-(datPF$min_y.x+90)-(datPF$min_y.y+90)
  datPF$change_max_x<-datPF$max_x.x-datPF$max_x.y
  datPF$change_min_x<-datPF$min_x.x-datPF$min_x.y
  
  #datPF$change_max_y[datPF$change_max_y<(-4.5)]<-NA
  
  datPF$Longitude<-datPF$x.x
  
  datPF$Disease<-"Disease"
  
  datPF$Disease[datPF$disease %in% unique(res6[res6$disease.y=="None reported" ,"disease" ])]<-"No Disease"
  
 datPF$disease2<-as.numeric(as.factor(datPF$disease))

 
 #positive plot
 datPF$posy<-"max_increase"
 datPF$posy[datPF$change_max_y<0]="max_decrease"
 
 
 datPF$miny<-"min_increase"
 datPF$miny[datPF$change_min_y<0]="min_decrease"
 
 
 datPF$poxy<-paste(datPF$posy,datPF$miny,sep="-")
 
 
 xt<-xtabs(data=datPF[datPF$Disease=="Disease",],dummy.x~Pole.x+poxy)
 chi1<-chisq.test(xt)
 chi1$observed
 xt/chi1$expected
 chi1$stdres
 
 
 ####no change
 
 aggregate(datPF[,c("change_max_y","change_min_y")],by=list(datPF$Region.x),median)
 
 
 
 #positive plot
 datPF$posy<-"increase"
 datPF$posy[round(datPF$change_max_y,10)<0]="decrease"
 #datPF$posy[round(datPF$change_max_y,1)>0]="increase"
 
 
 datPF$miny<-"increase"
 datPF$miny[round(datPF$change_min_y,10)<0]="decrease"
 #datPF$miny[round(datPF$change_min_y,1)>0]="increase"
 
 
 datPF$poxy<-paste(datPF$posy,datPF$miny,sep="-")
 
 
 
 ###how are expansions  experienced - shifts are static. Of the expansions what drives the expansion poleward or equatorial?
 
 
 #xt<-xtabs(data=datPF[,],
 chi1<-chisq.test(xt)
 chi1$observed
 xt/chi1$expected
 chi1$stdres
 
 
 ##change in max
 summary(lm(data=datPF[datPF$Disease=="Disease",],change_max_y~Region.x))
 aggregate(data=datPF[datPF$Disease=="Disease",],change_max_y~Region.x,mean)
 
 
 ##change in min
 summary(lm(data=datPF[datPF$Disease=="Disease",],change_min_y~Region.x))
 aggregate(data=datPF[datPF$Disease=="Disease",],change_min_y~Region.x,mean)
 
 #Region              0         1
 #Temperate North 0.4477124 1.1823085
 #Tropical North  1.6117647 0.7980583
 #Tropical South  1.2820856 0.9068844
 #Temperate South 2.1696833 0.6138910
 
 #X-squared = 17.405, df = 3, p-value = 0.0005834
 
 
 se <- function(x) sqrt(var(na.omit(x)) / length(na.omit(x)))
 
 # 0 host only, 1 vector and host, 2 vector only
 zx<-aggregate(datPF$change_max_y,by=list(datPF$Region.x),mean,na.rm=TRUE)
 
 zx<-aggregate(datPF$change_max_y,by=list(datPF$Region.x),se)

 zx$Higher=zx$x+zx2$x
 zx$Lower=zx$x-zx2$x
 zx
 zx*111
 
 # 0 host only, 1 vector and host, 2 vector only
 zx<-aggregate(datPF$change_min_y,by=list(datPF$Region.x),mean)
 
 zx<-aggregate(datPF$change_min_y,by=list(datPF$Region.x),se)
 
 zx$Higher=zx$x+zx2$x
 zx$Lower=zx$x-zx2$x
 zx
 zx*111
 
  ##make plot
 ## Figure 3
  gg <- ggplot(data=datPF[datPF$Disease=="Disease",])
 # gg <- gg + geom_point( aes(x=(change_x), y=change_y, group=disease,fill=Longitude), #
 #                      size=2,pch=21)
  gg <- gg + geom_point( aes(x=(Longitude), y=change_max_y),fill="seagreen",size=2,pch=21)
  gg <- gg + geom_point( aes(x=(Longitude), y=change_min_y),fill="orange",size=2,pch=21)
  #gg <- gg + geom_line( aes(x=(change_x), y=change_y, group=disease),lty=2)
  gg <- gg+theme_cowplot()
  gg <- gg + geom_hline(yintercept = 0)
  #gg <- gg + geom_vline(xintercept = 0)
  #gg <- gg + xlim(-5,5)
  #gg <- gg + ylim(-4.5,4.5)
  gg <- gg + xlab("Longitude")
  gg <- gg + ylab("Change in latitude")
  #gg <- gg + scale_fill_distiller(palette = "YlOrRd")
  gg <- gg + facet_wrap(~Region.x,ncol=1,scales="free_y")
  gg <- gg + theme(legend.position="right")
  #gg <- gg + scale_shape_manual(values=c(21,22))
  
  gg
  #set_palette(gg,"npg")
  
  
  length(unique(datPF[datPF$Disease=="Disease","disease"]))
  
  
  
 ##### Run resamples to check for disease impacts
#### control for geog and tax by ordering rarer combinations first
  
 
 ###what sensitivity groups 
  types<-c("Baseline","Accuracy","Regional","Vectored","RCP","SSP","Burden","Emerging","Hantaviridae","Kingdom","Tropical/Temperate","RCP/SSP","Threshold")
  
  res6z<-NULL
  for(tt in 1:length(types)){
  

    dones=0
    #ff=1
    ###run sensitivity
    for(ff in 1:100){
  
  #baseline
    #if(tt==1){
      samp1<-0.65#sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
      samp2<-NA#sample(c("Africa","Oceania","North America","South America","Asia","Europe"),6)
      samp3<-c(0,1)#sample(list(0,1,c(0,1)),1)[[1]]
      samp4<-c(2.6,4.5,6.0)#sample(c(2.6,4.5,6.0),1)
      samp5<-c(1:5)#sample(1:5,1)
      samp6<-c(0.25,0.5,0.75,1)#sample(c(0.25,0.5,0.75,1),1)
      samp7<-"None reported"
      samp8<-"XXX"
      samp9<-c("bacteria","parasite","fungi","virus")
      samp10<-c("Temperate","Tropical")
      samp11a<-c(2.6,4.5,6.0)#sample(c(2.6,4.5,6.0),1)
      samp11b<-c(1:5)#sample(1:5,1)
      samp12<-"PRESENT"#]c("MIN","MAX","PRESENT")
   # }
      rr=" "

  #test across accuracy    
    if(tt==2){
      
      samp1<-sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
      rr=samp1
      
    }
      
    #test across regions
    if(tt==3){
      
       samp2<-sort(sample(c("Africa","Oceania","North America","South America","Asia","Europe"),5))
       count1<-c("Africa","Oceania","North America","South America","Asia","Europe")[!c("Africa","Oceania","North America","South America","Asia","Europe") %in% samp2 ]
       rr=count1
    }
    
    #test across vector non vector
    if(tt==4){
        samp3<-sample(list(0,1),1)[[1]]
        rr=samp3
          }
    
      #test across RCP  
      if(tt==5){
          samp4<-sample(c(2.6,4.5,6.0),1)
          rr=samp4
                 }

      
      #test across SSP -cases only 
      if(tt==6){
          samp5<-sample(1:5,1)
          rr=samp5
      }
      
      
      #test across current case numbers
      if(tt==7){
          samp6<-sample(c(0.25,0.5,0.75,1),1)
          rr=samp6
      }
      
      #test non disease
      if(tt==8){
          samp7<-"XXX"
          rr=" "
      }  
        
      #test non disease
      if(tt==9){
        samp8<-"hantavirus"
        rr=" "
      }  
      
      #test non disease
      if(tt==10){
        samp9<-sample(c("bacteria","parasite","fungi","virus"),1)
        rr=samp9
      }
      
      #test non disease
      if(tt==11){
        samp10<-sample(c("Temperate","Tropical"),1)
        rr=samp10
      }
      
      
      #test non disease
      if(tt==12){
        
        samp11a<-sample(c(2.6,4.5,6.0),1)
        samp11b<-sample(1:5,1)
        
        rr=paste(samp11a,samp11b,sep="_")
      }
      
      #test non disease
      if(tt==13){
        samp12<-sample(c("MAX","PRESENT"),1)
        rr=samp12
      }
      
      
   # }else{
    
    #sample for sensitivity
   # samp1<-sample(c(0.55,0.6,0.65,0.7,0.75),1)
   # samp2<-sample(c("Africa","Oceania","North America","South America","Asia","Europe"),5)
   # samp3<-sample(list(0,1,c(0,1)),1)[[1]]
    
   # }  
  
    count1<-c("Africa","Oceania","North America","South America","Asia","Europe")[!c("Africa","Oceania","North America","South America","Asia","Europe") %in% samp2 ]
    
    res6$count1=rowSums(cbind(res6[,names(res6) %in% count1],res6[,names(res6) %in% count1]))
    
  ##remove 8.5 - opposite results
  res6x<-res6[res6$AUC>samp1 & res6$count1!=0 & res6$Vectored %in% samp3 & res6$RCP %in% samp4 & res6$SSP %in% samp5 & res6$spillover_rate2 %in% samp6 & res6$disease.y!=samp7 & res6$group2!=samp8 & res6$Kingdom %in% samp9 & res6$Tropical %in% samp10 & res6$RCP %in% samp11a & res6$SSP %in% samp11b & res6$threshold %in% samp12,]#  & 
  
  ##set data type
  res6x$results_type=types[tt]
  ##incase two values
  res6x$ff=rr
  
  ##subsample
  #if(ff==1){
  #  res6x$var1=rep("PRESENT",nrow(res6x))

   #     } else {
      
    res6x$var1=paste(paste(samp1,collapse=":"),paste(count1,collapse=":"),paste(samp3,collapse=":"),paste(samp4,collapse=":"),paste(samp5,collapse=":"),paste(samp6,collapse=":"),paste(samp7,collapse=":"),paste(samp8,collapse=":"),paste(samp9,collapse=":"),paste(samp10,collapse=":"),paste(samp11a,collapse=":"),paste(samp11b,collapse=":"),paste(samp12,collapse=":"),sep="_")
      
  #  }
  
  if(res6x$var1[1] %in% dones) {next}
    
  dones<-c(dones,paste(paste(samp1,collapse=":"),paste(count1,collapse=":"),paste(samp3,collapse=":"),paste(samp4,collapse=":"),paste(samp5,collapse=":"),paste(samp6,collapse=":"),paste(samp7,collapse=":"),paste(samp8,collapse=":"),paste(samp9,collapse=":"),paste(samp10,collapse=":"),paste(samp11a,collapse=":"),paste(samp11b,collapse=":"),paste(samp12,collapse=":"),sep="_"))

  e <- simpleError(1)
  
  
  res6y<-setDT(res6x)
  
  res6y[,cor1:=cor(richness,Year,method="spearman")[1]
,by=disease]
  
  res6y[,cor2:=cor(cases,Year,method="spearman")[1]
        ,by=disease]
  
  res6y[,cor3:=cor(poor_cases,Year,method="spearman")[1]
        ,by=disease]
  
  res6y[,lm1:=tryCatch(lm(scale(richness)~Year)$coefficients[2], error = function(e) e)[1]
        ,by=disease]
  
  res6y[,lm2:=tryCatch(lm(scale(cases)~Year)$coefficients[2], error = function(e) e)[1]
        ,by=disease]
  
  res6y[,lm3:=tryCatch(lm(scale(poor_cases)~Year)$coefficients[2], error = function(e) e)[1]
        ,by=disease]
  
  #res6y[complete.cases(res6y),lm1_fitted:=lm(scale(richness)~Year)$fitted.values,by=disease]
  
  #res6y[complete.cases(res6y),lm2_fitted:=lm(scale(cases)~Year)$fitted.values,by=disease]
  
  #res6y[complete.cases(res6y),lm3_fitted:=lm(scale(poor_cases)~Year)$fitted.values,by=disease]
  

  res6y[, tt:=tt]
  
  if(is.null(res6z)){res6z<-res6y; dones<-res6y$var1[1]} else {res6z<-rbindlist(list(res6z,res6y)); dones<-c(dones,res6y$var1[1])}

  
  print(ff)
  
    }

  }
  
  ###recreate unique identified
  res6z[,disease2:=paste0(RCP,SSP,disease,threshold)]
  
  
  ###Extended data figure 3 
  r6x<-res6z[res6z$results_type=="RCP/SSP",]
  setkey(r6x,disease2)
  
  r6z2<-res6z[res6z$results_type=="RCP/SSP" & res6z$Year==2030,]
  
  r6z2[,names(r6z2)[!names(r6z2) %in% c("disease2","lm1_fitted","lm2_fitted","lm3_fitted")]:=NULL]
  names(r6z2)[1]<-"X2030lm1"
  names(r6z2)[2]<-"X2030lm2"
  names(r6z2)[3]<-"X2030lm3"
  setkey(r6z2,disease2)
  
  r6x<-r6x[r6z2]
  
  r6x[,lm1_fitted2:=lm1_fitted-X2030lm1]
  r6x[,lm2_fitted2:=lm2_fitted-X2030lm2]
  r6x[,lm3_fitted2:=lm3_fitted-X2030lm3]
  
  
  #temp1<-r6x[r6x$lm1<0,]
  #temp1<-temp1[temp1$Year==2080 & temp1$lm1_fitted2<(-1),]
  #temp1<-temp1[temp1$lm1>(-0.035) & temp1$lm1<(-0.03) ]
  #temp1<-temp1[temp1$lm1>(-0.0325) & temp1$lm1<(-0.032) ]
  
  #length(unique(r6x$disease))
  
#  hist(temp1$lm1)
  
  ##results
 # r6x[,mean(lm1*abs(cor1)),by=.(Year,Kingdom)]
  
  r6x[,RCP:=paste("RCP ",RCP)]
  r6x[,SSP:=paste("SSP ",SSP)]
  
  ##remove very high mitigation cost scenarios
  #https://unfccc.int/sites/default/files/part1_iiasa_rogelj_ssp_poster.pdf
  
  
  p<-ggplot(r6x[ !r6x$ff %in% c("2.6_3","2.6_5","2.6_2") & r6x$threshold=="PRESENT",],aes(y=lm2_fitted2,x=jitter(Year),col=SSP)) +
  #p<-ggplot(r6x[ !r6x$ff %in% c("2.6_3","2.6_5","2.6_2") & r6x$threshold=="PRESENT",],aes(y=lm3_fitted2*abs(cor3),x=Year,col=SSP)) +
    #geom_histogram()+
    #geom_density(adjust=2.5)+
    #geom_density(data=res6z[res6z$var1=="PRESENT",],aes(x=cor1),lwd=1.5,col="black",adjust=2.5)+
    #res7$uni<-paste(res7$SSP,res7$rep)
    #ggplot(res7,aes(x=jitter(Year),y=richness,col=as.factor(SSP),group=uni)) +
    #stat_summary(fun =  mean, geom = "point",size=2,position = position_dodge(3)) + 
    #stat_summary(fun =  median, geom = "point",size=2,shape=1,position = position_dodge(3)) + 
    #stat_summary(fun =  function(x) quantile (x,0.25), geom = "line",size=1,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  function(x) quantile (x,0.75), geom = "line",size=1,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  median, geom = "point",size=2,position = position_dodge(3)) +
    #stat_summary(fun.data=stat_HPDI, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_point()+
    stat_smooth()+
    #ylim(-1,1)+
    #scale_y_continuous( expand = c(0, 0)) +
    #scale_x_continuous( expand = c(0, 0)) +
    #geom_line(stat="smooth",method = "lm",col="seagreen",
    #          size = 1.5,
    #          alpha = 0.1)+    #geom_point(size=3,alpha=0.5) +
    #guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    #geom_smooth(data=res7,aes(x=jitter(Year),y=richness,group=uni2),se=F,method="lm",col="black")+
    #geom_violin()+
    #scale_x_log10()+
  #xlim(2010,2080)+
  #ylim(-0.5,0.5)+
  theme(legend.position = "none")+
    theme_cowplot(16)+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("Year")+
    ylab("Mean Trend")+
    #geom_abline(intercept=0,slope=1,lty=2)+
    geom_hline(yintercept = 0,lty=2)+
    #xlim(2,9)+
    theme(axis.text.x = element_text(angle = 45,hjust=1))+
    #facet_grid(Vectored~RCP)#+
    #facet_wrap(.~Kingdom,scales="free")+
    theme(legend.position = "none")+
    facet_grid(RCP~SSP)
  
    set_palette(p,"aaas")
  


  res6a<-res6z[!duplicated(paste(res6z$results_type,res6z$ff,res6z$disease,sep="_")),]
  #res6a<-rbind(res6z,res6a)
  
  
  res6a[,cor_cat:="Expanding"]
  res6a$cor_cat[res6a$cor1<0.25]="Static"
  res6a$cor_cat[res6a$cor1<(-0.25)]="Contracting"
  
  res6a$Group<-paste(res6a$results_type,res6a$ff,sep=" ")
  
  res6a$Group<-gsub("Tropical/Temperate ","Region ",res6a$Group)
  res6a$Group<-gsub(" PRESENT"," Mean",res6a$Group)
  res6a$Group<-gsub(" MAX"," Max",res6a$Group)
  
  res6a$Group<-factor(res6a$Group,levels=c("Baseline  ","Emerging  ","Hantaviridae ",sort(unique(res6a$Group)[!unique(res6a$Group) %in% c("Baseline  ","Emerging  ","No Hantaviridae ")])))
  
  
  aggregate(res6[,c("richness")],by=list(res6$model),mean,na.rm=TRUE)
  aggregate(res6[,c("richness")],by=list(res6$model),sd,na.rm=TRUE)
  
  ##results
  res6a[,mean(lm1*abs(cor1)),by=Group]
  
  ggplot(res6a[!res6a$results_type %in% c("SSP","Regional","RCP/SSP") & res6a$threshold!="MIN" & res6a$Group!="Kingdom fungi",],aes(y=lm1*abs(cor1),x=Group,group=Group,col=results_type)) +
   # ggplot(res6a[,],aes(y=lm2*abs(cor2),x=Group,group=Group,col=results_type)) +
    #ggplot(res6a[!res6a$results_type %in% c("SSP","Regional","RCP") & res6a$threshold!="MIN" & res6a$Group!="Kingdom fungi",],aes(y=lm2*abs(cor2),x=Group,group=Group,col=results_type)) +
    
    #geom_histogram()+
    #geom_density(adjust=2.5)+
    #geom_density(data=res6z[res6z$var1=="PRESENT",],aes(x=cor1),lwd=1.5,col="black",adjust=2.5)+
    #res7$uni<-paste(res7$SSP,res7$rep)
    #ggplot(res7,aes(x=jitter(Year),y=richness,col=as.factor(SSP),group=uni)) +
    stat_summary(fun =  mean,na.rm=TRUE, geom = "point",size=2,position = position_dodge(3)) + 
    #stat_summary(fun =  median, geom = "point",size=2,shape=1,position = position_dodge(3),col="black") + 
    #stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  median, geom = "point",size=2,position = position_dodge(3)) +
    
    #stat_summary(fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3))+
    stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_point()+
    #ylim(-1,1)+
    #scale_y_continuous( expand = c(0, 0)) +
    #scale_x_continuous( expand = c(0, 0)) +
    #geom_line(stat="smooth",method = "lm",col="seagreen",
    #          size = 1.5,
    #          alpha = 0.1)+    #geom_point(size=3,alpha=0.5) +
    #guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    #geom_smooth(data=res7,aes(x=jitter(Year),y=richness,group=uni2),se=F,method="lm",col="black")+
    #geom_violin()+
    #scale_x_log10()+
    #xlim(2010,2080)+
    #ylim(0,20500)+
    theme(legend.position = "none")+
    theme_cowplot(14)+
    #geom_line(aes(group =  disease ),color="grey")+
    #xlab("Group")+
    ylab("Mean Trend")+
    #geom_abline(intercept=0,slope=1,lty=2)+
    geom_hline(yintercept = 0,lty=2)+
    #xlim(2,9)+
    theme(axis.text.x = element_text(angle = 45,hjust=1))+
    #facet_grid(Vectored~RCP)#+
  #facet_wrap(.~threshold,scales="free")+
  theme(legend.position = "none")
  #facet_grid(Vectored~SSP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")
  
  
  ############################################################### <------------figure 1
  
  ### figure 1 - general relation
  #Trend over years - with min mean max - accompanying corplot showing missing diseases at negative trajectories - main message is positive but with lots of variation - different subsets of size 10, 20, 30, 40
  #Need to show variation and general relationship
  
  res6$KingTrans<-paste(res6$Kingdom,res6$Vectored,sep="_")
  kngs<-unique(res6$KingTrans)
  
  ###what sensitivity groups 
  types<-c("Baseline","Accuracy","Regional","Vectored","RCP","SSP","Burden","Emerging","Hantaviridae","Kingdom","Tropical/Temperate","RCP/SSP","Threshold","KingTrans")
  
  typesX<-c("Baseline","Accuracy","Regional","Vectored","RCP","SSP","Burden","Emerging","Hantaviridae","Kingdom","Tropical/Temperate","RCP/SSP","Threshold","KingTrans")
  #"Baseline","Vectored","Emerging","Kingdom","Tropical/Temperate","KingTrans")
  
  
  res7<-NULL
  res7b<-NULL
  
  res6z<-NULL
  for(tt in 1:length(types)){
  
    if(!types[tt] %in% typesX){next}
    
    dones=0
    #ff=1
    ###run sensitivity
    for(ff in 1:100){
      
      #baseline
      #if(tt==1){
      samp1<-0.5#sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
      samp2<-"XXX"#c("Africa","Oceania","North America","South America","Asia","Europe")#,6)
      samp3<-c(0,1,2)#sample(list(0,1,c(0,1)),1)[[1]]
      samp4<-c(2.6,4.5,6.0)#sample(c(2.6,4.5,6.0),1)
      samp5<-c(1:5)#sample(1:5,1)
      samp6<-c(0.25,0.5,0.75,1)#sample(c(0.25,0.5,0.75,1),1)
      samp7<-"None reported"
      samp8<-"XXX"
      samp9<-c("bacteria","parasite","fungi","virus")
      samp10<-c("Temperate","Tropical")
      samp11a<-c(2.6,4.5,6.0)#sample(c(2.6,4.5,6.0),1)
      samp11b<-c(1:5)#sample(1:5,1)
      samp12<-"PRESENT"#c("MIN","MAX","PRESENT")
      samp13<-kngs
      # }
      rr=" "
      
      #test across accuracy    
      if(tt==2){
        
        samp1<-sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
        rr=samp1
        
      }
      
      #test across regions
      if(tt==3){
        
        samp2<-sort(sample(c("Africa","Oceania","North America","South America","Asia","Europe"),5))
        count1<-c("Africa","Oceania","North America","South America","Asia","Europe")[!c("Africa","Oceania","North America","South America","Asia","Europe") %in% samp2 ]
        rr=count1
      }
      
      #test across vector non vector
      if(tt==4){
        samp3<-sample(list(0,1,2),1)[[1]]
        rr=samp3
      }
      
      #test across RCP  
      if(tt==5){
        samp4<-sample(c(2.6,4.5,6.0),1)
        rr=samp4
      }
      
      
      #test across SSP -cases only 
      if(tt==6){
        samp5<-sample(1:5,1)
        rr=samp5
      }
      
      
      #test across current case numbers
      if(tt==7){
        samp6<-sample(c(0.25,0.5,0.75,1),1)
        rr=samp6
      }
      
      #test non disease
      if(tt==8){
        samp7<-"XXX"
        rr=" "
      }  
      
      #test non disease
      if(tt==9){
        samp8<-"hantavirus"
        rr=" "
      }  
      
      #test non disease
      if(tt==10){
        samp9<-sample(c("bacteria","parasite","fungi","virus"),1)
        rr=samp9
      }
      
      #test non disease
      if(tt==11){
        samp10<-sample(c("Temperate","Tropical"),1)
        rr=samp10
      }
      
      
      #test non disease
      if(tt==12){
        
        samp11a<-sample(c(2.6,4.5,6.0),1)
        samp11b<-sample(1:5,1)
        
        rr=paste(samp11a,samp11b,sep="_")
      }
      
      #test non disease
      if(tt==13){
        samp12<-sample(c("MAX","PRESENT","MIN"),1)
        rr=samp12
      }
      
      #test non disease
      if(tt==14){
        samp13<-sample(kngs,1)
        rr=samp13
      }
      
      # }else{
      
      #sample for sensitivity
      # samp1<-sample(c(0.55,0.6,0.65,0.7,0.75),1)
      # samp2<-sample(c("Africa","Oceania","North America","South America","Asia","Europe"),5)
      # samp3<-sample(list(0,1,c(0,1)),1)[[1]]
      
      # }  
      
      count1<-c("Africa","Oceania","North America","South America","Asia","Europe")[!c("Africa","Oceania","North America","South America","Asia","Europe") %in% samp2 ]
      
      res6$count1=rowSums(cbind(res6[,names(res6) %in% count1],res6[,names(res6) %in% count1]))
      
      ##remove 8.5 - opposite results
      res6x<-res6[res6$AUC>samp1 & res6$count1!=0 & res6$Vectored %in% samp3 & res6$RCP %in% samp4 & res6$SSP %in% samp5 & res6$spillover_rate2 %in% samp6 & res6$disease.y!=samp7 & res6$group2!=samp8 & res6$Kingdom %in% samp9 & res6$Tropical %in% samp10 & res6$RCP %in% samp11a & res6$SSP %in% samp11b & res6$threshold %in% samp12 & res6$KingTrans %in% samp13,]#  & 
  
      #if(nrow(res6x)==0){next}
      
      
      res6x$var1=paste(paste(samp1,collapse=":"),paste(count1,collapse=":"),paste(samp3,collapse=":"),paste(samp4,collapse=":"),paste(samp5,collapse=":"),paste(samp6,collapse=":"),paste(samp7,collapse=":"),paste(samp8,collapse=":"),paste(samp9,collapse=":"),paste(samp10,collapse=":"),paste(samp11a,collapse=":"),paste(samp11b,collapse=":"),paste(samp12,collapse=":"),paste(samp13,collapse=":"),sep="_")
      
      if(res6x$var1[1] %in% dones) {next}
      
      dones<-c(dones,res6x$var1[1])
      
      ##set data type
      res6x$results_type=types[tt]
      ##incase two values
      res6x$ff=rr
  #table(is.na(res6x[,"group2"]))
  
  #length(unique(res6x$disease))  
  
  #x=1
  ###foreach
  NDIS<-100#c(40,80,120,160) ###check sim p value
  
  ##do maps=1 don't=0
  #zzz=0

  for(qq in 1:length(NDIS)){
  ##make dummy
    res6x$dummy<-1
      
    resN2<-NULL
    resN2b<-NULL
    
    ##loop through stuff
    for(jj in 4){
        
    if(jj==1){res6x$sort<-res6x$dummy;xxx<-"All"}
      
    if(jj==2){res6x$sort<-res6x$range;xxx<-"Geographic"}
      
    if(jj==3){res6x$sort<-res6x$group;xxx<-"Taxonomic"}
      
    if(jj==4){res6x$sort<-paste0(res6x$range,res6x$group);xxx<-"Both"}
   
      
      #library(doParallel)
      cl <- makeCluster(10)
      registerDoParallel(cl)
      resN2<-foreach(x=1:50, .combine=rbind) %dopar% {
       ###take lots of bootstraps to detect impact of certain diseases
      #for (x in 1:100){
    
      ##make bootstap id
      res6x$rep=x
      
      ###find random sets
      inside<-ave(seq_along(res6x$sort),res6x$sort,FUN=function(x) sample(length(x)))
      outside<-ave(inside,inside,FUN=function(x) sample(seq_along(x)))
      res6x<-res6x[order(inside,outside),]  
    
      ##random so different disease chosen by duplicated
      #res6x<-res6x[sample(1:nrow(res6x)),]
    
      ##make uni
      #res6x$group_range<-paste(res6x$group2,res6x$range,sep="_")
    
      ##make temp to choose one diases
    
      #length(unique(res6x$group2[1:n]))
      #length(unique(res6x$disease[1:n]))
    
      #resNtemp<-res6x[!duplicated(res6x$group2),]
      #resNtemp<-res6x[res6x$group2 %in% DIS,]
      resNtemp<-res6x[1:NDIS[qq],]
    
      ##choose all options of that disease
      resN<-res6x[res6x$disease %in% resNtemp$disease,]
    
      resN$dummy<-1
      #colSums(xtabs(data=resN,dummy~group2+range))
    
      #### per disease loop
      slp<-rep(NA,length(unique(resN$disease)))
      sim_slp<-rep(NA,length(unique(resN$disease)))
      sim_sd<-rep(NA,length(unique(resN$disease)))
      
      
      slp2<-rep(NA,length(unique(resN$disease)))
      sim_slp2<-rep(NA,length(unique(resN$disease)))
      sim_sd2<-rep(NA,length(unique(resN$disease)))
      
     
      slp3<-rep(NA,length(unique(resN$disease)))
      sim_slp3<-rep(NA,length(unique(resN$disease)))
      sim_sd3<-rep(NA,length(unique(resN$disease)))
      
      
      ##run per disease loop
      for(zz in 1:length(unique(resN$disease))){
      
        dis1<-unique(resN$disease)[zz]
      
        r3<-resN[resN$disease==dis1 ,]
      
        if(length(unique(r3$Year))<3){next}
        
        #r4<-r3
        #coordinates(r4)<-~x+y
        #plot(r4)
        
        
        #r3$richness<-scale(r3$richness)
        r3<-r3[sample(1:nrow(r3)),]
        
        r3<-r3[!duplicated(r3$richness),]
      
        #r3$Year=r3$Year-2030
        
        #lm1<-lm(data=r3,formula="richness~Year")$coefficients
        lm1<-cor(r3$richness,r3$Year,method="spearman",use="pairwise.complete.obs")[c(1,1)]
        lm2<-cor(r3$cases,r3$Year,method="spearman",use="pairwise.complete.obs")[c(1,1)]
        lm3<-cor(r3$poor_cases,r3$Year,method="spearman",use="pairwise.complete.obs")[c(1,1)]
        
        slp[i]<-lm1[2]
        slp2[i]<-lm2[2]
        slp3[i]<-lm3[2]
        
        reps<-250
        cor1<-rep(NA,reps)
        cor2<-rep(NA,reps)
        cor3<-rep(NA,reps)
        
        for (hh in 1:reps){
          ##make simulated data
          r4<-data.frame(Year=r3$Year)
          r4$richness<-sample(r3$richness)
          r4$cases<-sample(r3$cases)
          r4$poor_cases<-sample(r3$poor_cases)
          if(nrow(r4)<3){next}
          
          
          #lm1<-lm(data=r3,formula="richness~Year")$coefficients
          lm1s<-cor(r4$richness,r4$Year,method="spearman",use="pairwise.complete.obs")[c(1,1)]
          lm2s<-cor(r4$cases,r4$Year,method="spearman",use="pairwise.complete.obs")[c(1,1)]
          lm3s<-cor(r4$poor_cases,r4$Year,method="spearman",use="pairwise.complete.obs")[c(1,1)]
          
          cor1[hh]<-lm1s[2]
          cor2[hh]<-lm2s[2]
          cor3[hh]<-lm3s[2]
          
          
           }## and of hh simulate
          
        sim_slp[zz]<-mean(cor1,na.rm=TRUE)
        sim_sd[zz]<-sd(cor1,na.rm=TRUE)
        sim_slp2[zz]<-mean(cor2,na.rm=TRUE)
        sim_sd2[zz]<-sd(cor2,na.rm=TRUE)
        sim_slp3[zz]<-mean(cor3,na.rm=TRUE)
        sim_sd3[zz]<-sd(cor3,na.rm=TRUE)
        #if(is.null(cor2)){cor2<-cor1} else {cor2<-c(cor2,cor1)}
        #print(i)
        }## for i in n diseases 
      
      
    
      #hist(slp,main=paste(i,"_",round(length(slp[slp>0])/length(slp),2)))
      #mean(slp)
    
      cordata1<-data.frame(disease=unique(resN$disease),correlation=slp,correlation2=slp2,correlation3=slp3,bootstrap=x,NDIS=NDIS[qq],ndis=length(unique(resN$disease)),panel=types[tt],panel2=rr,pval=sim_slp,sd=sim_sd,pval2=sim_slp2,sd2=sim_sd2,pval3=sim_slp3,sd3=sim_sd3,vars=rr)
    
      #cordata2<-data.frame(disease=unique(resN$disease),correlation=cor2,bootstrap=x,NDIS=NDIS[qq],ndis=length(unique(resN$disease)),panel=xxx,pval=NA)
      
      
      #resN<-merge(resN,cordata1,by="disease")
      #resNb<-merge(resN,cordata2,by="disease")
      #cordata1$results_type=types[tt]
      ##incase two values
      #cordata1$ff=rr
      

      
      ##make maps
      ## uniqe diseases
      #kk<-unique(resNtemp$disease)
    
      return(cordata1)
      #if(is.null(res7)) {res7<-cordata1} else {res7<-rbind(res7,cordata1)}
      
      #if(is.null(resN2b)) {resN2b<-resNb} else {resN2b<-rbind(resN2b,resNb)}
      
      
      print(x)
      } ##end of x
  
      stopCluster(cl)
      
      if(is.null(res7)) {res7<-resN2} else {res7<-rbind(res7,resN2)}

      
    } ### end of jjj
    
    
    #if(is.null(res7b)) {res7b<-resN2b} else {res7b<-rbind(res7b,resN2b)}
    
    
    }## end NDIS
  
    } ##end of combinations
    print(tt)
    } ### end of different groups
    
  
    
  res7$correlationS<-rnorm(n=rep(1,nrow(res7)),mean=res7$pval,sd=res7$sd)
  res7$correlationS[res7$correlation2<(-1)]<-(-1)
  res7$correlationS[res7$correlation2>(1)]<-(1)
  
  res7$correlationSb<-rnorm(n=rep(1,nrow(res7)),mean=res7$pval2,sd=res7$sd2)
  res7$correlationSb[res7$correlationSb<(-1)]<-(-1)
  res7$correlationSb[res7$correlationSb>(1)]<-(1)
  
  res7$correlationSc<-rnorm(n=rep(1,nrow(res7)),mean=res7$pval3,sd=res7$sd3)
  res7$correlationSc[res7$correlationSc<(-1)]<-(-1)
  res7$correlationSc[res7$correlationSc>(1)]<-(1)
    
  ##write to disck
  #write.csv(res7,file="C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\workings04052022.csv")
  res7<-read.csv(file="C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\workings06042022.csv",stringsAsFactors = F)


    #unique(res7$panel)  
    #[1] "Baseline"           "Accuracy"           "Vectored"          
    #[4] "RCP"                "Burden"             "Emerging"          
    #[7] "Hantaviridae"       "Kingdom"            "Tropical/Temperate"
    #[10] "RCP/SSP"  
    
    #res7b<-res7[!is.na(res7$correlation) & res7$panel2!="fungi_0" & res7$panel!="Emerging" & res7$panel=="Breadth",] #& res7$panel=="RCP"
  
    #res7b<-res7[!is.na(res7$correlation) & res7$panel2!="fungi_0" & res7$panel!="Emerging" & res7$panel=="RCP/SSP",] #& res7$panel=="RCP"
  
    res7b<-res7[!is.na(res7$correlation) & res7$panel2!="fungi_0" & res7$panel!="Emerging" & res7$panel!="KingTrans",] #& res7$panel=="RCP"
    
    
      length(unique(res7b$disease))
  
    #res7[1,]
    res7b$panel2[res7b$panel2==" "]<-""
    res7b$Panel<-paste(res7b$panel,res7b$panel2,sep=" ")
    #hist(res7$correlation,main=paste(i,"_",round(length(res7$correlation[res7$correlation>0])/length(res7$correlation),2)))
    #hist(res7$correlation2,add=TRUE,col="red")
    
    ##get rid of zeros - but why do we have them??? is it land-use dominated diseases?
    #res7<-res7[!(res7$richness==0 & res7$Year!=2010) ,]
  
    #res7$uni<-paste(res7$RCP,res7$rep)
    #res7$uni<-paste(res7$Vectored,res7$rep)
    #res7$uni<-paste(res7$disease,res7$rep)
    res7b$uni<-paste(res7b$bootstrap,res7b$NDIS,res7b$Panel)
    
    #res7b<-res7b#[!duplicated(paste0(res7$correlation,res7$bootstrap,res7$NDIS,res7$panel)),]
    
    res7b$uni2<-1
    
    
    ##choose what to show
    res7b$correlationX<-res7b$correlation
    res7b$correlationSX<-res7b$correlationS
    
        #res7b$Trend<-cut(res7b$correlation,breaks=c(-1,-0.75,-0.25,0.25,0.75,1))
    res7b$Trend<-cut(res7b$correlationX,breaks=c(-1,-0.25,0.25,1))
    #res7b$Trend<-cut(res7b$correlationX,breaks=c(-1,-0.1,0.1,1))
    #res7b$Trend<-gsub("0.1","0.25",res7b$Trend,fixed=TRUE)
    #res7b$Trend<-cut(res7b$correlationX,breaks=c(-1,-0.5,0.5,1))
   # res7b$Trend<-gsub("0.5","0.25",res7b$Trend,fixed=TRUE)
    
        
    #table(res7b$Trend)
    res7b$Trend<-gsub("(","",res7b$Trend,fixed=TRUE)
    res7b$Trend<-gsub("]","",res7b$Trend,fixed=TRUE)
    res7b$Trend<-gsub(","," to ",res7b$Trend,fixed=TRUE)
    
    #res7b$Trend<-cut(res7b$correlation,breaks=c(-1,-0.75,-0.25,0.25,0.75,1))
    res7b$Trend2<-cut(res7b$correlationSX,breaks=c(-1,-0.25,0.25,1))
    #res7b$Trend2<-cut(res7b$correlationSX,breaks=c(-1,-0.1,0.1,1))
    #res7b$Trend2<-gsub("0.1","0.25",res7b$Trend2,fixed=TRUE)
    #res7b$Trend2<-cut(res7b$correlationSX,breaks=c(-1,-0.5,0.5,1))
    #res7b$Trend2<-gsub("0.5","0.25",res7b$Trend2,fixed=TRUE)
    
    #table(res7b$Trend2)
    res7b$Trend2<-gsub("(","",res7b$Trend2,fixed=TRUE)
    res7b$Trend2<-gsub("]","",res7b$Trend2,fixed=TRUE)
    res7b$Trend2<-gsub(","," to ",res7b$Trend2,fixed=TRUE)
    
    #res7b$Trend<-factor(res7b$Trend, levels=c("-1 to -0.75","-0.75 to -0.25","-0.25 to 0.25","0.25 to 0.75","0.75 to 1"))
    
    #res7b$Trend2<-factor(res7b$Trend2, levels=c("-1 to -0.25","-0.25 to 0.25","0.25 to 1"))
    
        ##findout how often in categories
    res7c<-aggregate(res7b$uni2,by=list(res7b$bootstrap,res7b$NDIS,res7b$Panel,res7b$Trend),sum,na.rm=TRUE)
    
    names(res7c)<-c("bootstrap","NDIS","Panel","Trend","Count")
    
    res7c$Trend2<-"Contraction"
    res7c$Trend2[res7c$Trend=="-0.25 to 0.25"]<-"Static"
    res7c$Trend2[res7c$Trend=="0.25 to 1"]<-"Expansion"
    
    res7c$Trend<-factor(res7c$Trend2, levels=c("Contraction","Static","Expansion"))
    
    
    ##findout how often in categories
    ###simulated data
    res7c2<-aggregate(res7b$uni2,by=list(res7b$bootstrap,res7b$NDIS,res7b$Panel,res7b$Trend2),sum,na.rm=TRUE)
    
    names(res7c2)<-c("bootstrap","NDIS","Panel","Trend","Count")
    
    res7c2$Trend2<-"Contraction"
    res7c2$Trend2[res7c2$Trend=="-0.25 to 0.25"]<-"Static"
    res7c2$Trend2[res7c2$Trend=="0.25 to 1"]<-"Expansion"
    
    res7c2$Trend<-factor(res7c2$Trend2, levels=c("Contraction","Static","Expansion"))
    
    ##Change labels
    
    res7c$Panel<-gsub("Kingdom ","",res7c$Panel,fixed=TRUE)
    res7c$Panel<-gsub("Tropical/Temperate ","",res7c$Panel,fixed=TRUE)
    res7c$Panel<-gsub("Vectored 0","Non-vectored",res7c$Panel,fixed=TRUE)
    res7c$Panel<-gsub("Vectored 1","Host-Vectored",res7c$Panel,fixed=TRUE)
    res7c$Panel<-gsub("Vectored 2","Vectored",res7c$Panel,fixed=TRUE)
    res7c$Panel<-gsub("bacteria","Bacteria",res7c$Panel,fixed=TRUE)
    res7c$Panel<-gsub("virus","Virus",res7c$Panel,fixed=TRUE)
    res7c$Panel<-gsub("parasite","Parasite",res7c$Panel,fixed=TRUE)
    res7c$Panel<-gsub("Baseline","All",res7c$Panel,fixed=TRUE)
    
    
    res7c$Panel2<-gsub("KingTrans ","",res7c$Panel)
    rrr<-read.table(text=res7c$Panel2,sep="_",stringsAsFactors = F)
    names(rrr)<-c("Kingdom","Transmission")
    res7c<-cbind(res7c,rrr)
    
    res7c$Transmission<-gsub("0","Non-vectored",res7c$Transmission,fixed=TRUE)
    res7c$Transmission<-gsub("1","Host-Vectored",res7c$Transmission,fixed=TRUE)
    res7c$Transmission<-gsub("2","Vectored",res7c$Transmission,fixed=TRUE)
    
    
    #res7c$Panel2<-gsub("RCP/SSP ","",res7c$Panel)
    #rrr<-read.table(text=res7c$Panel2,sep="_",stringsAsFactors = F)
    #names(rrr)<-c("RCP","SSP")
    #res7c<-cbind(res7c,rrr)
    #res7c$Count=res7c$Count/sum(res7c$Count,na.rm=TRUE)
    
    #agg1<-aggregate(res7c$Count,by=list(res7c$Panel,res7c$Trend),mean)
    
    
    res7c2$Panel<-gsub("Kingdom ","",res7c2$Panel,fixed=TRUE)
    res7c2$Panel<-gsub("Tropical/Temperate ","",res7c2$Panel,fixed=TRUE)
    res7c2$Panel<-gsub("Vectored 0","Non-vectored",res7c2$Panel,fixed=TRUE)
    res7c2$Panel<-gsub("Vectored 1","Host-Vectored",res7c2$Panel,fixed=TRUE)
    res7c2$Panel<-gsub("Vectored 2","Vectored",res7c2$Panel,fixed=TRUE)
    res7c2$Panel<-gsub("bacteria","Bacteria",res7c2$Panel,fixed=TRUE)
    res7c2$Panel<-gsub("virus","Virus",res7c2$Panel,fixed=TRUE)
    res7c2$Panel<-gsub("parasite","Parasite",res7c2$Panel,fixed=TRUE)
    res7c2$Panel<-gsub("Baseline","All",res7c2$Panel,fixed=TRUE)
    
    #res7c2$Panel2<-gsub("RCP/SSP ","",res7c2$Panel)
    #rrr2<-read.table(text=res7c2$Panel2,sep="_",stringsAsFactors = F)
    #names(rrr2)<-c("RCP","SSP")
    #res7c2<-cbind(res7c2,rrr2)
    #res7c2$Count=res7c2$Count/sum(res7c2$Count,na.rm=TRUE)
    
    
    
    res7c2$Panel2<-gsub("KingTrans ","",res7c2$Panel)
    rrr<-read.table(text=res7c2$Panel2,sep="_",stringsAsFactors = F)
    names(rrr)<-c("Kingdom","Transmission")
    res7c2<-cbind(res7c2,rrr)
    
    res7c2$Transmission<-gsub("0","Non-vectored",res7c2$Transmission,fixed=TRUE)
    res7c2$Transmission<-gsub("1","Host-Vectored",res7c2$Transmission,fixed=TRUE)
    res7c2$Transmission<-gsub("2","Vectored",res7c2$Transmission,fixed=TRUE)
    
    
    #table(res7c$Panel,res7c$Trend)
    
    #res7c$Panel<-factor(res7c$Panel,levels=c("All ","Temperate","Tropical","Bacteria","Parasite","Virus","Non-vectored","Host-Vectored","Vectored"))  
    #res7c2$Panel<-factor(res7c2$Panel,levels=c("All ","Temperate","Tropical","Bacteria","Parasite","Virus","Non-vectored","Host-Vectored","Vectored"))  
    
    gg<-
      ggplot(data=res7c[,],aes(x=Trend,y=Count,col=Panel)) +
    #res7$uni<-paste(res7$SSP,res7$rep)]
    #geom_density(alpha=0.1,col=NA)+
    #ggplot(res7,aes(x=jitter(Year),y=richness,col=as.factor(SSP),group=uni)) +
    #stat_summary(fun =  mean, geom = "point",size=2) + 
    stat_summary(data=res7c2[,],aes(x=Trend,y=Count),col="light grey",alpha=0.75,fun =  median, geom = "point",size=2,position = position_dodge(3)) + 
    stat_summary(data=res7c2[,],aes(x=Trend,y=Count),col="light grey",alpha=0.75,fun.data=stat_HPDI, geom = "errorbar",width=0.2,position = position_dodge(3)) + 
    #stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    stat_summary(fun =  median, geom = "point",size=2,position = position_dodge(3)) + 
    #stat_summary(fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3))+
    stat_summary(fun.data = stat_HPDI, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_line()+
    #scale_y_continuous( expand = c(0, 0)) +
    #scale_x_continuous( expand = c(0, 0)) +
    #geom_line(stat="smooth",method = "lm",col="seagreen",
    #          size = 1.5,
    #          alpha = 0.1)+    #geom_point(size=3,alpha=0.5) +
    #guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    #geom_smooth(data=res7,aes(x=jitter(Year),y=richness,group=uni2),se=F,method="lm",col="black")+
    #geom_violin()+
    #scale_x_log10()+
    #xlim(2010,2080)+
    #ylim(0,20500)+
    theme_cowplot()+
    #scale_x_discrete(guide = guide_axis(n.dodge=3))+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("Trend")+
    ylab("Frequency")+
      theme(axis.text.x=element_text(angle=45, hjust=1),
            panel.border = element_rect(colour = "black", fill=NA, size=1))+
    #geom_abline(intercept=1,slope=0,lty=2)+
    #geom_hline(data=res7c,aes(yintercept=mean(Count,na.rm=TRUE),group=Panel),col="black",lty=2)+
      #xlim(2,9)+
    #facet_grid(RCP~SSP)+
    facet_grid(Kingdom~Transmission,scales="free_y")+
    #facet_wrap(.~Panel,ncol=3,scales="free_y")+
    theme(legend.position = "none")#+
      #ggtitle("C")
  #facet_grid(Vectored~SSP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")
  
    set_palette(gg,"aaas")
    
    
    
    
    
    table(res6$KingTrans)
    
    
    
    gg<-ggplot(data=res7c,aes(x=Trend,y=Count,col=as.factor(NDIS))) +
      #res7$uni<-paste(res7$SSP,res7$rep)]
      #geom_density(alpha=0.1,col=NA)+
      #ggplot(res7,aes(x=jitter(Year),y=richness,col=as.factor(SSP),group=uni)) +
      #stat_summary(fun =  mean, geom = "point",size=2) + 
      stat_summary(data=res7c2,aes(x=Trend,y=Count,group=as.factor(NDIS)),col="light grey",alpha=0.75,fun =  median, geom = "point",size=2,position = position_dodge(3)) + 
      stat_summary(data=res7c2,aes(x=Trend,y=Count,group=as.factor(NDIS)),col="light grey",alpha=0.75,fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3)) + 
      stat_summary(fun =  median, geom = "point",size=2,position = position_dodge(3)) + 
      #stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
      #stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
      #stat_summary(fun =  median, geom = "point",size=2,position = position_dodge(3)) + 
      stat_summary(fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3))+
      #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2)+
      #geom_line()+
      #scale_y_continuous( expand = c(0, 0)) +
      scale_x_discrete( expand = c(0, 0)) +
      #geom_line(stat="smooth",method = "lm",col="seagreen",
      #          size = 1.5,
      #          alpha = 0.1)+    #geom_point(size=3,alpha=0.5) +
      #guides(colour = guide_legend(override.aes = list(alpha = 1)))+
      #geom_smooth(data=res7,aes(x=jitter(Year),y=richness,group=uni2),se=F,method="lm",col="black")+
      #geom_violin()+
      #scale_x_log10()+
    #xlim(2010,2080)+
    #ylim(0,20500)+
    theme_cowplot()+
      #scale_x_discrete(guide = guide_axis(n.dodge=3))+
      #geom_line(aes(group =  disease ),color="grey")+
      xlab("Trend")+
      ylab("Frequency")+
      theme(axis.text.x=element_text(angle=45, hjust=1))+
      #geom_abline(intercept=1,slope=0,lty=2)+
      #geom_vline(xintercept = 0,lty=2)+
      #xlim(2,9)+
      #facet_grid(Vectored~RCP)#+
      facet_wrap(.~NDIS,scales="free")+
      theme(legend.position = "none")
    #facet_grid(Vectored~SSP,scales="free")
    #facet_grid(Vectored~SSP,scales="free")
    
    set_palette(gg,"aaas")
  
  
  ###maps
  
do_pres=FALSE
NDIS=200
qq=1
futureG2030<-NULL
futureG2050<-NULL
futureG2070<-NULL
futureG2080<-NULL
pp=2080
jj=6.0

for (xxx in 1:10){  
    ##only disease and good
    res6x<-res6[ res6$AUC>0.5 & res6$Disease=="Disease" & res6$threshold=="ALL" ,]
    
        ###choose subsection
        #mods4<-mods2[mods2$Year %in% pp & mods2$SSP=="ssp3" & mods2$RCP==jj,]# & mods2$SSP=="ssp2" ,]
        
        #  mods4<-mods2[mods2$RCP==6.0,]
        
        #ff2=0
        #sum1=NULL
        ##loop to create map
        #for (f in 1:250){
        
        #  res6q<-res6[ res6$AUC>0.5 & res6$Disease=="Disease" & res6$threshold=="ALL",]
    
    #unbiased
    res6x$sort<-paste0(res6x$range,res6x$group)
  
    
    ###find random sets
    inside<-ave(seq_along(res6x$sort),res6x$sort,FUN=function(x) sample(length(x)))
    outside<-ave(inside,inside,FUN=function(x) sample(seq_along(x)))
    res6x<-res6x[order(inside,outside),]  
    
    ##get subset
    resNtemp<-res6x[1:NDIS[qq],]
    
    ##choose all options of that disease
    #resN<-res6x[res6x$disease %in% resNtemp$disease,]
    
    
    kk<-unique(resNtemp$disease)
  
    #for(pp in c(2030,2050,2070,2080)){
     
    pp=2070 
      
      for (jj in c(2.6,4.5,6.0)){
        
    #jj=6.0
    
  # create loop for maps
 
     present2<-NULL
     futureG2<-NULL
     futureF2<-NULL
     j=1
     for (j in j:length(kk) ){
      
      ##check point data # point_data$LU  %in% d2c$disease
      pd<-point_data[(point_data$name_LU  %in% kk[j]) & point_data$Status!="Denied" & point_data$SumCases>0,]
      if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd);projection(pd)<-projection(raster())} else {rw1<-1000}
      
      ##disease cell id
      cells<-fread(ci1[ci2==kk[j]])
      
      cells<-cells[cells$year_RCP==paste(jj,pp,sep="_"),]
      if(nrow(cells)==0){next}
      cells[,dummy:=1]
      
      #times=length(unique(cells$year_RCP))
      
      ##get disease details
      cur1<-res6[res6$disease==kk[j],][1,]
      
      #countries
      countr1<-wrld_simpl3[wrld_simpl3$NAME!="Antartica" & wrld_simpl3$ISO2 %in% strsplit(cur1$countries.x,",")[[1]],]
      
      #regions
      reg1<-wrld_simpl3[wrld_simpl3$REGION %in% countr1$REGION,]
      reg2<-fasterize(st_as_sf(reg1),template)
      
      ##future
      cells2<-cells[time1!="present",]#.(meanval=sum(dummy)/times),by=cell.id]
      
      ##present
      cellsP<-cells[time1=="present",]#.(meanval=sum(dummy)/times),by=cell.id]
      
      ##gain
      cellsG<-cells2[!cells2$cell.id %in% cellsP$cell.id,]
      
      ##loss
      cellsL<-cellsP[!cellsP$cell.id %in% cells2$cell.id,]
      
      
      ######################## DO PRESENT??? #######################
      if(do_pres==TRUE){
      
      present1<-template2
      present1[cells$cell.id[cells$time1=="present"]]<-1
      
      ##make results
      dis_trans<-cur1
      
      ##if real points
      if(nrow(pd)>0){
        for (u in 1:5){
          real1<-raster::extract(present1,pd)#[,c("Longitude","Latitude")]
          rp<-randomPoints(reg2,5000)
          not_real<-raster::extract(present1,rp)
          
          #res6d<-data.frame(Obs=,Fit=c(real,not_real))
          
          e1<-evaluate(p=real1,a=not_real)
          thr1w<-data.frame(best.stat=e1@auc,cutoff=e1@t[which.max(e1@TPR + e1@TNR)], sensitivity=NA, specificity=NA)
          #thr1w<-Find.Optim.Stat(Stat="ROC",Fit=res6d$Fit,Obs=res6d$Obs)
          colnames(thr1w)<-paste("AUC_",colnames(thr1w),sep="")
          
          if(u==1){thr2w=thr1w} else {thr2w<-rbind(thr2w,thr1w)}
        }  
        realthr<-colMeans(thr2w,na.rm = TRUE)
        dis_trans$AUC_cutoff_real<-realthr[2]
        dis_trans$AUC_real<-realthr[1]
        dis_trans$rep=u
        
      } else {dis_trans$AUC_cutoff_real<-NA;dis_trans$AUC_real<-NA;dis_trans$rep=u} # end of if pb
      
      
      present2x<-present1
      present2x[present2x==0]<-NA
      
      if(!paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\present_maps\\present_",kk[j],".png",sep="") %in% list.files("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\present_maps\\",full.names=T)) {
        
        png(file=paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\present_maps\\present_",kk[j],".png",sep=""),width=1000,height=500)
        
        plot(mask1,col="#00000050",main=paste(kk[j],round(cur1$AUCmax,2),sep=" "),legend=FALSE)
        plot(present2x,add=TRUE,col="orange",legend=FALSE)
        plot(countr1,lwd=1,border="olivedrab",add=TRUE)
        if(nrow(pd)>0){points(pd)}
        dev.off()
      }
      
      if(is.null(present2)) {present2<-present1;dis_trans2<-dis_trans} else {present2<-present2+present1;dis_trans2<-rbind(dis_trans2,dis_trans)} 
      
      } #end of do_pres################
      
      
      ##future gain raster
      futureG<-template2
      futureG[cellsG$cell.id]<-cellsG$dummy
      
      ##future loss raster
      futureL<-template2
      futureL[cellsL$cell.id]<-cellsL$dummy
      
      ##loop though diseases
      if(is.null(futureG2)) {futureG2<-futureG;futureL2<-futureL} else {futureG2<-futureG2+futureG;futureL2<-futureL2+futureL} 
      
      #print(j)
      
      rm(cells,cells2,cellsG,cellsL);gc() 
    } ##end of j loop

     ##future limit    
     futureL2[futureL2>6]<-6
     futureG2[futureG2>6]<-6

     net1<-log(futureG2+1)-log(futureL2+1)
     
     
     
     
     ##loop though diseases
     if(pp==2030){
       if(is.null(futureG2030)) {futureG2030<-futureG2;futureL2030<-futureL2;net2030<-net1} else {futureG2030<-futureG2030+futureG2;futureL2030<-futureL2030+futureL2;net2030<-net2030+net1} 
     }
     
     
     ##loop though diseases
     if(pp==2050){
       if(is.null(futureG2050)) {futureG2050<-futureG2;futureL2050<-futureL2;net2050<-net1} else {futureG2050<-futureG2050+futureG2;futureL2050<-futureL2050+futureL2;net2050<-net2050+net1} 
     }
     
     
     ##loop though diseases
     if(pp==2070){
       if(is.null(futureG2070)) {futureG2070<-futureG2;futureL2070<-futureL2;net2070<-net1} else {futureG2070<-futureG2070+futureG2;futureL2070<-futureL2070+futureL2;net2070<-net2070+net1} 
     }
     
     
     ##loop though diseases
     if(pp==2080){
       if(is.null(futureG2080)) {futureG2080<-futureG2;futureL2080<-futureL2;net2080<-net1} else {futureG2080<-futureG2080+futureG2;futureL2080<-futureL2080+futureL2;net2080<-net2080+net1} 
     }
     
     
     
     #sort out 0s for plotting
     futureG2[futureG2==0]<-NA
     futureL2[futureL2==0]<-NA
     net1[net1==0]<-NA
     
      
     ##gain
     png(file=paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\future_maps\\AAfuture_gain_",paste(jj,pp,sample(1:1e6,1),sep="_"),".png",sep=""),width=1000,height=500)
     
     plot(mask1,col="#00000010",main=paste(jj,pp,"gain",sep=" "),legend=FALSE)
     plot(log(futureG2),add=TRUE,col=rev(viridis(50,option="C"))[10:50],legend=TRUE)
     #plot(countr1,lwd=1,border="olivedrab",add=TRUE)
     #if(nrow(pd)>0){points(pd)}
     dev.off()
     
     ##loss
     png(file=paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\future_maps\\AAfuture_loss_",paste(jj,pp,sample(1:1e6,1),sep="_"),".png",sep=""),width=1000,height=500)
     
     plot(mask1,col="#00000010",main=paste(jj,pp,"loss",sep=" "),legend=FALSE)
     plot(log(futureL2),add=TRUE,col=rev(viridis(50,option="C"))[10:50],legend=TRUE)
     #plot(countr1,lwd=1,border="olivedrab",add=TRUE)
     #if(nrow(pd)>0){points(pd)}
     dev.off()
     
     ##net
     png(file=paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\future_maps\\AAfuture_net_",paste(jj,pp,sample(1:1e6,1),sep="_"),".png",sep=""),width=1000,height=500)
     
     plot(mask1,col="#00000010",main=paste(jj,pp,"net change",sep=" "),legend=FALSE)
     plot(net1,add=TRUE,col=rev(viridis(50,option="C"))[8:50],legend=TRUE)
     #plot(countr1,lwd=1,border="olivedrab",add=TRUE)
     #if(nrow(pd)>0){points(pd)}
     dev.off()
     
     
 
   }##end of RCP
  
  
  #} #end of year
  
    print("END LOOP")
    print(xxx)
    
}### end of xxx reps
  

net2080b<-futureG2080-futureL2080
net2080b<-net2080b/10


#net2080b<-focal(net2080b,w=matrix(1/81,9,9),na.rm=TRUE)
library(terra)
#net2080b<-aggregate(net2080b,30,method='bilinear')
#net2080b<-disaggregate(net2080b,30,method='bilinear')

net2080b<-mask(net2080b,mask1)

net2080b[net2080b==0]<-NA

#net2080b[net2080b<(-1*max(values(net2080b),na.rm=TRUE))]=(-1*max(values(net2080b),na.rm=TRUE))

net2080b<-rast(net2080b)

net2080c<-abs((net2080b)/max(values(net2080b),na.rm=TRUE))*255
net2080c[is.na(net2080c)]<-0


net2080x<-futureG2080+futureL2080

#hist(net2080x)


net2080x[net2080x>quantile(values(net2080x),0.995)]<-quantile(values(net2080x),0.995)

net2080x<-((net2080x^2)/(max(values(net2080x^2),na.rm=TRUE)))#*255

#net2080x<-mask(net2080x,mask1)
net2080x<-rast(net2080x)

##net
#png(file=paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\future_maps\\BBfuture_net_",paste(jj,pp,sep="_"),".png",sep=""),width=1000,height=500)

plot(mask1,col="#00000010",main=paste(jj,pp,"net change",sep=" "),legend=FALSE)
plot(net2080b,add=TRUE,col=rev(viridis(50,option="B")),legend=TRUE,alpha=net2080x)
#plot(countr1,lwd=1,border="olivedrab",add=TRUE)
#if(nrow(pd)>0){points(pd)}
#dev.off()


