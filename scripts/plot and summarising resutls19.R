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

source("C:\\Users\\david.redding\\Dropbox\\R_scripts\\bivariate_plot_functions2.R")


##read in empress-i point data
point_data<-read.csv("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\all_empress_i_data4.csv",stringsAsFactors = FALSE)
point_data<-point_data[!is.na(point_data$SumCases),]
#point_data$LU<-paste(" ",point_data$name_LU,sep="")

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


###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl3),template)

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



res1<-list.files("D:\\more_results2\\",pattern="ALL.csv",full.names=TRUE)
res2<-list.files("D:\\more_results2\\",pattern="ALL.csv",full.names=FALSE)

#res1<-list.files("C:\\Users\\david.redding\\Documents\\more_results5\\",pattern="ALL.csv",full.names=TRUE)
#res2<-list.files("C:\\Users\\david.redding\\Documents\\more_results5\\",pattern="ALL.csv",full.names=FALSE)
res2<-gsub("_ALL","",res2)
res2<-gsub(".csv","",res2)
res2b<-read.table(text=res2,sep="_",stringsAsFactors = FALSE)

res3<-NULL
for (i in 1:length(res1)){
  res2x<-read.csv(res1[i],stringsAsFactors = F)
  if(length(names(res2x)[names(res2x) %in% c("disease")])==0){res2x$disease=res2b$V1[i]}
  if(is.null(res3)){res3<-res2x} else {res3<-rbind(res3,res2x)}
}

#res3a<-res3[res3$model=="present",]
#names(res3a)<-paste0(names(res3a),"ALL")
#res3$model<-gsub("present","presentX0_2010XALL",res3$model)
  
res3b<-res3[res3$model!="present",]
res3b2<-read.table(text=res3b$model,sep="X",stringsAsFactors = FALSE)
res3b3<-read.table(text=res3b2$V2,sep="_",stringsAsFactors = FALSE)
names(res3b2)[c(1,3)]<-c("SSP","DISEASE_GROUP")
names(res3b3)<-c("RCP","Year")

#res3c<-merge(res3b,res3a[,2:14],by.x="disease",by.y="diseaseALL")

res3z<-cbind(res3b,res3b2[c(1,3)],res3b3)

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
d1$Outbreak<-0
d1$Vectored[d1$Type %in% c("HOST->HUMAN->HUMAN","HOST->HUMAN")]<-0
d1$Outbreak[d1$Type %in% c("HOST->VECTOR->HUMAN->HUMAN","HOST->VECTOR->HUMAN->VECTOR->HUMAN->HUMAN","HOST->HUMAN->HUMAN","HUMAN->VECTOR->HUMAN->HUMAN")]<-1
d1<-d1[,c(1:20,46:48)]

d1[d1$name=="omsk","countries"]<-"RU"

for (w in 1:nrow(d1)){
  
  d2<-d1[w,]
  
  cc<-strsplit(d2$countries,",")[[1]]
  
  ws2<-wrld_simpl3[wrld_simpl3$ISO2 %in% cc, ]
  
  r2<-data.table(disease=d2$name,popn=sum(ws2$POP2005,na.rm=TRUE), ws2@data,stringsAsFactors = F)
  
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
  
  
  res6z<-NULL
  for(tt in 1:5){
  

  
    dones=0
    #ff=1
    ###run sensitivity
    for(ff in 1:100){
  

    if(tt==1){
      samp1<-0.5#sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
      samp2<-NA#sample(c("Africa","Oceania","North America","South America","Asia","Europe"),6)
      samp3<-c(0,1)#sample(list(0,1,c(0,1)),1)[[1]]
      samp4<-c(2.6,4.5,6.0)#sample(c(2.6,4.5,6.0),1)
    }
    
    if(tt==2){
      
      samp1<-sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
      samp2<-NA#sample(c("Africa","Oceania","North America","South America","Asia","Europe"),6)
      samp3<-c(0,1)#sample(list(0,1,c(0,1)),1)[[1]]
      samp4<-c(2.6,4.5,6.0)#sample(c(2.6,4.5,6.0),1)
      
    }
    if(tt==3){
      
      samp1<-0.5#sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
      samp2<-sample(c("Africa","Oceania","North America","South America","Asia","Europe"),5)
      samp3<-c(0,1)#sample(list(0,1,c(0,1)),1)[[1]]
      samp4<-c(2.6,4.5,6.0)#sample(c(2.6,4.5,6.0),1)
      
    }
    
    if(tt==4){
      
      samp1<-0.5#sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
      samp2<-NA#sample(c("Africa","Oceania","North America","South America","Asia","Europe"),6)
      samp3<-sample(list(0,1,c(0,1)),1)[[1]]
      samp4<-c(2.6,4.5,6.0)#sample(c(2.6,4.5,6.0),1)
      
    }
    
      if(tt==5){
        samp1<-0.5#sample(c(0.5,0.55,0.6,0.65,0.7,0.75),1)
        samp2<-NA#sample(c("Africa","Oceania","North America","South America","Asia","Europe"),6)
        samp3<-c(0,1)#sample(list(0,1,c(0,1)),1)[[1]]
        samp4<-sample(c(2.6,4.5,6.0),1)
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
  res6x<-res6[res6$AUC>samp1 & res6$count1==0 & res6$Vectored %in% samp3 & res6$RCP %in% samp4 & res6$disease.y!="None reported" ,]#  & 
  
  
  ##incase two values
  samp3<-paste(samp3,collapse=":")
  
  ##subsample
  #if(ff==1){
  #  res6x$var1=rep("ALL",nrow(res6x))

   #     } else {
      
    res6x$var1=paste(samp1,paste(count1,collapse=";"),samp3,sep="_")
      
  #  }
  
  if(res6x$var1[1] %in% dones) {next}
  
  e <- simpleError(1)
  
  
  res6y<-setDT(res6x)
  
  res6y[,cor1:=cor(richness,Year,method="spearman")[1]
,by=disease]
  
  res6y[,cor2:=cor(cases,Year,method="spearman")[1]
        ,by=disease]
  
  res6y[,lm1:=tryCatch(lm(scale(richness)~Year)$coefficients[2], error = function(e) e)[1]
        ,by=disease]
  
  res6y[,lm2:=tryCatch(lm(scale(cases)~Year)$coefficients[2], error = function(e) e)[1]
        ,by=disease]
  
  res6y[, tt:=tt]
  
  if(is.null(res6z)){res6z<-res6y; dones<-res6y$var1[1]} else {res6z<-rbindlist(list(res6z,res6y)); dones<-c(dones,res6y$var1[1])}

  
  print(ff)
  
    }

  }
  
 
  res6a<-res6z
  #res6a<-rbind(res6z,res6a)
  
  res6a$var1<-gsub("__","_ALL_",res6a$var1)
  
  res6a[ ,c("cutoff","region","vectored"):=tstrsplit(var1,"_") ]
   
  #install.packages("waffle", repos = "https://cinc.rud.is")
  
  
  res6a[ , cor1a:=cut(cor1,breaks =3,labels=FALSE)]
  
  res6a$dummy=1
  
  defs<-aggregate(res6a$dummy,by=list(res6a$cutoff,res6a$region,res6a$vectored,res6a$cor1a),sum)
  
  defs$Freq<-round(defs$x/100,0)
  
  ggplot(defs[defs$Group.1!=0.5 & defs$Group.2=="ALL" & defs$Group.3=="0:1" , ], aes(values = Freq, fill = as.factor(Group.4))) +
    geom_waffle(n_rows = 6, size = 0.33, colour = "white") +
    #geom_bar(position = "fill")+
    #geom_pie
    coord_equal() +
    theme_void()+
    facet_wrap(vars(Group.1),ncol=6)+
    #theme_cowplot(16)+
    theme(legend.position = "none")
  
  
  ggplot(defs[defs$Group.1==0.5 & defs$Group.2!="ALL" & defs$Group.3=="0:1" , ], aes(values = Freq, fill = as.factor(Group.4))) +
    geom_waffle(n_rows = 6, size = 0.33, colour = "white") +
    #geom_bar(position = "fill")+
    #geom_pie
    coord_equal() +
    theme_void()+
    facet_wrap(vars(Group.2),ncol=6)+
    #theme_cowplot(16)+
    theme(legend.position = "none")
  
  
  ggplot(defs[defs$Group.1==0.5 & defs$Group.2=="ALL" & defs$Group.3!="0:1" , ], aes(values = Freq, fill = as.factor(Group.4))) +
    geom_waffle(n_rows = 6, size = 0.33, colour = "white") +
    #geom_bar(position = "fill")+
    #geom_pie
    coord_equal() +
    theme_void()+
    facet_wrap(vars(Group.3),ncol=6)+
    #theme_cowplot(16)+
    theme(legend.position = "none")
  
  
  ggplot(defs[defs$Group.1==0.5 & defs$Group.2=="ALL" & defs$Group.3=="0:1" , ], aes(values = Freq, fill = as.factor(Group.4))) +
    geom_waffle(n_rows = 12, size = 0.33, colour = "white") +
    #geom_bar(position = "fill")+
    #geom_pie
    coord_equal() +
    theme_void()+
    facet_wrap(vars(Group.3),ncol=6)+
    #theme_cowplot(16)+
    theme(legend.position = "none")
  
  
  ggplot(res6z,aes(x=cor1,col=var1)) +
    geom_density(adjust=2.5)+
    geom_density(data=res6z[res6z$var1=="ALL",],aes(x=cor1),lwd=1.5,col="black",adjust=2.5)+
    #res7$uni<-paste(res7$SSP,res7$rep)
    #ggplot(res7,aes(x=jitter(Year),y=richness,col=as.factor(SSP),group=uni)) +
    #stat_summary(fun =  mean, geom = "point",size=2,position = position_dodge(3)) + 
    #stat_summary(fun =  median, geom = "line",size=2,position = position_dodge(3),col="red") + 
    #stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  median, geom = "point",size=2,position = position_dodge(3)) + 
    #stat_summary(fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_point()+
    scale_y_continuous( expand = c(0, 0)) +
    scale_x_continuous( expand = c(0, 0)) +
    #geom_line(stat="smooth",method = "lm",col="seagreen",
    #          size = 1.5,
    #          alpha = 0.1)+    #geom_point(size=3,alpha=0.5) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    #geom_smooth(data=res7,aes(x=jitter(Year),y=richness,group=uni2),se=F,method="lm",col="black")+
    #geom_violin()+
    #scale_x_log10()+
    #xlim(2010,2080)+
    #ylim(0,20500)+
    theme_cowplot(16)+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("Direction of Change")+
    ylab("Frequency")+
    #geom_abline(intercept=0,slope=1,lty=2)+
    geom_vline(xintercept = 0,lty=2)+
    #xlim(2,9)+
    #facet_grid(Vectored~RCP)#+
  #facet_wrap(.~Group.1,scales="free")+
  theme(legend.position = "none")
  #facet_grid(Vectored~SSP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")
  
  
  ##
  #length(unique(res6x$disease))
  
  #x=1
  ###foreach
  
  ##do maps=1 don't=0
  zzz=0
  
  ###take lots of bootstraps to detect impact of certain diseases
  for (x in 1:100){
    
    ##make bootstap id
    res6x$rep=x

    inside<-ave(seq_along(res6x$range),res6x$range,FUN=function(x) sample(length(x)))
    outside<-ave(inside,inside,FUN=function(x) sample(seq_along(x)))
    res6x<-res6x[order(inside,outside),]  
    
    ##random so different disease chosen by duplicated
    #res6x<-res6x[sample(1:nrow(res6x)),]
    
    ##make uni
    #res6x$group_range<-paste(res6x$group2,res6x$range,sep="_")
    
    ##make temp to choose one diases
    resNtemp<-res6x[!duplicated(res6x$group2),]
    
    ##choose all options of that disease
    resN<-res6x[res6x$disease %in% resNtemp$disease,]
    
    resN$dummy<-1
    #colSums(xtabs(data=resN,dummy~group2+range))
    
    
    slp<-rep(0,length(unique(resN$disease)))
    
    for(i in 1:length(unique(resN$disease))){
      
      dis1<-unique(resN$disease)[i]
      
      r3<-resN[resN$disease==dis1 ,]
      
      r3$richness<-scale(r3$richness)
      
      #r3$Year=r3$Year-2030
      
      #lm1<-lm(data=r3,formula="richness~Year")$coefficients
      lm1<-cor(r3$richness,r3$Year,method="spearman")[c(1,1)]
      
      slp[i]<-lm1[2]
      
    }  
    
    #hist(slp,main=paste(i,"_",round(length(slp[slp>0])/length(slp),2)))
    #mean(slp)
    
    cordata1<-data.frame(disease=unique(resN$disease),correlation=slp,bootstrap=x)
    
    resN<-merge(resN,cordata1,by="disease")
    
    
    ##make maps
    ## uniqe diseases
    kk<-unique(resNtemp$disease)
    
    if(zzz==1){
      for (j in 1:length(kk) ){
      
        
        ##check point data # point_data$LU  %in% d2c$disease
        pd<-point_data[(point_data$name_LU  %in% kk[j]) & point_data$Status!="Denied" & point_data$SumCases>0,]
        if(nrow(pd)>0){coordinates(pd)<-~Longitude+Latitude;rw1<-nrow(pd);projection(pd)<-projection(raster())} else {rw1<-1000}
        
   
        
      ##disease cell id
      cells<-fread(ci1[ci2==kk[j]])
      cells[,dummy:=1]
      times=length(unique(cells$year_RCP))
      
      ##get disease details
      cur1<-resN[resN$disease==kk[j],][1,]
      
      #countries
        countr1<-wrld_simpl3[wrld_simpl3$NAME!="Antartica" & wrld_simpl3$ISO2 %in% strsplit(cur1$countries,",")[[1]],]
        
      #regions
        reg1<-wrld_simpl3[wrld_simpl3$REGION %in% countr1$REGION,]
        reg2<-fasterize(st_as_sf(reg1),template)
        
        
      cells2<-cells[time1!="present",.(meanval=sum(dummy)/times),by=cell.id]
      
      present1<-template2
      present1[cells$cell.id[cells$time1=="present"]]<-1
      
      ##make results
      dis_trans<-cur1
      
      ##if real points
      if(nrow(pd)>0){
        for (u in 1:5){
          real1<-extract(present1,pd)
          rp<-randomPoints(reg2,5000)
          not_real<-extract(present1,rp)
          
          
          res6d<-data.frame(Obs=,Fit=c(real,not_real))
          
          e1<-evaluate(p=real1,a=not_real)
          thr1w<-data.frame(best.stat=e1@auc,cutoff=e1@t[which.max(e1@TPR + e1@TNR)], sensitivity=NA, specificity=NA)
          #thr1w<-Find.Optim.Stat(Stat="ROC",Fit=res6d$Fit,Obs=res6d$Obs)
          colnames(thr1w)<-paste("AUC_",colnames(thr1w),sep="")
          
          if(u==1){thr2w=thr1w} else {thr2w<-rbind(thr2w,thr1w)}
        }  
        realthr<-colMeans(thr2w,na.rm = TRUE)
        dis_trans$AUC_cutoff_real<-realthr[2]
        dis_trans$AUC_real<-realthr[1]
        dis_trans$rep=x
        
      } else {dis_trans$AUC_cutoff_real<-NA;dis_trans$AUC_real<-NA;dis_trans$rep=x} # end of if pb
      
      
      present2x<-present1
     present2x[present2x==0]<-NA
      
     if(!paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\present_maps\\present_",kk[j],".png",sep="") %in% list.files("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\present_maps\\",full.names=T)) {
      
          png(file=paste("C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\present_maps\\present_",kk[j],".png",sep=""),width=1000,height=500)
      
          plot(mask1,col="#00000050",main=paste(kk[j],round(cur1$AUCmax,2),sep=" "),legend=FALSE)
        plot(present2x,add=TRUE,col="orange",legend=FALSE)
         plot(countr1,lwd=1,border="olivedrab",add=TRUE)
      
        dev.off()
     }
      
      
        
      future1<-template2
      future1[cells2$cell.id]<-cells2$meanval
        
      ##loop though diseases
      if(j==1) {present2<-present1;future2<-future1;dis_trans2<-dis_trans} else {present2<-present2+present1;future2<-future2+future1;dis_trans2<-rbind(dis_trans2,dis_trans)} 
     
      rm(cells);gc() 
      } ##end of j loop
    
    ##work out change
    change1<-future2-present2
    
    ##make new res7 ## this is the loop with maps
        if(x==1) {res7<-resN;present3<-present2;future3<-future2;change2<-change1;dis_trans3<-dis_trans2} else {res7<-rbind(res7,resN);present3<-present3+present2;future3<-future3+future2;change2<-change2+change1;dis_trans3<-rbind(dis_trans3,dis_trans2)}

    
    ##make new res7 ## this is the loop without maps
    }else{
      
        if(x==1) {res7<-resN} else {res7<-rbind(res7,resN)}
      
      }##end of zzz if 
    
        
    print(x)
  }
  
  res7<-res7[!is.na(res7$correlation),]
  
  length(unique(res7$disease))
  
  hist(res7$correlation,main=paste(i,"_",round(length(res7$correlation[res7$correlation>0])/length(res7$correlation),2)))
  
  ##get rid of zeros - but why do we have them??? is it land-use dominated diseases?
  #res7<-res7[!(res7$richness==0 & res7$Year!=2010) ,]
  
  res7$uni<-paste(res7$RCP,res7$rep)
  res7$uni<-paste(res7$Vectored,res7$rep)
  res7$uni<-paste(res7$disease,res7$rep)
  
  res7$uni2<-1
  
  ggplot(res7,aes(x=jitter(Year),y=richness,group=rep)) +
   #res7$uni<-paste(res7$SSP,res7$rep)
  #ggplot(res7,aes(x=jitter(Year),y=richness,col=as.factor(SSP),group=uni)) +
    #stat_summary(fun =  mean, geom = "point",size=2,position = position_dodge(3)) + 
    #stat_summary(fun =  median, geom = "line",size=2,position = position_dodge(3),col="red") + 
    #stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  median, geom = "point",size=2,position = position_dodge(3)) + 
    #stat_summary(fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #geom_line()+
    scale_y_continuous( expand = c(0, 0)) +
    scale_x_continuous( expand = c(0, 0)) +
    geom_line(stat="smooth",method = "lm",col="seagreen",
              size = 1.5,
              alpha = 0.1)+    #geom_point(size=3,alpha=0.5) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    geom_smooth(data=res7,aes(x=jitter(Year),y=richness,group=uni2),se=F,method="lm",col="black")+
    #geom_violin()+
    #scale_x_log10()+
    #xlim(2010,2080)+
    #ylim(0,20500)+
    theme_cowplot()+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("Year")+
    ylab("Change in endemic area")+
    geom_abline(intercept=1,slope=0,lty=2)+
    #xlim(2,9)+
    facet_grid(Vectored~RCP)#+
    #facet_wrap(.~Group.1,scales="free")+
    #theme(legend.position = "none")
  #facet_grid(Vectored~SSP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")
  
  
  ggplot(res7,aes(x=jitter(Year),y=cases,group=rep)) +
    #stat_summary(fun =  mean, geom = "point",size=2,position = position_dodge(3)) + 
    #stat_summary(fun =  median, geom = "line",size=2,position = position_dodge(3),col="red") + 
    #stat_summary(fun =  function(x) quantile (x,0.45), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun =  function(x) quantile (x,0.55), geom = "line",size=2,position = position_dodge(3),col="grey") + 
    #stat_summary(fun.data=median_hilow, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position = position_dodge(3))+
    #scale_color_manual(values=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0'))+
    #scale_color_viridis(discrete=TRUE,option="A") +
    scale_y_continuous( expand = c(0, 0)) +
    scale_x_continuous( expand = c(0, 0)) +
    geom_line(stat="smooth",method = "lm",col="seagreen",
              size = 1.5,
              alpha = 0.1)+    #geom_point(size=3,alpha=0.5) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    geom_smooth(data=res7,aes(x=jitter(Year),y=cases,group=uni2),se=F,method="lm",col="black")+
    #geom_violin()+
    #scale_x_log10()+
    #xlim(2010,2080)+
    #ylim(0,20500)+
    theme_cowplot()+
    #geom_line(aes(group =  disease ),color="grey")+
    xlab("Year")+
    ylab("Change in expected cases")+
    geom_abline(intercept=1,slope=0,lty=2)+
    #xlim(2,9)+
    facet_grid(.~SSP)+
    theme(legend.position = "none")
  #facet_grid(Vectored~SSP,scales="free")
  #facet_grid(Vectored~SSP,scales="free")
  
  
  ###maps
  
  
  gain<-change1
  loss<-change1
  gain[gain<=0]<-NA
  loss[loss>=0]<-NA
  loss=loss*-1
  
  gain2<-aggregate(gain,150,sum,na.rm=TRUE)
  loss2<-aggregate(loss,150,sum,na.rm=TRUE)
  gain2[gain2==0]<-NA
  loss2[loss2==0]<-NA
  
  # But let's use this simple code. You can change "nquantiles" to generate color matrices with different color schemes. For example, change it to 4 to produce a 4x4 color scheme.
  # You can specify the number of quantiles, colors and labels of your color matrix. Example:
  #col.matrix<-colmat(nquantiles=10, upperleft="blue", upperright="yellow", bottomleft="green", bottomright="red", xlab="Richness", ylab="Diversitication Rate")
  col.matrix<-colmat(nquantiles=10)
  
  bivmap<-bivariate.map((gain2),(loss2),colormatrix=col.matrix, nquantiles=10)
  
  plot(bivmap,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
  map(interior=T,add=T,col="#00000050") 