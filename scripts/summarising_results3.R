

###run disease validation
library(dismo)
library(maptools)
library(data.table)
library(sp)
library(stringr)

data(wrld_simpl)

library(raster)
library(fasterize)
library(dismo)
library(sf)
library(rgdal)

library(viridis)
library(sp)

library(ggplot2)
library(velox)
library(cowplot)

####new summary file

##read in admins to move away from GRID
ad1<-readOGR("C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\admin1\\ne_10m_admin_1_states_provinces.shp","ne_10m_admin_1_states_provinces")

###read final results
res11<-fread(file="C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\final_risk1.csv")
res11<-res11[!is.na(res11$year_RCP),]
res11<-res11[!is.na(res11$diss_me),]

##number of diseases
res11[,n:=1]

##make broad grouping
res11[,group3:=group,]
res11[str_detect(res11[["group"]],"vir"),"group3"]<-"virus"
#set(res11, i=which(str_detect(res11[["group"]],"vir")), j=group3, value="virus")

###agregate
#res12<-	res11[auc>0.7 & cases_per_year>10  ,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]]

### compare global patterns for future
for(k in 1:3){
##potential

  #if(k==1){res12<-	res11[auc>0.5 & cases_per_year<1  ,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]];type="potential"}

##emerging
  #if(k==2){res12<-	res11[auc>0.5 & cases_per_year>1 & cases_per_year<=1000,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]];type="emerging"}

##endemic
  #if(k==3){res12<-	res11[auc>0.5 & cases_per_year>1000 & cases_per_year<=100000,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]];type="endemic"}

##hyper endemic
  #if(k==4){res12<-	res11[auc>0.5 & cases_per_year>100000,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]];type="hyperendemic"}


  ###types non vector
  if(k==1){res12<-	res11[group3=="virus" & auc>0.5 & Type %in% c("HOST->HUMAN","HOST->HUMAN->HUMAN") & cases_per_year>10 ,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]];type="non_vector_virus"}
  
  ##vector
  if(k==2){res12<-	res11[group3=="virus" & auc>0.5 & Type %in% unique(res11$Type)[c(2,4,5,6)] & cases_per_year>10,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]];type="vector_virus"}
  
  ###types non vector
  if(k==3){res12<-	res11[group3=="bacteria" & auc>0.5 & Type %in% c("HOST->HUMAN","HOST->HUMAN->HUMAN") & cases_per_year>10 ,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]];type="non_vector_bacteria"}
  
  ##vector
  if(k==4){res12<-	res11[group3=="bacteria" & auc>0.5 & Type %in% unique(res11$Type)[c(2,4,5,6)] & cases_per_year>10,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]];type="vector_bacteria"}
  
  
  
  
  year_RCP<-c("2.6_2030","4.5_2030","6_2030","8.5_2030","2.6_2050","4.5_2050","6_2050","8.5_2050","2.6_2070","4.5_2070","6_2070","8.5_2070","2.6_2080","4.5_2080","6_2080","8.5_2080")
  year_RCP<-sort(year_RCP)

  #for(h in 1:length(year_RCP)){

  for(h in c(4,length(year_RCP))){
    
    ##choose one scenario
    res13<-res12[res12$year_RCP==year_RCP[h],]

    ## plot result
    ad1b<-ad1#[ad1$diss_me %in% res13$diss_me,]
    res13<-res13[match(ad1$diss_me,res13$diss_me),]
    res13[is.na(res13)]<-0
    ad1b@data<-cbind(ad1b@data,as.data.frame(res13))

    ###just convert to sf and then can used natively
    ad1c<-st_as_sf(ad1b)
    ad1c<-ad1c[ad1c$name!="Antarctica",]
    for (j in 2){
      SSPx<-c("SSP1","SSP2","SSP3")[j]
      if(SSPx=="SSP1"){delta="delta_riskSSP1"}
      if(SSPx=="SSP2"){delta="delta_riskSSP2"}
      if(SSPx=="SSP3"){delta="delta_riskSSP3"}

      ##when best breaks are found 
      jpeg(file=paste("C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\Summary_",year_RCP[h],"_",SSPx,"_",type,"_",sample(1:1000,1),".jpg",sep=""),width=2500,height=1000)
      print(
        ggplot(data = ad1c) +
          geom_sf(aes(fill = (!!ensym(delta))/n),colour=NA) +
  #      scale_fill_viridis_c(option = "plasma")+
          scale_fill_gradientn(colors=viridis_pal(option = "plasma")(25),limits=c(-1.5,2.5),na.value="light grey")+
                  labs(title=paste("RCP",year_RCP[h],SSPx,type,"n",round(mean(ad1c$n,na.rm=TRUE),0),sep="-"))+
          theme(plot.title=element_text(size=32),axis.text = element_text(size=28))
    #+	facet_wrap(.~year_RCP)
        )
      dev.off()

    }##end of j loop

    if(h==1){res13<-res12}else{res13<-rbind(res13,res12)}
    print(h)
  }##end of h loop
  ##create final line plots (Hazard and Risk)

} ##end of k loop

##agregate
res12<-	res11[res11$auc>0.5 & cases_per_year>100 & cases_per_year<500000 ,lapply(.SD,function (x) sum(x,na.rm=TRUE)),by=.(year_RCP,diss_me),.SDcols=names(res11)[c(22:36,84:88)]]

res12[,delta_risk:=rowMeans(res12[,.(delta_riskSSP1,delta_riskSSP2,delta_riskSSP3)],na.rm=TRUE)]

year_RCP<-c("2.6_2030","4.5_2030","6_2030","8.5_2030","2.6_2050","4.5_2050","6_2050","8.5_2050","2.6_2070","4.5_2070","6_2070","8.5_2070","2.6_2080","4.5_2080","6_2080","8.5_2080")
#year_RCP<-sort(year_RCP)

years<-c(2030,2050,2070,2080)

for(h in 1:length(years)){

   ##choose one scenario
  res13b<-res12[res12$year_RCP %in% year_RCP[((h*4)-3):(h*4)],]
  res13b[,c("RCP","year"):=tstrsplit(year_RCP,"_")]
  
  for (j in 1:4){
    res13<-res13b[RCP==unique(res13b$RCP)[j],]
    
    ## plot result
  ad1b<-ad1#[ad1$diss_me %in% res13$diss_me,]
  res13<-res13[match(ad1$diss_me,res13$diss_me),]
  res13[is.na(res13)]<-0
  ad1b@data<-cbind(ad1b@data,as.data.frame(res13))
  
  ###just convert to sf and then can used natively
  ad1c<-st_as_sf(ad1b)
  ad1c<-ad1c[ad1c$name!="Antarctica",]
    
  assign(paste("ad1c",j,sep=""),ad1c)
  ##when best breaks are found 
    p<-ggplot(data = get(paste("ad1c",j,sep=""))) +
        geom_sf(aes(fill = delta_riskSSP2/n),colour=NA) +
        #      scale_fill_viridis_c(option = "plasma")+
        scale_fill_gradientn(colors=viridis_pal(option = "plasma")(25),limits=c(-1.5,2.5),na.value="light grey")+
        labs(title=paste("Year",ad1c$year[1],"- RCP", ad1c$RCP[1],sep=" "))+
        theme(plot.title=element_text(size=32),axis.text = element_text(size=28))
    assign(paste("p",j,sep=""),p)
    
  }##end of j loop
  
  jpeg(file=paste("C:\\Users\\David\\Dropbox\\New_Global_MAXENT\\RCP_grid_Summary_",years[h],"_",sample(1:1000,1),".jpg",sep=""),width=2500,height=1000)
  
  print(
  plot_grid(p1,p2,p3,p4)
  )
  dev.off()
  
  }##end of h loop




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



