library(data.table)
library(raster)
library(reshape2)
#library(velox)
#library(ggplot2)
#library(ggpubr)
#library(cowplot)
#library(ggthemes)
library(sp)
library(rgdal)
library(sf)
library(fasterize)
library(RColorBrewer)


##base raster template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")



##change projecion??
##change projecion??
#load(file="D:/Users/xxxx/Documents/wrld_simpl3.r")
#wrld_simpl3<-spTransform(wrld_simpl3,CRS=projection(template))
#writeOGR(wrld_simpl3,driver="ESRI Shapefile",dsn="D:/Users/xxxx/Documents/wrld_simpl3.shp","wrld_simpl3")
wrld_simpl3<-readOGR("D:/Users/xxxx/Documents/wrld_simpl3.shp","wrld_simpl3")
#wrld_simpl3[wrld_simpl3$NAME %in% c("Russia","France"),]
wrld_simpl3@data[wrld_simpl3$NAME=="Russia","SUBREGION"]<-wrld_simpl3@data[wrld_simpl3$NAME=="Ukraine","SUBREGION"]
#wrld_simpl3[wrld_simpl3$NAME %in% c("Siberia","Mongolia"),]
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","SUBREGION"]<-30
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","REGION"]<-142
wrld_simpl3@data$ISO2<-as.character(wrld_simpl3@data$ISO2)
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","ISO2"]<-"XX"
wrld_simpl3@data$FIPS<-as.character(wrld_simpl3@data$FIPS)
wrld_simpl3@data[wrld_simpl3$NAME=="Siberia","FIPS"]<-"XX"
wrld_simpl3@data[wrld_simpl3$NAME=="Chile","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Paraguay","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Uruguay","SUBREGION"]<-4
wrld_simpl3@data[wrld_simpl3$NAME=="Argentina","SUBREGION"]<-4
###find neighbours

#neigh1a<-read.csv(file="D:/Users/xxxx/Documents/country_borders.csv",stringsAsFactors = FALSE)
#neigh1a<-neigh1a[!is.na(neigh1a$country_border_code),]

#wrld_simpl3[wrld_simpl3$NAME=="Western Sahara","AREA"]<-999
#wrld_simpl2[wrld_simpl2$NAME=="Luxembourg","AREA"]
#wrld_simpl3<-wrld_simpl3[wrld_simpl2$NAME!="Antarctica" & wrld_simpl2$AREA>10,]

###rasterize world simple
###rasterize world simple
ws2<-fasterize(st_as_sf(wrld_simpl3),template)


### plot and summarising resutls
### areas that have increased due to climate versus those that increased due to
### land use. Size of dots - burden estimate?
### panels for each year - facet across year and RCP?
### if over/under threshold due climate versus 

#gini<-stack("C:\\Users\\Public\\Documents\\GINI_STACK.tif")
#names(gini)<-read.csv("C:\\Users\\Public\\Documents\\GINI_STACK_NAMES.csv")$x

#gini3<-extract(gini,xyFromCell(template,1:ncell(template)))
#fwrite(gini3,file="C:\\Users\\Public\\Documents\\GINI_STACK.csv")
#names(gini3)<-read.csv("C:\\Users\\Public\\Documents\\GINI_STACK_NAMES.csv")$x
#gini4<-fread(file="C:\\Users\\Public\\Documents\\GINI_STACK.csv")
#gini4[,cell.id:=1:ncell(template)]

#gini5<-melt(gini4,id.vars="cell.id")
#setkey(gini5,value)
#gc()
#gini5<-gini5[!is.na(value),]
#gc()
#gini5[,c("year","scrap","SSP","scrap3"):= tstrsplit(variable,"_")]
#gc()
#gini5[,c("variable", "scrap","scrap3"):=NULL]
#gc()
#gini5[,year:=as.numeric(gsub("X","",year))]
#gc()
#gini5[,cell_by_year:=paste(cell.id,year,sep="_")]


#gini6<-gini5[year==2010,]
#gini6<-gini6[!duplicated(cell.id),]
#gini6[,c("year","SSP","cell_by_year"):=NULL]
#names(gini6)[2]<-"present_value"
#setkey(gini6,cell.id)
#gini5<-gini5[year!=2010,]
#gini5[, future_value2:=value]

#setkey(gini5,cell.id)
#gc()
#gini7<-gini6[gini5]
#gc()
#gini7[,future_value2:=NULL]
#fwrite(gini7,file="C:\\Users\\Public\\Documents\\GINI_by_cell.csv")

#gini5<-fread(file="D:/Users/xxxx/Documents/data/GINI_by_cell1.csv")
#setkey(gini5,cell_by_year)
#gini5<-dcast(gini4,year_dissme+year+diss_me~SSP,value.var = "value")
#setkey(gini5,year_dissme)


ad1<-list.files("D:/Users/xxxx/Documents/datasets1/",pattern="all_data3",full.names=TRUE)
ad2<-list.files("D:/Users/xxxx/Documents/datasets1/",pattern="all_data3",full.names=FALSE)
ad2<-gsub("_all_data4.csv","",ad2)
ad2<-gsub("_X_",";",ad2)
ad3<-read.table(text=ad2,sep=";",stringsAsFactors = F)$V1

gl1<-list.files("D:/Users/xxxx/Documents/gainlossresults/",pattern="_2.csv",full.names=TRUE)
gl2<-list.files("D:/Users/xxxx/Documents/gainlossresults/",pattern="_2.csv",full.names=FALSE)
gl2<-gsub("_2.csv","",gl2)
gl2<-gsub("_X_",";",gl2)
gl3<-read.table(text=gl2,sep=";",stringsAsFactors = F)$V1

##popdgp
#popgdp<-fread("C:\\Users\\Public\\Documents\\future_gdppop_maskedb.csv")
#setkey(popgdp,cell_by_year)

##popdgp
#pop<-fread("D:/Users/xxxx/Documents/data/future_pop_maskedb.csv")
#names(pop)[9]<-"humans_present"
#setkey(pop,cell_by_year)

##gini model (called lm1)
#load(file="D:/Users/xxxx/Documents/data/fitted_gini_10th.r")

err = simpleError("Error in read1")
err2 = simpleError("Error in read2")

#res5<-NULL
for (j in sample(1:length(ad1))){
  res11<-NULL
  
  
  ##load up each disease
  resF<-tryCatch(fread(ad1[j]), error=function(err) err)
  if(class(resF)[1]=="simpleError"){ next }
  
  ##add in name
  resF[,name:=strsplit(ad3[j],"_")[[1]][1]]
  
  if(paste(resF$name[1],"_",2,"_cellids.csv",sep="") %in% list.files("D:/Users/xxxx/Documents/cellids/",full.names=FALSE)){next}
  
  nnn<-resF$name[1]
  ##make RCP
  #resF[,c("RCP","scrap"):=tstrsplit(year_RCP,"_")]
  
  #remove unneeded columns
  #resF[,c("i.cell.id","i.year","i.cell.id.1","time","scrap"):=NULL]
  
  ### large ggplot for each disease with difference as once columns and scenario
  ### compare to combined scenario
  unis<-unique(resF$year_RCP)
  
  ##what columns?
  cols1<-c("people_in_povertySSP5","people_in_povertySSP5","people_in_povertySSP5","people_in_povertySSP4","people_in_povertySSP4","people_in_povertySSP4","people_in_povertySSP3","people_in_povertySSP3","people_in_povertySSP3","people_in_povertySSP2","people_in_povertySSP2","people_in_povertySSP2","people_in_povertySSP1","people_in_povertySSP1","people_in_povertySSP1","peopleSSP5","peopleSSP5","peopleSSP5","peopleSSP4","peopleSSP4","peopleSSP4","peopleSSP3","peopleSSP3","peopleSSP3","peopleSSP2","peopleSSP2","peopleSSP2","peopleSSP1","peopleSSP1","peopleSSP1","dummy","dummy","dummy","clim_up","clim_up","clim_up","clim_down","clim_down","clim_down","lu_up","lu_up","lu_up","lu_down","lu_down","lu_down")
  cols1<-unique(cols1)
  
  #par(mfrow=c(3,3))
  for(yy in 1:length(unis)){
    
    
    res9<-resF[resF$year_RCP==unis[yy] ,]
    
    res9<-as.data.frame(res9)
    #titls<-paste(c("HazardT","HazardL","HazardFLSSP3","Hazard Change"),round(mean(thr2[,1]),2),sep=" ")
    #cols1<-c("clim_present_mean2","Hazard","HazardT","HazardTD","Hazard_future","Hazard_futureTD","host_range","host_range2","RiskL")
    
    xx=11
    #for(xx in 1:length(cols1)){
    
    #plot(wrld_simpl)
    pl1<-ws2
    pl1[pl1==1]<-0
    
    ##future
    pl2<-ws2
    pl2[pl2==1]<-0
    
    
    ##remove cells outside intersection of countries and host ranges
    #res9b<-res9[ !is.na(res9$countr2),]
    res9b<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding)  & !is.na(res9$present2010) & res9$HazardT==1 ,]
    res9b$time1="present"
    res9c<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding)  & !is.na(res9$presentfuture) & res9$Hazard_futureT==1,]
    if(nrow(res9c)==0){next}
    res9c$time1="future"
    
    
    ##populate with HazardT present
    pl1[res9b$cell.id]<-res9b[ ,names(res9b)[names(res9b)==cols1[xx] ],drop=TRUE]
    #future
    pl2[res9c$cell.id]<-res9c[ ,names(res9c)[names(res9c)==cols1[xx]],drop=TRUE]
    
    #if(zz==9){xxx<-(log(res9$Hazard+1))-(log(res9$Hazard_future+1));xxx[xxx==0]=NA;pl1[res9$cell.id]<-xxx}
    #crop to general area
    #pl1<-crop(pl1,extent( min(res9b$x)-5, max(res9b$x)+5,min(res9b$y)-5, max(res9b$y)+5))
    #pl1[pl1==0]<-NA
    names(pl1)<-paste(cols1[xx],unis[yy],"present",sep="-")
    #pl2[pl2==0]<-NA
    names(pl2)<-paste(cols1[xx],unis[yy],"future",sep="-")
    
    #nnn<-res9b$name[1]     
    
    pdf(file=paste("D:/Users/xxxx/Documents/future_maps1/",nnn,"_",sample(1:1000,1),"all_plots.pdf",sep=""),width=8,height=8)
    # pdf(file=paste("C:\\Users\\david.redding\\Documents\\future_maps1\\",res9b$name[1],"_",sample(1:1000,1),"-all_plots.pdf",sep=""),width=16,height=8)
    
    plot(pl2-pl1,col=brewer.pal(3, "RdBu"),main=paste("In future - ",nnn," ",unis[yy],sep=""),colNA="seagreen")
    
    dev.off()
    
    ##choose cells in present and future
    res10<-rbind(res9b[,c("name","year_RCP","cell.id","time1")],res9c[,c("name","year_RCP","cell.id","time1")])
    
    
    if(is.null(res11)){res11<-res10} else {res11<-rbind(res11,res10)}
    
    ##
    rm(res9,res9b,res9c,res10)
    #plot(pl1)
    
    
  }
  write.csv(res11,file=paste("D:/Users/xxxx/Documents/cellids/",nnn,"_",2,"_cellids.csv",sep=""))
  
}

#write.csv(res11,file="D:/Users/xxxx/Documents/res11.csv")

