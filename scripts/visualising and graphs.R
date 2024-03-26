

library(ggplot2)

template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
values(template)<-1:ncell(template)


load("X:\\DRtemp\\per_disease\\ tanganya full.r")


res7[,rcp_year_ssp:=paste(RCP,year,ssp,sep="_")]

res7[,x:=(xyFromCell(template,res7$cell.id))[,"x"]]

res7[,y:=(xyFromCell(template,res7$cell.id))[,"y"]]



loops<-unique(res7$rcp_year_ssp)


for (i in 1:length(loops)){
  

  res8<-res7[res7$rcp_year_ssp==loops[i],]
    
  t2<-template
  
  values(t2)[!is.na(values(t2))]<-0
  
  t2[res8$cell.id]<-res8$future_realised_clim_upper
  
  plot(t2)
  
  
}
