library(raster)
library(data.table)
library(sp)
library(sf)
library(maptools)
library(fasterize)
library(rgdal)
data(wrld_simpl)


##base raster template
template<-raster(nrow=3600, ncol=8640,ext=extent(-180, 180, -60, 90),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


##get gini data
## move cases in poor country
## worse reporting in poor countries
## per grid cell poverty - so number of people in poverty
gini<-fread("C:\\Users\\xxxx\\Dropbox\\New_Global_MAXENT/ssp_gini_countries2.csv",header=TRUE)
cdata<-merge(wrld_simpl@data,gini,by.x="ISO3",by.y="ISO3",all.x=TRUE)
cdata<-cdata[order(match(cdata$NAME,wrld_simpl$NAME)),]
cdata[,(15:ncol(cdata))] <- lapply(cdata[,(15:ncol(cdata))], function(x) ifelse(is.na(x), mean(x, trim=0.1,na.rm = TRUE), x))
identical(cdata$NAME,wrld_simpl$NAME)
wrld_simpl@data<-cdata


for (i in c(21,27:46)){
  
#if(i %in% c(22:26,30,31,35,36,40,41,45,46)){next}
  
  gini3<-fasterize(st_as_sf(wrld_simpl),template,field=names(cdata)[i])
  names(gini3)<-names(cdata)[i]
  if(i==21){gin<-gini3}else{gin<-stack(gin,gini3)}
}


writeRaster(gin,"C:\\Users\\Public\\Documents\\GINI_STACK1.tif",format="GTiff",overwrite=TRUE)

write.csv(names(gin),"C:\\Users\\Public\\Documents\\GINI_STACK_NAMES1.csv")

