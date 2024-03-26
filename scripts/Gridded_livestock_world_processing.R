
library(raster)

template<-raster(nrow=3432, ncol=8640,ext=extent(-180, 180, -58, 85),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

livestock<-stack(list.files("E:\\Downloads\\GLWV2\\", pattern=".tif",full.names=T)[c(1,2,6,10,13,14)])
names(livestock)<-c("Cattle","Chickens","Ducks","Goats","Pigs","Sheep")

#xys<-xyFromCell(template,1:ncell(template))
#xys2<-as.data.frame(xys)
#coordinates(xys2)<-~x+y
#lstk<-extract(livestock,xys2)
#class(lstk)
#lt2<-as.data.frame(lstk)
#lt2$cell.id<-1:nrow(lt2)
#names(lt2)[1:6]<-c("Cattle","Chickens","Ducks","Goats","Pigs","Sheep")
#save(lt2,compress="xz",file="E:\\Dropbox\\legion_script\\livestock.r")

load(file="E:\\Dropbox\\legion_script\\livestock.r")


###### FUTURE

xxx1<-list.files("E:\\Documents\\countries\\",full.names=TRUE)

require(XLConnect)

for (i in 1:length(xxx1)){

wb = loadWorkbook(xxx1[i])
df = readWorksheet(wb, 1, header = TRUE)

iso3<-df$ISO3
x2000<-df[,7,drop=F]
x2000[is.na(x2000[,1]),1]<-median(x2000[,1],na.rm=T)
type=gsub("00","",names(x2000))

x2030<-df[,8,drop=F]
x2030[is.na(x2030[,1]),1]<-median(x2030[,1],na.rm=T)

res1<-data.frame(ISO3=iso3,x2000=x2000,x2030=x2030,type=type)
names(res1)[2:3]<-c("x2000","x2030")

if(i==1){res2<-res1} else { res2<-rbind(res2,res1)}

}

load(file="E:\\Dropbox\\Disease_analyses\\livestock_temp.r")
load(file="E:\\Dropbox\\Disease_analyses\\livestock_values_country.r")

library(maptools)
data(wrld_simpl)

ws1<-wrld_simpl
ws2<-wrld_simpl@data
ws2$ISO3<-as.vector(ws2$ISO3)

type1<-as.vector(unique(res2$type))

for (y in 1:length(type1)){

res3<-res2[res2$type==type1[y],]
res3$ISO3<-as.vector(res3$ISO3)
ws3<-merge(ws2,res3,all.x=T)
ws3$type<-as.vector(ws3$type)
ws4<-ws3[match(ws2$ISO3,ws3$ISO3),]
ws1@data<-ws4

x2000<-rasterize(ws1,template,field="x2000")
x2030<-rasterize(ws1,template,field="x2030")
s1<-stack(x2000,x2030)
names(s1)<-paste(type1[y],c(2000,2030),sep="_")


if(y==1) {livestock2<-s1} else {livestock2<-stack(livestock2,s1)}
#livestock3<-stack(livestock2,s1)
print(y)

}


save(livestock2,format="raster",file="livestock_temp.r")
save(res2,file="livestock_values_country.r")



####extract to template

xys<-xyFromCell(template,1:ncell(template))
xys2<-as.data.frame(xys)


coordinates(xys2)<-~x+y
lstk<-extract(livestock2,xys2)


lt2<-as.data.frame(lstk)
lt2$cell.id<-1:nrow(lt2)

#names(lt2)[1:6]<-c("Cattle","Chickens","Ducks","Goats","Pigs","Sheep")

save(lt2,compress="xz",file="E:\\Dropbox\\legion_script\\livestock_future.r")


setwd("E:\\Dropbox\\legion_script\\")


load(".//livestock_future.r")
lt3<-lt2
load(".//livestock.r")
cut_extreme<-function(x,threshold=0.95) {x[x>quantile(x,threshold,na.rm=TRUE)]<-quantile(x,threshold,na.rm=TRUE):return(x)}
lt2$Cattle<-cut_extreme(lt2$Cattle)
lt2$Chickens<-cut_extreme(lt2$Chickens)
lt2$Ducks<-cut_extreme(lt2$Ducks)
lt2$Goats<-cut_extreme(lt2$Goats)
lt2$Pigs<-cut_extreme(lt2$Pigs)
lt2$Sheep<-cut_extreme(lt2$Sheep)
lt2$Ducks[lt2$Ducks<0]<-0

lt2$Cattle2070=lt2$Cattle*((lt3$Beef_Prod_2030/lt3$Beef_Prod_2000)*2)
lt2$Chickens2070=lt2$Chickens*((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*2)
lt2$Ducks2070=lt2$Ducks*((((lt3$Egg_Prod_2030/lt3$Egg_Prod_2000)*2)+((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*2))/2)
lt2$Goats2070=lt2$Goats*((((lt3$Milk_Prod_2030/lt3$Milk_Prod_2000)*2)+((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*2))/2)
lt2$Pigs2070=lt2$Pigs*((lt3$Pork_Prod_2030/lt3$Pork_Prod_2000)*2)
lt2$Sheep2070=lt2$Sheep*((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*2)
lt2$Sheep2070[lt2$Sheep2070==Inf]<-0

lt2$Cattle2070<-cut_extreme(lt2$Cattle2070)
lt2$Chickens2070<-cut_extreme(lt2$Chickens2070)
lt2$Ducks2070<-cut_extreme(lt2$Ducks2070)
lt2$Goats2070<-cut_extreme(lt2$Goats2070)
lt2$Pigs2070<-cut_extreme(lt2$Pigs2070)
lt2$Sheep2070<-cut_extreme(lt2$Sheep2070)


save(lt2,compress="xz",file="E:\\Dropbox\\legion_script\\livestock_future2.r")

####### landcover

lc7<-raster("E:\\Dropbox\\Disease_analyses\\modis\\land2070_med_modis_486891.tif")
lc5<-raster("E:\\Dropbox\\Disease_analyses\\modis\\land2050_med_modis_486891.tif")

lc7b<-raster("E:\\Dropbox\\Disease_analyses\\modis\\land2070_med_modis_4778046.tif")
lc5b<-raster("E:\\Dropbox\\Disease_analyses\\modis\\land2050_med_modis_4778046.tif")


now<-raster("E:\\Dropbox\\Disease_analyses\\2005_modis_reduced.tif")

fins<-stack(lc5,lc5b,lc7,lc7b)

fins2<-resample(fins,template,method='ngb')

now2<-resample(now,template,method='ngb')

lc<-extract(stack(now2,fins2),xys2)

lc2<-data.frame(cell.id=1:nrow(lc),lc)
names(lc2)<-c("cell.id","landc","fland2050","fland2070","fland2050b","fland2050b")
save(lc2,compress="xz",file="E:\\Dropbox\\legion_script\\landcover_future3.r")

load(file="E:\\Dropbox\\legion_script\\landcover_future3.r")

lc2<-round(lc2,0)

save(lc2,compress="xz",file="E:\\Dropbox\\legion_script\\landcover_future4.r")


#### daily walking distance

names(WD)<-"daily_walking_distance"





