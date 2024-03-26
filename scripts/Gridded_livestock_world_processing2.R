
library(raster)

template<-raster(nrow=3432, ncol=8640,ext=extent(-180, 180, -58, 85),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

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

lt2$Cattle2030=lt2$Cattle*((lt3$Beef_Prod_2030/lt3$Beef_Prod_2000)*1)
lt2$Chickens2030=lt2$Chickens*((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*1)
lt2$Ducks2030=lt2$Ducks*((((lt3$Egg_Prod_2030/lt3$Egg_Prod_2000)*1)+((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*1))/2)
lt2$Goats2030=lt2$Goats*((((lt3$Milk_Prod_2030/lt3$Milk_Prod_2000)*1)+((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*1))/2)
lt2$Pigs2030=lt2$Pigs*((lt3$Pork_Prod_2030/lt3$Pork_Prod_2000)*1)
lt2$Sheep2030=lt2$Sheep*((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*1)
lt2$Sheep2030[lt2$Sheep2030==Inf]<-0

lt2$Cattle2030<-cut_extreme(lt2$Cattle2030)
lt2$Chickens2030<-cut_extreme(lt2$Chickens2030)
lt2$Ducks2030<-cut_extreme(lt2$Ducks2030)
lt2$Goats2030<-cut_extreme(lt2$Goats2030)
lt2$Pigs2030<-cut_extreme(lt2$Pigs2030)
lt2$Sheep2030<-cut_extreme(lt2$Sheep2030)

lt2$Cattle2050=lt2$Cattle*((lt3$Beef_Prod_2030/lt3$Beef_Prod_2000)*1.66666667)
lt2$Chickens2050=lt2$Chickens*((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*1.66666667)
lt2$Ducks2050=lt2$Ducks*((((lt3$Egg_Prod_2030/lt3$Egg_Prod_2000)*1.66666667)+((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*1.66666667))/2)
lt2$Goats2050=lt2$Goats*((((lt3$Milk_Prod_2030/lt3$Milk_Prod_2000)*1.66666667)+((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*1.66666667))/2)
lt2$Pigs2050=lt2$Pigs*((lt3$Pork_Prod_2030/lt3$Pork_Prod_2000)*1.66666667)
lt2$Sheep2050=lt2$Sheep*((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*1.66666667)
lt2$Sheep2050[lt2$Sheep2050==Inf]<-0

lt2$Cattle2050<-cut_extreme(lt2$Cattle2050)
lt2$Chickens2050<-cut_extreme(lt2$Chickens2050)
lt2$Ducks2050<-cut_extreme(lt2$Ducks2050)
lt2$Goats2050<-cut_extreme(lt2$Goats2050)
lt2$Pigs2050<-cut_extreme(lt2$Pigs2050)
lt2$Sheep2050<-cut_extreme(lt2$Sheep2050)

lt2$Cattle2070=lt2$Cattle*((lt3$Beef_Prod_2030/lt3$Beef_Prod_2000)*2.333333)
lt2$Chickens2070=lt2$Chickens*((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*2.333333)
lt2$Ducks2070=lt2$Ducks*((((lt3$Egg_Prod_2030/lt3$Egg_Prod_2000)*2.333333)+((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*2.333333))/2)
lt2$Goats2070=lt2$Goats*((((lt3$Milk_Prod_2030/lt3$Milk_Prod_2000)*2.333333)+((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*2.333333))/2)
lt2$Pigs2070=lt2$Pigs*((lt3$Pork_Prod_2030/lt3$Pork_Prod_2000)*2.333333)
lt2$Sheep2070=lt2$Sheep*((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*2.333333)
lt2$Sheep2070[lt2$Sheep2070==Inf]<-0

lt2$Cattle2070<-cut_extreme(lt2$Cattle2070)
lt2$Chickens2070<-cut_extreme(lt2$Chickens2070)
lt2$Ducks2070<-cut_extreme(lt2$Ducks2070)
lt2$Goats2070<-cut_extreme(lt2$Goats2070)
lt2$Pigs2070<-cut_extreme(lt2$Pigs2070)
lt2$Sheep2070<-cut_extreme(lt2$Sheep2070)

lt2$Cattle2080=lt2$Cattle*((lt3$Beef_Prod_2030/lt3$Beef_Prod_2000)*2.66666667)
lt2$Chickens2080=lt2$Chickens*((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*2.66666667)
lt2$Ducks2080=lt2$Ducks*((((lt3$Egg_Prod_2030/lt3$Egg_Prod_2000)*2.66666667)+((lt3$Poul_Prod_2030/lt3$Poul_Prod_2000)*2.66666667))/2)
lt2$Goats2080=lt2$Goats*((((lt3$Milk_Prod_2030/lt3$Milk_Prod_2000)*2.66666667)+((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*2.66666667))/2)
lt2$Pigs2080=lt2$Pigs*((lt3$Pork_Prod_2030/lt3$Pork_Prod_2000)*2.66666667)
lt2$Sheep2080=lt2$Sheep*((lt3$Mut_Prod_2030/lt3$Mut_Prod_2000)*2.66666667)
lt2$Sheep2080[lt2$Sheep2080==Inf]<-0

lt2$Cattle2080<-cut_extreme(lt2$Cattle2080)
lt2$Chickens2080<-cut_extreme(lt2$Chickens2080)
lt2$Ducks2080<-cut_extreme(lt2$Ducks2080)
lt2$Goats2080<-cut_extreme(lt2$Goats2080)
lt2$Pigs2080<-cut_extreme(lt2$Pigs2080)
lt2$Sheep2080<-cut_extreme(lt2$Sheep2080)


lt2$cell.id=1:ncell(template)

sortx<-rowSums(lt2,na.rm=TRUE)
lt2<-lt2[lt2$cell.id!=sortx,]

save(lt2,compress="xz",file="E:\\Dropbox\\legion_script\\livestock_future_2030_2050_2070_2080.r")

