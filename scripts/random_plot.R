
pl1<-ws2
pl1[pl1==1]<-0

##remove cells outside intersection of countries and host ranges
#res9b<-res9[ !is.na(res9$countr2),]
res9b<-res9[ !is.na(res9$host_range) & !is.na(res9$surrounding),]

##populate with HazardT
pl1[res9b$cell.id]<-res9b[ ,names(res9b)[names(res9b)=="Hazard"],drop=TRUE]

#if(zz==9){xxx<-(log(res9$Hazard+1))-(log(res9$Hazard_future+1));xxx[xxx==0]=NA;pl1[res9$cell.id]<-xxx}
#crop to general area
pl1<-crop(pl1,extent( min(res9b$x)-5, max(res9b$x)+5,min(res9b$y)-5, max(res9b$y)+5))
pl1[pl1==0]<-NA


plot(pl1)
plot(ad2,add=TRUE,border="red")
