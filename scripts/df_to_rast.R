
df_to_rast<-function(template,data,column){

##plot to see where the data are 
#plot(wrld_simpl)
pl1<-template
pl1[pl1==1]<-0

##remove cells outside intersection of countries and host ranges
#res9b<-res9[ !is.na(res9$countr2),]
#res9b<-pres2

##populate with HazardT
pl1[data$cell.id]<-data[[ (1:length(data))[names(data)==column] ]]
names(pl1)<-column

return(pl1)

}


#if(zz==9){xxx<-(log(res9$Hazard+1))-(log(res9$Hazard_future+1));xxx[xxx==0]=NA;pl1[res9$cell.id]<-xxx}
#crop to general area
#pl1<-crop(pl1,extent( min(res9b$lon)-5, max(res9b$lon)+5,min(res9b$lat)-5, max(res9b$lat)+5))
#pl1[pl1==0]<-NA

