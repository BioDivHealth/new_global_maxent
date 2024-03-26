

ws3<-ws2

ws3[!is.na(values(ws3))]<-0
rr<-res6[res6$year==2080,]

rr<-as.data.frame(rr)

ws3[rr$cell.id]<-rr$Hazard_future

plot(ws3)

ad1$AREA<-st

###test lot admins

res10<-res9[res9$year==2080,]
res10<-res10[match(ad1$diss_me,res10$diss_me),]

## plot result
ad1b<-ad1
ad1b@data<-cbind(ad1@data,as.data.frame(res10))

###just convert to sf and then can used natively
ad1c<-st_as_sf(ad1b)
ad1c$AREA<-as.numeric(st_area(ad1c))
ad1c<-ad1c[ad1c$name!="Antarctica",]

ggplot(data = ad1c) +
  geom_sf(aes(fill = (RiskSSP2)),colour=NA) +
  scale_fill_viridis_c(option = "plasma")

