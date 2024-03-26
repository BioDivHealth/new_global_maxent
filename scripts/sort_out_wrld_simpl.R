library(smoothr)
library(sp)
library(rgeos)

data(wrld_simpl)
wrld_simpl[wrld_simpl$SUBREGION==5,"REGION"]<-20

wrld_simpl<-wrld_simpl[wrld_simpl$NAME!="Antarctica",]
wrld_simpl$FIPS<-as.character(wrld_simpl$FIPS)
wrld_simpl$ISO2<-as.character(wrld_simpl$ISO2)
wrld_simpl$ISO3<-as.character(wrld_simpl$ISO3)
wrld_simpl$NAME<-as.character(wrld_simpl$NAME)



# make SpatialPoints
points <- sp::SpatialPoints(cbind(c(65.74730195113234, 51.45324979300407),c(69.70728168243, 46.375529205168085)))

#81.91376746912066, 80.87150494464777
#39.80455725744677, 42.814670025058525
#69.70728168243, 65.74730195113234
#46.375529205168085, 51.45324979300407
# use as to convert to line
sp_line <- as(points,"SpatialLines")
projection(sp_line)<-projection(wrld_simpl)

# make SpatialPoints
points2 <- sp::SpatialPoints(cbind(c(180, 180),c(90, -90)))

# use as to convert to line
sp_line2 <- as(points2,"SpatialLines")
projection(sp_line2)<-projection(wrld_simpl)
sp_lines<-rbind(sp_line,sp_line2)

russia<-wrld_simpl[wrld_simpl$NAME=="Russia",]

projection(russia)<-projection(wrld_simpl)

lpi <- gIntersection(russia, sp_lines)               
# intersect your line with the polygon
blpi <- gBuffer(lpi, width = 0.000001)  
# create a very thin polygon buffer of the intersected line
dpi <- gDifference(russia, blpi)                # split using gDifference

dpi2 <- disaggregate(dpi)

rus_sib<-dpi2[c(3,2),]

new_data1<-rbind(wrld_simpl@data[wrld_simpl$NAME=="Russia",],wrld_simpl@data[wrld_simpl$NAME=="Russia",])
new_data1$NAME[2]<-"Siberia"
new_data1$FIPS[2]<-"SR"
new_data1$ISO2[2]<-"SS"
new_data1$ISO3[2]<-"SIB"
new_data1$LON[2]<-(102)
new_data1$LAT[2]<-(68.5)

rownames(new_data1)<-getSpPPolygonsIDSlots(rus_sib)

rus_sib2<-SpatialPolygonsDataFrame(rus_sib,data=new_data1)

wrld_simpl2<-rbind(wrld_simpl[wrld_simpl$NAME!="Russia",],rus_sib2)

#US

# make SpatialPoints
points <- sp::SpatialPoints(cbind(c(-168, -168),c(73, 51)))


#73.13852782859281, -168.01135557405593
#51.26370721808503, -169.41760545546762

# use as to convert to line
sp_line <- as(points,"SpatialLines")
projection(sp_line)<-projection(wrld_simpl)

# make SpatialPoints
points2 <- sp::SpatialPoints(cbind(c(-141, -141),c(73, 41)))

# use as to convert to line
sp_line2 <- as(points2,"SpatialLines")
projection(sp_line2)<-projection(wrld_simpl)
sp_lines<-rbind(sp_line)#,sp_line2)

US<-wrld_simpl[wrld_simpl$NAME=="United States",]

projection(US)<-projection(wrld_simpl)

lpi <- gIntersection(US, sp_lines)               
# intersect your line with the polygon
blpi <- gBuffer(lpi, width = 0.000001)  
# create a very thin polygon buffer of the intersected line
dpi <- gDifference(US, blpi)                # split using gDifference

dpi2 <- disaggregate(dpi)

us_als<-dpi2[c(1,2),]

new_data1<-rbind(wrld_simpl@data[wrld_simpl$NAME=="United States",],wrld_simpl@data[wrld_simpl$NAME=="United States",])
new_data1$NAME[2]<-"Alaska"
new_data1$FIPS[2]<-"AK"
new_data1$ISO2[2]<-"AK"
new_data1$ISO3[2]<-"ALK"
new_data1$LON[2]<-(-151)
new_data1$LAT[2]<-(65)

rownames(new_data1)<-getSpPPolygonsIDSlots(us_als)

us_als2<-SpatialPolygonsDataFrame(us_als,data=new_data1)

##add together
wrld_simpl2<-rbind(wrld_simpl[wrld_simpl$NAME!="United States" & wrld_simpl$NAME!="Russia", ],us_als2,rus_sib2)

wrld_simpl3<-drop_crumbs(wrld_simpl2, 1000000000, drop_empty = TRUE)

save(wrld_simpl3,file="C:\\Users\\david.redding\\Dropbox\\New_Global_MAXENT\\wrld_simpl3.r")

