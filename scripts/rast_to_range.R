
require(sp)
require(rgeos)
require(smoothr)

rast_to_range<-function(x,crumb_size=1e9,present){

  pl2<-rasterToPolygons(x,dissolve=TRUE,na.rm=TRUE)

  pl3<-drop_crumbs(pl2, crumb_size, drop_empty = TRUE)

  pl4<-disaggregate(pl3)

  if(class(present)=="RasterLayer"){drops<-raster::extract(present,pl4,fun=mean,na.rm=TRUE,method='simple')}else{drops<-over(pl4,present)}
  
  if(length(na.omit(drops[,1]))==0){return(NA)}
  
  pl5<-pl4[!is.na(drops[,1]),]

  pl6<-fill_holes(pl5,threshold=1e32)

  pl6$type=cols1[zz]

  return(pl6)

}

##make buffer suitable for year
## https://science.sciencemag.org/content/333/6045/1024
## says mean 16.9km a year


##different between threshold poly for present pl7 and future pl6
#if(zz==2){pl8<-gDifference(pl6,pl7)}

#save present hazard
#if(zz==3){pl9<-pl6}

##same for hazard ##waht about shinking - need both ways 
#if(zz==4){pl10<-gDifference(pl6,pl9)}


#if(zz==1){pl7<-pl6} else {pl7<-rbind(pl7,pl6)}

##nop


#rm(pl1,pl2,pl3,pl4,pl5)#