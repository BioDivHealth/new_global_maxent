
#template<-raster(nrow=3432, ncol=8640,ext=extent(-180, 180, -58, 85),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# wrld_simpl[wrld_simpl$NAME=="Russia","SUBREGION"]<-161
# wrld_simpl[wrld_simpl$NAME=="Russia","REGION"]<-180
# wrld_simpl[wrld_simpl$SUBREGION==13,"REGION"]<-18
# wrld_simpl[wrld_simpl$SUBREGION==34,"REGION"]<-141
# wrld_simpl[wrld_simpl$SUBREGION==5,"REGION"]<-40
# wrld_simpl[wrld_simpl$SUBREGION==35,"REGION"]<-125
# wrld_simpl[wrld_simpl$NAME=="Greenland","REGION"]<-15
# wrld_simpl[wrld_simpl$NAME=="Madagascar","REGION"]<-7
# wrld_simpl2<-wrld_simpl#[wrld_simpl$AREA>999,]
#save(wrld_simpl2,file="E:\\Dropbox\\Disease_analyses\\wrld_simpl2.r")
#regions<-rasterize(wrld_simpl,template,field="REGION",fun='last')
#save(regions,file="E:\\Dropbox\\Disease_analyses\\regions.r")
#subregions<-rasterize(wrld_simpl,template,field="SUBREGION",fun='last')
#save(subregions,file="E:\\Dropbox\\Disease_analyses\\subregions.r")


remove_small_islands<-function(polys,min_val=5){
temp3<-NULL
for(i in 1:length(polys)){
temp1x<-polys@polygons[[i]]
#print(length(temp1x@Polygons))
for(z in 1:length(temp1x@Polygons)){
temp2x<-temp1x@Polygons[[z]]
#print(slot(temp2x,"area"))}
if(slot(temp2x,"area")>min_val){
if(is.null(temp3)){temp3<-list(temp2x)} else {temp3<-append(temp3,temp2x)}#
}}}
temp4<-Polygons(temp3,"1")
temp4a<-Polygons(temp3,"2")
temp5<-SpatialPolygons(list(temp4,temp4a),1:2)

return(extent(temp5))
}



geog<-function(countries,wrld_simpl2,regions,subregions){
 #conts with ranges
countries<-gsub(" ","",fixed=T,countries)
list_count<-strsplit(countries,",")
 ws2<-wrld_simpl2[wrld_simpl2$ISO2 %in% unlist(list_count),]
 ws3<-wrld_simpl2[wrld_simpl2$SUBREGION %in% ws2$SUBREGION,]
 ws4<-wrld_simpl2[wrld_simpl2$REGION %in% ws2$REGION,]
ext1<-remove_small_islands(ws2)
if(ext1@xmin<(-165)){ext1@xmin<-(-165)}
ext2<-remove_small_islands(ws3)
if(ext2@xmin<(-165)){ext2@xmin<-(-165)}
ext3<-remove_small_islands(ws4)
if(ext3@xmin<(-165)){ext3@xmin<-(-165)}
ws2a<-suppressWarnings(crop(ws2,ext1))
ws3a<-suppressWarnings(crop(ws3,ext2))
ws4a<-suppressWarnings(crop(ws4,ext3))
ext1a<-extent(ws2a)
ext2a<-extent(ws3a)
ext3a<-extent(ws4a)

sr1<-NULL#rasterize(ws3a,regions)
r1<-rasterize(ws4a,regions)

return(list(ext1=ext1a,ext2=ext2a,ext3=ext3a,sr1=sr1,r1=r1))
}



geog2<-function(data1,geogX,wrld_simpl2,regions,subregions){

data2<-data1$data1[!is.na(extract(geogX$r1,data1$data1)),]
countries<-as.vector(over(data2,wrld_simpl2)$ISO2)
list_count<-unique(countries[!is.na(countries)])
 ws2<-wrld_simpl2[wrld_simpl2$ISO2 %in% unlist(list_count),]
 ws3<-wrld_simpl2[wrld_simpl2$SUBREGION %in% ws2$SUBREGION,]
 ws4<-wrld_simpl2[wrld_simpl2$REGION %in% ws2$REGION,]
ext1<-remove_small_islands(ws2)
if(ext1@xmin<(-165)){ext1@xmin<-(-165)}
ext2<-remove_small_islands(ws3)
if(ext2@xmin<(-165)){ext2@xmin<-(-165)}
ext3<-remove_small_islands(ws4)
if(ext3@xmin<(-165)){ext3@xmin<-(-165)}
ws2a<-suppressWarnings(crop(ws2,ext1))
ws3a<-suppressWarnings(crop(ws3,ext2))
ws4a<-suppressWarnings(crop(ws4,ext3))
ext1a<-extent(ws2a)
ext2a<-extent(ws3a)
ext3a<-extent(ws4a)

sr1<-NULL#rasterize(ws3a,regions)
r1<-rasterize(ws4a,regions)

return(list(ext1=ext1a,ext2=ext2a,ext3=ext3a,sr1=sr1,r1=r1))
}



process_species<-function(key,keyX,geog1,template){


locations3<-NULL
ext1<-geog1$ext2
locations1<-tryCatch(occ2df(occ(query=key, limit=50000,geometry=c(ext1[1],ext1[3],ext1[2],ext1[4]))),error=function(e) e)
if(class(locations1)[1]=="simpleError"){locations1<-data.frame(a=NULL,b=NULL)}
if(nrow(locations1)>0){
locations2<-locations1[!is.na(locations1$longitude),]
sortx<-paste(locations2$longitude,locations2$latitude,sep="")
locations3<-locations2[!duplicated(sortx),]
names(locations3)[2:3]<-c("lon","lat")
coordinates(locations3)<-~lon+lat
ext2<-ext1;sens=1}

if(is.null(locations3)){zzz=TRUE}else{if(nrow(locations3)<25){zzz=TRUE}else{zzz=FALSE}}

if(zzz==TRUE){
print("switched to region")
ext1<-geog1$ext3
locations1<-tryCatch(occ2df(occ(query=key, limit=50000,geometry=c(ext1[1],ext1[3],ext1[2],ext1[4]))),error=function(e) e)
if(class(locations1)[1]=="simpleError"){locations1<-data.frame(a=NULL,b=NULL)}
if(nrow(locations1)>0){
locations2<-locations1[!is.na(locations1$longitude),]
sortx<-paste(locations2$longitude,locations2$latitude,sep="")
locations3<-locations2[!duplicated(sortx),]
names(locations3)[2:3]<-c("lon","lat")
coordinates(locations3)<-~lon+lat
ext2<-extent(template);sens=2
}}

if(is.null(locations3)){zzz=TRUE}else{if(nrow(locations3)<25){zzz=TRUE}else{zzz=FALSE}}

#if(zzz==TRUE){
#print("switched to globe")
#locations1<-tryCatch(occ2df(occ(query=key, limit=50000)),error=function(e) e)
#if(class(locations1)[1]=="simpleError"){locations1<-data.frame(a=NULL,b=NULL)}
#if(nrow(locations1)>0){
#locations2<-locations1[!is.na(locations1$longitude),]
#sortx<-paste(locations2$longitude,locations2$latitude,sep="")
#locations3<-locations2[!duplicated(sortx),]
#names(locations3)[2:3]<-c("lon","lat")
#coordinates(locations3)<-~lon+lat
#ext2<-extent(template);sens=3
#}}
#if(is.null(locations3)){zzz=TRUE}else{if(nrow(locations3)<25){zzz=TRUE}else{zzz=FALSE}}

if(zzz==TRUE){
key<-keyX
print("switched to genus")
ext1<-geog1$ext2
locations1<-tryCatch(occ2df(occ(query=key, limit=50000,geometry=c(ext1[1],ext1[3],ext1[2],ext1[4]))),error=function(e) e)
if(class(locations1)[1]=="simpleError"){locations1<-data.frame(a=NULL,b=NULL)}
if(nrow(locations1)>0){
locations2<-locations1[!is.na(locations1$longitude),]
sortx<-paste(locations2$longitude,locations2$latitude,sep="")
locations3<-locations2[!duplicated(sortx),]
names(locations3)[2:3]<-c("lon","lat")
coordinates(locations3)<-~lon+lat
ext2<-extent(template);sens=2
}}

if(is.null(locations3)){zzz=TRUE}else{if(nrow(locations3)<25){zzz=TRUE}else{zzz=FALSE}}

if(zzz==TRUE){
print("switched to genus_region")
ext1<-geog1$ext3
locations1<-tryCatch(occ2df(occ(query=key, limit=50000,geometry=c(ext1[1],ext1[3],ext1[2],ext1[4]))),error=function(e) e)
if(class(locations1)[1]=="simpleError"){locations1<-data.frame(a=NULL,b=NULL)}
if(nrow(locations1)>0){
locations2<-locations1[!is.na(locations1$longitude),]
sortx<-paste(locations2$longitude,locations2$latitude,sep="")
locations3<-locations2[!duplicated(sortx),]
names(locations3)[2:3]<-c("lon","lat")
coordinates(locations3)<-~lon+lat
ext2<-extent(template);sens=2
}}

#(ids <- get_ids(names="Lutzomyia", db = 'gbif'))
#occ(ids = ids, from='gbif', limit=20)
#get_gbifid("Satyrium", rank = "genus")

if(is.null(locations3)){return(NA)}

##sort points
projection(locations3)<-projection(regions)
locations2a<-raster::extract(geog1$r1,locations3) #overlay GBIF locations on wrld_simpl
data1<-locations3[!is.na(locations2a),] # get rid of in sea or miles away
print(length(data1))

#### not enought GBIF records ## return just continents with value of 1??<-CHECK OUT
if(length(data1)<15){return(NA)}

#### IDEALS FOR DATA ###
###crop by subregion within continent if only within one sub region
#data2a<-extract(geog1$sr1,data1)
#data2<-data1[!is.na(data2a),]
#if(!is.null(data2)){if(nrow(data2)>=20){data1<-data2}}

#work out precision
decimalnumcount<-function(x){stopifnot(class(x)=="character"); x<-gsub("(.*)(\\.)|([0]*$)","",x); nchar(x) } 
data1$mincert<-pmin(decimalnumcount(as.character(data1$lon)),decimalnumcount(as.character(data1$lat)))
data2<-data1[data1$mincert>2,]
if(!is.null(data2)){if(nrow(data2)>=100){data1<-data2}}

#work out precision
data2<-data1[data1$mincert>3,]
if(!is.null(data2)){if(nrow(data2)>=100){data1<-data2}}

#work out precision
data2<-data1[data1$mincert>4,]
if(!is.null(data2)){if(nrow(data2)>=100){data1<-data2}}

# reduce over sampling
r1<-raster(data1);res(r1)<-1
data1$cell<-cellFromXY(r1,coordinates(data1))
data2<-data1[!duplicated(data1$cell),]
if(!is.null(data2)){if(nrow(data2)>=100){data1<-data2}}

return(list(data1=data1,key=key,ext2=ext2,sens=sens))

}## end of process species


crop_predictors<-function(geog1,regions,predicts){

if(!geog1$ext3==extent(regions)){
	cl <- makeCluster(6)
	registerDoParallel(cl)
	predicts2<-foreach(j=0:7,.combine=raster::stack,.packages=c("rgdal","raster")) %dopar% {
	### CROP PREDICTORS masked predictors & hostR
	k=((j*5)+1):((j*5)+5);k=k[!k==40]
	ppp<-raster::mask(raster::subset(predicts,k),geog1$r1)
	return(raster::crop(ppp,geog1$ext3))
	}
stopCluster(cl)
gc()
	} else {predicts2<-predicts}

####process cropped predictors
len<-length(names(predicts2))
predictors<-subset(predicts2,1:(len/3))
predictors2005<-subset(predicts2,1:(len/3))
predictors2050<-subset(predicts2,((len/3)+1):(len*(2/3)))
predictors2070<-subset(predicts2,((len*(2/3))+1):len)
namesP<-c("soil","landcover","altitude","Mean_Diurnal_Range","Max_Temperature_of_Warmest_Month","Min_Temperature_of_Coldest_Month","Temperature_Annual_Range","Mean_Temperature_of_Warmest_Quarter","Mean_Temperature_of_Coldest_Quarter","Annual_Precipitation","Precipitation_of_Wettest_Month","Precipitation_of_Driest_Month","water_dist")
names(predictors)<-namesP
names(predictors2005)<-namesP
names(predictors2050)<-namesP
names(predictors2070)<-namesP
predictors$water_dist[predictors$water_dist<0]<-0
predictors$water_dist<-log(predictors$water_dist+1)
predictors2005$water_dist[predictors2005$water_dist<0]<-0
predictors2005$water_dist<-log(predictors2005$water_dist+1)
predictors2050$water_dist<-predictors2005$water_dist
predictors2070$water_dist<-predictors2005$water_dist

rm(predicts2);
gc()

return(list(predictors=predictors,predictors2005=predictors2005,predictors2050=predictors2050,predictors2070=predictors2070))
}#end of crop predictors


############## process species #################
do_MAX<-function(datax,n,predictors,background_n,vary_back=FALSE,rnd_mask,cluster_n=3){

bg_n<-background_n

cl <- makeCluster(cluster_n)
registerDoParallel(cl)

  
xms<-foreach(j=1:n,.combine=list,.packages=c("sp","dismo","rgdal","raster")) %dopar% {

	## test number of background points
	if(vary_back==TRUE){ background_n<-round(abs(sample(rnorm(1000,bg_n,bg_n/2),1)),0) +1}

	##set number of points
	background_n2=nrow(datax)
	if(background_n2<background_n){background_n2<-background_n}


	##dvide data into random groups
	group <- dismo::kfold(datax, 5)
	pres_train <- cbind(datax@data,sp::coordinates(datax))[group != 1, ]
	pres_test <- cbind(datax@data,sp::coordinates(datax))[group == 1, ]

	sp::coordinates(pres_train)<-~lon+lat
	sp::coordinates(pres_test)<-~lon+lat

	#### create background set with 500 or more random points
	backg2<-as.data.frame(dismo::randomPoints(rnd_mask,p=datax,n=background_n2,tryf=15,excludep=T))
	#if(is.null)
	coordinates (backg2)<-~x+y	
	backg<-backg2[!is.na(raster::extract(raster::subset(predictors,1),backg2)),]

	#colnames(backg) = c('lon', 'lat')
	group <- dismo::kfold(backg, 5)
	backg_train <- backg[group != 1, ]
	backg_test <- backg[group == 1, ]

	### ADD IN DIAGNOSTICS AUC etc.

	argsX<-c("nolinear","noquadratic","noproduct","nothreshold")
	compl1<-sample(1:(length(argsX)-1),1)
	args1X<-sample(argsX,compl1,replace=F)
	args1X<-args1X[order(args1X)]
	args1X2<-c("noautofeature",args1X)
	beta1<-sample(1:10,1)
	betas<-paste("betamultiplier=",beta1,sep="")
	args1<-c(paste("betamultiplier=",beta1,sep=""),args1X2)

	if(nrow(pres_train)>250){
		xm <- dismo::maxent(predictors, removeDuplicates=T, p=coordinates(pres_train),a=backg_train,args=args1)
		}else{
		xm <- dismo::maxent(predictors, removeDuplicates=F, p=coordinates(pres_train),a=backg_train,args=args1)
		}
	#xms[z]<-xm

	e3 <- dismo::evaluate(pres_test, backg_test, xm, predictors)
	#auc1[z]<-e3@auc
	return(list(xm=xm,auc=e3@auc))
	#print(paste(z," - ",auc1[z],sep=""))
	}

stopCluster(cl)
#auc1<-unlist(unlist(unlist(xms))[seq(from=2,to=length(unlist(unlist(xms))),by=2)])
#return((unlist(unlist(xms)))[(seq(from=1,to=length(unlist(unlist(xms))),by=2))[auc1==max(auc1,na.rm=TRUE)]]   )
return(xms)
} ##end of domax


#seq1<-rep(1:nlayers(predictors),numpred)
#for (x in 1:numpred){
#pb2<-foreach(j=1:length(xms),.combine=raster::stack,.packages=c("gbm","dismo","rgdal","raster")) %dopar% {



