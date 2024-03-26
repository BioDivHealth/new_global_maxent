

files<-list.files("X:\\DRtemp\\temp\\points\\", pattern="points", full.names=TRUE)

points2<-NULL
for(i in 1:length(files)){
  
  load(files[i])
  dt1<-data.frame(coordinates(data1$data1),data1$data1@data,stringsAsFactors = FALSE)
  dt1$filen<-files[i]
  if(!"date" %in% names(dt1)){dt1<-data.frame(dt1[,c("lon","lat","name","prov")],date=NA,dt1[,c("key","mincert","cell","filen")],stringsAsFactors = FALSE)}
  if(!"key" %in% names(dt1)){dt1<-data.frame(dt1[,c("lon","lat","name","prov","date")],key=NA,dt1[,c("mincert","cell","filen")],stringsAsFactors = FALSE)}
  
  if(is.null(points2)){points2<-dt1}else{points2<-rbind(points2,dt1)}
  print(i)
}

points2$speciesX<-points2$filen
points2$speciesX<-gsub("X:\\DRtemp\\temp\\points\\","",points2$speciesX, fixed=TRUE)
points2$speciesX<-gsub("_points.r","",points2$speciesX, fixed=TRUE)
points2$speciesX<-gsub("_1_subs_",";",points2$speciesX, fixed=TRUE)
points2$speciesX<-gsub("_2_subs_",";",points2$speciesX, fixed=TRUE)
points2$speciesX<-gsub("_3_subs_",";",points2$speciesX, fixed=TRUE)
points2$speciesX<-gsub("_4_subs_",";",points2$speciesX, fixed=TRUE)
#others1<-read.table(text=points2$speciesX,sep=";",stringsAsFactors=FALSE)
#names(others1)<-c("species2","region")
points2<-setDT(points2)
points2[,c("species2","region"):=tstrsplit(speciesX,split=";",fixed=TRUE)]

fwrite(points2,file="X:\\DRtemp\\All_points_combined.csv")
