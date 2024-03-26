

yy<-list.files("./resultsY/",pattern="XXX.tif",full.names=TRUE)

yy<-gsub("XXX.tif","XXX.r",yy)
yy<-gsub("resultsY","resultsX",yy)


for(i in 1:length(yy)){

	tt<-NULL
	save(tt,file=yy[i])
print(i)
}
