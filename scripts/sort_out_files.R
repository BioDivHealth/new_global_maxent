
library(raster)

#"V:/resultsY2/aedes_ochlerotatus_2_subs_5_11_14_15_17_18_29_35_XXX_2030_2.6.tif"
"V:/resultsY2/aedes_ochlerotatus_2_subs_5_11_14_15_17_18_29_35_XXX_2030_2.6.tif"

load(file="D:\\fin1.r")

fin1<-fin1[sample(1:nrow(fin1),nrow(fin1)),]

for(i in nrow(fin1):1){
  
  print(i)
  
  filen<-fin1[i,"filen"]
  filen<-gsub("V:","D:",filen,fixed=TRUE)
  
  if(filen %in% c(list.files("D:/resultsY2",full.names=TRUE),list.files("D:/resultsY",full.names=TRUE))){next}
  
  e <- simpleError("test error")
  rrr<-tryCatch(raster(fin1[i,"filen"]),error = function(e) e)
  
  if(class(rrr)[1]=="simpleError"){print("problem reading");next}
 
 rrr<-tryCatch(writeRaster(rrr,file=filen,overwrite=FALSE),error = function(e) e)
  if(class(rrr)[1]=="simpleError"){print("problem writing");next}
  
  
}
