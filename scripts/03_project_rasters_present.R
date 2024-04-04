

#setwd('/scratch/scratch/ucbtdw0/sdms')

suppressMessages(library(dismo))
suppressMessages(library(raster))
suppressMessages(library(foreach))
#suppressMessages(library(doMC))
suppressMessages(library(rJava))


spdata<-stack("V://temp//predictorsX.tif")
names(spdata)<-read.csv("V://temp//names1.csv",stringsAsFactors=FALSE)$x


files1<-list.files("V://temp//",pattern="maxent",recursive=TRUE,full.names=TRUE)
files2<-list.files("V://temp//",pattern="maxent",recursive=TRUE,full.names=FALSE)

files2<-gsub("disease_analyses/","",files2,fixed=TRUE)
files2<-gsub("disease_analyses2/","",files2,fixed=TRUE)


files1<-files1[!duplicated(files2)]
files2<-files2[!duplicated(files2)]

predictors<-spdata

#files1<-files1[files2==paste0("apodemus peninsulae_apodemus peninsulae_1_subs_30_161_","maxent.r")]
#files2<-files2[files2==paste0("apodemus peninsulae_apodemus peninsulae_1_subs_30_161_","maxent.r")]


files1<-files1[files2 %in% files2x]
files2<-files2[files2 %in% files2x]


############## PROCESS SPECIES #################

endCluster()


link3b<-link3[link3$year==2010,]

link3b<-link3b[!duplicated(link3b$filen),]

probs1<-NULL
jj<-1
for (ii in 1:nrow(link3b)){
  
  ###future climate data
  clims<-raster(link3b$filen[ii])
  clims2<-tryCatch(crop(clims,template2),error=function(e) e)
  if((class(clims2)[1]=="simpleError")==TRUE){probs1[jj]<-link3b$filen[ii];jj=jj+1}
  print(ii)
}

files2b<-gsub("X:/DRtemp/resultsY/","",probs1,fixed=TRUE)
files2b<-gsub("_present_XXX.tif","_maxent.r",files2b,fixed=TRUE)

files1<-files1[files2 %in% files2b]
files2<-files2[files2 %in% files2b]
