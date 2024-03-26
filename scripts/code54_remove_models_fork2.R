
require(raster)
require(rgdal)
require(dismo)
library(data.table)
library(sp)
library(taxize)
library(maptools)
library(doParallel)


e <- simpleError("test error")

####Get all future climate niches
link2b<-list.files("X:/DRtemp/resultsY", pattern=".tif",full.names=TRUE)
link2<-list.files("X:/DRtemp/resultsY", pattern=".tif",full.names=FALSE)
link2<-gsub("_ncc_",";",link2, fixed=TRUE)
link2<-gsub("_mri_",";",link2, fixed=TRUE)
link2<-gsub("_bcc_",";",link2, fixed=TRUE)
link2<-gsub("_nimr_",";",link2, fixed=TRUE)
link2<-gsub("_cccma_",";",link2, fixed=TRUE)
link2<-gsub("_cesm1_",";",link2, fixed=TRUE)
link2<-gsub("_csiro_",";",link2, fixed=TRUE)
link2<-gsub("_ec_",";",link2, fixed=TRUE)
link2<-gsub("_fio_",";",link2, fixed=TRUE)
link2<-gsub("_gfdl_",";",link2, fixed=TRUE)
link2<-gsub("_bnu_",";",link2, fixed=TRUE)
link2<-gsub("_giss_",";",link2, fixed=TRUE)
link2<-gsub("_ipsl_",";",link2, fixed=TRUE)
link2<-gsub("_inm_",";",link2, fixed=TRUE)
link2<-gsub("_miroc_",";",link2, fixed=TRUE)
link2<-gsub("_lasg_",";",link2, fixed=TRUE)
link2<-gsub("_ncar_",";",link2, fixed=TRUE)
link2<-gsub("_mpi_",";",link2, fixed=TRUE)
link2<-gsub("_mohr_",";",link2, fixed=TRUE)
link2<-gsub("_mri_",";",link2, fixed=TRUE)
link2<-gsub("_mohc_",";",link2, fixed=TRUE)
link2<-gsub("_present_",";present;2010;",link2, fixed=TRUE) ##check when in
link2<-gsub("_XXX.tif","",link2, fixed=TRUE)
link2<-gsub("_2.6",";2.6",link2, fixed=TRUE)
link2<-gsub("_4.5",";4.5",link2, fixed=TRUE)
link2<-gsub("_6.0",";6.0",link2, fixed=TRUE)
link2<-gsub("_8.5",";8.5",link2, fixed=TRUE)
link2<-gsub("_2030",";2030",link2, fixed=TRUE)
link2<-gsub("_2050",";2050",link2, fixed=TRUE)
link2<-gsub("_2070",";2070",link2, fixed=TRUE)
link2<-gsub("_2080",";2080",link2, fixed=TRUE)
link2<-gsub("-astOJIRBNlQ0","",link2, fixed=TRUE)
link2<-gsub("_points.r","",link2, fixed=TRUE)
link3<-read.table(text=link2,sep=";",stringsAsFactors=FALSE)
link3$filen<-link2b
names(link3)[1:4]<-c("species","model","year","RCP")
link3$RCP<-as.numeric(gsub("XXX.tif","999",link3$RCP))
###make other species regions available
link3$speciesX<-link3$species
link3$speciesX<-gsub("_1_subs_",";",link3$speciesX, fixed=TRUE)
link3$speciesX<-gsub("_2_subs_",";",link3$speciesX, fixed=TRUE)
link3$speciesX<-gsub("_3_subs_",";",link3$speciesX, fixed=TRUE)
link3$speciesX<-gsub("_4_subs_",";",link3$speciesX, fixed=TRUE)
others1<-read.table(text=link3$speciesX,sep=";",stringsAsFactors=FALSE)
names(others1)<-c("species2","region")
link3<-cbind(link3,others1)
rm(link2,link2b)

link3$uni1<-paste(link3$year,link3$RCP,sep="_")

sss<-unique(link3$species)

#beginCluster(8)

##randomise
sss<-sample(sss,length(sss),replace=FALSE)

for (hh in 1:length(sss)){
  
  l3<-link3[link3$year>2010 & link3$species==sss[hh],]
  sss2<-unique(l3$uni1)
  
           for (jj in 1:length(sss2)) {
            
             l4<-l3[l3$uni1==sss2[jj],]
             
             if(paste("V:\\ResultsY2\\",l4$species[1],"_XXX_",sss2[jj],".tif",sep="") %in% list.files("V:\\ResultsY2\\",full.names=TRUE)){next}
             
              
          tt<-stack(l4[,"filen"])
          
          if(nrow(l4)==1){
            
            names(tt)<-paste(l4$species[1],"_XXX_",sss2[jj],sep="")
            writeRaster(tt,format="GTiff",file=paste("V:\\ResultsY2\\",l4$species[1],"_XXX_",sss2[jj],".tif",sep=""))
            next
           }        
          
          #ras.mean <- clusterR(tt, calc, args=list(mean, na.rm=T))
          ras.mean <- tryCatch(calc(tt,mean, na.rm=TRUE),error=function(e) e)
          
          if(class(ras.mean)[1]=="simpleError"){
            
            
            for(ik in 1:nlayers(tt)){
              
              
              
              rasmean <- tryCatch(plot(subset(tt,ik)),error=function(e) e)
              
              if(class(rasmean)[1]=="simpleError"){break}
                
                  }
            
            tt<-subset(tt,(1:nlayers(tt))[(1:nlayers(tt)) !=ik])
            
            ras.mean <- calc(tt,mean, na.rm=TRUE)
            
          }
          
          
          names(ras.mean)<-paste(l4$species[1],"_XXX_",sss2[jj],sep="")
          
          #ras.sd <- clusterR(tt, calc, args=list(sd, na.rm=T))
  
          writeRaster(ras.mean,format="GTiff",file=paste("V:\\ResultsY2\\",l4$species[1],"_XXX_",sss2[jj],".tif",sep=""))
          rm(l4)
           }
  print(paste(hh," of ",length(sss),sep=""))
  
}

#endCluster()

