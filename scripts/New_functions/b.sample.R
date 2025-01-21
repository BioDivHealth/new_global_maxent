#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Selection of background data for the processing of SDM  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
b_sample<-function(x,start=20,plot.l=FALSE){
  
  # Get the possible range of ovservations----
  warning("Variables with homogeneous distributions might influence the results of the background sample selection!")
  
  if(ncol(x)<2){
    stop("At least two variables are needed to compare distributions!")
  }
  
  # Create the sequence of testing values
  n_bk_test <- seq(from=min,to=max,length.out=breaks.r)
  n_bk_test <- round(n_bk_test,digits=0)
  
  # Run the loop over the different objects
  Ov.predictors<-list()
  
  # Use here in while() to increase the number of background data until we reach a certain 
  # level of density distribution overlapping. In our case we are going to set this to 75%
  
  ov_val<-0
  ov_par<-list()
  
  s.size<-start
  
  for(ff in 1:ncol(x)){
    
    var1 <- 
    
    while(ov_val<75 & s.size<30000){
     
       ov_val <- D.overlap(x=list(Sutyd_area=var1,
                       Sample=var1[sample(1:length(var1),size=s.size)]),
                mean.d = F, plot.d = F)$OV*100
    
        Ov.predictors[[ff]] <- do.call("c",overlap)
        
        s.size <- s.size + 2000
      }
  print(ff)
  
  }  
    
  for(ff in 1:ncol(x)){
      var1 <- x[,ff]
      # var1 <- var1[!is.na(var1)] #; var1<-sqrt(var1)
    
      overlap<-lapply(n_bk_test, function(v){D.overlap(x=list(Sutyd_area=var1,Sample=var1[sample(1:length(var1),size=v)]),
                                                     mean.d = F, plot.d = F)$OV*100})
    
      Ov.predictors[[ff]] <- do.call("c",overlap)
      print(ff)
  }
  
  ov.vals<-do.call("cbind",Ov.predictors)
  colnames(ov.vals)<-names(x)
  
  if(plot.l==T){
  par(bg="white")
  
  plot(x=n_bk_test,y=ov.vals[,1],xlab="Number sampled points",ylab="Distribution overlap",
       axes=T,col=NA,ylim=range(ov.vals),bty="n")
  
  apply(ov.vals,2,function(p) lines(n_bk_test,p,type="b",pch=19))
  }
  # As a rule of thumb we are going to select a threshold of 70% of coverage as sufficient for the modelling
  ov.vals <-cbind(n_bk=n_bk_test,mean=rowMeans(ov.vals),ov.vals)
  
  if(TRUE %in% c(ov.vals[,"mean"] > 75)){
    n_bk <- ov.vals[ov.vals[,"mean"] > 75,"n_bk"][1]
  } else {
    
    n_bk <- ov.vals[order(ov.vals[,"mean"]),"n_bk"][1]
  }
  return(list(sample.size.bk=n_bk,ov.data=ov.vals))  
  
  }
  
#
# x <- matrix(rnorm(1000*1000,1,sd=3.5),ncol=100,nrow=100,byrow=TRUE)
# y <- matrix(dpois(1:c(1000*1000),1,lambda=3),ncol=100,nrow=100)
# p <- matrix(dgamma((1:10000)/100000,shape=6),ncol=100,nrow=100,byrow=TRUE)
# 
# # x <- lapply(list(x,y,z,p),terra::rast)
# # x <- do.call("c",x)
# 
# x <- list(x,y,p)
# b_sample(x=x,n.rec = 1000,max.cap = 0.75,breaks=100,plot.l = T)
# 
# plot(rast(lapply(x,rast)))

