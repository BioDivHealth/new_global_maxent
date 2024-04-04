
require(biomod2)

do_auc<-function(data1,test_column,group_column,timesP=1000,timesA=10000, replications=25,summarise=TRUE,weights=NULL){

###start

  data1=as.data.frame(data1)
  
  if(is.null(weights)){data1$dummyx1<-1;weights="dummyx1"}
  
  if(replications==FALSE){replications==1}

  for (i in 1:replications){
    
  ##random absence points
  rand2=data1[,test_column][sample((1:nrow(data1))[is.na(data1[,group_column])],ceiling(timesP),prob=data1[is.na(data1[,group_column]),weights],replace=TRUE)]
  real3=data1[,test_column][sample((1:nrow(data1))[!is.na(data1[,group_column])],ceiling(timesA),prob=data1[!is.na(data1[,group_column]),weights],replace=TRUE)]

  ##comparison by 'group_column'y/admin areas
  #res6d<-rbind(data.frame(ea2=rep("0",length(rand2)),value=rand2),data.frame(ea2=rep("1",length(real3)),value=real3))
  #res6d<-res6d[!is.na(res6d$value),]

  e1<-dismo::evaluate(p=real3,a=rand2)
  thr1<-data.frame(best.stat=e1@auc,cutoff=e1@t[which.max(e1@TPR + e1@TNR)], sensitivity=NA, specificity=NA)
  #thr1<-Find.Optim.Stat(Stat="ROC",Fit=res6d$value,Obs=res6d$ea2)
  colnames(thr1)<-paste("AUC_",colnames(thr1),sep="")

  if(i==1) {thr2<-thr1} else {thr2<-rbind(thr2,thr1)}
  
  }
  
  if(replications==1|replications==FALSE) {return(thr2)}
  
  if(summarise==TRUE) {return(colMeans(thr2,na.rm = TRUE))} else {return(thr2)}

}

#do_auc(data1=res6a,test_column="Hazard",group_column="countr",summarise=FALSE)
