

###This removes the columns with greater then the threshold level of within dataframe autocorrelation
###Removes the variable with the highest average correlation to other variables or random if the same
###Does not work for factors

autocor<-function(dataf,threshold=0.75,cor.method="spearman"){

corm<-cor(dataf,method=cor.method,use="pairwise.complete.obs")

corm<-abs(corm)

diag(corm)<-NA

mv<-max(corm,na.rm=T)

while(mv>threshold & ncol(corm)>1){
	
	cnam<-colnames(corm)

	topp<-cnam[apply(corm,2,function(x) max(x,na.rm=T))==mv]

	cm<-colMeans(corm[,topp],na.rm=T)

	lose1<-topp[cm==max(cm)]

	if(length(lose1)>1){lose1<-lose1[sample(1:length(lose1),1)]}

	corm<-corm[,-((1:length(cnam))[cnam==lose1])]	

	corm<-corm[-((1:length(cnam))[cnam==lose1]),]	

	mv<-max(corm,na.rm=T)

}



return(dataf[,colnames(corm)])

}