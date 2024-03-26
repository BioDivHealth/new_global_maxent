 install.packages("ineq")
 
 
 library("ineq")
 library(splines)
 
 gini=function(x){
    n=length(x)
    mu=mean(x)
    g=2/(n*(n-1)*mu)*sum((1:n)*sort(x))-(n+1)/(n-1)
    return(g)}
 
 
plot(Lc.lognorm,param=1.5,col="red")
tt<-lines(Lc.lognorm,param=1.2,col="blue")
lines(Lc.lognorm,param=.8,col="green")

library(ggplot2)



nums<-seq(from=0,to=2,by=0.01)

GDP=50000
pop1=1000

res2<-NULL
for(i in nums){

  if(i==0){next}
  
size1=GDP/pop1
  
##simualate and test
incomes1<-Lc.lognorm(seq(from=0, to=1,by=1/pop1),i)

income1<-incomes1*GDP

incomes2<-c(income1[2:length(income1)],income1[length(income1)])-income1

incomes2<-incomes2[1:(length(income1)-1)]

g1<-gini(incomes2)

#length(incomes2[incomes2<=size1])/length(incomes2)
#plot(incomes2,main=round(g1,2))
#abline(h=size1)

#plot(seq(from=0, to=1,by=size1),incomes1,main=round(g1,2))
#lines(seq(from=0, to=1,by=size1),seq(from=0, to=1,by=size1))

  seq2<-size1/2#seq(from=0, to=size1,by=size1/100)
  jj<-rep(0,1)
  for(j in 1:1){
    
    jj[j]<-length(incomes2[incomes2<=seq2[j]])/length(incomes2)
    
    
  }

  res1<-data.frame(gini=g1,income=seq2,prop=jj,GDP=GDP,pop=pop1)

  if(is.null(res2)){res2<-res1} else { res2<-rbind(res2,res1)}
  print(i)
} 
    
#head(res2)

ggplot(res2[,],aes(x=gini,y=prop))+
  geom_point()

plot(x=res2$gini,y=res2$prop)

#lm1<-(lm(data=res2,prop~-1+(gini+I(gini^2)+I(gini^3)+gini:I(gini^3)+I(gini^2):I(gini^3))))
lm1<-(lm(data=res2,prop~-1+(gini+I(gini^2)+gini:I(gini^2))))

#lm1<-(lm(data=res2,prop~gini))

summary(lm1)
points(x=res2$gini,y=lm1$fitted,col="red")

save(lm1, file="C:/Users/david.redding/Dropbox/New_Global_MAXENT/fitted_gini_50th.r")


stats::predict.lm(newdata=data.frame(gini=0.99),lm1)

