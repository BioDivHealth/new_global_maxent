


##need to redo name so powassan is the same and make decisions accross all of them
coefs_area<-res14d[AUC_best.stat>0.65 & thresh==1,list(intercept=coef(lm(log(area+1)~year))[1], coef=coef(lm(log(area+1)~year))[2],prob=(summary(lm(area~year)))$coefficients[2,4]),by=name]

hist(coefs_area$coef[coefs_area$prob<0.05])

coefs_hazard<-res14d[AUC_best.stat>0.65,list(intercept=coef(lm(Hazard_future~year))[1], coef=coef(lm(Hazard_future~year))[2],prob=(summary(lm(Hazard_future~year)))$coefficients[2,4]),by=name]

setorder(coefs_hazard,coef)

hist(coefs_hazard$coef[coefs_hazard$prob<0.05])

coefs_risk<-res14d[AUC_best.stat>0.65,list(intercept=coef(lm(RiskF~year))[1], coef=coef(lm(RiskF~year))[2],prob=(summary(lm(RiskF~year)))$coefficients[2,4]),by=name]

hist(coefs_risk$coef[coefs_risk$prob<0.05])

setorder(coefs_risk,coef)




