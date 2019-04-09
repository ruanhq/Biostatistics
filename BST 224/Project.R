##############################################################
########################
#Project for STA 224(Longitudinal data analysis): Here we explore the and do corresponding diagnosis
#as well as comparing different possible models. 
#We compare GLMM and GEE on the skin cancer data
###############
#load data and required packages:
library(foreign)
library(geepack)
library(lme4)
library(MASS)
library(lattice)
library(VIM)
library(mi)
skin_cancer<-read.dta("skin.dta")

#Here we restrict the analysis in the first group:
skin<-skin_cancer[skin_cancer$center==1,]
##Exploratory analysis:
#Group mean and population(global) mean:
#Transform the data into wide format:
skin_wide<-reshape(skin,v.names="y",idvar="id",timevar="year",direction="wide")
response_placebo<-skin_wide[skin_wide[,7]==0,8:12]
response_beta_caterone<-skin_wide[skin_wide[,7]==1,8:12]
#calculate the group mean after removing the outliers:
population_mean<-apply(skin_wide[,8:12],2,mean,na.rm=TRUE)
placebo_mean<-apply(response_placebo,2,mean,na.rm=TRUE)
beta_caterone_mean<-apply(response_beta_caterone,2,mean,na.rm=TRUE)
#plot the mean:
png('Group_Population_Means.png',width=600,height=600)
plot(population_mean,type="l",ylim=c(0,0.5),lty=1,lwd=0.89,xlab="Year",ylab="Mean",main="Group and Population Means",)
lines(placebo_mean,col="blue",lty=2,lwd=1.20)
lines(beta_caterone_mean,col="red",lty=3,lwd=1.10)
legend("topright",legend=c("Population Mean","Placebo Mean","Beta_Caterone Mean"),col=c("black","blue","red"),lty=c(1,2,3),cex=1.2)
dev.off()
################
#GEE:
model_gee1<-geeglm(y~year+trt+I(trt*year)+age+skin+gender+exposure+I(trt*age)+I(trt*skin)+I(trt*exposure)+I(trt*gender)
	+I(year^2)+I(year^3),id=id,data=skin,family=poisson("log"),corstr="unstructured")
model_gee2<-geeglm(y~trt+exposure+I(trt*age)+I(trt*gender),id=id,data=skin,family=poisson("log"),corstr="unstructured")
model_gee3<-geeglm(y~trt+exposure,id=id,data=skin,family=poisson("log"),corstr="unstructured")
model_gee4<-geeglm(y~trt+exposure+I(trt*age),id=id,data=skin,family=poisson("log"),corstr="unstructured")
anova(model_gee1,model_gee2)
anova(model_gee2,model_gee4)
anova(model_gee4,model_gee3)
model1_gee<-geeglm(y~trt+exposure+I(trt*age)+I(trt*gender),id=id,data=skin,family=poisson("log"),corstr="unstructured")
model2_gee<-geeglm(y~trt+exposure+I(trt*age)+I(trt*gender),id=id,data=skin,family=poisson("log"),corstr="ar1")
model3_gee<-geeglm(y~trt+exposure+I(trt*age)+I(trt*gender),id=id,data=skin,family=poisson("log"),corstr="exchangeable")
model4_gee<-geeglm(y~trt+exposure+I(trt*age)+I(trt*gender),id=id,data=skin,family=poisson("log"),corstr="independence")
#select the best model by the difference between the naive  Correlation estimation and the robust correlation estimation:
dif1<-sum(abs(model1_gee[["geese"]]$vbeta-model1_gee[["geese"]]$vbeta.naiv))
dif2<-sum(abs(model2_gee[["geese"]]$vbeta-model2_gee[["geese"]]$vbeta.naiv))
dif3<-sum(abs(model3_gee[["geese"]]$vbeta-model3_gee[["geese"]]$vbeta.naiv))
dif4<-sum(abs(model4_gee[["geese"]]$vbeta-model4_gee[["geese"]]$vbeta.naiv))
#So here we use the unstructured covariance structure:



#Diagnostic w.r.t GEE and identification of outliers:
#Calculate the pearson residuals:
fitted_value<-fitted(model1_gee)
#Pearson residuals:
res_0<-skin$y-fitted_value
resid<-(skin$y-fitted_value)/sqrt(1.21*fitted_value)
#Perform chisquare test:
sum(resid^2)
#Residual plot:
png('Residual_Diagnostics.png',width=1000,height=700)
par(mfrow=c(1,2))
plot(fitted_value,resid,xlab="Fitted Value",ylab="Pearson Residuals")
lines(smooth.spline(fitted_value,resid,spar=1.3))
scatter.smooth(skin$year,resid,xlab="Year",ylab="Pearson Residuals")
dev.off()
#Calculate and plot the mahalanobis distance for each one:
#transform it to wide:
skin_wide<-reshape(skin,v.names="y",idvar="id",timevar="year",direction="wide")
respons<-skin_wide[,8:12]
nn<-apply((1-is.na(respons)),1,sum)
m<-length(nn)
pval<-rep(0,m)
d<-rep(0,m)
cormat<-VarCorr(model1_gee)
for(i in 1:m){
	end<-sum(nn[1:i])
	start<-end-nn[i]+1
	#C is the estimated correlation structure
	na_index<-which(is.na(skin_wide[i,8:12]))
	C<-cormat[is.na(pmatch(c(1:5),na_index)),is.na(pmatch(c(1:5),na_index))]
	if(start!=end){
		A<-diag(sqrt(1.27*fitted_value[start:end]))
	}
	else{
		A<-sqrt(1.27*fitted_value[start:end])
	}
	V=A%*%C%*%A
	d[i]<-t(as.matrix(res_0[start:end]))%*%solve(V)%*%as.matrix(res_0[start:end])
	pval[i]<-pchisq(d[i],nn[i],lower.tail=FALSE)
}
png('Mahalanobis.png',height=600,width=600)
plot(d,cex=0.9,main="Mahalanobis Distance for GEE",xlab="index",ylab="Mahalanobis Distance",pch=5)
points(order(pval)[1],max(d),col="red")
dev.off()
#Identify the largest mahalanobis distance:
index_outlier<-which.max(d) #241
fitted_gee<-fitted(model1_gee)
end_point<-sum(nn[1:index_outlier])
start_point<-end_point-nn[index_outlier]+1
#Plot the trajectory of this object with the global mean and beta-carotene group mean:
png('Outlier_Visualization.png',height=600,width=600)
plot(population_mean,type="l",ylim=c(0,7.1),lty=1,lwd=0.89,xlab="Year",ylab="Mean",main="Outlier Visualization")
lines(beta_caterone_mean,col="red",lty=3,lwd=1.10)
lines(c(1:5),skin_wide[index_outlier,8:12],col="green",lty=4,lwd=1.85)
lines(fitted_gee[start_point:end_point],col="blue",lty=5,lwd=1.41)
legend("topright",legend=c("Population Mean","Beta_Caterone Mean","Outlier Trajectory","Fitted Trajectory"),col=c("black","red","green","blue"),lty=c(1,3,4,5))
dev.off()
#Then we remove the outlier object and fit the GEE again:
skin_wide<-skin_wide[-241,]
skin<-skin[-(start_point:end_point),]
model5_gee<-geeglm(y~year+trt+I(trt*year)+age+skin+gender+exposure+I(trt*age)+I(trt*skin)+I(trt*exposure)+I(trt*gender)
	+I(year^2)+I(year^3),id=id,data=skin,family=poisson("log"),corstr="unstructured")
summary(model5_gee)
model6_gee<-geeglm(y~trt+exposure+I(trt*gender)+I(trt*age),id=id,data=skin,family=poisson("log"),corstr="unstructured")
anova(model6_gee,model5_gee)

model7_gee<-geeglm(y~trt+exposure+I(trt*gender)+I(trt*age),id=id,data=skin,family=poisson("log"),corstr="ar1")
model8_gee<-geeglm(y~trt+exposure+I(trt*gender)+I(trt*age),id=id,data=skin,family=poisson("log"),corstr="exchangeable")
model9_gee<-geeglm(y~trt+exposure+I(trt*gender)+I(trt*age),id=id,data=skin,family=poisson("log"),corstr="independent")
dif1<-sum(abs(model6_gee[["geese"]]$vbeta-model1_gee[["geese"]]$vbeta.naiv))
dif2<-sum(abs(model7_gee[["geese"]]$vbeta-model2_gee[["geese"]]$vbeta.naiv))
dif3<-sum(abs(model8_gee[["geese"]]$vbeta-model3_gee[["geese"]]$vbeta.naiv))
dif4<-sum(abs(model9_gee[["geese"]]$vbeta-model4_gee[["geese"]]$vbeta.naiv))
#Choose exchangeable covariance structure:
summary(model8_gee)

################
#GLMM:
model_glmm1<-glmer(y~year+trt+age+skin+gender+exposure+I(trt*year)+(year|id),data=skin,family=poisson("log"))
model_glmm2<-glmer(y~age+gender+exposure+(year|id),data=skin,family=poisson("log"))
model_glmm3<-glmer(y~year+age+gender+exposure+(1|id),data=skin,family=poisson("log"))
#apply the likelihood ratio test to compare different models:
anova(model_glmm1,model_glmm2,test="Chi")
anova(model_glmm2,model_glmm3,test="Chi")
#So model_glmm2 is our final model for GLMM. 
#Diagnostics:
empi_blup<-ranef(model_glmm2)$id
png('Histogram_EBLUP.png',width=800,height=400)
par(mfrow=c(1,2))
hist(empi_blup[,1],main="Histogram for Empirical BLUP for Random Intercept")
hist(empi_blup[,2],main="Histogram for Empirical BLUP for Random Slope")
dev.off()

#First calculate the degree of freedom:
m_random<-VarCorr(model_glmm2)
m_fixed<-fixef(model_glmm2)
df<-3+length(m_fixed)
test_stat<-sum(residuals(model_glmm2,type="pearson")^2)
pvalue_glmm<-pchisq(test_stat,nrow(skin)-df,lower.tail=FALSE)
#Here we have identified that the 52th object is the outlier by the histogram of EBLUP:
outlier<-which.max(empi_blup[,1])
plot(1:5,skin_wide[52,8:12],type="l",xlab="year",ylab="Response(# of new skin cancer)",col="red")
lines(apply(skin_wide[skin_wide$trt==0,8:12],2,mean,na.rm=TRUE))

#####################
#Trajectory Visualization and comparison:
#Here we randomly choose 16 objects randomly and display their trajectory and 
#the fitted curve by GLMM and GEE.
samp<-sample(1:m,16)
fitted_gee<-fitted(model_gee1)
fitted_glmm<-fitted(model_glmm2)
yy<-skin_wide[,8:12]
png('GEE_trajectory_Visualization.png',width=1600,height=1600)
par(mfrow=c(4,4))
for(i in 1:16){
	end<-sum(nn[1:samp[i]])
	start<-end-nn[samp[i]]+1
	ylow<-min(yy[samp[i],(1:nn[samp[i]])])-0.5
	yhigh<-max(yy[samp[i],(1:nn[samp[i]])])+0.5
	plot(fitted_gee[start:end],ylim=c(ylow,yhigh),type="l",xlab="Index of sample",ylab="Response")
	points(1:nn[samp[i]],yy[samp[i],(1:nn[samp[i]])],type="p")
}
dev.off()
png('GLMM_trajectory_Visualization.png',width=800,height=800)
par(mfrow=c(4,4))
for(i in 1:16){
	end<-sum(nn[1:samp[i]])
	start<-end-nn[samp[i]]+1
	ylow<-min(yy[samp[i],(1:nn[samp[i]])])-0.5
	yhigh<-max(yy[samp[i],(1:nn[samp[i]])])+0.5
	plot(fitted_glmm[start:end],ylim=c(ylow,yhigh),type="l",xlab="Index of sample",ylab="Response")
	points(1:nn[samp[i]],yy[samp[i],(1:nn[samp[i]])],type="p")
}
dev.off()
#####################

