#Project for Survival Analysis?
library(survminer)
library(spatstat)
library(survival)
library(ggplot2)
library(rms)
library(survMisc)
library(zoo)
library(KMsurv)
library(MASS)
library(cluster)
library(muhaz)
library(fitdistrplus)
library(logspline)
library(aftgee)
library(factoextra)
library(magrittr)
library(lss)
library(timereg)



###############
#exploratory analysis and processing
data(larynx)
head(larynx)
names(larynx)
n<-nrow(larynx)
pairs(larynx[,2:4],main="Scatterplot Matrix")
par(mfrow=c(2,2))
hist(larynx$stage,main="Histogram for stages",freq=FALSE)
lines(density(larynx$stage))
hist(larynx$time,main="Histogram for TIME",freq=FALSE)
lines(density(larynx$time))
hist(larynx$age,main="Histogram for Age",freq=FALSE)
lines(density(larynx$age))
hist(larynx$diagyr,main="Histogram for diagnosis Year",freq=FALSE)
lines(density(larynx$diagyr))
treat1<-larynx[1:33,]
treat2<-larynx[34:50,]
treat3<-larynx[51:77,]
treat4<-larynx[78:90,] #Four stages'
par(mfrow=c(2,2))
boxplot(larynx$delta~larynx$stage,main="Delta vs Stage")
boxplot(larynx$time~larynx$stage,main="Time vs Stage")
boxplot(larynx$age~larynx$stage,main="Age vs Stage")
boxplot(larynx$diagyr~larynx$stage,main="diagyr vs stage")
par(mfrow=c(2,2))
boxplot(larynx$age~larynx$diagyr,main="Age vs diagyr")
boxplot(larynx$time~larynx$diagyr,main="Time vs diagyr")
boxplot(larynx$time~larynx$age,main="Time vs Age")
boxplot(larynx$age~larynx$time,main="Age vs Time")
interaction.plot(larynx$age,larynx$stage,larynx$time,main="Interaction plot of Time with Age/stage")
interaction.plot(larynx$diagyr,larynx$stage,larynx$time,main="Interaction plot of Time with diagyr/stage")
####
##Quantile-quantile plot:
par(mfrow=c(1,3))
qqnorm(larynx$age,main="Q-Q plot of Age")
qqline(larynx$age)
qqnorm(larynx$diagyr,main="Q-Q plot of diagyr")
qqline(larynx$diagyr)#
qqnorm(larynx$time,main="Q-Q plot for time")
qqline(larynx$time)
pairs(larynx[,2:4],main="Scatterplot Matrix")
#######
#grouping the continuous covariates(age,diagyr):
quantile(larynx$diagyr)
quantile(larynx$age)
age_cate<-c(1:n)
for(i in 1:n){
	if(larynx[i,3]<=55){
		age_cate[i]=1
	}
	else if(55<larynx[i,3] && larynx[i,3]<=65){
		age_cate[i]=2
	}
	else if(65<larynx[i,3] && larynx[i,3]<=75){
		age_cate[i]=3
	}
	else{
		age_cate[i]=4
	}
}
new_data<-cbind(larynx,factor(age_cate))
diagyr_cate<-c(1:n)
for(i in 1:n){
	if(larynx[i,4]<=72){
		diagyr_cate[i]<-1
	}
	else if(larynx[i,4]<=74 && larynx[i,4]>72){
		diagyr_cate[i]<-2
	}
	else if(larynx[i,4]<=76 && larynx[i,4]>74){
		diagyr_cate[i]<-3
	}
	else{
		diagyr_cate[i]<-4
	}
}
new_data<-cbind(new_data,factor(diagyr_cate))
############
#Preliminary Modeling:
time<-sort(unique(larynx[,2]),decreasing=FALSE)
nn<-length(time)
n_risk<-rep(0,nn)
n_event<-rep(0,nn)
n_risk[1]<-nrow(larynx)
for(i in 2:nn){
	n_risk[i]<-n_risk[i-1]-length(which(larynx[,2]==time[i-1]))
}
for(i in 1:nn){
	index<-larynx[which(larynx[,2]==time[i]),5]
	n_event[i]<-length(which(index==1))
}
N_A<-cumsum(n_event/n_risk)
KM<-rep(1,nn)
for(i in 1:nn){
	for(j in 1:i){
		if(n_event[j]!=0){
			KM[i]<-KM[i]*(1-n_event[j]/n_risk[j])
		}
	}
}
#Interval Estimation via log-log approach:
#Construct the interval estimation by hand
var_KM_L<-rep(0,nn)
for(i in 1:nn){
	di<-0
	for(j in 1:i){
		di<-di+(n_event[j])/((n_risk[j]-n_event[j])*n_risk[j])
	}
	var_KM_L[i]<-di/(log(KM[i])^2)
}
lowci_KM<-KM^(exp(1.96*sqrt(var_KM_L)))
highci_KM<-KM^(exp(-1.96*sqrt(var_KM_L)))
plot(time,KM,type="l",lty=2,main="KM estimator and 95% interval estimation")
lines(time,lowci_KM,type="l",lty=3)
lines(time,highci_KM,type="l",lty=4)
var_NA_L<-cumsum(((n_risk-n_event)*n_event)/((n_risk-1)*n_risk^2))
lowci_NA<-N_A*exp(-1.96*var_NA_L/N_A)
highci_NA<-N_A*exp(1.96*var_NA_L/N_A)
plot(time,N_A,type="l",lty=1,main="Nelson Aalen estimator and 95% interval estimation")
lines(time,lowci_NA,type="l",lty=2)
lines(time,highci_NA,type="l",lty=3)

##################
#Plot the survival curve among different categories:
fit_age<-survfit(Surv(time,delta)~factor(age_cate),data=new_data)
ggsurv<-ggsurvplot(fit_age,data=larynx,risk.table=TRUE,pval=TRUE,conf.int=TRUE,palette=rainbow(4),xlim=c(0,15),xlab="Times",break.time.by=2,ggtheme=theme_light(),risk.table.y.text.col=T,risk.table.height=0.25,risk.table.y.text=FALSE,ncensor.plot=TRUE,ncensor.plot.height=0.25,surv.median.line="hv",legend.labs=c("Age<55","55<Age<=65","65<Age<=75","75<Age<=86"))
ggsurvplot(fit_age,conf.int=TRUE,risk.table.col="strata",ggtheme=theme_bw(),palette=rainbow(4),xlim=c(0,12),fun="event",break.time.by=1,risk.table.y.text=FALSE,legend.labs=c("Age<55","55<Age<=65","65<Age<=75","75<Age<=86"),risk.table.height=0.15,ncensor.plot.height=0.25)
fit_stage<-survfit(Surv(time,delta)~factor(stage),data=new_data)
ggsurv<-ggsurvplot(fit_stage,data=larynx,risk.table=TRUE,pval=TRUE,conf.int=TRUE,palette=rainbow(4),xlim=c(0,15),xlab="Times",break.time.by=2,ggtheme=theme_light(),risk.table.y.text.col=T,risk.table.height=0.25,risk.table.y.text=FALSE,ncensor.plot=TRUE,ncensor.plot.height=0.25,surv.median.line="hv",legend.labs=c("Stage1","Stage2","Stage3","Stage4"))
#Cumulative hazard plot
ggsurvplot(fit_stage,conf.int=TRUE,risk.table.col="strata",ggtheme=theme_bw(),palette=rainbow(4),xlim=c(0,12),fun="event",break.time.by=1,risk.table.y.text=FALSE,legend.labs=c("Stage1","Stage2","Stage3","Stage4"),risk.table.height=0.15,ncensor.plot.height=0.25)
#For diagnosed year
fit_diagyr<-survfit(Surv(time,delta)~factor(diagyr_cate),data=new_data)
ggsurv<-ggsurvplot(fit_diagyr,data=larynx,risk.table=TRUE,pval=TRUE,conf.int=TRUE,palette=rainbow(4),xlim=c(0,15),xlab="Times",break.time.by=2,ggtheme=theme_light(),risk.table.y.text.col=T,risk.table.height=0.25,risk.table.y.text=FALSE,ncensor.plot=TRUE,ncensor.plot.height=0.25,surv.median.line="hv",legend.labs=c("diagyr<=72","72<diagyr<=74","74<diagyr<=76","76<diagyr"))
ggsurvplot(fit_diagyr,conf.int=TRUE,risk.table.col="strata",ggtheme=theme_bw(),palette=rainbow(4),xlim=c(0,12),fun="event",break.time.by=1,risk.table.y.text=FALSE,legend.labs=c("diagyr<=72","72<diagyr<=74","74<diagyr<=76","76<diagyr"),risk.table.height=0.15,ncensor.plot.height=0.25)
#fit the full model:
fit_full<-survfit(Surv(time,delta)~1,data=larynx)
ggsurv<-ggsurvplot(fit_full,data=larynx,risk.table=TRUE,pval=TRUE,conf.int=TRUE,palette=rainbow(4),xlim=c(0,12),xlab="Times",break.time.by=2,ggtheme=theme_light(),risk.table.y.text.col=T,risk.table.height=0.25,risk.table.y.text=FALSE,ncensor.plot=TRUE,ncensor.plot.height=0.25,surv.median.line="hv",legend.labs="Full Model")
#Fit complex survival curves:
fit_interact<-survfit(Surv(time,delta)~factor(stage)+factor(diagyr_cate)+factor(age_cate),data=new_data)
ggsurv<-ggsurvplot(fit_interact,fun="event",conf.int=TRUE,ggtheme=theme_bw())
ggsurv$plot+theme(legend.position="right")+facet_grid(factor(age_cate)~factor(diagyr_cate))
#################
#Estimate the mean survival time:
EST_EXP_T<-time[1]*(KM[1]+1)/2
for(i in 1:(nn-1)){
	EST_EXP_T=EST_EXP_T+(time[i+1]-time[i])*(KM[i]+KM[i+1])/2
}
EST_EXP_T
integrate_interval<-function(i){
	INTE<-0
	for(j in i:(nn-1)){
		st<-(KM[j]+KM[j+1])/2
		INTE<-INTE+(time[j+1]-time[j])*st
	}
	INTE<-INTE*INTE*(n_event[i]/(n_risk[i]*(n_risk[i]-n_event[i])))
	return(INTE)
}
var_miu<-0
for(i in 1:(nn-1)){
	var_miu=var_miu+integrate_interval(i)
}
#The 95% confidence interval:
c(EST_EXP_T-qnorm(0.975)*sqrt(var_miu),EST_EXP_T+qnorm(0.975)*sqrt(var_miu))
EST_EXP_T
mean(larynx[,2])
#Estimate the median:
medi_point<-length(KM[KM>0.5]) #36
#So the median survival time is:
EST_MEDIAN<-(time[medi_point]+time[medi_point+1])/2 #Smooth
EST_MEDIAN
density<-(KM[medi_point]-KM[medi_point+1])/(time[medi_point+1]-time[medi_point])
VAR_MEDIAN<-0
for(i in 1:33){
	VAR_MEDIAN<-VAR_MEDIAN+n_event[i]/(n_risk[i]*(n_risk[i]-n_event[i]))
}
VAR_MEDIAN<-VAR_MEDIAN/(4*density^2)
VAR_MEDIAN
#95% interval of median:
c(EST_MEDIAN-qnorm(0.975)*sqrt(VAR_MEDIAN),EST_MEDIAN+qnorm(0.975)*sqrt(VAR_MEDIAN))
EST_MEDIAN
#plot the histogram to estimate the distribution:
m<-length(breaks)
KM_break<-c(1,KM[breaks])
KM_interval<-rep(0,m)
for(i in 1:m){
	KM_interval[i]<-KM_break[i]-KM_break[i+1]
}
plot(n,type="n",xlim=c(0,11),ylim=c(0,0.2),xlab="Time",ylab="hazard rate",main="Histogram for hazard rate")
segments(-100,0,5000,0)
segments(0,0,0,KM_interval[1])
segments(0,KM_interval[1],time[breaks][1],KM_interval[1])
segments(time[breaks][1],KM_interval[1],time[breaks][1],0)
for(i in 2:m){
	segments(time[breaks][i-1],0,time[breaks][i-1],KM_interval[i])
	segments(time[breaks][i-1],KM_interval[i],time[breaks][i],KM_interval[i])
	segments(time[breaks][i],KM_interval[i],time[breaks][i],0)
}
time_diff<-c(0,m)
time_break_diff<-c(0,time[breaks])
for(i in 1:m){
	time_diff[i]<-time_break_diff[i+1]-time_break_diff[i]
}
time_midd<-rep(0,m)
for(i in 1:m){
	time_midd[i]<-(time_break_diff[i]+time_break_diff[i+1])/2
}
lines(time_midd,KM_interval)
#
breaks<-seq(2,54,2)
m<-length(breaks)
KM_break<-c(1,KM[breaks])
KM_interval<-rep(0,m)
for(i in 1:m){
	KM_interval[i]<-KM_break[i]-KM_break[i+1]
}
plot(n,type="n",xlim=c(0,11),ylim=c(0,0.2),xlab="Time",ylab="hazard rate",main="Histogram for hazard rate")
segments(-100,0,5000,0)
segments(0,0,0,KM_interval[1])
segments(0,KM_interval[1],time[breaks][1],KM_interval[1])
segments(time[breaks][1],KM_interval[1],time[breaks][1],0)
for(i in 2:m){
	segments(time[breaks][i-1],0,time[breaks][i-1],KM_interval[i])
	segments(time[breaks][i-1],KM_interval[i],time[breaks][i],KM_interval[i])
	segments(time[breaks][i],KM_interval[i],time[breaks][i],0)
}
time_diff<-c(0,m)
time_break_diff<-c(0,time[breaks])
for(i in 1:m){
	time_diff[i]<-time_break_diff[i+1]-time_break_diff[i]
}
time_midd<-rep(0,m)
for(i in 1:m){
	time_midd[i]<-(time_break_diff[i]+time_break_diff[i+1])/2
}
lines(time_midd,KM_interval)
#Try Minimum chi-square technique for the survival data.
plot(time,type="l",main="Survival Time")
lines(sort(rexp(90,0.09)),lty=3) #Nearly like exp(0.09)?
CHI<-NULL
larynxx<-larynx[sort.list(larynx[,2]),]
statuss<-larynxx[,5]
#Randomly generate a break to construct a possibly gapped histogram
breaks<-break_specify
#breaks<-c(0,breaks)
for(i in 1:1000){
	alpha<-0.001*i
	chi_stat<-0
	m<-length(breaks)
	Z_1<-sqrt(n)*(-log(KM[breaks[2]]))
	DENO<-0
	for(j in 1:breaks[2]){
	DENO=DENO+n*status[j]/((n-j)*(n-j+1))
	}
	chi_stat=chi_stat+Z_1*Z_1/DENO
	for(i in 2:(m-2)){
		Z<-sqrt(n)*(-log(KM[breaks[i+1]])+log(KM[breaks[i]])-alpha*(time[breaks[i+1]]-time[breaks[i]]))
		DEOO<-0
	for(j in (breaks[i]):(breaks[i+1]-1)){
		DEOO=DEOO+n*status[j]/((n-j)*(n-j+1))
		}
	if(DEOO>0){
		chi_stat=chi_stat+Z*Z/DEOO
		}
}
CHI<-c(CHI,chi_stat)
}
sequence<-seq(0.001,1,0.001)
plot(sequence,CHI,type="l",xlab="alpha",ylab="Chi-Square Statistic",main="CHI-Square Statistic for the goodness of fit test
")
min(CHI) #133 Awful fit.


#plot the histogram of survival time:
hist(time,breaks=seq(0,11,1),freq=FALSE)
number<-c(1:11)
for(i in 1:11){
	number[i]<-length(which(time<=i && time>(i-1)))/length(time)
}
xdex<-seq(0.5,10.5,1)
lines(number,xdex)
#Estimating the survival function:
descdist(KM,discrete=FALSE)
#Non-parametric Tests:
MH_stage<-survdiff(formula=Surv(time,delta)~stage,data=larynx)
MH_age<-survdiff(formula=Surv(time,delta)~factor(age_cate),data=larynx)
MH_diagyr<-survdiff(formula=Surv(time,delta)~factor(diagyr_cate),data=larynx)
#Cox Model fitting:
cox.fit.full<-coxph(formula=Surv(time,delta)~.^2,data=larynx)
cox.fit.null<-coxph(formula=Surv(time,delta)~1,data=larynx)
cox.fit1<-stepAIC(cox.fit.full,direction="backward",k=2)
#cox.fit1(model1) denotes the model with stage,age,diagyr,age:diagyr
ggcoxzph(cox.zph(cox.fit1))
ggcoxdiagnostics(cox.fit.full,type="dfbeta",linear.predictions=FALSE,ggtheme=theme_bw())
#cox.fit2(model2) denotes the model with age,stage,diagyr
cox.fit2<-coxph(formula=Surv(time,delta)~age+factor(stage)+diagyr,data=larynx)
cox.zph(cox.fit2)
summary(cox.fit2)
#Model fitting for age,stage(model3)
cox.fit3<-coxph(formula=Surv(time,delta)~age+factor(stage),data=larynx)
ggcoxzph(cox.zph(cox.fit3))
#Compare model fit2 and fit3:
lrt_stat<-cox.fit2$loglik[2]-cox.fit3$loglik[2]
lrt_stat
lrt_stat>qchisq(0.95,3)
C<-rbind(c(0,0,0,0,1))
Ccov<-C%*%cox.fit2$var%*%t(C)
X.w=t(C%*%cox.fit2$coef)%*%ginv(Ccov)%*%(C%*%cox.fit2$coef)
X.w
X.w>qchisq(0.95,3)
#So we can ignore the diagyr term temporarily in this model fitting
#model fitting for stage(model4):
cox.fit4<-coxph(formula=Surv(time,delta)~factor(stage),data=larynx)
summary(cox.fit4)
#Residuals:
coxsnell_fit3<-larynx$time-resid(cox.fit3,type="martingale")
fit3_<-survfit(Surv(coxsnell_fit3,larynx$delta)~1)
htilde<-fit3_$n.event/fit3_$n.risk
plot(fit3_$time,htilde,type="s",col="blue")
abline(0,1,col="green",lty=2)
#Local Wald test whether the coefficients of stages different:
A<-rbind(c(0,1,-1,0),c(0,0,1,-1))
local_cov<-A%*%cox.fit3$var %*% t(A)
W_stat<-t(A%*%coef(cox.fit3))%*%solve(local_cov) %*%(A%*% coef(cox.fit3)) 
W_stat>qchisq(0.95,3) #True
W_stat_whole<-t(coef(cox.fit3))%*%solve(cox.fit3$var) %*%(coef(cox.fit3)) 
W_stat_whole>qchisq(0.95,3) #True
#Then we choose cox.fit3 and cox.fit4
#Cox-Snell Residuals(Indeed the alternative for Goodness of fit test)
par(mfrow=c(1,2))
coxsnellres3<-larynx$delta-resid(cox.fit3,type="martingale")
fitres3<-survfit(coxph(Surv(coxsnellres3,larynx$delta)~1,method="breslow"),type='aalen')
plot(fitres3$time,-log(fitres3$surv),type="s",xlab="Cox-Snell Residuals",main="Cox Model Goodness of Fit test for model3")
abline(0,1,col='red',lty=2)
coxsnellres4<-larynx$delta-resid(cox.fit4,type="martingale")
fitres4<-survfit(coxph(Surv(coxsnellres4,larynx$delta)~1,method="breslow"),type='aalen')
plot(fitres4$time,-log(fitres4$surv),type="s",xlab="Cox-Snell Residuals",main="Cox Model Goodness of Fit test for model4")
abline(0,1,col='red',lty=2)
#Do diagnostics(influential observations) via deviance residual:
ggcoxdiagnostics(cox.fit3,type="dfbeta")
ggcoxdiagnostics(cox.fit3)
#Prediction on age,stage:
blc<-data.frame(stage=1)
blc$stage<-factor(blc$stage)
newdata1<-data.frame(age=55,stage=1:4)
newdata2<-data.frame(age=65,stage=1:4)
newdata3<-data.frame(age=75,stage=1:4)
newdata4<-data.frame(age=85,stage=1:4)
model.new1<-survfit(cox.fit3,newdata=newdata1,type="kaplan-meier")
model.new2<-survfit(cox.fit3,newdata=newdata2,type="kaplan-meier")
model.new3<-survfit(cox.fit3,newdata=newdata3,type="kaplan-meier")
model.new4<-survfit(cox.fit3,newdata=newdata4,type="kaplan-meier")
par(mfrow=c(2,2))
plot(model.new1,col=rainbow(4),xlab="Age",ylab="Estimated Survival Function S(t)",main="Estimated Survival Time for Age 55")
legend("topright",c("Stage I","Stage II","Stage III","Stage IV"),col=rainbow(4),lty=1,cex=0.5)
plot(model.new2,col=rainbow(4),xlab="Age",ylab="Estimated Survival Function S(t)",main="Estimated Survival Time for Age 65")
legend("topright",c("Stage I","Stage II","Stage III","Stage IV"),col=rainbow(4),lty=1,cex=0.5)
plot(model.new3,col=rainbow(4),xlab="Age",ylab="Estimated Survival Function S(t)",main="Estimated Survival Time for Age 75")
legend("topright",c("Stage I","Stage II","Stage III","Stage IV"),col=rainbow(4),lty=1,cex=0.5)
plot(model.new4,col=rainbow(4),xlab="Age",ylab="Estimated Survival Function S(t)",main="Estimated Survival Time for Age 85")
legend("topright",c("Stage I","Stage II","Stage III","Stage IV"),col=rainbow(4),lty=1,cex=0.5)
#Cox model with factor(age) and factor(stage):
cox.fit5<-coxph(Surv(time,delta)~factor(age_cate)+factor(stage),data=new_data)
summary(cox.fit5)
ggcoxdiagnostics(cox.fit5)
cox.fit6<-coxph(Surv(time,delta)~factor(age_cate),data=new_data)
summary(cox.fit6)
ggcoxzph(cox.zph(cox.fit6))
cox.fit7<-coxph(Surv(time,delta)~factor(age_cate)+factor(stage)+factor(age_cate)*factor(stage),data=new_data)
summary(cox.fit7) #age group2:stage4 is significant
cox.fit8<-coxph(Surv(time,delta)~factor(age_cate)+factor(stage)+factor(age_cate==2)*factor(stage==4),data=new_data)


########################
#AFT model:
#Buckley James Estimator:
buckley_james<-bj(Surv(time,delta)~diagyr+stage,x=TRUE,y=TRUE,data=larynx) #not converge
f1<-bj(Surv(time,delta)~diagyr,data=larynx,x=TRUE,y=TRUE)
f2<-bj(Surv(time,delta)~age,data=larynx,x=TRUE,y=TRUE)
f3<-bj(Surv(time,delta)~age+stage,data=larynx,x=TRUE,y=TRUE)
f1
f2
f3
validate(f2,B=15)
validate(f3,B=15)
validate(f1,B=15)
#Semi-Parametric model by aftgee:
fit.ind.aft_1<-aftgee(formula=Surv(time,delta)~factor(stage),data=larynx)
fit.ind.aft_2<-aftsrr(formula=Surv(time,delta)~age+factor(stage),data=larynx)
#Fit with respect to interaction terms:
fit_aft_2<-aftsrr(formula=Surv(time,delta)~age+factor(stage)+age:factor(stage),data=larynx)
summary(fit.ind.aft) # |age| coefficients pretty small comparing with coefficients of stages.
summary(fit.ind.aft.diagyr)# diagyr coefficient pretty small compared with stage coefficients
#Using LSS function to fit semi-parametric model:
fit_aft_1<-lss(Surv(log(time),delta)~factor(stage),data=larynx,trace=TRUE,maxiter=20)
fit_aft_2<-lss(Surv(time,delta)~age+factor(stage),data=larynx,trace=TRUE,maxiter=20)
#Compare two methods:
summary(fit.ind.aft_1)
summary(fit.ind.aft_2)
fit_aft_1
fit_aft_2
#Estimating the acceleration factor
acce_factor<-exp(-coef(fit.ind.aft_1)[2:4])
acce_factor
#####
##Goodness of fit test on the residual from aftgee function in the two models:
res_1<-resid(fit.ind.aft_1)
U_1<-exp(sort(res_1))
res_2<-resid(fit.ind.aft_2)
U_2<-exp(sort(res_2))
breaks<-seq(10,90,10)
#For U_1:
n<-nrow(larynx)
k<-length(breaks)
Obj<-rep(10,k)
E<-c(1:k)
E[1]<-n*(pexp(U_1[breaks[1]],1))
for(i in 1:8){
	E[i+1]<-n*(pexp(U_1[breaks[i+1]],1)-pexp(U_1[breaks[i]],1))
}
Chi_Stat<-sum(((E-Obj)/E)^2) #1.431 very well.
Chi_Stat<qchisq(0.99,8) #True
#For U_2:
n<-nrow(larynx)
k<-length(breaks)
Obj<-rep(10,k)
E<-c(1:k)
E[1]<-n*(pexp(U_2[breaks[1]],1))
for(i in 1:7){
	E[i+1]<-n*(pexp(U_2[breaks[i+1]],1)-pexp(U_2[breaks[i]],1))
}
Chi_Stat<-sum(((E-Obj)/E)^2) #It dosen't work 
Chi_Stat<qchisq(0.99,8)  #False So we should apply other methods:
