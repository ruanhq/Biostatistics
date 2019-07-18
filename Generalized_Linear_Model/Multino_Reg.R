######################################################
#Project for STA 223:Generalized Linear Model
#Here we use the dataset
library(MASS)
library(mgcv)
library(caret)
library(lawstat)
stu<-read.table("student-mat.csv",sep=";",header=TRUE)
head(stu)
stude<-stu[,-c(31,32)]
head(stude)
n<-nrow(stu)
#G3 is the final grade. Here we aim to explore the functional relationship between the final grade of students and the predictors. There are 30 predictors which is a mixture of categorical variables and continuous variables. Some variables are binary while some others are nominal(multinomial) and others are numeric.
n<-nrow(stu)
#histogram of the students' grade:
hist(stude$G3,main="Histogram of G3")
table(stude$G3)
A_B<-16
B_C<-14
C_D<-12
D_F<-10
#Here we may choose various parametrization of the grades:
#(1): We choose the A,B,C,D,F parametrization: 16-20 means(A), 14-15 means good(B), 12-13 means (C),10-11 means(D),below 10 means (F)
#(2):We split the grade into pass(over 10) and fail(below 10)
passing<-c(1:n)
for(i in 1:n){
	if(stude[i,31]<D_F){
		passing[i]<-0
	}
	else{
		passing[i]<-1
	}
}
level<-c(1:n)
for(i in 1:n){
	if(stude[i,31]<10){
		level[i]<-"F"
	}
	else if((9<stude[i,31]) && (stude[i,31]<C_D)){
		level[i]<-"D"
	}
	else if((11<stude[i,31])&&(stude[i,31]<B_C)){
		level[i]<-"C"
	}
	else if((13<stude[i,31]) &&(stude[i,31]<A_B)){
		level[i]<-"B"
	}
	else{
		level[i]<-"A"
		}
}
stupas<-data.frame(cbind(stude[,1:30],passing))
#Passing 265 cases, passing 130 cases.
stulevel<-data.frame(cbind(stude[,1:30],level))
#  A    B   C    D    F 
# 40  60  62 103 130 
 #Perform logistic regression model:
yes<-rep(0,n)
no<-rep(0,n)
for(i in 1:n){
	if(stude[i,31]<D_F){
		no[i]<-1
	}
	else{
		yes[i]<-1
	}
}
fit_logi1<-glm(cbind(yes,no)~.+failures*schoolsup,data=stude[,1:30],family=binomial())
fitted_logi1<-fitted(fit_logi1)
fit<-rep(0,n)
for(i in 1:n){
	if(fitted_logi1[i]<0.5){
		fit[i]<-0
	}
	else{
		fit[i]<-1
	}
}
tab<-table(fit,passing)
1-sum(diag(tab))/sum(tab)# 24.0% classification error
tab
coeftable<-coef(summary(fit_logi1))
which(coeftable[,4]<0.05) #age,failure,schoolsup_yes,famsupyes,goout
coeftable[coeftable[,4]<0.05,]
#Checking Interaction effect:
fit_logi_inter<-glm(cbind(yes,no)~(age+failures+schoolsup+famsup+goout)*(age+failures+schoolsup+famsup+goout),data=stude[,1:30],family=binomial())
#failures:schoolsup, final logistic model:
fit_logi_int<-glm(cbind(yes,no)~age+failures+schoolsup+famsup+goout+failures*schoolsup,data=stude[,1:30],family=binomial())
summary(fit_logi_int)
fitttt<-stepAIC(fit_logi1,trace=FALSE)
fitted_logi1<-fitted(fitttt)
fit_logistic<-rep(0,n)
for(i in 1:n){
  if(fitted_logi1[i]<0.5){
    fit_logistic[i]<-0
  }
  else{
    fit_logistic[i]<-1
  }
}
tab<-table(fit_logistic,passing)
1-sum(diag(tab))/sum(tab)
#10 fold Cross Validation to assess the classifier:
errr<-rep(0,10)
sample<-createFolds(1:n,10,list=TRUE,returnTrain=FALSE)
for(i in 1:10){
	fit_logistic_i<-glm(cbind(yes[-as.numeric(unlist(sample[i]))],no[-as.numeric(unlist(sample[i]))])~.,data=stude[-as.numeric(unlist(sample[i])),1:30],family=binomial)
	fitte<-fitted(fit_logistic_i)
	fit_l<-predict(fit_logistic_i,stude[as.numeric(unlist(sample[i])),1:30])
	nn<-length(fit_l)
	for(j in 1:nn){
	  if(fit_l[j]<0.5){
	    fit_l[j]<-0
	  }
	  else{
	    fit_l[j]<-1
	  }
	}
	tab<-table(fit_l,passing[as.numeric(unlist(sample[i]))])
	errr[i]<-1-tab[2]/sum(tab)
}
mean(errr) 
plot(x=c(1:10),y=sort(errr),type="l",main='')
#Try to add some potential interactions identified from previous running:
fit_logi3<-glm(cbind(yes,no)~age+failures+schoolsup+famsup+goout,data=stude[,1:30],family=binomial())#Model1
par(mfrow=c(2,2))
res.D<-resid(fit_logi1,type="deviance")
res.P<-resid(fit_logi1,type="pearson")
boxplot(cbind(res.D,res.P),label=c("Deviance Residual","Pearson Residual"))
plot(fitted(fit_logi1),res.D,main="Deviance Residual vs Fitted Value")
lines(smooth.spline(fitted(fit_logi1),res.D,spar=0.97))
plot(fitted(fit_logi1),res.P,main="Pearson Residual vs Fitted Value")
lines(smooth.spline(fitted(fit_logi1),res.P,spar=0.9))
plot(fitted(fit_logi3),res.P,main="Pearson Residual vs Fitted Value ")
fitted_logi1<-fitted(fit_logi_final)
fit_logistic<-rep(0,n)
for(i in 1:n){
  if(fitted_logi1[i]<0.5){
    fit_logistic[i]<-0
  }
  else{
    fit_logistic[i]<-1
  }
}
tab<-table(fit_logistic,passing)
1-sum(diag(tab))/sum(tab)
#Diagnostics:Leverages and Cook's Distances:
par(mfrow=c(1,2))
leverage<-hatvalues(fit_logi1)
plot(names(leverage),leverage,xlab="Index",type="h",main="Leverage Plot")
points(names(leverage),leverage,cex=0.4,type="p",pch=3,col="red")
abline(h=2*mean(leverage),col=2,lwd=2,lty=2)
cookk<-cooks.distance(fit_logi1)
plot(cookk,ylab="Cook's Distance",pch=4,cex=0.2,main="COOK DISTANCE PLOT")
index1<-which(leverage>2*mean(leverage))
points(index,cookk[index],col="red",type="h",cex=0.6)
abline(h=4/(n-length(coef(fit_logi1))))
index2<-which(cookk>4/(n-length(coef(fit_logi1))))
#######################################################

#Multinomial models:
#Proportional Odds Model:
fitplr_1<-polr(factor(level)~.,data=stulevel)
pred_labl_prob<-predict(fitplr_1,stulevel)
table1<-table(pred_labl_prob,level)
1-sum(diag(table1))/sum(table1) #misclassification error 56.2% Awful!
coeftable_plr<-coef(summary(fitplr_1))
which(coeftable_plr[,3]>qt(0.95,n-length(coef(fitplr_1))))

#age,Mjobteacher,failures,schoolsupyes    famsupyes,goout,health
#There maybe some interaction terms in this model. Here we may explore whether school/family support is helpful after previous failure.
fitplr_2<-polr(factor(level)~age+Mjob+failures+goout+health+schoolsup+famsup+failures+failures*schoolsup,data=stulevel)#Model 2
pred_labl_prob_plr<-predict(fitplr_2,stulevel)
table_plr<-table(pred_labl_prob_plr,level)
#Misclassification error 59.5% What a awful fit!
##############
#Baseline Odds Model:
library(nnet)
level<-relevel(factor(level),ref="F")
fitbo<-nnet::multinom(level~.,data=stulevel)
prd_prob_bo<-predict(fitbo,stulevel)
table2<-table(prd_prob_bo,level)
1-sum(diag(table2))/sum(table2) #49.9%
#The predicting error is pretty high.
std_err_bo<-summary(fitbo)$standard.errors
coef_bo<-coef(summary(fitbo))
coef_bo/std_err_bo>qt(0.95,n-length(coef(fitbo)))
exp(coef_bo)
#Baseline Odds model may not be appropriate here.

#We have prior information that G1 and G2 are highly informative to predict G3.
#So we try to perform PCA on the G1,G2,G3 to extract as much information as possible from the data matrix.
grade_matrix<-as.matrix(stu[,31:33])
pc_1<-prcomp(grade_matrix,scale=FALSE)
summary(pc_1) #The first princomp component explain 90% of total variation.
combined_grade<-grade_matrix%*%pc_1$rotation[,1]
studen<-cbind(stu[,1:30],combined_grade)
#Then we transform the response to the whole number axis by first transforming them into percentage of highest possible score
highest_possible_score<-20*sum(pc_1$rotation[,1])
scaled_score<-logit(combined_grade/highest_possible_score)
hist(scaled_score)
#The scaled score is generally following normal distribution. So we try general linear regression.
fit_lm<-lm(scaled_score~.,data=stu[,1:30])
#studytime,sexM are positively related to score
#failures,schoolsup,famsup,goout,romanticyes(explainable) are negatively related to the score. 
fit_transformed_model<-lm(scaled_score~sex+studytime+failures+schoolsup+famsup+romantic+goout,data=stu[,1:30])
summary(fit_transformed_model)
fit_transformed_final<-stepAIC(fit_transformed_model,trace=FALSE)
png('fitted_PC.png',width=800,height=800)
par(mfrow=c(2,2))
plot(fit_transformed_final,1)
plot(fit_transformed_final,2)
plot(fit_transformed_final,3)
plot(fit_transformed_final,4)
dev.off()
################################################
#We merge the 5 categories(ABCDF) to three categories(High, Medium, Low) 
#to alleviate the highly imbalance between categories
level2<-rep(0,n)
for(i in 1:n){
	if(stude[i,31]<D_F){
		level2[i]<-1
	}
	else if((stude[i,31]>D_F-1 &&(stude[i,31]<A_B))){
		level2[i]<-2
	}
	else{
		level2[i]<-3
	}
}
stude_three<-cbind(stude[,1:30],level2)
fitplr_3_cate<-polr(factor(level2)~.,data=stude_three)
pred<-predict(fitplr_3_cate,stude_three[,1:30])
table<-table(pred,level2)
1-sum(diag(table))/sum(table)#36.2%
table
#####Performance Evaluation:
#Calculate Multi-label F score:
pre1<-table[1,1]/sum(table[1,])
pre2<-table[2,2]/sum(table[2,])
pre3<-table[3,3]/sum(table[3,])
pre<-mean(pre1,pre2,pre3)
rec<-mean(table[1,1]/sum(diag(table)-table[1,1]),table[2,2]/sum(diag(table)-table[2,2]),table[3,3]/sum(diag(table)-table[3,3]))
accuracy<-sum(diag(table))/sum(table)
F_score<-2*pre*rec/(pre+rec)
F_score_1.2<-2.44*pre*rec/(1.44*pre+rec)
Fscore<-rep(0,100)
for(i in 1:100){
  Fscore[i]<-(1+(1+0.01*i)^2)*pre*rec/((1+0.01*i)^2*pre+rec)
}
plot(c(1:100),Fscore,type="l",cex=0.7,col="blue")


coeftable<-coef(summary(fitplr_3_cate))
tval<-2*(1-pt(abs(coeftable[,3]),n-length(coef(fitplr_3_cate))))
coeftable<-cbind(coeftable,tval)
abs(coeftable[,3])>qt(0.975,n-length(coef(fitplr_3_cate)))
#age,Fjobteacher,failures,schoolsup,goout are significant.
fitted_plr<-fitted(fitplr_3_cate)
fi<-polr(factor(level2)~famsup,data=stude_three)
ffi<-polr(factor(level2)~schoolsup,data=stude_three)
#####
#plot(x=fitted_plr[,3],rdi_high,xlab="Probability of Medium",ylab="Resid",main="Deviance Residual_High")
#lines(smooth.spline(x=fitted_plr[,3],y=rdi_medium,spar=1.25),col="blue")
boxplot(cbind(rpi_medium,rdi_medium),label=c("Pearson","Deviance"),main="Medium")
boxplot(cbind(rpi_high,rdi_high),label=c("Pearson","Deviance"),main="High")
runs.test(rdi_high)
runs.test(rpi_high)
runs.test(rdi_medium)
runs.test(rpi_medium)

################
#Evaluate the prediction performance here after merging categories:
fitbo<-multinom(factor(level2)~.+failures*schoolsup,data=stude_three)
pred<-predict(fitbo,stude_three[,1:30])
table<-table(pred,level2)
1-sum(diag(table))/sum(table)#31.9%
std_err_bo<-summary(fitbo)$standard.errors
coef_bo<-coef(summary(fitbo))
coef_medium<-coef_bo[1,]
coef_high<-coef_bo[2,]
pva_medium<-2*(1-pt(abs(coef_medium),n-length(fitbo)))
pva_high<-2*(1-pt(abs(coef_high),n-length(fitbo)))
taa<-cbind(coef_medium,pva_medium,coef_high,pva_high)
taa

#Check the fit of the basel
fitted_bo<-fitted(fitbo)
#For class medium:
rpi_medium<-rep(0,n)
rdi_medium<-rep(0,n)
for(i in 1:n){
  rpi_medium[i]<-(as.numeric(level2[i]==2)-fitted_bo[i,2])/sqrt(fitted_bo[i,2]*(1-fitted_bo[i,2]))
}
for(i in 1:n){
  di<-(as.numeric(level2[i]==2))*log(fitted_bo[i,2]/(1-fitted_bo[i,2]))+log(1-fitted_bo[i,2])
  rdi_medium[i]<-sign(di)*sqrt(2*abs(di))
}
par(mfrow=c(2,2))
plot(x=fitted_plr[,2],rpi_medium,xlab="Probability of Medium",ylab="Resid",main="Pearson Residual_Medium",cex=0.7)
lines(smooth.spline(x=fitted_bo[,2],y=rpi_medium,spar=0.9),col="blue")
#plot(x=fitted_plr[,2],rdi_medium,xlab="Probability of Medium",ylab="Resid",main="Deviance Residual_Medium")
#lines(smooth.spline(x=fitted_bo[,2],y=rdi_medium,spar=1.2),col="blue")
boxplot(cbind(rpi_medium,rdi_medium),label=c("Pearson","Deviance"))
#For class high:
rpi_high<-rep(0,n)
rdi_high<-rep(0,n)
for(i in 1:n){
  rpi_high[i]<-(as.numeric(level2[i]==3)-fitted_bo[i,3])/sqrt(fitted_bo[i,3]*(1-fitted_bo[i,3]))
}
for(i in 1:n){
  di<-(as.numeric(level2[i]==3))*log(fitted_bo[i,3]/(1-fitted_bo[i,3]))+log(1-fitted_bo[i,3])
  rdi_high[i]<-sign(di)*sqrt(2*abs(di))
}
plot(x=fitted_plr[,3],rpi_high,xlab="Probability of Medium",ylab="Resid",main="Pearson Residual_High",cex=0.7)
lines(smooth.spline(x=fitted_bo[,3],y=rpi_medium,spar=1.6),col="blue")
boxplot(cbind(rpi_high,rdi_high),label=c("Pearson","Deviance"))
runs.test(rpi_high)
runs.test(rdi_high)
runs.test(rpi_medium)
runs.test(rpi_medium)