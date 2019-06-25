
####################
#Code for Project: 
#   Investigation of crime in US
#    Heqiao Ruan
# Department of Statistics, 
# Genome Center
# University of California, Davis
################


library(MASS)
library(glmnet)
library(pls)
library(DAAG)
library(e1071)
crime<-read.csv("crime.csv", header=TRUE)
n<-nrow(crime)

#summary statistics:
apply(crime, 2, summary)

#Preprocess:
fit0<-lm(Crime~., data=crime)
png('prefit.png', width = 1000, height = 1000)
par(mfrow = c(2,2))
plot(fit0, which=1)
plot(fit0, which=2)
plot(fit0, which=3)
plot(fit0, which=4)
dev.off()

#Pairwise scatterplot:
plot(crime)

#Histogram and transformation:
par(mfrow = c(3,4))
crime<-data.frame(crime)
sapply(colnames(crime[,2:4]),function(x){
  hist(crime[[x]], main=sprintf("Hist of %s",x), xlab=x)
  hist(log(crime[[x]]), main=sprintf("Hist of log(%s)",x), xlab=sprintf("%s",x))
  hist(crime[[x]]^2, main=sprintf("Hist of square(%s)",x), xlab=sprintf("%s",x))
  hist(sqrt(crime[[x]]), main=sprintf("Hist of squareroot(%s)",x), xlab=sprintf("%s",x))
})
par(mfrow=c(3,4))
sapply(colnames(crime[,5:7]),function(x){
  hist(crime[[x]], main=sprintf("Hist of %s",x), xlab=x)
  hist(log(crime[[x]]), main=sprintf("Hist of log(%s)",x), xlab=sprintf("%s",x))
  hist(crime[[x]]^2, main=sprintf("Hist of square(%s)",x), xlab=sprintf("%s",x))
  hist(sqrt(crime[[x]]), main=sprintf("Hist of squareroot(%s)",x), xlab=sprintf("%s",x))
})
par(mfrow=c(3,4))
sapply(colnames(crime[,8:10]), function(x){
  hist(crime[[x]], main=sprintf("Hist of %s",x), xlab=x)
  hist(log(crime[[x]]), main=sprintf("Hist of log(%s)",x), xlab=sprintf("%s",x))
  hist(crime[[x]]^2, main=sprintf("Hist of square(%s)",x), xlab=sprintf("%s",x))
  hist(sqrt(crime[[x]]), main=sprintf("Hist of squareroot(%s)",x), xlab=sprintf("%s",x))
})
par(mfrow = c(3,4))
sapply(colnames(crime[,11:13]),function(x){
  hist(crime[[x]], main=sprintf("Hist of %s",x), xlab=x)
  hist(log(crime[[x]]), main=sprintf("Hist of log(%s)",x), xlab=sprintf("%s",x))
  hist(crime[[x]]^2, main=sprintf("Hist of square(%s)",x), xlab=sprintf("%s",x))
  hist(sqrt(crime[[x]]), main=sprintf("Hist of squareroot(%s)",x), xlab=sprintf("%s",x))
})
par(mfrow = c(1,2))
hist(crime$Crime)
hist(log(crime$Crime))

#Then we do transformations: log(UE2),log(PE.1),squareroot(Pop),log(NW),square(UE1)
crime$UE2 <- log(crime$UE2)
crime$PE.1 <- log(crime$PE.1)
crime$Pop <- sqrt(crime$Pop)
crime$NW <- log(crime$NW)
crime$UE1 <- (crime$UE1)^2
crime$Crime <- log(crime$Crime)
crime_std <- as.data.frame(scale(crime))
#pairwise scatterplot:
pairs(crime_std)

#Draw the correlation plot after transformation:
library(corrplot)
par(mfrow=c(1,1), mar=c(1,1,1,1))
corrplot(cor(crime_std), method="circle", pch=1, tl.cex=0.75)
corrplot(cor(crime_std), method="number", pch=1, tl.cex=0.75)


#Calculate the Variance Inflation factor:
fit0 <- lm(Crime~., data=crime_std)
ols_vif_tol(fit0)
crime_std <- as.data.frame(as.matrix(crime_std[,-4]))
fit0 <- lm(Crime~., data=crime_std)
ols_vif_tol(fit0)

#Implement the Stepwise regression and perform diagnostics:
fit1<-lm(Crime~.,data=crime_std)
fit1_AIC<-stepAIC(fit1, trace=FALSE)
par(mfrow=c(2,2))
plot(fit1_AIC, which=1)
plot(fit1_AIC, which=2)
plot(fit1_AIC, which=6)
plot(fit1_AIC, which=4)
ols_vif_tol(fit1_AIC)# Evaluate the VIF after model selection.


#Ridge regression:
X <- cbind(rep(1,n), as.matrix(crime_std[,-1]))
ridge1 <- lm.ridge(Crime~., data=crime_std, lambda=seq(0,2,0.0001))
plot(seq(0,2,0.0001), ridge1$GCV, main="Plot of GCV", xlab="k", ylab="GCV", cex=0.2)
k_select <- 1.2548
fit_ridge <- lm.ridge(Crime~., data=crime_std, lambda=k_select)
beta <- matrix(coef(fit_ridge))
fv <- X%*%beta
res <- crime_std$Crime-fv


#Model diagnostics:
par(mfrow=c(2,2))
#Fitted value ~ observation plot
plot(crime_std$Crime, fv, main="Observed vs Fitted", xlab="Fitted", ylab="Observed Values")
#Residual plot
plot(fv, res, main="residuals vs fitted values", xlab="Fitted values", ylab="residuals")
lines(smooth.spline(fv,res,spar=1.06))
#Histgram of residual
hist(res, main="Histogram of residuals", xlab="residuals")
#Q-Q plot
qqnorm(res)
qqline(res)


#Coefficients estimation for ridge regression:
Y <- as.matrix(crime_std$Crime)
D <- t(X)%*%X
H <- X%*%solve(D+k_select*diag(ncol(X)))%*%t(X)
beta_e <- solve(D+k_select*diag(ncol(X)))%*%t(X)%*%Y
residual <- Y-X%*%beta_e
dematrix <- diag(ncol(H))-H
sigma2 <- sum(residual^2)/sum(diag(dematrix%*%dematrix))
cov_beta <- sigma2*solve(D+k_select*diag(ncol(X)))%*%D%*%solve(D+k_select*diag(ncol(X)))
std_error <- sqrt(diag(cov_beta))
estimation <- beta_e
data.frame(estimation,std_error)

#VIF values for ridge regression:
vif_cov <- solve(D+k_select*diag(ncol(X)))%*%D%*%solve(D+k_select*diag(ncol(X)))
vif_ridge <- diag(D)*diag(vif_cov)
as.matrix(vif_ridge)
coef_ridge <- as.matrix(coef(fit_ridge))
sort(coef_ridge, decreasing=TRUE)

#Lasso:
fit_lasso <- cv.glmnet(X, Y, intercept=FALSE)
plot(fit_lasso)
fit_lasso$lambda.min
coef(fit_lasso)

#lasso model diagnostic:
Ylasso <- predict(fit_lasso, newx=X)
par(mfrow=c(2,2))
plot(Ylasso, Y, xlab="Fitted Value", ylab="Observed value", main="Obs vs Fitted")
plot(Ylasso, Y-Ylasso, xlab="Fitted Value", ylab="Residuals", main="Residual vs Fitted")
lines(smooth.spline(Ylasso,Y-Ylasso,spar=1.3))
hist(Y-Ylasso, main="Histogram of Residuals", xlab="Residuals")
qqnorm(Y-Ylasso)
qqline(Y-Ylasso)



#Partial Least Square:
set.seeds(1999)
fit_pls <- plsr(Crime~.,10,data=crime_std,validation="CV")
summary(fit_pls)
sse <- rep(0,10)
for(i in 1:10){
  sse[i] <- sum((residuals(fit_pls)[(47*i-46):(47*i)])^2)
}
ssto <- sum((crime_std$Crime)^2)
R_square <- rep(0,10)
for(i in 1:10){
  R_square[i] <- 1-sse[i]/ssto
}
#Select the order:
k <- c(1:10)
data.frame(k,sse,R_square)
F_k <- rep(0,10)
F_k[1] <- (nrow(crime_std)-1-1)*R_square[1]/(1-R_square[1])
for(i in 2:10){
  F_k[i] <- (nrow(crime_std)-i-1)*(R_square[i]-R_square[i-1])/(1-R_square[i])
}
qF <- rep(0,10)
for(i in 1:10){
  qF[i] <- qf(0.95,1,nrow(crime_std)-i-1)
}
data.frame(k,sse, R_square, F_k, qF)
plot(fit_pls, plottype="scores",comp=1:7)
loadings(fit_pls)[,1:7]
fit_pls_final <- plsr(Crime~., 7, data=crime_std, validation="CV")


#################################
#Model selection by sum of 10-fold cross validation error:
#################################
cv <- function(traindata,testdata,method){
  if(method=="stepregr"){
    fit <- lm(Crime~.,data=traindata)
    stepfit <- stepAIC(fit,trace=FALSE)
    res <- testdata$Crime-predict(stepfit,testdata)
    return(sum(res^2)/nrow(testdata))
  }
  else if(method=="ridge"){
    rid <- lm.ridge(Crime~.,data=traindata,lambda=seq(0,30,0.001))
    k_select <- 0.001*which.min(rid$GCV)
    fit <- lm.ridge(Crime~.,data=traindata,lambda=k_select)
    beta <- matrix(coef(fit))
    X <- cbind(rep(1,nrow(testdata)),testdata[,-1])
    fv <- X%*%beta
    res <- testdata$Crime-fv
    return(sum(res^2)/nrow(testdata))
  }
  else if(method=="lasso"){
    X <- cbind(rep(1,nrow(traindata),traindata[,-1]))
    Y <- as.matrix(traindata[,1])
    fit <- cv.glmnet(X,Y,intercept=FALSE)
    X_pred <- cbind(rep(1,nrow(testdata)),testdata[,-1])
    Y_pred <- as.matrix(testdata[,1])
    yfit <- predict(fit,X_pred)
    res <- Y_pred-yfit
    return(sum(res^2)/nrow(testdata))
  }
  else if(method=="pls"){
    fit_pls = plsr(Crime~.,10,data=traindata, validation="CV")
    k_select = which.min(fit_plr$validation$PRESS)
    fit = plsr(Crime~.,k_select,data=traindata, validation="CV")
    n <- nrow(testdata)
    fv <- predict(fit,newdata=testdata)[((k_select-1)*n+1):(k_select*n)]
    res <- testdata$Crime-fv
    return(sum(res^2)/nrow(testdata))
  }
}
getCV<-function(method,datanew,nfolder){
  n <- nrow(crime)
  cv_value <- rep(0,nfolder)
  n = floor(nrow(datanew)/nfolder)
  for(i in 1:nfolder){
    testdata = datanew[((i-1)*n+1):(i*n),]
    traindata = datanew[((i-1)*n+1):(i*n),]
    cv_value[i] <- CV(traindata,testdata,method)
  }
  print(method)
  return(sum(cv_value)/nfolder)
}
sampl <- sample(1:n,n,relplace=FALSE)
datanew <- crime_std[sampl,]
method = c("stepreg","ridge","lasso","pls")
datanew <- sample

#Calculate the 10-fold CV error for each of the method:
CV_result <- sapply(method, function(x) getCV(x,datanew,10))
CV_result
