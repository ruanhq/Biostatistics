##############

#Logistic_Regression Group-Lasso for student dataset:
library(grplasso)
library(bootstrap)
library(fastDummies)
library(ROCR)
library(PRROC)

stu<-read.table("student-mat.csv",sep=";",header=TRUE)
#Define 12(60%) as the passing grade for the final:
passing <- as.numeric(stu$G3 > 12)
stu$pass <- passing
stu_final <- cbind(stu[,1:30], pass = stu$pass)
n <- nrow(stu_final)

#transform the features into the factors:
for(i in 1:30){
  if(class(stu_final[,i]) != 'factor'){
    stu_final[,i] <- as.factor(stu_final[,1])
  }
}
X_stu <- stu_final[,1:30]

#Note that 'age' and 'absences' are numerical variables and others are variables.
Age <- X_stu$age
Absence <- X_stu$absences
health <- X_stu$health
famrel <- X_stu$famrel
freetime <- X_stu$freetime
goout <- X_stu$goout
Medu <- X_stu$Medu
Fedu <- X_stu$Fedu

X_stu_categorical <- subset(X_stu, select = -c(age, absences, health, famrel, freetime, goout, Medu, Fedu))
for(i in 1:ncol(X_stu_categorical)){
  X_stu_categorical[,i] <- as.factor(X_stu_categorical[,i])
}
n_cate <- ncol(X_stu_categorical)

#Then create the dummy representation for the categorical parts:
X_dummy_cate <- dummy_cols(X_stu_categorical, remove_first_dummy = TRUE)
n_tot <- ncol(X_dummy_cate)
X_dummy <- X_dummy_cate[,(n_cate + 1):n_tot]

#Then construct the continuous parts:
X_continuous <- as.data.frame(cbind(rep(1, nrow(stu)),Age = Age, Absence = Absence, health = health,
                                    famrel = famrel, freetime = freetime, goout = goout,
                                    Medu = Medu, Fedu = Fedu))

#Then construct the group information in preparation for the group lasso:
X_train <- as.data.frame(cbind(X_continuous, X_dummy))

#Then construct the group index:
groups_index <- c(NA, 1, 2, 3, 4, 5, 6, 7, 8)
for(i in 1:n_cate){
  n_group <- len(unique(X_stu_categorical[,i])) - 1
  groups_index <- c(groups_index, rep(i + 8, n_group))
}

##########
#Apply the 10-fold Cross-Validation to choose the optimal model for downstream interpretation:
Pass <- stu$pass

theta_predict <- function(model, X){
  output <- predict(model, X, type = 'response')
}
####
#Estimate the CV-error:
CV_error <- function(train_df, train_label, lambda, k = 10){
  index <- groups_index
  list_val <- crossval(train_df, train_label, grplasso, theta_predict,
                       index = index, model = LogReg(), lambda = lambda, center = TRUE, standardize = TRUE, ngroup = k)
  val_error <- mean(as.numeric(list_val$cv.fit>0.5)!= train_label )
  return(val_error)
}
#Perform the leave_one_out cross-validation:
lambda_max <- 50
m_lambda <- 50
lambda_sequence <- rep(0, 45)
for(i in 1:m_lambda){
  lambda_sequence[i] <- lambda_max * (0.9)^(i-1)
}

#Fit the regularization path selected by the LooCV:
oob_pred <- rep(0, n)
oob_score <- rep(0, n)
mcc_seq <- rep(0, n)
aupr_score <- rep(0, m_lambda)
accuracy <- rep(0, m_lambda)
for(k in 1:m_lambda){
  #Perform the LooCV and select optimal lambda by MCC and AUPR
  cur_lambda <- lambda_sequence[k]
  for(i in 1:n){
    X_tra <- as.matrix(X_train[-i,])
    Y_tra <- Pass[-i]
    model0 <- grplasso(X_tra, Y_tra, index = groups_index, lambda = cur_lambda)
    oob_score[i] <- predict(model0, as.matrix(X_train[i,]), type = 'response')
    oob_pred[i] <- as.numeric(cur_pred > 0.5)
  }
  aupr_score[k] <- pr.curve(oob_score, Pass)$auc.integral
  accuracy[k] <- mean(Pass == oob_pred)
}
##########

##########
#Visualize the important predictors:
#optimal_lambda <- lambda_sequence[which.max(aupr_score)]
final_model <- grplasso(x = as.matrix(X_train), y = Pass, index = groups_index, lambda = 7)
coeff <- final_model$norm.pen
names(coeff) <- c('Intercept',
                  'Age', 'Absence', 'health', 'famrel', 'freetime', 'goout', 'Medu', 'Fedu', '')
#Return the variables by its importance:
names(sort(coeff[coeff!=0], decreasing = TRUE))







