##
library(sae)
data("cornsoybean")
data("cornsoybeanmeans")
# delete 33rd sample point
cornsoybean <- cornsoybean[-33,]
m <- 12
PX1 <- cornsoybeanmeans$MeanCornPixPerSeg
PX2 <- cornsoybeanmeans$MeanSoyBeansPixPerSeg
PX1sq <- PX1^2
PX2sq <- PX2^2
PX1X2 <- PX1*PX2
PX = cbind(PX1,PX2,PX1sq,PX2sq,PX1X2)
ni <- c(1,1,1,2,3,3,3,3,4,5,5,5)
X1 <- cornsoybean$CornPix
X2 <- cornsoybean$SoyBeansPix
X1sq <- X1^2
X2sq <- X2^2
X1X2 <- X1*X2
Y <- cornsoybean$CornHec
county = as.factor(rep(1:12,ni,each=T))
fit_F = lmer(Y~X1+X2+X1sq+X2sq+X1X2+(1|county))
#####

### model_selection

nested_model_selection<-function(response,
                                 data_matrix,
                                 candidate_variable=
                                   c("CornPix","SoyBeansPix","I(CornPix^2)","I(SoyBeansPix^2)","CornPix:SoyBeansPix"),
                                 AIC= TRUE){
  #candidate_variable<-c("CornPix","SoyBeansPix","I(CornPix^2)","I(SoyBeansPix^2)","CornPix:SoyBeansPix")
  n<-length(response)
  data_matrix<-cbind(response,data_matrix)
  colnames(data_matrix)[1]<-'response'
  County<-data_matrix$County
  n_county<-max(County)
  n_var<-length(candidate_variable)
  n_possible<-32
  candidate_encoding<-expand.grid(replicate(5,0:1,simplify=FALSE))[2:n_possible,]
  list_criterion<-rep(0,nrow(candidate_encoding))
  #Traverse all the possible conditions:
  for(i in 1:nrow(candidate_encoding)){
    involved_terms<-which(candidate_encoding[i,]!=0)
    #current_data<-data_matrix[involved_variables,]
    formula_1<-""
    for(k in 1:length(involved_terms)){
      formula_1<-paste(formula_1,"+",candidate_variable[involved_terms[k]],sep="")
    }
    current_formula<-paste("response~ (1|County) ",formula_1,sep="")
    current_model<-lmer(formula=current_formula,data=data_matrix,REML=FALSE)#,REML=FALSE)
    k1<-ncol(model.matrix(current_model))+2
    if(AIC==TRUE){
      list_criterion[i]<-(-2*logLik(current_model))+2*k1
    }
    else{
      #For BIC we use the alternative BIC criterion:
      list_criterion[i]<-(-2*logLik(current_model)+log(n)*(k1-2)+log(n_county)*2)
    }
  }
  #####
  #select the indice of what we select:
  variables_indice<-which(candidate_encoding[which.min(list_criterion),]!=0)
  formula_1<-""
  for(k in 1:length(variables_indice)){
    formula_1<-paste(formula_1,"+",candidate_variable[variables_indice[k]],sep="")
  }
  selected_formula<-paste("response~ (1|County) ",formula_1,sep="")
  selected_model<-lmer(formula=selected_formula,data=data_matrix,REML=FALSE)
  variable_selected<-candidate_variable[variables_indice]
  return(list(selected_model=selected_model,variables_indice=variables_indice,variable_selected=variable_selected))
}





#fit_F2<-lmer()
fit_F2<-lmer(CornHec~CornPix+SoyBeansPix+(1|County),data=cornsoybean)
##### theta_hat;; here out is the out of nested_model_selection
theta= function(out,fit_F2)
{
  candidate = out$variables_indice[out$variables_indice<3]
  #Check it whether there are the different ones:
  fit = out$selected_model
  beta = fixef(fit)
  betanew = beta[(1:(length(candidate)+1))]
  ref = unlist(ranef(fit))
  if(length(betanew==1)){
    result<-betanew+ref
  }
  else{
    result = betanew[1]+as.matrix(PX[,candidate])%*%beta[2:length(betanew)]+ref
  }
  result
}

###### variance fit is the full model
variance = function(fit)
{
  sigma2_v = unlist(VarCorr(fit))
  sigma2_e = sigma(fit)^2
  sigma2_v*sigma2_e/(ni*sigma2_v+sigma2_e)
}


##### 
cond_exp <- function(Y, lmerF){
  f_beta <- fixef(lmerF)
  f_sig_v <- sqrt(unlist(VarCorr(lmerF)))
  f_sig_e <- sigma(lmerF)
  PX_new = cbind(1,PX)
  
  diag <- list()
  for (ns in ni) {
    diag <- c(diag, list(rep(1/ns,ns)))
  }
  Z <- bdiag(diag)
  ybar <- t(Z)%*%Y
  f_xbar <- t(Z)%*%cbind(X1,X2,X1sq,X2sq,X1X2)
  
  exp_result <- PX_new %*% f_beta + ni*f_sig_v^2/(ni*f_sig_v^2+f_sig_e^2)*
    (ybar-f_xbar%*%f_beta[2:6]-f_beta[1])
  
  return(exp_result)
}



#compute a 
#cond_exp: out 1, variance: out 2, theta: out 3.
a<-function(out1,out2,out3){
  a<-(out3-out1)^2+out2
  return(a)
}

### 
d_k = function(fit_F)
{
  Y_new  = simulate(fit_F)
  out_new = nested_model_selection(response = Y_new,cornsoybean,AIC)
  names(Y_new) = "Y_new"
  fit_F_new = lmer(Y_new~X1+X2+X1sq+X2sq+X1X2+(1|county),data = cbind(Y_new,X1,X2,X1sq,X2sq,X1X2,county))
  out1 = cond_exp(unlist(Y_new),fit_F)
  out2 = variance(fit_F)
  out3 = theta(out_new)
  a_former = a(out1,out2,out3)
  out12 = cond_exp(unlist(Y_new),fit_F_new)
  out22 = variance(fit_F_new)
  out32 = theta(out_new)
  a_latter = a(out12,out22,out32)
  a_former-a_latter
}



#####
fit_selected = nested_model_selection(Y,cornsoybean)
out1 = cond_exp(Y,fit_F)
out2 = variance(fit_F)
out3 = theta(fit_selected)



K = 12
d = rep(0,12)
for(i in 1:K)
{
  d = d+d_k(fit_F)
}


d/12+a(out1,out2,out3)


