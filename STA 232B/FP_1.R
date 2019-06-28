######Model Selection for mixed model:



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
  #Expand all the candidate variables as binary encoding(involved or not involved)
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

    #Note that for the mixed model selection, we should use the ML instead of REML.
    current_model<-lmer(formula=current_formula,data=data_matrix,REML=FALSE)#,REML=FALSE)

    #Calculate the AIC and BIC, here we apply the alternative BIC criterion
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
