####
#MCEM(Monte-Carlo EM algorithm to find the MLE of the linear mixed effects model):
####

#####
#First create the data inputs:
#// Reference: Lab7
nlist <- list()
nlist[[1]] <- c(15,1)
nlist[[2]] <- c(6,1,2)
nlist[[3]] <- c(6,6)
nlist[[4]] <- c(7,2,3,NA,2)
nlist[[5]] <- c(16,9,3,3,1)
nlist[[6]] <- c(57,38,17,2,2)
nlist[[7]] <- c(119,81,45,6,1,NA,NA,1)
nlist[[8]] <- c(173,118,57,16,3,NA,NA,NA,1)
nlist[[9]] <- c(136,103,50,13,6,1,1)
nlist[[10]] <- c(54,51,32,5,1,NA,NA,NA,NA,1)
nlist[[11]] <- c(13,15,12,3,1)
nlist[[12]] <- c(NA,4,3,1)
nlist[[13]] <- c(NA,NA,1,NA,NA,NA,NA,1)
myobs <- vector()
count <- 0
for (k in 1:13) {
  for (l in 1:length(nlist[[k]])) { 
    if(!is.na(nlist[[k]][l])){
      for(m in 1:nlist[[k]][l]){
        count <- count + 1
        myobs <- rbind(myobs, cbind( c(rep(0,k-l+1),rep(1,l-1) ) , count))
      }}
  }}
myobs <- data.frame(y = myobs[,1], litter = factor(myobs[,2]))


########
library(lme4)
library(MASS)
library(nlme)
library(HLMdiag)


#Get the initial fitting by glmer:
n0 <- glmer(y ~ (1|litter), data = myobs, family = binomial, nAGQ = 20)
mu_n0 <- fixef(n0)
sigma_n0 <- sqrt(as.numeric(VarCorr(n0)))

#Get the initial fitting for glmmPQL:
pql <- glmmPQL(y ~ 1, random = ~1 | litter, family = binomial, data = myobs)
mu_pql <- fixef(pql)
sigma_pql <- sqrt(getVarCov(pql)[1])

mcem_main_func <- function(
  mcmc_burnin = 1200, mcmc_iter = 2000,
  iter = 1000, epsi = 0.00001, mu_starter = mu_pql, sigma_starter = sigma_pql
){
  #Get the summary statistic:
  y_i <- split(myobs$y, myobs$litter)
  sum_y_i <- sapply(y_i, sum)
  n_i <- sapply(y_i, length)
  m = length(unique(myobs$litter))
  
  #Initialization:
  mu_hist = rep(0, iter)
  mu_hist[1] = mu_starter
  sigma_hist = rep(0, iter)
  sigma_hist[1] = sigma_starter
  Alppha = matrix(0, mcmc_iter, m)
  
  #Main loop:
  for(j in 2:iter){
    
    #MCMC to sample from Q|y,\theta(E-step)
    for(d in 2:mcmc_iter){
      current_draw = rnorm(m, 0, sigma_hist[j-1])
      Az = pmin(1, (((1 + exp(mu_hist[j-1] + Alppha[d-1, ]))/(1 + exp(mu_hist[j-1] + current_draw)))^n_i) * exp((current_draw - Alppha[d-1, ]) * sum_y_i))
      accepted = runif(m, 0, 1) < Az
      Alppha[d, ] = ifelse(accepted, current_draw, Alppha[d-1, ])
    }
    
    #Apply the burnin for MCMC
    alpha_current = Alppha[(mcmc_burnin+1): mcmc_iter,]
    
    #Calculate the objective function to optimize:
    Qvalue = function(Mu){
      expit = exp(Mu + alpha_current) / (1 + exp(Mu + alpha_current))
      q = (sum_y_i %*% colSums(log(expit)) + (n_i - sum_y_i) %*% colSums(log(1 - expit))) / nrow(alpha_current)
      return(q)
    }
    
    #M-step:
    Alppha[1,] = Alppha[mcmc_iter, ]
    mu_hist[j] = optimize(Qvalue, interval = c(-5, 0), maximum = TRUE)$maximum
    sigma_hist[j] = sqrt(mean(alpha_current ^ 2))
    
    #Print the current result for the estimation:
    print(paste('Iter', j,': Current_Mu_Estimation: ', mu_hist[j], ', Current_Sigma_Estimation: ', sigma_hist[j], sep = ''))
    
    #Checking convergence:
    dif = max(abs(mu_hist[j] - mu_hist[j-1]), abs(sigma_hist[j] - sigma_hist[j-1]))
    
    #If the difference smaller than the tolerance, stop fitting.
    if(dif < epsi){
      mu_final = mu_hist[j]
      sigma_final = sigma_hist[j]
      break
    }
  }
  
  #Not escape in the middle:
  if(j == iter){
    mu_final = mu_hist[iter]
    sigma_final = sigma_hist[iter]
  }
  
  png(paste('MCMC', mcmc_iter, 'Iter', iter, 'Init_mu', mu_starter,'Mu_Trajectory.png', sep = ''), width = 500, height = 500)
  plot(mu_hist[1:j], type = 'l', lwd = 1.5, xlab = 'Iter', ylab = 'Mu', main = 'Mu_Trajectory')
  dev.off()
  png(paste('MCMC', mcmc_iter, 'Iter', iter, 'Init_sigma', sigma_starter, 'Sigma_Trajectory.png', sep = ''), width = 400, height = 400)
  plot(sigma_hist[1:j], type = 'l', lwd = 1.5, xlab = 'Iter', ylab = 'Sigma', main = 'Sigma_Trajectory')
  dev.off()
  return(list(sigma_final, mu_final))
}

result_pql <- mcem_main_func()
result_glmer <- mcem_main_func(mu_starter = mu_n0, sigma_starter = sigma_n0)
