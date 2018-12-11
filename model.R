###################################################
##### JAGS (Just Another Gibbs Sampler) model #####
###################################################
model_string<-"
model{
  for(i in 1:n.countries){
    for(j in 1:N.states[i]){
        #################
        #Observation model
        #################
        w_hat[1:ts.length[i,j], j,i] ~ dmnorm(w_true[i,j, 1:ts.length[i,j]], log_rr_prec_all[1:ts.length[i,j], 1:ts.length[i,j], j,i])
        #################
        #Model of 'true' time series data
        #################
        w_true[i,j, 1:ts.length[i,j]] ~ dmnorm(reg_mean[i,j, 1:ts.length[i,j]], w_true_cov_inv[i,j, 1:ts.length[i,j], 1:ts.length[i,j]])
      for(t in 1:ts.length[i,j] ){    
      ##################### 
      #CHANGE POINT MODEL #
      #####################
      reg_mean[i,j, t]<-(beta[i,j, 1] 
            + step(t-cp1[i,j])*(1-step(t-cp2[i,j]))*beta[i,j, 2]*(t-cp1[i,j]) #Slope for period between CP1 and CP2
            + step(t-cp2[i,j])*beta[i,j, 2]*(cp2[i,j]-cp1[i,j])
            )
       }
      for(k1 in 1:ts.length[i,j]){
        for(k2 in 1:ts.length[i,j]){
        w_true_cov_inv[i,j,k1,k2]<-ifelse(k1==k2, w_true_var_inv[i,j], 0)
        }
      }
    w_true_var_inv[i,j]<-1/(w_true_sd[i,j]*w_true_sd[i,j])
    w_true_sd[i,j] ~ dunif(0, 1000)

    #In beta matrix, beta1=intercept, beta2=slope, beta3=changepoint, beta4=
    #This ensures CP[2] is after CP1
      cp1[i,j]<-exp(beta[i,j, 3])
      cp2.add[i,j]<-exp(beta[i,j, 4])
      cp2[i,j]<-cp1[i,j] +cp2.add[i,j]
    
    ##############################################################
    #Second Stage Statistical Model
    ##############################################################
    beta[i,j, 1:p] ~ dmnorm(mu1[i,j, 1:p], Sigma_inv[i, 1:p, 1:p])
for(k in 1:p){
    mu1[i,j,k] <- z[i,j, 1:q]%*%gamma[i,k, 1:q]
}
}
for(k in 1:p){
###############################################################
#Third Stage Statistical Model
###############################################################
gamma[i,k, 1:q] ~ dmnorm(mu2[i,k, 1:q], Omega_inv[k, 1:q, 1:q])
for(l in 1:q){
mu2[i,k,l] <- w[i, 1:m]%*%theta[k,l, 1:m]
}
} 
Sigma_inv[i, 1:p, 1:p] ~ dwish(I_Sigma[1:p, 1:p], (p + 1))
Sigma[i, 1:p, 1:p] <- inverse(Sigma_inv[i, 1:p, 1:p])
}
#############################################################
#Remaining Prior Distributions
#############################################################
##IF q=1 (intercept only for predictors of slopes and intercepts)
for(k in 1:p){
Omega_inv[k, 1:q, 1:q] <- 1/(Omega[k, 1:q, 1:q]*Omega[k, 1:q, 1:q])
Omega[k, 1:q, 1:q] ~ dunif(0, 1000)
for(l in 1:q){
for(r in 1:m){
theta[k,l,r] ~ dnorm(0, 0.0001)
  }
 }
 }
}
"
