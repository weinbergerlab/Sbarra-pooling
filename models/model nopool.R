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
        reg_mean[i,j, 1:ts.length[i,j]]<-spl.t.std[1:ts.length[i,j],]%*%beta[i,j, 1:p]
      for(k1 in 1:ts.length[i,j]){
        for(k2 in 1:ts.length[i,j]){
        w_true_cov_inv[i,j,k1,k2]<-ifelse(k1==k2, w_true_var_inv[i,j], 0)
        }
      }
    w_true_var_inv[i,j]<-1/(w_true_sd[i,j]*w_true_sd[i,j])
    w_true_sd[i,j] ~ dunif(0, 1000)
    ##############################################################
    #Second Stage Statistical Model
    ##############################################################
    beta[i,j, 1:p] ~ dmnorm(mu1[1:p], Omega_inv[1:p, 1:p])
    
    }
}
#############################################################
#Remaining Prior Distributions
#############################################################
##IF q=1 (intercept only for predictors of slopes and intercepts)
for(k1 in  1:p){   
  mu1[k1]<-0
  for(k2 in 1:p){
  Omega_inv[k1,k2]<-ifelse(k1==k2, (1/1000), 0)
  }
}

}
"