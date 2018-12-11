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
            + step(t-cp[1,i,j])*(1-step(t-cp[1,i,j]))*beta[i,j, 2]*(t-cp[1,i,j])
            + step(t-cp[2,i,j])*beta[i,j, 2]*(cp[2,i,j])
            )
       }
      for(k1 in 1:ts.length[i,j]){
        for(k2 in 1:ts.length[i,j]){
        w_true_cov_inv[i,j,k1,k2]<-ifelse(k1==k2, w_true_var_inv[i,j], 0)
        }
      }
    w_true_var_inv[i,j]<-1/(w_true_sd[i,j]*w_true_sd[i,j])
    w_true_sd[i,j] ~ dunif(0, 1000)

    #This ensures CP[2] is after CP[11..see JAGS manual
      for(v in c(1:2)){
       cp0[v,i,j]~dgamma(0.001,0.001)
      }
      cp[1:2,i,j]<-sort(cp0[,i,j])

    ##############################################################
    #Second Stage Statistical Model
    ##############################################################
    beta[i,j, 1] ~ dnorm(0, 1e-4) #Prior on intercept
    beta[i,j, 2] ~ dnorm(0, 1e-4)T(,0) #Truncate slope to be 0 or negative
    }
  }
}
"
