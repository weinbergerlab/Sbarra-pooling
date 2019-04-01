###################################################
##### JAGS (Just Another Gibbs Sampler) model #####
###################################################
model_string<-"

model{

for(i in 1:n.countries){
for(j in 1:n.states[i]){
for(v in 1:ts.length[i,j]){       

##################
#Observation model
##################
w_hat[v,j,i] ~ dnorm(w_true[i,j, v], w_hat_prec[ v, j,i])

#################################
#Model of 'true' time series data
#################################
w_true[i,j, v] ~ dnorm(reg_mean[i,j,v], w_true_prec[i,j]) 

##################### 
#CHANGE POINT MODEL #
#####################
reg_mean[i,j,v]<-(beta[i,j,1] 
+ step(time.index[v] - cp1[i,j])*(1 - step(time.index[v] - cp2[i,j]))*beta[i,j,2]*(time.index[v] - cp1[i,j]) 
+ step(time.index[v] - cp2[i,j])*beta[i,j,2]*(cp2[i,j] - cp1[i,j])
)
reg_unbias[i,j,v]<-(
step(time.index[v] - cp1[i,j])*(1 - step(time.index[v] - cp2[i,j]))*beta[i,j,2]*(time.index[v] - cp1[i,j]) 
+ step(time.index[v] - cp2[i,j])*beta[i,j,2]*(cp2[i,j] - cp1[i,j])
)
}

w_true_prec[i,j]<-1/(w_true_sd[i,j]*w_true_sd[i,j])
w_true_sd[i,j] ~ dunif(0, 100)

cp1[i,j]<-exp(beta[i,j,3])
cp2.add[i,j]<-exp(beta[i,j,4])
cp2[i,j]<-cp1[i,j] +cp2.add[i,j]  + 1/max.time.points   #ensure Cp2 is at least 1 unit after CP1
 
##############################################################
    #Second Stage Statistical Model
##############################################################
beta[i,j, 1:4] ~ dmnorm(mu1[i,j, 1:4], Omega_inv[i, 1:4, 1:4])
for(k in 1:4){
mu1[i,j,k] <- z[i,j, 1:q]%*%gamma[i,k, 1:q]
}
########################################################
#Third Stage Statistical Model
########################################################
gamma[i, 1:4] ~ dmnorm(lambda[1:4], Sigma_inv[1:4, 1:4])
}

#######################################################
#Remaining Prior Distributions
#######################################################
Omega_inv[1:4, 1:4] ~ dwish(I_Omega[1:4, 1:4], (4 + 1))
Omega[1:4, 1:4]<-inverse(Omega_inv[1:4, 1:4])
Sigma_inv[1:4, 1:4] ~ dwish(I_Sigma[1:4, 1:4], (4 + 1))
Sigma[1:4, 1:4]<-inverse(Sigma_inv[1:4, 1:4])
for(j in 1:4){
lambda[j] ~ dnorm(0, 1e-4) 
}

}
"