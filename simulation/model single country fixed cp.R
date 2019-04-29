###################################################
##### JAGS (Just Another Gibbs Sampler) model #####
###################################################
model_string<-"

model{


for(j in 1:n.states){
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

cp1[i,j]<- 0
cp2[i,j]<- 24/max.time.points   #3xactly 24 m after cp1
 
##############################################################
    #Second Stage Statistical Model
##############################################################
beta[i,j, 1:2] ~ dmnorm(lambda[1:2], Sigma_inv[1:2, 1:2])
}


Sigma_inv[1:2, 1:2] ~ dwish(I_Sigma[1:2, 1:2], (2 + 1))
Sigma[1:2, 1:2] <- inverse(Sigma_inv[1:2, 1:2])

for(j in 1:2){
lambda[j] ~ dnorm(0, 1e-4) 
}
}
"