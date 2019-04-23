
library (RCurl)
library(boot)
library(reshape2)
library(RColorBrewer)
library(matlib)
library(parallel)
library(rstanarm)
library(loo)  #iMPORTANT: Need to install development version of Loo (SEE ABOVE)
library(splines, quietly = TRUE)


#############################################
data.start<-1 #Determines which year to include; set to 1 to include full TS, to 61 to just 2008+
#############################################

#Import the data and select the right age group, get rid of columns with NAs


n_cores <- detectCores() - 1

t=1:120
intro.date=85
n=length(t)

xx1<- t/100
xx2 <- t/100 + t/100*t/100
xx3= -0.02 *t
  
ncovars=5
n.states=27
results<- vector("list", n.states) #combine models into a list
cp.true<-rep(NA, n.states)
log.rr.true<-rep(NA, n.states)
log_rr.q <- array( dim=c(n.states,3,2000) ) 
log_rr.sd <- array( dim=c(n.states,2000) ) 
log_rr_full_t_samples.prec <-  array(dim=c(n.states,n,n) ) 

set.seed(123)
for (k in 1:n.states){
  print(k)
b1<-rnorm(n=ncovars  )
b2<-rnorm(n=ncovars)
b3<-rnorm(n=ncovars)
int<- rnorm(n=1, mean=3, sd=0.75)
#hist(int<- rnorm(n=100, mean=3, sd=0.75))
covars<-as.data.frame(matrix(nrow=length(t), ncol=ncovars))
  
  for(i in 1:ncovars){
covars[,i]<- rpois(n=length(t), exp(int+ b1[i]*xx1 + b2[i]*xx2  + b3[i]*xx3) )
names(covars)[i]<-paste0('covar',i)
  }
  matplot(log(covars+0.5))

  spl<-1:120
  cp=sample(85:108, 1)
  spl<-spl-cp
  spl[spl<0]=0
  log.rr.true[k]<- rnorm(n=1, mean=log(0.8), sd=0.05)   #allow for some true variation in vaccine effect between locations
  # hist(rnorm(n=1000, mean=log(0.8), sd=0.051))
  vax.effect=  log.rr.true[k]/max(spl) #20% decline by end of time series
  J12_18= rpois( n=length(t),exp( int+ b1[i]*xx1 + b2[i]*xx2 + spl*vax.effect))

y.fit=J12_18
y.fit[t>=85]<-NA
mod1<-glm(y.fit~ ., data=log(covars+0.5), family="poisson")
pred1=predict(mod1, newdata=log(covars+0.5), type='response')
rr<-J12_18/pred1
plot(rr)
sim.data<-cbind.data.frame(k, J12_18,covars)
state=k

####################
#STAN MODEL
CORES=n_cores
CHAINS=2
SEED=123
data.fit<- sim.data
outcome.pre<-data.fit$J12_18
outcome.pre[t>intro.date]=NA
data.fit[,3:7]<-log(data.fit[,3:7]+0.5)
stan_glm.mod <- stan_glm(J12_18 ~ covar1 +covar2+covar3+covar4 +covar5,
                         data = data.fit, family = poisson()   , QR=FALSE, 
                         prior = normal(0,5), prior_intercept = normal(0,1000), prior_aux=NULL , 
                         chains = CHAINS, cores = CORES, seed = SEED ,iter=2000)
preds   <-  posterior_predict(stan_glm.mod, newdata=data.fit)
log_rr<-log(t(data.fit$J12_18/preds))
log_rr.q[k,,]<-apply(log_rr,2 ,quantile, probs=c(0.025,0.5,0.975), na.rm=TRUE)
log_rr.sd[k,]<-apply(log_rr,2 ,sd, na.rm=TRUE)
log_rr_full_t_samples.covar<-cov(t(log_rr))
log_rr_full_t_samples.prec[k,,]<-solve(log_rr_full_t_samples.covar)

cp.true[k]<cp


}
results<- list( cp.true, log.rr.true,  log_rr.q,log_rr.sd, log_rr_full_t_samples.prec )
names(results[[k]])<-c('cp.true','log_rr.true', 'log_rr.q', 'log_rr.sd', 'log_rr_full_t_samples.prec')


saveRDS(results,'C:\\Users\\dmw63\\Desktop\\My documents h\\GATES\\CausalImpact code\\Meta analysis\\SIM\\log_rr_sim.rds' )

