N.countries=length(countries)
log.rr.compile <-vector("list", length(countries)) 
log.rr.prec.mat.all<-vector("list", length(countries)) 
log.rr.prec.diag.mat.all<-vector("list", length(countries)) 
#N.states<-rep(NA, N.countries)
ts.length=tot_time
c=1
country<-countries[c]
  print(country)
  
  #####################################################################
    eval_period <- c(as.Date('2010-01-01'), as.Date('2012-12-01'))
    keep.grp='09'
 
  ####################################################################
  
  #Import the data for each country and extract the relevant sub groups and time periods
  input_directory <- paste0(getwd(),'/first stage estimates/',country, '/')
  sim1.ds<-readRDS(file=paste0('./simulation/', 'log_rr_sim cp.rds'))
  log_rr_q<-sim1.ds$log_rr.q
  #log_rr_prec<-sim1.ds$log_rr_full_t_samples.prec
  log_rr_sd<-sim1.ds$log_rr.sd
  time<-seq.Date(from=as.Date('2003-01-01'), by='month', length.out=120)
  
  index.post<-which(time>=(eval_period[1] %m-% months(pre.vax.time)) & time<=eval_period[2])  
  keep.index<- c((min(index.post)-12):(min(index.post)-1),index.post)
  if (length(index.post)>tot_time){ index.post<-index.post[1:tot_time] } 
  
  log_rr_q<-log_rr_q[,2,, drop=F]
  log_rr_q_all <- aperm(log_rr_q, c(3,1,2))
  log_rr_q_all<-log_rr_q_all[keep.index,,]
  
  log_rr_prec<-array(1/log_rr_sd^2, dim=c(nrow(log_rr_sd), ncol(log_rr_sd),1))
  log_rr_prec_diag_all <- aperm(log_rr_prec, c(2,1,3))
  log_rr_prec_diag_all<-log_rr_prec_diag_all[keep.index,,]



#PREPARE SPLINE BASIS using gam()
set.seed(15785)
N.knots <- 4  #N-1 internal knots
t.index<-1:max.time.points
out<-rnorm(n=max.time.points)
ds.sp<-cbind.data.frame(out,t.index)

spl.t.std<-bs(1:max.time.points, degree=3,df=N.knots, intercept=FALSE) 
zeros.pre<-matrix(0, nrow=pre.vax.time, ncol=ncol(spl.t.std))
spl.t.std<-rbind(zeros.pre,spl.t.std )
rowSums(spl.t.std)
matplot(spl.t.std, type='l')
spl.t.std<-cbind(rep(1, nrow(spl.t.std)), spl.t.std)

ts.length_mat<-matrix(NA, nrow=N.countries, ncol=max(N.states))
for(i in 1:N.countries){
  ts.length_mat[i,1:N.states[i]]<-ts.length[i]
}
identity<-diag(nrow(spl.t.std))
ts.length.vec<-as.vector(t(ts.length_mat))[!is.na(as.vector(t(ts.length_mat)))]

#Time index
time.index<-c(rep(0,times=pre.vax.time), seq.int(from=1, to=max.time.points, by=1))/max.time.points
#time.index<-c(rep(0,times=pre.vax.time), seq.int(from=1, to=max.time.points, by=1))

##Matrices to input to JAGS
#I_Sigma<-replicate( N.countries, diag(p) )
p=4 #intercept and slope for each time series
q=4  #2nd level predictors of intercept and slope q=1 for intercept only
I_Sigma<-diag(p)  #For p=N vcovariates of intercept and slope
I_Omega<-diag(q)  # q= number of predictors of the p slopes 
m<-1  #m=number of country-level covariates
w<-matrix(NA, nrow=N.countries, ncol=m) 
for(i in 1:N.countries){
  for(k in 1:m){
    w[i,k]<-1
  }
}
#Z needs to be an array i,j,q ; with assignment of covariate based on hdi.index df
# hdi<-t(matrix(word(state.labels,3),nrow=nrow(state.labels), ncol=ncol(state.labels)))
# hdi[hdi=='A']<-"Low"  #if only category available is "A", assume it is low (ie nc)
# hdi.low<-matrix(0, nrow=nrow(hdi),ncol=ncol(hdi))
# hdi.low[hdi=="Low"]<-1
# hdi.med<-matrix(0, nrow=nrow(hdi),ncol=ncol(hdi))
# hdi.med[hdi=="Med"]<-1
hdi1<-read.csv('./covariate data/hdi.csv')
hdi1<-hdi1[,c('brazil','mexico', 'ecuador','chile')]
hdi1.spl<-split(t(hdi1) ,1:ncol(hdi1))
library(dummies)
hdi1<-sapply(hdi1.spl, function(x) dummies::dummy(x))

##THIS NEEDS TO BE FIXED!!
#z<-array(0, dim=c(dim(log_rr_q_all)[3], dim(log_rr_q_all)[2],q))
# z[,,1]<-1 #intercept for z 
 
# z[,,2]<-hdi.low #intercept for z 
# z[,,3]<-hdi.med #intercept for z 
