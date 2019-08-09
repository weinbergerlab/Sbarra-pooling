library(lubridate)
library(reshape2)
library(splines)
pre.vax.time=12
max.time.points.post=48
#rmarkdown::render('./first stage estimates/Chile_state/chile_state.Rmd')
#rmarkdown::render('./first stage estimates/Ecuador_state/ecuador_state.Rmd')
#rmarkdown::render('./first stage estimates/Mexico_state/mexico_state.Rmd')
#rmarkdown::render('./first stage estimates/Brazil_state/brazil_state.Rmd')

 chi1<-readRDS('./first stage estimates/Chile_state/Results_Chile_state.Rds')
 ec1<-readRDS('./first stage estimates/Ecuador_state/Results_Ecuador_state.Rds')
 mx1<-readRDS('./first stage estimates/Mexico_state/Results_Mexico_state.Rds')
 br1<-readRDS('./first stage estimates/Brazil_state/Results_Brazil_state.Rds')


#filter out national and regional results for Brazil
br.names<-dimnames(br1$results$impact$best$log_rr_quantiles)[[3]]
aggregate.regions<-grep('A',br.names)
if(length(aggregate.regions)>0){
br1$results$impact$best$log_rr_quantiles<-br1$results$impact$best$log_rr_quantiles[,,-aggregate.regions]
br1$results$impact$best$log_rr_sd<-br1$results$impact$best$log_rr_sd[,,-aggregate.regions, drop=F]
}

all.countries<-list('Chile_state'=chi1, 'Ecuador_state'=ec1, 'Mexico_state'=mx1, 'Brazil_state'=br1)
subnational=c(1,1,1,1) #Vector indicating whether dataset contains subnational estmatesmax.time.points=48
countries<-names(all.countries)
all.dates<-lapply(all.countries, function(x) unique(x$input$data$date))
all.post.start<-lapply(all.countries, function(x) unique(x$input$post_period[1]))
all.post.end<-lapply(all.countries, function(x) unique(x$input$post_period[2]))

all.log.rr<-lapply(all.countries, function(x) x$results$impact$best$log_rr_quantiles[,'50%',])
all.log.rr.prec<-lapply(all.countries, function(x) (1/x$results$impact$best$log_rr_sd[1,,]^2))

#Trim the time series to pre.vax.time to max.time.points
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}
post.index.func<-function(ds.date, ds.start){
  index.post<-elapsed_months(ds.date,ds.start)  
  #index.post[index.post<0]<-0
  return(index.post)
}
all.post.index<-mapply(post.index.func,all.dates,all.post.start)
trim.func<-function(ds.time,ds1){
  ds2<-ds1[ds.time>=(-pre.vax.time) & ds.time < max.time.points.post,]
  return(ds2)
}
all.log.rr<-mapply( trim.func, all.post.index,all.log.rr)
all.log.rr.prec<-mapply( trim.func, all.post.index,all.log.rr.prec)

N.states<-sapply(all.log.rr, function(x)  dim(x)[2])
N.time.series <- sum(N.states)
ts.length<-sapply(all.log.rr, function(x)  dim(x)[1])
N.countries<-length(all.log.rr)
state.labels.list<-sapply(all.log.rr, function(x)  dimnames(x)[[2]])
max.time.points<- max.time.points.post+pre.vax.time

##Combine together estimates from each country into a single array
tot_time<-pre.vax.time+ max.time.points
log_rr_q_all<-array(NA, dim=c(tot_time,max(N.states),N.countries))
log_rr_prec_diag_all<-array(NA, dim=c(tot_time,max(N.states),N.countries))
state.labels<-array(NA, dim=dim(log_rr_q_all)[c(2:3)]   )

for(i in 1:N.countries){
  for(j in 1:N.states[i]){
    log_rr_q_all[1:ts.length[i],j,i]<-all.log.rr[[i]][,j] #Extract the median
    log_rr_prec_diag_all[1:ts.length[i],j,i]<-all.log.rr.prec[[i]][,j]
  }
  if(N.states[i]>1){
    state.labels[1:N.states[i] ,i]<-dimnames(all.log.rr[[i]])[[2]]
  }else{
    state.labels[1,i]<-'national'
  }
}
N.time.series <- sum(N.states)

##Spline basis

#PREPARE SPLINE BASIS using gam()
set.seed(15785)
N.knots <- 4  #N-1 internal knots
log.rr.q.extr<-log_rr_q_all[,,1]
sp.df<-as.data.frame(log.rr.q.extr)
names(sp.df)<-paste0('state', 1:ncol(sp.df))
sp.df$time=1:nrow(sp.df)
sp.m<-melt(sp.df ,id='time'  )
sp.c<- acast(sp.m, variable~time)
t.index<-1:max.time.points
out<-rnorm(n=max.time.points)
ds.sp<-cbind.data.frame(out,t.index)

spl.t.std<-bs(1:max.time.points.post, degree=3,df=N.knots, intercept=FALSE) 
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
# hdi1<-read.csv('./covariate data/hdi.csv')
# hdi1<-hdi1[,c('brazil','mexico', 'ecuador','chile')]
# hdi1.spl<-split(t(hdi1) ,1:ncol(hdi1))
# library(dummies)
# hdi1<-sapply(hdi1.spl, function(x) dummies::dummy(x))

##THIS NEEDS TO BE FIXED!!
#z<-array(0, dim=c(dim(log_rr_q_all)[3], dim(log_rr_q_all)[2],q))
#z[,,1]<-1 #intercept for z 

# z[,,2]<-hdi.low #intercept for z 
# z[,,3]<-hdi.med #intercept for z 



#Need to remove aggregate state and national for brazil (A and AA)


