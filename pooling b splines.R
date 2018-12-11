
rm(list=ls(all=TRUE))
require(reshape)
require(reshape2)
require(abind)
library(mgcv)   
library(R2jags) 
library(splines)
library(lubridate)

output_directory<- paste0(dirname(getwd()), "/Results/")
ifelse(!dir.exists(output_directory), dir.create(output_directory), FALSE)

####SET INPUTS PARAMETERS#############################################################################################################
#setwd('C:/Users/dmw63/Dropbox (Personal)/meta_analysis_results/stack all states/')
countries<-c('Brazil_state', 'Chile_state')
subnational=c(1,1) #Vector indicating whether dataset contains subnational estmates
max.time.points=48
pre.vax.time<-12 #how many to use when anchoring at t=0
tot_time<-max.time.points+pre.vax.time
#####################################################################################################################################
#format data
source('compile stage1 results.R') #Read in data
source('model.R')  #Read in model

#Run Model
model_jags<-jags.model(textConnection(model_string),
                      data=list('n.countries' = N.countries, 
                                'N.states' = c(5,5), 
                                'w_hat' = log_rr_q_all[,1:5,, drop=FALSE],
                                'log_rr_prec_all' = log_rr_prec_all[,,1:5,, drop=FALSE],  
                                'ts.length' = ts.length_mat,
                                'p'=4,
                                'q'=q,
                                'm'=m,
                                'w'=w,
                                'z'=z,
                                'I_Omega'= I_Omega,
                                'I_Sigma'=I_Sigma), n.chains=2, n.adapt=1000) 

#Posterior Sampling
update(model_jags, 
       n.iter=5000)  #Burnin for 10,000 Samples
# dic.mod1a <- dic.samples(model_jags, 5000,"pD")  #DIC=2161 with Brazil state, Mexico state, US
# dic.mod1a

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("reg_mean", 'cp',"beta", "sigma_regression", "w_true", "alpha.C", "beta_k_q",'beta_prec_ts'),
                                thin=10,
                                n.iter=10000)
#plot(posterior_samples, ask=TRUE)


##########################################################################################
##### test for Convergence - trace plots (specify which to test in "variable.names") #####
##########################################################################################
#par(ask=TRUE)
#dev.off()
#pdf(file = "Trace Plots.pdf")
#par(mfrow = c(3,2))
#plot(posterior_samples)
#dev.off()

##################################################################
##### model fit (DIC) - must change n.chains=2 in model_jags #####
##################################################################
# dic.samples(model_jags, 500) #with pooling josh's fix (2 runs w 500 each): deviance: -717, -715; penalty=269,266, DIC=-447, -449
 #dic.samples(model_jags, 500)

# par(mfrow=c(1,1))
# caterplot(posterior_samples)

#############################################################################
##### create posterior estimates into data frames (w.true and reg.mean) #####
#############################################################################
betas.m1<-posterior_samples[[1]][,grep("^beta.*,1]",dimnames(posterior_samples[[1]])[[2]])]
betas.m2<-posterior_samples[[1]][,grep("^beta.*,2]",dimnames(posterior_samples[[1]])[[2]])]
betas.m3<-posterior_samples[[1]][,grep("^beta.*,3]",dimnames(posterior_samples[[1]])[[2]])]
betas.m4<-posterior_samples[[1]][,grep("^beta.*,4]",dimnames(posterior_samples[[1]])[[2]])]
betas.m5<-posterior_samples[[1]][,grep("^beta.*,5]",dimnames(posterior_samples[[1]])[[2]])]
betas.m6<-posterior_samples[[1]][,grep("^beta.*,6]",dimnames(posterior_samples[[1]])[[2]])]
beta.labs<-colnames(betas.m1) # labels from matrix
beta.labs.extract<- sub(".*\\[(.*)\\].*", "\\1", beta.labs, perl=TRUE)    #extract index numbers from the labels
beta.labs.extract2<-matrix(as.numeric(unlist(strsplit(beta.labs.extract, ","))),ncol=3, byrow=TRUE) #format into a matrix


preds<-array(NA,dim=c(nrow(betas.m1), nrow(spl.t.std), ncol(betas.m1)))
dimnames(preds)[[3]]<-beta.labs.extract

for(i in 1:nrow(betas.m1)){ #by iteration
  #ds.select 
  for(j in 1:ncol(betas.m1[i,, drop=FALSE] )){
   # print(i)
    #print(j)
    spl.sub<-as.matrix(spl.t.std[1:ts.length.vec[j], ,drop=FALSE])
    preds[i,1:ts.length.vec[j],j]<-(spl.sub[,1, drop=FALSE] %*% betas.m1[i,j, drop=FALSE] 
                                    +spl.sub[,2, drop=FALSE] %*% betas.m2[i,j, drop=FALSE] 
                                    +spl.sub[,3, drop=FALSE] %*% betas.m3[i,j, drop=FALSE] 
                                    +spl.sub[,4, drop=FALSE] %*% betas.m4[i,j, drop=FALSE] 
                                    +spl.sub[,5, drop=FALSE] %*% betas.m5[i,j, drop=FALSE] )
      }
}
preds.q<-as.data.frame(t(apply(preds,c(2,3),quantile, probs=c(0.5),na.rm=TRUE)))
preds.ucl<-as.data.frame(apply(preds,c(2,3),quantile, probs=c(0.975),na.rm=TRUE))
preds.lcl<-as.data.frame(apply(preds,c(2,3),quantile, probs=c(0.025),na.rm=TRUE))
saveRDS(preds, file=paste0(output_directory,"reg_mean_with_pooling bsplines.rds"))
saveRDS(state.labels, file=paste0(output_directory,"state labels.rds"))

###Remove bias term (intercept)
preds.nobias<-array(NA,dim=c(nrow(betas.m1), nrow(spl.t.std), ncol(betas.m1)))
dimnames(preds.nobias)[[3]]<-beta.labs.extract
for(i in 1:nrow(betas.m1)){ #by iteration
  #ds.select 
  for(j in 1:ncol(betas.m1[i,, drop=FALSE] )){
    # print(i)
    #print(j)
    spl.sub<-as.matrix(spl.t.std[1:ts.length.vec[j], ,drop=FALSE])
    preds.nobias[i,1:ts.length.vec[j],j]<-(spl.sub[,2, drop=FALSE] %*% betas.m2[i,j, drop=FALSE] 
                                    +spl.sub[,3, drop=FALSE] %*% betas.m3[i,j, drop=FALSE] 
                                    +spl.sub[,4, drop=FALSE] %*% betas.m4[i,j, drop=FALSE] 
                                    +spl.sub[,5, drop=FALSE] %*% betas.m5[i,j, drop=FALSE]  )
  }
}
preds.nobias.q<-as.data.frame(t(apply(preds.nobias,c(2,3),quantile, probs=c(0.5),na.rm=TRUE)))
preds.nobias.ucl<-as.data.frame(apply(preds.nobias,c(2,3),quantile, probs=c(0.975),na.rm=TRUE))
preds.nobias.lcl<-as.data.frame(apply(preds.nobias,c(2,3),quantile, probs=c(0.025),na.rm=TRUE))

saveRDS(preds.nobias, file=paste0(output_directory,"reg_mean_with_pooling bsplines nobias.rds"))



##small multiples plot
country=beta.labs.extract2[,1]
state=beta.labs.extract2[,2]
t=beta.labs.extract2[,3]
preds.nobias.q.alt<-apply(preds.nobias,c(2,3),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
preds.nobias.q.alt<-preds.nobias.q.alt[,,order(country,state )]  #sort by country, then state
pred.labs<-dimnames(preds.nobias.q.alt)[[3]] # labels from matrix
pred.labs.extract<- sub(".*\\[(.*)\\].*", "\\1", pred.labs, perl=TRUE)    #extract index numbers from the labels
pred.labs.extract2<-matrix(as.numeric(unlist(strsplit(pred.labs.extract, ","))),ncol=3, byrow=TRUE) #format into a matrix
country2<-pred.labs.extract2[,1] 
state2<-pred.labs.extract2[,2]
country.labs2<-c('US', 'Fiji','Denmark' ,'Brazil', 'Mexico', 'Ecuador', 'Chile')
cols.plot<- c('#e41a1c','#377eb8','#4daf4a','#984ea3', '#ff7f00','#a65628','#f781bf')
shiny.ds<-list(preds.nobias.q.alt,country2 ,state2,cols.plot,country.labs2)
names(shiny.ds)<-c('preds.nobias.q.alt','country2' ,'state2','cols.plot','country.labs2')
saveRDS(shiny.ds,paste0(output_directory,'shiny.small.multiples.rds'))
tiff(paste0(output_directory,'small multiple bsplines.tiff'), width = 12, height = 9, units = "in",res=200)
par(mfrow=c(10,10) , mar=c(0.5,1,0.6,0))
for (i in 1:dim(preds.nobias.q.alt )[3]){
  y.all<- t(exp(preds.nobias.q.alt[,,i]))
  keep.obs<-as.data.frame(complete.cases(y.all))[,1]
  y.all<-y.all[keep.obs,]
  n.times<- nrow(y.all)
  xx<- c(1:n.times, rev(1:n.times))
  yy<- c(y.all[,1], rev(y.all[,3] ))
  matplot(t(exp(preds.nobias.q.alt[,,i])), bty='l', type='l', col='white', lty=c(2,1,2), xlim=c(0,48), ylim=c(0.5,2), axes=F)
  polygon(xx, yy, col='gray90', border=NA)
  points( y.all[,2], col=cols.plot[country2[i]], type='l')
  abline(h=c(0.8,1,1.2), lty=c(3,2,3), col=c('gray','black','gray'), lwd=c(0.5,1,0.5)) 
  Axis(side=1, labels=FALSE, at=c(0,12,24,36,48), tck=-0.005,  col='gray')
  Axis(side=2, labels=FALSE, at=c(0.5,1,2.0), tck=-0.005, col='gray')
  if(state2[i]==1){
    title(country.labs2[country2[i]], col.main=cols.plot[country2[i]] )
    }
}
dev.off()



#########################
##PLOTS FOR FITTED VALUES, INCLUDING INTERCEPT
tiff(paste0(output_directory,'ucl and LCL smooth with pooling josh fix all states bsplines.tiff'), width = 7, height = 4, units = "in",res=200)
    par(mfrow=c(1,2))
    matplot(preds.ucl, type='l', bty='l',col='gray')
    title("Upper CrI")
    abline(h=0, col='red',lty=3)
    matplot(preds.lcl, type='l', bty='l',col='gray')
    title("Lower CrI")
    abline(h=0, col='red',lty=3)
dev.off()

country=beta.labs.extract2[,1]
state=beta.labs.extract2[,2]
t=beta.labs.extract2[,3]
tiff(paste0(output_directory, 'fix all states bspline.tiff'), width = 10.5, height = 4, units = "in",res=200)
    par(mfrow=c(1,N.countries))
    for(i in 1: N.countries){
      print(i)
      ds.select<-t(preds.q[country==i,,drop=FALSE])
      
      N.states<-ncol(ds.select)
      log_rr_prec_state_mean<-rep(NA, times=N.states)
      for (m in 1: N.states){
        log_rr_prec_state_mean[m]<-mean(diag(log_rr_prec_all[,,m,i]), na.rm=TRUE)
      }
      inv.var.ave.scale<-log_rr_prec_state_mean/max(log_rr_prec_state_mean, na.rm=TRUE)
      inv.var.ave.scale[is.na(inv.var.ave.scale)]<-0
      col.plot<-rgb(0,0,0, alpha=inv.var.ave.scale)
      matplot( 1:ts.length[i],ds.select[1:ts.length[i],], col=col.plot,lty=1, ylim=c(-1.0,1.0), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
      abline(h=0)
      title(countries[i])
    }
dev.off()   
      
tiff(paste0(output_directory,'ucl and LCL smooth with pooling josh fix all states bsplines nobias.tiff'), width = 7, height = 4, units = "in",res=200)
      par(mfrow=c(1,2)) 
      matplot(preds.nobias.ucl, type='l', bty='l',col='gray')
      title("Upper CrI")
      abline(h=0, col='red',lty=3)
      matplot(preds.nobias.lcl, type='l', bty='l',col='gray')
      title("Lower CrI")
      abline(h=0, col='red',lty=3)
      
dev.off()
  
tiff(paste0(output_directory,'/smooth with pooling josh fix all states bspline nobias.tiff'), width = 10.5, height = 4, units = "in",res=200)
  country=beta.labs.extract2[,1]
  state=beta.labs.extract2[,2]
  t=beta.labs.extract2[,3]
    par(mfrow=c(1,N.countries))
  for(i in 1: N.countries){
    print(i)
    ds.select<-t(preds.nobias.q[country==i,,drop=FALSE])
    
    N.states<-ncol(ds.select)
    log_rr_prec_state_mean<-rep(NA, times=N.states)
    for (m in 1: N.states){
      log_rr_prec_state_mean[m]<-mean(diag(log_rr_prec_all[,,m,i]), na.rm=TRUE)
    }  
        inv.var.ave.scale<-log_rr_prec_state_mean/max(log_rr_prec_state_mean, na.rm=TRUE)
        inv.var.ave.scale[is.na(inv.var.ave.scale)]<-0
        col.plot<-rgb(0,0,0, alpha=inv.var.ave.scale)
        matplot( 1:ts.length[i],ds.select[1:ts.length[i],], col=col.plot,lty=1, ylim=c(-1.0,1.0), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
        abline(h=0)
        title(countries[i])
         }
dev.off()


# 
# tiff('C:/Users/dmw63/Dropbox (Personal)/meta_analysis_results/stack all states/smooth with pooling josh fix all states bspline nobias plus cis.tiff', width = 10.5, height = 4, units = "in",res=200)
# country=beta.labs.extract2[,1]
# state=beta.labs.extract2[,2]
# t=beta.labs.extract2[,3]
# par(mfrow=c(1,N.countries))
# for(i in 1: N.countries){
#   print(i)
#   ds.select<-t(preds.nobias.q[country==i,,drop=FALSE])
#   ds.select.ucl<- t(preds.nobias.ucl[country==i,,drop=FALSE]) 
#   ds.select.lcl<- t(preds.nobias.lcl[country==i,,drop=FALSE]) 
#   N.states<-ncol(ds.select)
#   log_rr_prec_state_mean<-rep(NA, times=N.states)
#   for (m in 1: N.states){
#     log_rr_prec_state_mean[m]<-mean(diag(log_rr_prec_all[,,m,i]), na.rm=TRUE)
#   }
#   inv.var.ave.scale<-log_rr_prec_state_mean/max(log_rr_prec_state_mean, na.rm=TRUE)
#   inv.var.ave.scale[is.na(inv.var.ave.scale)]<-0
#   col.plot<-rgb(0,0,0, alpha=inv.var.ave.scale)
#   matplot( 1:ts.length[i],ds.select[1:ts.length[i],], col='black',lty=1, ylim=c(-1.0,1.0), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
#   matplot( 1:ts.length[i],ds.select.ucl[1:ts.length[i],], col='gray',lty=2, ylim=c(-1.0,1.0), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
#   matplot( 1:ts.length[i],ds.select.lcl[1:ts.length[i],], col='gray',lty=2, ylim=c(-1.0,1.0), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
#   
#     abline(h=0)
#   title(countries[i])
# }
# dev.off()

save(list = ls(all.names = TRUE),file=paste0(output_directory,"pooling no covars all states b splines.RData"),envir = .GlobalEnv)



#Project estimates back onto original time series
##Extract Fiji
fiji.dates<-seq.Date(    from=as.Date('2012-09-01'),to= as.Date('2015-12-01'), by='month')
fiji.preds<-t(preds[,1:length(fiji.dates),which(countries=="Fiji")])
fiji.raw<-read.csv('C:\\Users\\DMW63\\Desktop\\My documents h\\GATES\\CausalImpact code\\Input data\\Fiji\\prelog_Fiji_processed_data4.csv')
fiji.raw2<-fiji.raw[fiji.raw$age_group=='ac_pneumonia_Age1_0',]
fill.length<- nrow(fiji.raw2) - length(fiji.dates)
pre.fill<-matrix(0, nrow=fill.length, ncol=ncol(fiji.preds))
fiji.preds.fill<-exp(rbind(pre.fill,fiji.preds ))
counterfact<- 1/fiji.preds.fill * fiji.raw2$pneuvar
countfact.ci<-t(apply(counterfact,1,quantile, probs=c(0.025,0.5,0.975)))
countfact.ci[1:fill.length,] <-NA
#Plot the observed vs counterfactual cases
matplot(countfact.ci, type='l', col='gray', lty=1)
points(fiji.raw2$pneuvar)
dev.off()

#Plot the RR with CIs
preds.ci<-t(apply(fiji.preds,1,quantile, probs=c(0.025,0.5,0.975)))
matplot(exp(preds.ci), type='l', col='gray', lty=1, bty='l')
abline(h=1, lty=2, col='gray')

#Aggregate the counterfactual and observed by year
 fiji.years<-year(fiji.dates)
 fiji.preds.year.split<-lapply(split(counterfact[-c(1:fill.length),], fiji.years), matrix , ncol=ncol(fiji.preds)  )
 fiji.preds.year.sum<-sapply(fiji.preds.year.split , function(x) apply(x,2, sum) ) 
 fiji.preds.year.sum.length<-sapply(fiji.preds.year.split , function(x) nrow(x) ) 
 fiji.preds.year.sum.q<-t(apply(fiji.preds.year.sum , 2, quantile, probs=c(0.025,0.5,0.975)) )*12/fiji.preds.year.sum.length
obs.year<- aggregate( fiji.raw2$pneuvar[-c(1:fill.length)], by=list(fiji.years), sum)[,2]*12/fiji.preds.year.sum.length
tiff(paste0(output_directory,'/project onto original aggregated time series fiji.tiff'), width = 10.5, height = 4, units = "in",res=200)
 par(mfrow=c(1,1))
  matplot(unique(fiji.years), fiji.preds.year.sum.q, type='l', col='gray', bty='l', ylim=c(0, max(fiji.preds.year.sum.q)))
  points(unique(fiji.years), obs.year, pch=16, col='black')
  dev.off()
  #