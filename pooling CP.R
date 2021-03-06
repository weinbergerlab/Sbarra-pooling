
rm(list=ls(all=TRUE))
require(reshape)
require(reshape2)
require(abind)
library(mgcv)   
library(R2jags) 
library(splines)
library(lubridate)

output_directory<- paste0(dirname(getwd()), "/Results CP no hdi/")
ifelse(!dir.exists(output_directory), dir.create(output_directory), FALSE)

####SET INPUTS PARAMETERS#############################################################################################################
#setwd('C:/Users/dmw63/Dropbox (Personal)/meta_analysis_results/stack all states/')
#countries<-c('US', 'Fiji','Denmark' ,'Brazil_state', 'Mexico_state', 'Ecuador_state', 'Chile_state')
#####################################################################################################################################
#format data
source('./first stage estimates/run all markdown files.R') #Read in data
source('models/model.R')  #Read in model

#Run Model
model_jags<-jags.model(textConnection(model_string),
                      data=list('n.countries' = N.countries, 
                                'n.states' = N.states, 
                                 'w_hat' = log_rr_q_all,
                                 'w_hat_prec' = log_rr_prec_diag_all, 
                                'ts.length' = ts.length_mat,
                                'I_Omega'= I_Omega,
                                'max.time.points'=max.time.points,
                                'time.index'=time.index,
                                'I_Sigma'=I_Sigma), n.chains=1, n.adapt=10000) 

#Posterior Sampling
update(model_jags, 
       n.iter=10000)  #Burnin for 10,000 Samples
# dic.mod1a <- dic.samples(model_jags, 5000,"pD")  #DIC=2161 with Brazil state, Mexico state, US
# dic.mod1a

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("reg_mean",'reg_unbias' 
                                                 ,'cp1','cp2',"beta", 
                                                 'gamma','lambda'),
                                thin=10,
                                n.iter=50000)
#plot(posterior_samples, ask=TRUE)
saveRDS(posterior_samples,paste0(output_directory,'posterior_samples_CP.rds'))
#posterior_samples<-readRDS(paste0(output_directory,'posterior_samples_CP.rds'))


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


x.func <- function(x){ sub(".*\\[(.*)\\].*", "\\1", x, perl=TRUE)} 

#############################################################################
##### create posterior estimates into data frames (w.true and reg.mean) #####
#############################################################################
beta1<-posterior_samples[[1]][,grep("^beta.*,1]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
beta1.lab<-x.func(dimnames(beta1)[[2]]) 
beta1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(beta1.lab.spl)<-c('country','state')

#Extract first change point time:
beta3<-posterior_samples[[1]][,grep("^beta.*,3]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
cp1<-max.time.points*exp(beta3) #Intercept
beta3.lab<-x.func(dimnames(beta3)[[2]]) 
beta3.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta3.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(beta3.lab.spl)<-c('country','state')
for(i in c(1:9)){hist(beta3[,i])}
quant.cp1<-t(max.time.points*exp(apply(beta3,2,quantile, probs=c(0.025,0.5,0.975))))
var.cp1<-apply(beta3,2,var)
cp1.lab<-x.func(dimnames(quant.cp1)[[1]]) 
cp1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(cp1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(cp1.lab.spl)<-c('country','state')
par(mfrow=c(2,2))
plot(y=1:nrow(quant.cp1), x=quant.cp1[,'50%'], bty='l')
library(ggplot2)
plot.data<-cbind.data.frame('strata'=1:nrow(quant.cp1),'median.cp'=quant.cp1[,'50%'], 'inv.var.cp'=1/var.cp1,beta3.lab.spl)
plot.data<-plot.data[order(plot.data$country),]
plot.data$order2<-1:nrow(plot.data)
plot.data$country2<-NA
for(i in 1:length(countries)){plot.data$country2[plot.data$country==i] <-countries[i] }
tiff(paste0(output_directory,"_cp.tiff"))
plot1<-ggplot(data=plot.data, aes(x=median.cp, y=order2, color=country2)) +
  geom_point(aes(size=inv.var.cp)) +
  scale_size_continuous(range=c(1,15)) +
  theme_bw()+
  guides( size = FALSE)+
  # theme(legend.position = "none")+
  scale_color_manual(values=c('#7fc97f','#beaed4', '#fdc086','#ffff99','#386cb0', '#f0027f','#bf5b17'))
print(plot1)
dev.off()
############################
#Extract slope
beta2<-posterior_samples[[1]][,grep("^beta.*,2]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
slp1<-beta2 #Intercept
beta2.lab<-x.func(dimnames(beta2)[[2]]) 
beta2.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta2.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(beta2.lab.spl)<-c('country','state')
for(i in c(1:9)){hist(beta2[,i])}
quant.slp1<-t(apply(beta2,2,quantile, probs=c(0.025,0.5,0.975)))
var.slp1<-apply(beta2,2,var)
slp1.lab<-x.func(dimnames(quant.slp1)[[1]]) 
slp1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(slp1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(slp1.lab.spl)<-c('country','state')
par(mfrow=c(2,2))
plot(y=1:nrow(quant.slp1), x=quant.slp1[,'50%'], bty='l')
library(ggplot2)
plot.data<-cbind.data.frame('strata'=1:nrow(quant.slp1),'median.slp'=quant.slp1[,'50%'], 'inv.var.slp'=1/var.slp1,beta2.lab.spl)
plot.data<-plot.data[order(plot.data$country),]
plot.data$order2<-1:nrow(plot.data)
plot.data$country2<-NA
for(i in 1:length(countries)){plot.data$country2[plot.data$country==i] <-countries[i] }
tiff(paste0(output_directory,"_slope.tiff"))
plot2<-ggplot(data=plot.data, aes(x=median.slp, y=order2, color=country2)) +
  geom_point(aes(size=inv.var.slp)) +
  scale_size_continuous(range=c(1,15)) +
  theme_bw()+
  guides( size = FALSE)+
  # theme(legend.position = "none")+
  scale_color_manual(values=c('#7fc97f','#beaed4', '#fdc086','#ffff99','#386cb0', '#f0027f','#bf5b17'))
print(plot2)
dev.off()
##########################################

#Relationship between change point location and slope
tiff(paste0(output_directory,"_slope vs cp.tiff"))
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(quant.slp1[,'50%'], quant.cp1[,'50%'], pch=16, col=beta2.lab.spl$country,bty='l', ylab="Change point", xlab="Slope")
dev.off()

##melt and cast predicted values into 4D array N,t,i,j array
reg_mean<-posterior_samples[[1]][,grep("reg_mean",dimnames(posterior_samples[[1]])[[2]])]
pred1<-reg_mean
pred.indices<- x.func(dimnames(pred1)[[2]]) 
pred.indices.spl<-matrix(unlist(strsplit(pred.indices, ',',fixed=TRUE)), ncol=3, byrow=TRUE)
pred.indices.spl<-as.data.frame(pred.indices.spl)
names(pred.indices.spl)<- c('country','state','time')
reg_mean2<-cbind.data.frame(pred.indices.spl,t(reg_mean))
reg_mean_m<-melt(reg_mean2, id=c('time','country','state'))
reg_mean_c<-acast(reg_mean_m, variable~time~country~state)
preds.q<-apply(reg_mean_c,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
dimnames(preds.q)[[2]]<- as.numeric(as.character(dimnames(preds.q)[[2]]))

#Unbias
reg_unbias<-posterior_samples[[1]][,grep("reg_unbias",dimnames(posterior_samples[[1]])[[2]])]
pred1<-reg_unbias
pred.indices<- x.func(dimnames(pred1)[[2]]) 
pred.indices.spl<-matrix(unlist(strsplit(pred.indices, ',',fixed=TRUE)), ncol=3, byrow=TRUE)
pred.indices.spl<-as.data.frame(pred.indices.spl)
names(pred.indices.spl)<- c('country','state','time')
reg_unbias2<-cbind.data.frame(pred.indices.spl,t(reg_unbias))
reg_unbias_m<-melt(reg_unbias2, id=c('time','country','state'))
reg_unbias_c<-acast(reg_unbias_m, variable~time~country~state)
reg_unbias_c<-reg_unbias_c[,order(as.numeric(dimnames(reg_unbias_c)[[2]])),order(as.numeric(dimnames(reg_unbias_c)[[3]])),order(as.numeric(dimnames(reg_unbias_c)[[4]]))]
preds.unbias.q<-apply(reg_unbias_c,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
dimnames(preds.unbias.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))

######################
for(i in c(1:length(countries))){
  tiff(paste0(output_directory,countries[i],"RR trajectory by state.tiff"))
  par(mfrow=c(5,6), mar=c(4,2,1,1))
    for(j in 1:N.states[i]){
        plot.data<-t(preds.unbias.q[,,i,j])
        plot.data2<-plot.data[complete.cases(plot.data),]
        final.rr<-paste0(round(exp(plot.data2[nrow(plot.data2),'50%', drop=F]),2),
                         ' (' ,round(exp(plot.data2[nrow(plot.data2),'2.5%', drop=F]),2),',',
                         round(exp(plot.data2[nrow(plot.data2),'97.5%', drop=F]),2),")")
        start.post<-which( ((1:tot_time)-pre.vax.time)==0)
        matplot( ((1:tot_time)-pre.vax.time), plot.data,type='l',yaxt='n', xlab='Months post-PCV',xlim=c(0, max.time.points),  
                    ylim=c(-1,1), col='gray', lty=c(2,1,2), bty='l')
        abline(h=0)
        axis(side=2, at=c(-1,-0.5,0,5,1), las=1,labels=round(exp(c(-1,-0.5,0,5,1)),1 ))
       # abline(v=0)
       # title(countries[i])
        text(0, 1, final.rr, xpd=NA, offset=0, pos=4)
      
    }
  dev.off()
}
saveRDS(preds.unbias.q, file=paste0(output_directory,"reg_mean_with_pooling cp nobias.rds"))

# 
# ##small multiples plot
# country=pred.indices.spl[,1]
# state=pred.indices.spl[,2]
# t=pred.indices.spl[,3]
# preds.nobias.q.alt<-apply(preds.nobias,c(2,3),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
# preds.nobias.q.alt<-preds.nobias.q.alt[,,order(country,state )]  #sort by country, then state
# pred.labs<-dimnames(preds.nobias.q.alt)[[3]] # labels from matrix
# pred.labs.extract<- sub(".*\\[(.*)\\].*", "\\1", pred.labs, perl=TRUE)    #extract index numbers from the labels
# pred.labs.extract2<-matrix(as.numeric(unlist(strsplit(pred.labs.extract, ","))),ncol=3, byrow=TRUE) #format into a matrix
# country2<-pred.labs.extract2[,1] 
# state2<-pred.labs.extract2[,2]
# country.labs2<-c('US', 'Fiji','Denmark' ,'Brazil', 'Mexico', 'Ecuador', 'Chile')
# cols.plot<- c('#e41a1c','#377eb8','#4daf4a','#984ea3', '#ff7f00','#a65628','#f781bf')
# shiny.ds<-list(preds.nobias.q.alt,country2 ,state2,cols.plot,country.labs2)
# names(shiny.ds)<-c('preds.nobias.q.alt','country2' ,'state2','cols.plot','country.labs2')
# saveRDS(shiny.ds,paste0(output_directory,'shiny.small.multiples.rds'))
# tiff(paste0(output_directory,'small multiple CP.tiff'), width = 12, height = 9, units = "in",res=200)
# par(mfrow=c(10,10) , mar=c(0.5,1,0.6,0))
# for (i in 1:dim(preds.nobias.q.alt )[3]){
#   y.all<- t(exp(preds.nobias.q.alt[,,i]))
#   keep.obs<-as.data.frame(complete.cases(y.all))[,1]
#   y.all<-y.all[keep.obs,]
#   n.times<- nrow(y.all)
#   xx<- c(1:n.times, rev(1:n.times))
#   yy<- c(y.all[,1], rev(y.all[,3] ))
#   matplot(t(exp(preds.nobias.q.alt[,,i])), bty='l', type='l', col='white', lty=c(2,1,2), xlim=c(0,48), ylim=c(0.5,2), axes=F)
#   polygon(xx, yy, col='gray90', border=NA)
#   points( y.all[,2], col=cols.plot[country2[i]], type='l')
#   abline(h=c(0.8,1,1.2), lty=c(3,2,3), col=c('gray','black','gray'), lwd=c(0.5,1,0.5)) 
#   Axis(side=1, labels=FALSE, at=c(0,12,24,36,48), tck=-0.005,  col='gray')
#   Axis(side=2, labels=FALSE, at=c(0.5,1,2.0), tck=-0.005, col='gray')
#   if(state2[i]==1){
#     title(country.labs2[country2[i]], col.main=cols.plot[country2[i]] )
#     }
# }
# dev.off()


# 
# #########################
# ##PLOTS FOR FITTED VALUES, INCLUDING INTERCEPT
# tiff(paste0(output_directory,'ucl and LCL smooth with pooling josh fix all states CP.tiff'), width = 7, height = 4, units = "in",res=200)
#     par(mfrow=c(1,2))
#     matplot(preds.ucl, type='l', bty='l',col='gray')
#     title("Upper CrI")
#     abline(h=0, col='red',lty=3)
#     matplot(preds.lcl, type='l', bty='l',col='gray')
#     title("Lower CrI")
#     abline(h=0, col='red',lty=3)
# dev.off()
# 
# country=beta.labs.extract2[,1]
# state=beta.labs.extract2[,2]
# t=beta.labs.extract2[,3]
# tiff(paste0(output_directory, 'fix all states bspline.tiff'), width = 10.5, height = 4, units = "in",res=200)
#     par(mfrow=c(1,N.countries))
#     for(i in 1: N.countries){
#       print(i)
#       ds.select<-t(preds.q[country==i,,drop=FALSE])
#       
#       N.states<-ncol(ds.select)
#       log_rr_prec_state_mean<-rep(NA, times=N.states)
#       for (m in 1: N.states){
#         log_rr_prec_state_mean[m]<-mean(diag(log_rr_prec_all[,,m,i]), na.rm=TRUE)
#       }
#       inv.var.ave.scale<-log_rr_prec_state_mean/max(log_rr_prec_state_mean, na.rm=TRUE)
#       inv.var.ave.scale[is.na(inv.var.ave.scale)]<-0
#       col.plot<-rgb(0,0,0, alpha=inv.var.ave.scale)
#       matplot( 1:ts.length[i],ds.select[1:ts.length[i],], col=col.plot,lty=1, ylim=c(-1.0,1.0), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
#       abline(h=0)
#       title(countries[i])
#     }
# dev.off()   
#       
# tiff(paste0(output_directory,'ucl and LCL smooth with pooling josh fix all states CP nobias.tiff'), width = 7, height = 4, units = "in",res=200)
#       par(mfrow=c(1,2)) 
#       matplot(preds.nobias.ucl, type='l', bty='l',col='gray')
#       title("Upper CrI")
#       abline(h=0, col='red',lty=3)
#       matplot(preds.nobias.lcl, type='l', bty='l',col='gray')
#       title("Lower CrI")
#       abline(h=0, col='red',lty=3)
#       
# dev.off()
#   
# tiff(paste0(output_directory,'/smooth with pooling josh fix all states bspline nobias.tiff'), width = 10.5, height = 4, units = "in",res=200)
#   country=beta.labs.extract2[,1]
#   state=beta.labs.extract2[,2]
#   t=beta.labs.extract2[,3]
#     par(mfrow=c(1,N.countries))
#   for(i in 1: N.countries){
#     print(i)
#     ds.select<-t(preds.nobias.q[country==i,,drop=FALSE])
#     
#     N.states<-ncol(ds.select)
#     log_rr_prec_state_mean<-rep(NA, times=N.states)
#     for (m in 1: N.states){
#       log_rr_prec_state_mean[m]<-mean(diag(log_rr_prec_all[,,m,i]), na.rm=TRUE)
#     }  
#         inv.var.ave.scale<-log_rr_prec_state_mean/max(log_rr_prec_state_mean, na.rm=TRUE)
#         inv.var.ave.scale[is.na(inv.var.ave.scale)]<-0
#         col.plot<-rgb(0,0,0, alpha=inv.var.ave.scale)
#         matplot( 1:ts.length[i],ds.select[1:ts.length[i],], col=col.plot,lty=1, ylim=c(-1.0,1.0), type='l', bty='l', ylab="Log(Rate ratio)", xlab="Months post-PCV")
#         abline(h=0)
#         title(countries[i])
#          }
# dev.off()

