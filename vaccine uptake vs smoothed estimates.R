library(plyr)

output_directory<- paste0(dirname(getwd()), "/Results CP/")
output_directory_nopool<- paste0(dirname(getwd()), "/Results nopool/")

countries<-c('US', 'Fiji','Denmark' ,'Brazil_state', 'Mexico_state', 'Ecuador_state', 'Chile_state')
preds.pool<-readRDS( file=paste0(output_directory,"reg_mean_with_pooling cp nobias.rds"))
preds.no.pool<-readRDS( file=paste0(output_directory_nopool,"reg_mean_with_pooling bsplines nobias.rds"))
preds.no.pool.bias<-readRDS( file=paste0(output_directory_nopool,"reg_mean_with_pooling bsplines.rds"))
 
state.labels<-readRDS( file=paste0(output_directory_nopool,"state labels.rds"))

cov1.br<-read.csv("C:/Users/dmw63/Dropbox (Personal)/meta_analysis_results/stack all states/PCVCoverage_Brazil_paho.csv")
cov1.br$state.id<-as.numeric(substr(as.character(cov1.br$state),1,2))
cov1.br<-cov1.br[complete.cases(cov1.br),]

#Sort the coverage data to be in same order as the hospitalization data, get rid of states not analyzed 
cov1.ec<-read.csv("C:/Users/dmw63/Dropbox (Personal)/meta_analysis_results/stack all states/PCVCoverage_Ecuador_paho.csv")
state.code.ec<-state.labels[,countries=='Ecuador_state'][!is.na(state.labels[,countries=='Ecuador_state'] )]
state.code.ec<-substring(state.code.ec,first=3)
cov1.ec.2013<-cov1.ec[cov1.ec$year==2013,]
order.ec<- pmatch(as.character(cov1.ec.2013$province.code) , state.code.ec)
cov1.ec.2013<-cov1.ec.2013[order(order.ec),]
cov1.ec.2013<-cov1.ec.2013[as.character(cov1.ec.2013$province.code) %in% state.code.ec,]

#Mexico coverage
cov1.mx<-read.csv("C:/Users/dmw63/Dropbox (Personal)/meta_analysis_results/stack all states/PCVCoverage_Mexico_paho.csv")
state.code.mx<-state.labels[,countries=='Mexico_state'][!is.na(state.labels[,countries=='Mexico_state'] )]
state.code.mx<-substring(state.code.mx,first=3)
cov1.mx.2013<-cov1.mx[cov1.mx$Year==2010,]
order.mx<- pmatch(as.character(cov1.mx.2013$State.Code) , state.code.mx)
cov1.mx.2013<-cov1.mx.2013[order(order.mx),]
cov1.mx.2013<-cov1.mx.2013[as.character(cov1.mx.2013$State.Code) %in% state.code.mx,]

beta.labs.extract2<-matrix(as.numeric(unlist(strsplit(dimnames(preds.pool)[[3]], ","))),ncol=3, byrow=TRUE) #format into a matrix
country=countries[beta.labs.extract2[,1]]
state=beta.labs.extract2[,2]

preds.pool.q<-aaply(preds.pool, c(2,3,4), .fun=quantile, na.rm=TRUE, probs=c(0.025,0.5,0.975))
preds.no.pool.q<-aaply(preds.no.pool, c(2,3), .fun=quantile, na.rm=TRUE, probs=c(0.025,0.5,0.975))
preds.no.pool.bias.q<-aaply(preds.no.pool.bias, c(2,3), .fun=quantile, na.rm=TRUE, probs=c(0.025,0.5,0.975))


par(mfrow=c(1,1))
plot(preds.no.pool.q[49,,2],preds.pool.q[49,,2], bty='l', ylim=c(-15,15), xlim=c(-15,15))
abline(a=0, b=1, col='grey', lty=2)
abline(h=0, v=0, col='grey', lty=2)

plot(preds.no.pool.q[49,,2],preds.pool.q[49,,2], bty='l', ylim=c(-1,1), xlim=c(-1,1))
abline(a=0, b=1, col='grey', lty=2)
abline(h=0, v=0, col='grey', lty=2)

########################
#Brazil only
#########################
preds.pool.q.Br<-preds.pool.q[,country=="Brazil_state",]
matplot(exp(preds.pool.q.Br[,,2]), type='l', col='gray', lty=1, bty='l')
abline(h=c(0.8,1,1.2), lty=3)
br.state<-state[country=="Brazil_state"]
preds.no.pool.q.Br<-preds.no.pool.q[,country=="Brazil_state",]
preds.no.pool.bias.Br<-preds.no.pool.bias.q[,country=="Brazil_state",]
#Merge with vax cov--some states don't have data for all 4 years post
#Need to chekc this--assumes that states analyzed in numerical order n 1st stage model)
#POOLING VS COVRAGE
par(mfrow=c(1,3))
plot(cov1.br$X2013 ,exp(preds.pool.q.Br[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
title("pooled estimate vs coverage")
#summary(lm(preds.pool.q.Br[37,,2]~ cov1.br$X2013 ))
##NO POOLING VS COVERAGE
plot(cov1.br$X2013 ,exp(preds.no.pool.q.Br[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
title("Unpooled estimate vs coverage")
##NO POOLING, NO UNBIASING VS COVERAGE
plot(cov1.br$X2013 ,exp(preds.no.pool.bias.Br[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
title("Unpooled estimate (without unbiasing) vs coverage")

#summary(lm(preds.no.pool.bias.Br[37,,2]~ cov1.br$X2013 ))
cor( cbind( cov1.br$X2013,preds.pool.q.Br[37,,2],preds.no.pool.q.Br[37,,2],preds.no.pool.bias.Br[37,,2]), method='spearman')

########################
#Ecuador only
#########################
preds.pool.q.ec<-preds.pool.q[,country=="Ecuador_state",]
matplot(preds.pool.q.ec[,,2], type='l', col='gray', lty=1, bty='l')
abline(h=0, lty=3)
ec.state<-state[country=="Ecuador_state"]
preds.no.pool.q.ec<-preds.no.pool.q[,country=="Ecuador_state",]

preds.no.pool.bias.ec<-preds.no.pool.bias.q[,country=="Ecuador_state",]
#Merge with vax cov--some states don't have data for all 4 years post
#Need to chekc this--assumes that states analyzed in numerical order n 1st stage model)
#POOLING VS COVRAGE
plot(cov1.ec.2013$coverage ,exp(preds.pool.q.ec[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
#summary(lm(preds.pool.q.ec[37,,2]~ cov1.ec$X2013 ))
##NO POOLING VS COVERAGE
plot(cov1.ec.2013$coverage ,exp(preds.no.pool.q.ec[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
#summary(lm(preds.no.pool.q.ec[37,,2]~ cov1.ec$X2013 ))

##NO POOLING, NO UNBIASING VS COVERAGE
plot(cov1.ec.2013$coverage ,exp(preds.no.pool.bias.ec[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
#summary(lm(preds.no.pool.bias.ec[37,,2]~ cov1.ec$X2013 ))
cor( cbind( cov1.ec.2013$coverage,preds.pool.q.ec[37,,2],preds.no.pool.q.ec[37,,2],preds.no.pool.bias.ec[37,,2]), method='spearman')


########################
#Mexico only
#########################
preds.pool.q.mx<-preds.pool.q[,country=="Mexico_state",]
matplot(preds.pool.q.mx[,,2], type='l', col='gray', lty=1, bty='l', ylim=c(-2,1))
abline(h=0, lty=3)
mx.state<-state[country=="Mexico_state"]
preds.no.pool.q.mx<-preds.no.pool.q[,country=="Mexico_state",]
matplot(preds.no.pool.q.mx[,,2], type='l', col='gray', lty=1, bty='l', ylim=c(-2,1))
abline(h=0, lty=3)

preds.no.pool.bias.mx<-preds.no.pool.bias.q[,country=="Mexico_state",]
#Merge with vax cov--some states don't have data for all 4 years post
#Need to chekc this--assumes that states analyzed in numerical order n 1st stage model)
#POOLING VS COVRAGE
plot(cov1.mx.2013$Coverage..total.3rd.dose.total.eligible.children. ,exp(preds.pool.q.mx[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
#summary(lm(preds.pool.q.mx[37,,2]~ cov1.mx$X2013 ))
##NO POOLING VS COVERAGE
plot(cov1.mx.2013$Coverage..total.3rd.dose.total.eligible.children. ,exp(preds.no.pool.q.mx[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
#summary(lm(preds.no.pool.q.mx[37,,2]~ cov1.mx$X2013 ))

##NO POOLING, NO UNBIASING VS COVERAGE
plot(cov1.mx.2013$Coverage..total.3rd.dose.total.eligible.children. ,exp(preds.no.pool.bias.mx[37,,2]), type='p', col='black', lty=1, bty='l', ylab="Rate Ratio", xlab="Uptake (%)")
abline(h=1, col='gray', lty=2)
#summary(lm(preds.no.pool.bias.mx[37,,2]~ cov1.mx$X2013 ))
cor( cbind( cov1.mx.2013$Coverage..total.3rd.dose.total.eligible.children.,preds.pool.q.mx[37,,2],preds.no.pool.q.mx[37,,2],preds.no.pool.bias.mx[37,,2]), method='spearman')

