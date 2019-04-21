#first we assume k fixed at 10
rm(list = ls())
library(boot)


lambda=3
kmax=30
alpha=1
beta=200
c=1
k=10
coal1=(coal-1851)*365
L=(1962-1851)*365

pk=function(k,lambda,kmax){
  pk=dpois(k,lambda)/ppois(kmax,lambda)
  return(pk)
}
s_sampler=function(k,L){
  samp=runif(2*k+1,min=0,max=L)
  samp=samp[order(samp)][seq(2,2*k+1,2)]
  return(samp)
}
loglik=function(s,h){
  loglikhood=sum(log(h))-sum(diff(c(0,s_init,L))*h_init)
  return(loglikhood)
}

s_init=s_sampler(k,L)
h_init=rgamma(k+1,shape=alpha,rate=beta)

bk=c*min(1,pk(k+1,lambda,kmax)/pk(k,lambda,kmax))
dk_plus1=c*min(1,pk(k,lambda,kmax)/pk(k+1,lambda,kmax))

like=loglik(s=s_init,h_init)

#Initializing Chain Data Structures
nsamp<-5000; nburn<-2000; nthin<-1
nits<-nburn+nsamp*nthin
ico<-0
s.samp<-array(0,c(nsamp,1,k+1))
hsamp<-array(0,c(nsamp,1,k))



