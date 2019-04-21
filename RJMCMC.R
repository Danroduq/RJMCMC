#first we assume k fixed at 10
rm(list = ls())
library(boot)#library  that contains the coal miner data
library(lubridate)#library that changes decimal dates


lambda=3
kmax=30
alpha=1
beta=200
c=1
k=10
tempdate=lapply(coal,FUN=date_decimal)
tempdate=sapply(tempdate,FUN=as_date)
coal1=difftime(as.Date(tempdate,format="%Y-%m-%d",origin="1970-01-01"),"1851-01-01", units = c("days"))
L=as.numeric(difftime("1962/12/31" ,"1851/01/01", units = c("days")))

pk=function(k,lambda,kmax){
  pk=dpois(k,lambda)/ppois(kmax,lambda)
  return(pk)
}
s_sampler=function(k,L){
  samp=runif(2*k+1,min=0,max=L)
  samp=samp[order(samp)][seq(2,2*k+1,2)]
  return(samp)
}
evaluate=function(x,s,h){
  ind=min(which(x<=c(0,s,L)))
  return(h[ind-1])
}
evaluatev=Vectorize(evaluate,vectorize.args=c("x"))
loglik=function(s,h){
  sapply(coal1,FUN=evaluatev,s=s_init,h=h_init)
  loglikhood=sum(log(h))-sum(diff(c(0,s_init,L))*h_init)
  return(loglikhood)
}

s_init=s_sampler(k,L)
h_init=rgamma(k+1,shape=alpha,rate=beta)

#finding the c value as defined in paper
  c=seq(0.01,1,0.01)
  min=rep(-1,31)
  for (k in 0:30){
    bk=c*min(1,pk(k+1,lambda,kmax)/pk(k,lambda,kmax))
    dk=c*min(1,pk(k-1,lambda,kmax)/pk(k,lambda,kmax))
    #note that the relationship between c and the total is linear
    tot=bk+dk
    c_index=max(which(tot<=0.9))
    min[k+1]=c[c_index]
  }
  c=min(min)
  k=10
plot(c,tot)

like=loglik(s=s_init,h_init)

#Initializing Chain Data Structures
nsamp<-5000; nburn<-2000; nthin<-1
nits<-nburn+nsamp*nthin
ico<-0
s.samp<-array(0,c(nsamp,1,k+1))
hsamp<-array(0,c(nsamp,1,k))
k.samp<-matrix(0,nrow=nsamp,ncol=1)

check=matrix(0,nits,2)
mh.sig<-c(0.1,0.05)
for(iter in 1:nits){
   new.ssq<-new.prior<-0
  for(k in 1:2){
    new.theta<-old.theta
    new.theta[k]<-old.theta[k]+rnorm(1)*mh.sig[k]
    
    new.ssq<-sum((y-mu.func(x,new.theta[1],new.theta[2]))^2)
    new.prior<-dmvn(new.theta,mu,Sigma,log=T)
    
    if(log(runif(1)) < (-0.5*new.ssq/old.sigsq+new.prior)+
       (0.5*old.ssq/old.sigsq-old.prior)){
      check[iter,k]=1
      old.theta[k]<-new.theta[k]
      old.ssq<-new.ssq
      old.prior<-new.prior
    }
  }
  
  
  #Update sigmas
  a.n<-a+12
  b.n<-old.ssq+b
  old.sigsq<-1/rgamma(1,a.n/2,b.n/2)
  
  
  if(iter > nburn ){
    ico<-ico+1
    theta.samp[ico,,]<-old.theta
    sigsq.samp[ico,]<-old.sigsq
  }
}





