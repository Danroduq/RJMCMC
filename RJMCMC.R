#first we assume k fixed at 10
rm(list = ls())
library(boot)#library  that contains the coal miner data
library(lubridate)#library that changes decimal dates


lambda=3
kmax=30
alpha=1
beta=200
c=1
K=10
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
loglik=function(s,h,L){
  sapply(coal1,FUN=evaluatev,s=s,h=h)
  loglikhood=sum(log(h))-sum(diff(c(0,s,L))*h)
  return(loglikhood)
}
bkf=function(c,k,kmax,lambda){
  tr=c*min(1,pk(k+1,lambda,kmax)/pk(k,lambda,kmax))
  return(tr)
}
dkf=function(c,k,kmax,lambda){
  tr=c*min(1,pk(k-1,lambda,kmax)/pk(k,lambda,kmax))
  return(tr)
}

old.s=s_sampler(k,L)
oldh=rgamma(k+1,shape=alpha,rate=beta)

#finding the c value as defined in paper
  c=seq(0.01,1,0.01)
  min=rep(-1,31)
  for (k in 0:30){
    bk=bkf(c,k,kmax,lambda)
    dk=dkf(c,k,kmax,lambda)
    #note that the relationship between c and the total is linear
    tot=bk+dk
    c_index=max(which(tot<=0.9))
    min[k+1]=c[c_index]
  }
  c=min(min)

eta.density<-function(evec,Kv,xm){
  dvec<-diff(c(0,evec,xm))
  return(lgamma(2*Kv+2)+sum(log(dvec))-(2*Kv+1)*log(xm))
}

y<-bootstrap::cholost$y[order(cholost$z)]
x<-bootstrap::cholost$z[order(cholost$z)]
n<-length(y)
xmax<-100

prior.lam<-0.01
prior.a<-10
prior.b<-100
prior.gam<-10

old.K<-5  #This is the number of breaks, number of segments is K+1
eta<-sort(runif(2*old.K+1,0,100))
old.eta<-eta[2*(1:old.K)]

old.Xmat<-rep(1,n)
for(k in 1:old.K){
  old.Xmat<-cbind(old.Xmat,(x-old.eta[k])^r*(x>old.eta[k]))
}

old.np<-ncol(old.Xmat)
old.prior.L<-prior.lam*diag(1,old.np)
old.prior.be0<-rep(0,old.np)

XTX<-t(old.Xmat)%*%old.Xmat
be.hat<-solve(XTX) %*% (t(old.Xmat) %*% y)
y.fit<-old.Xmat %*% be.hat
sigsq.hat<-sum((y-y.fit)^2)/(n-old.np)
sig.hat<-sqrt(sigsq.hat)

old.be<-as.numeric(be.hat)
old.yfit<-old.Xmat %*% old.be
old.sig<-sig.hat

ysq<-t(y) %*% y
old.c<-(ysq + t(old.prior.be0) %*% old.prior.L %*% old.prior.be0 -
          t(post.mean) %*% post.L %*% post.mean)[1,1]

old.marg.like<-0.5*log(det(old.prior.L))-0.5*log(det(XTX+old.prior.L))-
  0.5*(n+prior.a)*log(0.5*(prior.b+old.c))

old.prior.K<-dpois(old.K,prior.gam,log=T)

old.prior.eta<-eta.density(old.eta,old.K,xmax)

nburn<-1000
nsamp<-10000
nthin<-20
nits<-nburn+nthin*nsamp

sig.post<-rep(0,nsamp)
be.post<-eta.post<-matrix(0,nrow=nsamp,ncol=100)
K.post<-rep(0,nsamp)
ico<-0

par(mar=c(4,4,2,0))
plot(x,y,type='p',pch=19,cex=0.5)
for(iter in 1:nits){
 
  bk=bkf(c,k,kmax,lambda)
  dk=dkf(c,k,kmax,lambda)
  split=1-bk-dk
  pik=split/2
  etak=split/2
  
  coin=runif(1)
  if(coin<bk){
    #birth 
    
    ss=runif(0,L)
    s_ind=min(which(ss<old.s))
    new.s=rep(0,old.K+1)
    new.s[1:s_ind]=old.s[1:s_ind]
    new.s[s_ind+1]=ss
    new.s[s_ind+2:old.K+1]=old.s[s_ind+1:old.K]
    u=runif(0,1)
   
    
    c=(old.s[s_ind+1]-old.s)*log(old.h[s_ind])
    hp1=exp(old.s-old.s)*loglog((1-u)/u)) 
    new.like=loglik(s,h,L)
    if(log(runif(1)) < (new.like+new.prior+new.proposal + new.Jacobian) -
       (old.like+old.prior+old.proposal)){
      
      old.K<-new.K
      old.s=new.s
      old.h=new.h
      
      old.like=new.like
      old.prior=new.prior
      old.proposal=new.proposal
      old.Jacobian=new.Jacobian
    }
    
  }else if (coin<=bk+dk){
    #Death
    if(old.K == 0) break #Prior prob on fewer than 0 breaks is zero
    
    new.K<-old.K-1
    k<-sample(1:old.K,size=1)
    new.eta<-old.eta[-k]
    if(new.K == 0){
      new.Xmat<-matrix(1,ncol=1,nrow=n)
    }else{
      new.Xmat<-rep(1,n)
      for(k in 1:new.K){
        new.Xmat<-cbind(new.Xmat,(x-new.eta[k])^r*(x>new.eta[k]))
      }
    }	
    new.np<-ncol(new.Xmat)
    new.prior.L<-prior.lam*diag(1,new.np)
    new.prior.be0<-rep(0,new.np)
    XTX<-t(new.Xmat)%*%new.Xmat
    post.L<-XTX+new.prior.L
    post.var<-solve(post.L)
    post.mean<-post.var %*% (t(new.Xmat) %*% y + new.prior.L %*% new.prior.be0)
    new.c<-(ysq + t(new.prior.be0) %*% new.prior.L %*% new.prior.be0 -
              t(post.mean) %*% post.L %*% post.mean)[1,1]
    
    new.marg.like<-0.5*log(det(new.prior.L))-0.5*log(det(XTX+new.prior.L))-
      0.5*(n+prior.a)*log(0.5*(prior.b+new.c))
    
    new.prior.K<-dpois(new.K,prior.gam,log=T)
    new.prior.eta<-eta.density(new.eta,new.K,xmax)
    
    evec<-c(0,old.eta,xmax)
    new.q<--log(evec[k+2]-evec[k])
    
    if(log(runif(1)) < (new.like+new.prior+new.proposal + new.Jacobian) -
                          (old.like+old.prior+old.proposal)){
      old.K<-new.K
      old.s=new.s
      old.h=new.h
      
      old.like=new.like
      old.prior=new.prior
      old.proposal=new.proposal
      old.Jacobian=new.Jacobian
    }
  }else if(coin<bk+dk+etak){
    #change in height
    h_ind=sample(seq(0:old.K),1)
    new.h=old.h
    new.h[h_ind+1]=old.h[h_ind+1]*exp(runif(-1/2,1/2))
    new.like=loglik(old.s,new.h,L)
    
    if(log(runif(1)) < (new.like-old.like)+alpha*log(new.h[h_ind+1]/new.h[h_ind+1])-beta*(new.h[h_ind+1]-new.h[h_ind+1]))
      old.h=new.h
      
      old.like=new.like
    }
    
  }else{
    #changing location of knot
    s_ind=sample(seq(1:old.K),1)
    new.s=old.s
    new.s[s_ind]=runif(c(0,old.s,L)[s_ind-1],c(0,old.s,L)[s_ind+1])
    new.like=loglike(s=new.s,h=old.h,L=L)
    if(log(runif(1)) < (new.like-old.like)+log(c(0,old.s,L)[s_ind+1]-new.s[s_ind])+log(new.s[s_ind]-c(0,old.s,L)[s_ind-1])
                                          -log(c(0,old.s,L)[s_ind+1]-old.s[s_ind])-log(old.s[s_ind]-c(0,old.s,L)[s_ind-1])){
      old.s=new.s
      old.like=new.like
    }
  }
  
  #Shift one of the knots
  if(old.K >0){
    new.eta<-old.eta
    evec<-c(0,old.eta,xmax)
    k<-sample(1:old.K,size=1)
    new.eta[k]<-runif(1,evec[k],evec[k+2])
    new.Xmat<-rep(1,n)
    for(k in 1:old.K){
      new.Xmat<-cbind(new.Xmat,(x-new.eta[k])^r*(x>new.eta[k]))
    }
    XTX<-t(new.Xmat)%*%new.Xmat
    post.L<-XTX+old.prior.L
    post.var<-solve(post.L)
    post.mean<-post.var %*% (t(new.Xmat) %*% y + old.prior.L %*% old.prior.be0)
    new.c<-(ysq + t(old.prior.be0) %*% old.prior.L %*% old.prior.be0 -
              t(post.mean) %*% post.L %*% post.mean)[1,1]
    new.marg.like<-0.5*log(det(old.prior.L))-0.5*log(det(XTX+old.prior.L))-
      0.5*(n+prior.a)*log(0.5*(prior.b+new.c))
    new.prior.eta<-eta.density(new.eta,old.K,xmax)
    if(log(runif(1)) < (new.marg.like+new.prior.eta) -
       (old.marg.like+old.prior.eta)){
      old.eta<-new.eta
      old.Xmat<-new.Xmat
      old.c<-new.c
      old.marg.like<-new.marg.like
      old.prior.eta<-new.prior.eta
    }
  }
  
  if(iter > nburn & iter %% nthin == 0){
    ico<-ico+1
    #Sample the parameters given old.eta
    XTX<-t(old.Xmat)%*%old.Xmat
    post.L<-XTX+old.prior.L
    post.var<-solve(post.L)
    post.mean<-post.var %*% (t(old.Xmat) %*% y + old.prior.L %*% old.prior.be0)
    old.be<-rmvn(1,mu=post.mean,sigma=old.sig^2*post.var)
    
    old.yfit<-old.Xmat %*% t(old.be)
    old.ssq<-sum((y-old.yfit)^2)
    post.a<-prior.a+n
    post.b<-prior.b+old.ssq
    
    old.sig<-1/sqrt(rgamma(1,post.a/2,post.b/2))
    
    old.c<-(ysq + t(old.prior.be0) %*% old.prior.L %*% old.prior.be0 -
              t(post.mean) %*% post.L %*% post.mean)[1,1]
    
    old.marg.like<-0.5*log(det(old.prior.L))-0.5*log(det(XTX+old.prior.L))-
      0.5*(n+prior.a)*log(0.5*(prior.b+old.c))
    
    be.post[ico,1:(old.K+1)]<-old.be
    sig.post[ico]<-old.sig
    K.post[ico]<-old.K
    eta.post[ico,1:old.K]<-old.eta
  }
  
  if(iter %% 1000 ==0){
    Xm<-rep(1,nx)
    for(k in 1:old.K){
      Xm<-cbind(Xm,(xvec-old.eta[k])^r*(xvec>old.eta[k]))
    }
    XTX<-t(old.Xmat)%*%old.Xmat
    post.L<-XTX+old.prior.L
    post.var<-solve(post.L)
    post.mean<-post.var %*% (t(old.Xmat) %*% y + old.prior.L %*% old.prior.be0)
    old.be<-rmvn(1,mu=post.mean,sigma=old.sig^2*post.var)
    y.fit.b<-Xm %*% t(old.be)
    lines(xvec,y.fit.b,col='red')
  }
}


