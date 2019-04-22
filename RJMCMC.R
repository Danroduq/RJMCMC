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
#processing the coal mining dates to be number of days since origin 1981/01/01
tempdate=lapply(coal,FUN=date_decimal)
tempdate=sapply(tempdate,FUN=as_date)
coal1=difftime(as.Date(tempdate,format="%Y-%m-%d",origin="1970-01-01"),"1851-01-01", units = c("days"))
L=as.numeric(difftime("1962/12/31" ,"1851/01/01", units = c("days")))

pk=function(k,lambda,kmax){
  pk=dpois(k,lambda)/ppois(kmax,lambda)
  return(pk)
}
evaluate=function(x,s,h){
  ind=min(which(x<=c(0,s,L)))
  return(h[ind-1])
}
#vectorizing the evaluate function
evaluatev=Vectorize(evaluate,vectorize.args=c("x"))
loglik=function(s,h,L){
  sapply(coal1,FUN=evaluatev,s=s,h=h)
  loglikhood=sum(log(h))-sum(diff(c(0,s,L))*h)
  return(loglikhood)
}
bkf=function(c,k,kmax,lambda){
  if(k!=kmax){tr=c*min(1,pk(k+1,lambda,kmax)/pk(k,lambda,kmax))}
  else{tr=0}
  return(tr)
}
dkf=function(c,k,kmax,lambda){
  if(k!=0){ tr=c*min(1,pk(k-1,lambda,kmax)/pk(k,lambda,kmax))}
  else{tr=0}
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

####
# Initializing Values
####
  old.K=3
  old.s=c(13000,28000,36000)#from visual inspection
  h1=sum(coal1<=13000)/13000
  h2=sum(coal1>13000 & coal1<=28000)/(28000-13000)
  h3=sum(coal1>28000 & coal1<=36000)/(36000-28000)
  h4=sum(coal1>36000 )/(L-28000)
  old.h=c(h1,h2,h3,h4)
  
nburn<-1000
nsamp<-1000
nthin<-20
nits<-nburn+nthin*nsamp

sig.post<-rep(0,nsamp)
be.post<-eta.post<-matrix(0,nrow=nsamp,ncol=100)
K.post<-rep(0,nsamp)



for(iter in 1:nits){
 
  bk=bkf(c,k=old.K,kmax,lambda)
  dk=dkf(c,k=old.K,kmax,lambda)
  split=1-bk-dk
  if(old.K!=0){
  pik=split/2
  etak=split/2
  }else{
    pik=0
    etak=split
  }
  
  coin=runif(1)
  if(coin<bk){
    #birth 
    new.K=old.K+1
    ss=runif(0,L)
    s_ind=min(which(ss<old.s))
    new.s=rep(0,old.K+1)
    new.s[1:s_ind]=old.s[1:s_ind]
    new.s[s_ind+1]=ss
    new.s[s_ind+2:old.K+1]=old.s[s_ind+1:old.K]
    u=runif(0,1)
   
    a=ss-old.s[s_ind]
    b=old.s[s_ind+1]-ss
    c=(old.s[s_ind+1]-old.s[s_ind])*log(old.h[s_ind])
    d=(1-u)/u
    
    new.h=rep(0,old.K+2)
    new.h[1:s_ind-1]=old.h[1:s_ind]
    new.h[s_ind+1]=exp(1/(a+b)(c+a*log(d)))
    new.h[s_ind]=new.h[s_ind+2]/d
    new.h[s_ind+2:old.K+2]=old.h[s_ind+2:old.K+2]
    
    new.like=loglik(new.s,new.h,L)
    prior_ratio=log(pk(k+1,lambda,kmax)/pk(k,lambda,kmax))
                +log(2*(k+1)*(2*k+3)/L^2)
                +log(a*b/old.s[s_ind+1]-old.s[s_ind])
                +log(beta^alpha/gamma(alpha))
                +(alpha-1)*log(new.h[s_ind+1]*new.h[s_ind]/old.h[s_ind])
                -beta*(new.h[s_ind+1]+new.h[s_ind]-old.h[s_ind])
    proposal_ratio=log(dkf(c,old.k+1,kmax,lambda)*L/(bkf(c,oldk,kmax,lambda)*(k+1)))
    Jacobian=log((new.h[s_ind+1]+new.h[s_ind])^2/old.h[s_ind])
    
    if(log(runif(1)) < (new.like-old.like)+prior_ratio+proposal_ratio+Jacobian){
      
      old.K<-new.K
      old.s=new.s
      old.h=new.h
      
      old.like=new.like
    }
  }else if (coin<=bk+dk){
    #Death
    new.K=old.K-1
    s_ind=sample(seq(1:old.K),1)
    new.s=rep(0,new.K)
    new.s[1:s_ind-1]=old.s[1:s_ind-1]
    new.s[s_ind:new.K]=old.s[s_ind+1:old.K]
    
    new.h=rep(0,new.K+1)
    new.h[1:s_ind-1]=old.h[1:s_ind-1]
    a=(new.s[s_ind]-new.s[s_ind])
    b=(new.s[s_ind+1]-new.s[s_ind]) *log(old.h[s_ind]) # this is a good example of the indexing
    c=(new.s[s_ind]-new.s[s_ind-1]) *log(old.h[s_ind-1])
    new.h[ind]=exp(b+c/a)
    new.h[s_ind:new.K+1]=old.h[s_ind+1:old.K+1]
      
    new.like=loglik(new.s,new.h,L)
    
    prior_ratio=-log(pk(k+1,lambda,kmax)/pk(k,lambda,kmax))
    -log(2*(k+1)*(2*k+3)/L^2)
    -log((old.s[s_ind]-old.s[s_ind-1])*(old.s[s_ind+1]-old.s[s_ind])/(new.s[s_ind]-new.s[s_ind-1]))
    -log(beta^alpha/gamma(alpha))
    -(alpha-1)*log(old.h[s_ind+1]*old.h[s_ind]/new.h[s_ind])
    +beta*(old.h[s_ind+1]+old.h[s_ind]-new.h[s_ind])
    
    proposal_ratio=-log(dkf(c,old.k+1,kmax,lambda)*L/(bkf(c,old.k,kmax,lambda)*(k+1)))
    Jacobian=-log((old.h[s_ind+1]+old.h[s_ind])^2/new.h[s_ind])
    
    if(log(runif(1)) < ((new.like-old.like)+prior_ratio+proposal_ratio+Jacobian)){
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


