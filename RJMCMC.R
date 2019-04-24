#first we assume k fixed at 10
rm(list = ls())
library(boot)#library  that contains the coal miner data
library(lubridate)#library that changes decimal dates


lambda=3
kmax=30
alpha=1
beta=200
cc=1
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
loglike=function(s,h,L){
  sapply(coal1,FUN=evaluatev,s=s,h=h)
  loglikhood=sum(log(h))-sum(diff(s)*h)
  return(loglikhood)
}
bkf=function(c,k,kmax,lambda){
  if(k!=kmax){tr=cc*min(1,pk(k+1,lambda,kmax)/pk(k,lambda,kmax))}
  else{tr=0}
  return(tr)
}
dkf=function(c,k,kmax,lambda){
  if(k!=0){ tr=cc*min(1,pk(k-1,lambda,kmax)/pk(k,lambda,kmax))}
  else{tr=0}
  return(tr)
}

#finding the c value as defined in paper
  cc=seq(0.01,1,0.01)
  min=rep(-1,31)
  for (k in 0:30){
    bk=bkf(cc,k,kmax,lambda)
    dk=dkf(cc,k,kmax,lambda)
    #note that the relationship between c and the total is linear
    tot=bk+dk
    c_index=max(which(tot<=0.9))
    min[k+1]=cc[c_index]
  }
  cc=min(min);cc

####
# Initializing Values
####
  old.K=3
  old.s=c(0,13000,28000,36000,L)#from visual inspection
  h1=sum(coal1<=13000)/13000
  h2=sum(coal1>13000 & coal1<=28000)/(28000-13000)
  h3=sum(coal1>28000 & coal1<=36000)/(36000-28000)
  h4=sum(coal1>36000 )/(L-28000)
  old.h=c(h1,h2,h3,h4)
  old.like=loglike(old.s,old.h,L)
###
#Initializing Markov chain characteristics
###
nburn=0
nsamp=100000
nthin=1
nits=nburn+nthin*nsamp
ico=0


s.post=matrix(0,nrow=nsamp,ncol=kmax+2)
h.post=matrix(0,nrow=nsamp,ncol=kmax+1)
K.post=rep(0,nsamp)

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
    ss=runif(1,min=0,max=L)
    #updating the s
      new.s=rep(0,new.K+2)
      s_ind=min(which(old.s>ss)) 
      new.s[1:(s_ind-1)]=old.s[1:(s_ind-1)]
      new.s[s_ind]=ss
      new.s[(s_ind+1):(new.K+2)]=old.s[s_ind:(old.K+2)]
    #updating the h
    #first must generate the new h
      u=runif(1,min=0,max=1)
      a=ss-old.s[s_ind-1]
      b=old.s[s_ind]-ss
      c=(old.s[s_ind]-old.s[s_ind-1])*log(old.h[s_ind-1])
      d=(1-u)/u
    
    #now updating the h
      new.h=rep(0,old.K+1+1)
      new.h[1:(s_ind-2)]=old.h[1:(s_ind-2)]
      new.h[s_ind]=exp(1/(a+b)*(c+a*log(d)))
      new.h[s_ind-1]=new.h[s_ind]/d
      if(s_ind!=(old.K+2)){new.h[(s_ind+1):(new.K+1)]=old.h[s_ind:(old.K+1)]}

      
    new.like=loglike(new.s,new.h,L)
    
    prior_ratio=log(pk(old.K+1,lambda,kmax)/pk(old.K,lambda,kmax)) +
                log(2*(old.K+1)*(2*old.K+3)/L^2) +
                log(a*b/(old.s[s_ind]-old.s[s_ind-1])) +
                log((beta^alpha)/gamma(alpha)) +
                (alpha-1)*log(new.h[s_ind]*new.h[s_ind-1]/old.h[s_ind-1])-
                beta*(new.h[s_ind-1]+new.h[s_ind]-old.h[s_ind-1])
    
    proposal_ratio=log(dkf(c,old.K+1,kmax,lambda)*L/(bkf(c,old.K,kmax,lambda)*(k+1)))
    Jacobian=log((new.h[s_ind]+new.h[s_ind-1])^2/old.h[s_ind-1])
    
    if(log(runif(1)) < ((new.like-old.like)+prior_ratio+proposal_ratio+Jacobian)){
      
      old.K=new.K
      old.s=new.s
      old.h=new.h
      
      old.like=new.like
    }
  }else if (coin<=bk+dk){
    #Death
    new.K=old.K-1
    
    #first deleting the new s
    if(old.K>1){
      s_ind=sample(seq(2,(old.K+1)),1)
      }else {s_ind=2}
    
    new.s=rep(0,new.K+2)
    new.s[1:(s_ind-1)]=old.s[1:(s_ind-1)]
    new.s[s_ind:(new.K+2)]=old.s[(s_ind+1):(old.K+2)]
    
    #creating one h from two
    new.h=rep(0,new.K+1)
    if(s_ind>2){new.h[1:(s_ind-2)]=old.h[1:(s_ind-2)]} #if not first gap
    c=(new.s[s_ind]-new.s[s_ind-1])
    b=(old.s[s_ind+1]-old.s[s_ind])
    b1=b*log(old.h[s_ind]) # this is a good example of the indexing
    a=(old.s[s_ind]-old.s[s_ind-1])
    a1=a*log(old.h[s_ind-1])
    new.h[s_ind-1]=exp((a1+b1)/c)
    if(s_ind!=(old.K+1)){new.h[s_ind:(new.K+1)]=old.h[(s_ind+1):(old.K+1)]} #if not the last gap
      
    new.like=loglike(new.s,new.h,L)
    
    prior_ratio=-log(pk(k+1,lambda,kmax)/pk(k,lambda,kmax))-
    log(2*(k+1)*(2*k+3)/L^2)-
    log((a*b)/c)-
    log(beta^alpha/gamma(alpha))-
    (alpha-1)*log((old.h[s_ind]*old.h[s_ind-1])/new.h[s_ind-1])+
    beta*(old.h[s_ind]+old.h[s_ind-1]-new.h[s_ind-1])
    
    proposal_ratio=-log(dkf(c,old.K+1,kmax,lambda)*L/(bkf(c,old.K,kmax,lambda)*(k+1)))
    
    Jacobian=-log((old.h[s_ind]+old.h[s_ind-1])^2/new.h[s_ind-1])
    
    if(log(runif(1)) < ((new.like-old.like)+prior_ratio+proposal_ratio+Jacobian)){#note added negatives to compensate
      old.K=new.K
      old.s=new.s
      old.h=new.h
      
      old.like=new.like
    }
  }else if(coin<bk+dk+etak){
    #change in height
    if(old.K>0){h_ind=sample(seq(1,old.K+1),1)}
    else{h_ind=1}
    new.h=old.h
    new.h[h_ind]=old.h[h_ind]*exp(runif(1,min=-1/2,max=1/2))
    new.like=loglike(old.s,new.h,L)
    
    if(log(runif(1)) < (new.like-old.like)+alpha*log(new.h[h_ind]/new.h[h_ind])-beta*(new.h[h_ind]-new.h[h_ind])){
      old.h=new.h
      old.like=new.like
    }
    
  }else{
    #changing location of knot
    if (old.K>1){
      s_ind=sample(seq(2,old.K+1),size=1)
    }else {s_ind=2}
    new.s=old.s
    new.s[s_ind]=runif(1,min=old.s[s_ind-1],max=old.s[s_ind+1])
    new.like=loglike(s=new.s,h=old.h,L=L)
    if(log(runif(1)) < (new.like-old.like)+log(old.s[s_ind+1]-new.s[s_ind])+log(new.s[s_ind]-old.s[s_ind-1])
                                          -log(old.s[s_ind+1]-old.s[s_ind])-log(old.s[s_ind]-old.s[s_ind-1])){
      old.s=new.s
      old.like=new.like
    }
  }
  
  if(iter > nburn & iter %% nthin == 0){
    ico<-ico+1
    h.post[ico,1:(old.K+1)]=old.h
    s.post[ico,1:(old.K+2)]=old.s
    K.post[ico]=old.K
  }
}


