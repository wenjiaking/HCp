library(ggplot2)
library(cowplot)
library(rlist)
library(Rmpfr)
library(Rfast)
library(GHC)
#install.packages("SetTest")
library(SetTest)

#######Prepare Data for Appendix Figure I #######
phc.mod <- function(q, M, k0, k1, LS = F, ZW = F, onesided=FALSE,MHC=FALSE){
  if(LS == F && ZW == F){
    return(pphi.mod(q=q, M=M, k0=k0, k1=k1, s=2, onesided=onesided,MHC = MHC)) #ZW analytic method only depends on pphi function, s=2 means HC2004 test statistic
  }else{
    n = length(M[1,])
    beta = k1/n
    if(LS == T){
      return(1-hcls(q, n, beta))
    }else{
      return(1-hczw(q, n, beta))
    }
  }
}


pphi.mod <- function(q, M, k0, k1, s=2, t=30, onesided=FALSE,MHC=FALSE)
{
  qtemp = q
  q = max(0.01, q) #why this step? generally observed statistic should be larger than 0.01
  NODES=seq(-4,4,length.out=40) #for independent studies, NODES is not required
  n = length(M[1,]) # the number os studies
  nrep = n
  nthr = 5
  lower = q*0.9#max(0.5,q*0.9)
  upper = q*1.1
  thr = seq(lower,upper,length.out=nthr) #actually calculate nthr points around the observed statistic
  rho = sample(M[row(M)!=col(M)], nrep) #sample from off-diagnal elements of correlation matrix
  rho = abs(rho)
  digi = rep(0.02, nrep)
  rho = round(rho/digi)*digi
  rho[rho==1] = 0.99
  tbl = table(rho) #for independent studies, all rho should be 0
  rho_unique = as.numeric(names(tbl)) # unique rho
  rho_freq = as.numeric(tbl) # frequency of unique rho
  
  nrep_uniq = length(rho_unique) #for independent studies, the length should be 1
  pvalue_cal = c()
  for(i in 1:nrep_uniq){
    freq = rho_freq[i] # for independent studies, the freq should be n
    #  y = sapply(1:nthr, function(x)pphi.rho(threshold=thr[x], n=n, rho=rho_unique[i], k0=k0, k1=k1, NODES=NODES, s=s, t=t, onesided=onesided))
    y = log(sapply(1:nthr, function(x)pphi.rho.mod(threshold=thr[x], n=n, rho=rho_unique[i], k0=k0, k1=k1, NODES=NODES, s=s, t=t, onesided=onesided,MHC = MHC)))
    #print(y)
    pvalue_cal = append(pvalue_cal,rep(y, freq)) #only depends on pphi.rho function
  }
  # data = data.frame(y=pvalue_cal, x=rep(thr, nrep))
  data = data.frame(y=pvalue_cal, x=rep(log(thr), nrep))
  #print(data[1:5,])
  fit = loess(y~x, data = data)
  # res = predict(fit, q)
  res = exp(predict(fit, log(q)))
  if(qtemp>=0.01){
    return(res)
  }else{
    warning(paste0("Left-tail prob. < ",round(res,3), " or p-value > ", 1-round(res,3)))
    return(res)
  }
}


phi.f.inv.i.mod <- function(q, i, n, s){
  #ep = 10^(-15)
  ep = 10^(-18)
  if(i!=n){
    if(q>0){
      if (s!=1&&s!=0){
        f <- function(x)(1-(i/n)^s*x^(1-s)-(1-i/n)^(s)*(1-x)^(1-s))*2*n/(s-s^2) - q^2 #when s=2 the phi divergence is HC2004, f is the HC(x)^2-q^2 where x is p_(i)
        if (f(ep)<0){
          CC=0 #if the root<10^(-15) then set the root as 0, which may push mu[1],mu[2] to be 0 when q is very large (e.g. HC, b=10^7 K=50)
        }else{
          CC=uniroot(f,c(ep,i/n),tol=1e-100)$root #find the root at range (10^(-15),0). Maybe consider to set ep much smaller
        }
      }else if (s==1){
        f <- function(x)2*n*((i/n)*log(i/n/x)+(1-i/n)*log((1-i/n)/(1-x)))-q^2
        if (f(ep)<0){
          CC=0
        }else{
          CC=uniroot(f,c(ep,i/n),tol=1e-100)$root
        }
      }else{
        f <- function(x)2*n*((x)*log(x/(i/n))+(1-x)*log((1-x)/(1-i/n)))-q^2
        if (f(ep)<0){
          CC=0
        }
        else{
          CC=uniroot(f,c(ep,i/n),tol=1e-100)$root
        }
      }
    }else if(q<0){
      if (s!=1&&s!=0){
        f <- function(x)(1-(i/n)^s*x^(1-s)-(1-i/n)^(s)*(1-x)^(1-s))*2*n/(s-s^2) - q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
        }
      }else if (s==1){
        f <- function(x)2*n*((i/n)*log(i/n/x)+(1-i/n)*log((1-i/n)/(1-x)))-q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
        }
      }else{
        f <- function(x)2*n*((x)*log(x/(i/n))+(1-x)*log((1-x)/(1-i/n)))-q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
        }
      }
    }else{
      CC = i/n
    }
  }else{
    if(q<=0||s<=0){
      CC = 1
    }else{
      if(s==1){
        f <- function(x)-2*log(x)*n - q^2
        CC=uniroot(f,c(ep,1),tol=1e-100)$root
      }else{
        f <- function(x)(1-x^(1-s))*2*n/(s-s^2) - q^2
        CC=uniroot(f,c(ep,1),tol=1e-100)$root
      }
    }
  }
  return(CC)
}

pphi.rho.mod <- function(threshold, n, rho, k0, k1, NODES=seq(-4,4,length.out=40),
                         s=2, t=30, onesided=T,MHC=TRUE){
  #t = min(t,k1-k0+1) #why minimum of t=30 and k1-k0+1. It should be k1-k0+1
  t = k1-k0+1
  #if t<k1-k0+1, then this is not the exact result, approximation exist (what if without this approximation?)
  u = sapply(1:k1,function(x)phi.f.inv.i.mod(threshold, x, n, s))
  #print(range(u))
  um = u[length(u)]
  #print(u)
  #u = phiInverse(threshold, s, k1=beta*n, n, t=beta*n)$bound
  if(rho!=0){
    SUM = sum(dnorm(NODES))
    result = 0
    if(onesided==T){
      for(z in NODES){
        d = pnorm((qnorm(1-u)-sqrt(rho)*z)/sqrt(1-rho),lower.tail = F)
        dm = pnorm((qnorm(1-um)-sqrt(rho)*z)/sqrt(1-rho),lower.tail = F)
        result = result + dnorm(z)*UnifCross_v1.mod(d,dm,t=t,n=n,k0=k0,k1=k1)/SUM
      }
    }else{
      for(z in NODES){
        d = pnorm((qnorm(1-u/2)-sqrt(rho)*z)/sqrt(1-rho),lower.tail = F) +
          pnorm((-qnorm(1-u/2)-sqrt(rho)*z)/sqrt(1-rho))
        dm = pnorm((qnorm(1-um/2)-sqrt(rho)*z)/sqrt(1-rho),lower.tail = F) +
          pnorm((-qnorm(1-um/2)-sqrt(rho)*z)/sqrt(1-rho))
        result = result + dnorm(z)*UnifCross_v1.mod(d,dm,t=t,n=n,k0=k0,k1=k1)/SUM
      }
    }
    return(result)
  }else{
    if (MHC) return(UnifCross_v1_MHC(u,um,t=t,n=n,k0=k0,k1=k1))
    else return(UnifCross_v1.mod(u,um,t=t,n=n,k0=k0,k1=k1)) #for independent studies, only depends on UnifCross_v1 function
  }
}


#Crossing probablity of uniform order statistics
UnifCross_v1.mod <- function(u, um, t, n, k0, k1){
  m = floor(k1)
  pp = rep(NA, t)
  pp[1] = pbeta(um,m,n-m+1,lower.tail=F)
  pp[2:t] = (lfactorial(n) - lfactorial(n-(k0:(k0+t-2))))+log(pbeta(um,m-(k0:(k0+t-2)),n-m+1,lower.tail=F))
  a = rep(1,t)
  a[2] = -dpois(k0,u[k0])*exp(u[k0])
  if(t>1){
    for (i in 2:(t-1)){
      d = dpois(c(1:(i-1), i+k0-1),u[i+k0-1])*exp(u[i+k0-1])
      
      a[i+1]=-a[i:1]%*%d
    }
    a_pos=which(a[2:t]>0)+1
    a_neg=which(a[2:t]<0)+1
    p_val=a[1]*pp[1]+sum(exp(pp[a_pos]+log(a[a_pos])))-sum(exp(pp[a_neg]+log(-a[a_neg])))
    return(p_val)
  }else{
    p = pbeta(u[k0], k0, n-k0+1)
    return(1-drop(p))
  }
}


UnifCross_v1_MHC <- function(u, um, n, t,k0, k1){
  alpha=1/n 
  m = floor(k1)
  u=sapply(u,function(x) max(x,alpha))
  um=max(um,alpha)
  pp = mpfr(rep(NA, t),80)
  pp[1]=mpfr(sum(pbinom(0:(k0-1),size = n,prob=alpha)*pbeta((um-alpha)/(1-alpha),m:(m-k0+1),n-m+1,lower.tail=F)),80)
  pp[2:t] = mpfr(exp(lfactorial(n) - lfactorial(n-(k0:(k0+t-2))))*pbeta(um,m-(k0:(k0+t-2)),n-m+1,lower.tail=F),80)
  a = mpfr(rep(1,t),80)
  a[2] =mpfr((alpha^k0-pbinom(k0-1,k0,alpha/u[k0])*u[k0]^k0)/factorial(k0),80)
  if(t>1){
    for (i in 2:(t-1)){
      d = mpfr(c((-dpois(1:(i-1),u[i+k0-1])*exp(u[i+k0-1])),
                 (alpha^(i+k0-1)-pbinom(k0-1,i+k0-1,alpha/u[i+k0-1])*u[i+k0-1]^(i+k0-1))/factorial(i+k0-1)),80)
      
      a[i+1]=a[i:1]%*%d
    }
    
    return(asNumeric(drop(pp[1:t]%*%a[1:t])))
    
  }else{
    p = pbeta(u[k0], k0, n-k0+1)
    return(1-drop(p))
  }
}

HC=function(p.values) {
  K=length(p.values)
  sort.p = sort(p.values)
  GHC.mesh=((1:K)-K*sort.p)/sqrt(K*sort.p*(1-sort.p))
  GHCstat = max(GHC.mesh)
  return(stat=GHCstat)
}

ImportSamp=function(k,lambda,family="Gaussian",mu=1,s=5) {
  if (family=="Gaussian") {
    x=sapply(lambda, function(x) rnorm(1,mean = 0,sd=sqrt(x)))
    p=2*pnorm(abs(x),mean=0,sd=sqrt(mu),lower.tail = FALSE)
    w=prod(sqrt(lambda))/mu^(k/2)*exp(sum((1/lambda-1/mu)*x^2/2)) 
    w.ext=w*x^2
    stat.val=HC(p.values = p)
  }
  if (family=="Beta") {
    x=sapply(lambda, function(x) rbeta(1,shape1 = x,shape2 = 1))
    w=1/prod(lambda*x^(lambda-1))
    w.ext=w*log(x)
    stat.val=HC(p.values = x)
  }
  
  if (family=="Exponential") {
    x=sapply(lambda, function(x) rexp(1,rate =1/x))
    w=prod(lambda)/mu^k*exp(sum((1/lambda-1/mu)*x))
    w.ext=w*x
    p=pexp(x,rate=1/mu,lower.tail = FALSE)
    stat.val=HC(p.values = p)
  }
  
  if (family=="Weibull") {
    # the smaller the shape s, the heavier tail of weibull distribution
    x=sapply(lambda, function(x) rweibull(1,shape = s,scale = x^(1/s)))
    w=prod(lambda)/mu^(k)*exp(sum((1/lambda-1/mu)*x^s)) 
    w.ext=w*x^s
    p=pweibull(x,shape=s,scale=mu^(1/s),lower.tail = FALSE)
    stat.val=HC(p.values = p)
  }
  
  return(list(stat.val=stat.val,weight=w,weight.ext=w.ext))
}


par.update=function(object,ro, family="Gaussian") {
  stat.val=sapply(object,"[[",1)
  gamma=quantile(stat.val,1-ro)
  weight=sapply(object,"[[",2)
  weight.ext=sapply(object,"[[",3)
  if (family=="Gaussian") {
    lambda.den=apply(weight.ext,1,function(x) sum(x[stat.val>=gamma]))
    lambda.num=sum(weight[stat.val>=gamma])
    lambda=-lambda.num/lambda.den
  }
  else {
    lambda.num=apply(weight.ext,1,function(x) sum(x[stat.val>=gamma]))
    lambda_den=sum(weight[stat.val>=gamma])
    lambda=lambda.num/lambda_den
  }
  
  return(list(gamma=gamma,lambda=lambda))
}

CE=function(k,ro,N=10^5,q,family="Gaussian",mu=1,s=5) {
  lambda=rep(mu,k)
  gamma=-Inf
  t=0
  while (gamma<q) {
    t=t+1
    obj=mclapply(1:N,function(i) ImportSamp(k=k,lambda=lambda,family=family,mu=mu, s=s),mc.cores = 30)
    par=par.update(object = obj,ro=ro,family=family)
    lambda=par$lambda
    gamma=par$gamma
    message(paste0(gamma,", ",lambda))
  }
  message("stop")
  obj.final=mclapply(1:N,function(i) ImportSamp(k=k,lambda=lambda,family=family,mu=mu, s=s),mc.cores = 30)
  stat.val.final=sapply(obj.final, "[[",1)
  weight.final=sapply(obj.final, "[[",2)
  p.val=1/N*sum(weight.final[stat.val.final>=q])
  return(list(p.val=p.val,iter=t,lambda=lambda))
}

libs_HC.CE=function(q.val.set,ro=0.01,N=10^4,K,family="Gaussian",mu=1,s=5) {
  temp.lib=lapply(q.val.set, function(x) CE(k=K,ro=ro,N=N,q=x,family=family,mu=mu,s=s))
  temp.p=sapply(temp.lib,"[[",1)
  return(temp.p)
}

#CE.beta(k=50,ro=0.01,N=10^4,method="HC",q=q.vals50[50])
set.seed(1025)
q.vals50=exp(seq(log(10^3),log(10^7),length.out = 200))
libs.FHC50.Exp=lapply(1:10,function(x) libs_HC.CE(q.val.set=q.vals50[c(50,80,100,150)],ro=0.01,N=10^4,K=50,family="Exponential"))
list.save(libs.FHC50.Exp,"libs.FHC50.Exp.RData")
libs.FHC50.Beta=lapply(1:10,function(x) libs_HC.CE(q.val.set=q.vals50[c(50,80,100,150)],ro=0.01,N=10^4,K=50,family="Beta"))
list.save(libs.FHC50.Beta,"libs.FHC50.Beta.RData")
libs.FHC50.Gauss=lapply(1:10,function(x) libs_HC.CE(q.val.set=q.vals50[c(50,80,100,150)],ro=0.01,N=10^4,K=50,family="Gaussian"))
list.save(libs.FHC50.Gauss,"libs.FHC50.Gauss.RData")
libs.FHC50.Weibull_s5=lapply(1:10,function(x) libs_HC.CE(q.val.set=q.vals50[c(50,80,100,150)],ro=0.01,N=10^4,K=50,family="Weibull",s=5))
list.save(libs.FHC50.Weibull_s5,"libs.FHC50.Weibull_s5.RData")
libs.FHC50.Weibull_s1=lapply(1:10,function(x) libs_HC.CE(q.val.set=q.vals50[c(50,80,100,150)],ro=0.01,N=10^4,K=50,family="Weibull",s=1))
list.save(libs.FHC50.Weibull_s1,"libs.FHC50.Weibull_s1.RData")
libs.FHC50.mixGrad=lapply(1:10, function(i) libs_HC.CE.mixGrad(q.val.set=q.vals50[c(50,80,100,150)],K=50,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=50, idx=1,thre=FALSE))
list.save(libs.FHC50.mixGrad,"libs.FHC50.mixGrad.RData")
saveRDS(sapply(q.vals50[c(50,80,100,150)],function(x) return(1-phc.mod(x, M=diag(rep(1,50)), k0=1, k1=50))),"FHC50.analyticP.rds")

HC.stat=function(x, k0=1, k1,thre=FALSE) {
  t_mesh = Rfast::rowSort(abs(x),descending = TRUE)
  K=dim(t_mesh)[2]
  n=dim(t_mesh)[1]
  tailprob = 2*pnorm(t_mesh,lower.tail=F)
  tailprob=tailprob[,k0:k1]
  if (thre) {
    tailprob[which(tailprob<1/K)]=k1/K
  }
  GHC.mesh=(matrix(rep(k0:k1,times=n),nrow = n,byrow = T)-K*tailprob)/sqrt(K*tailprob*(1-tailprob))
  GHCstat = rowMaxs(GHC.mesh,value = T)
  return(stat=GHCstat)
}

CE.HC.mixed.grad=function(K,theta,k0=1,k1,N1=1000,idx=1,thre=FALSE) {
  grad=1/((1:K)^idx+1)
  mix.prop=matrix(rbinom(K*N1,size=1,prob=grad),nrow=K)
  mix.prop=transpose(mix.prop)
  
  x=matrix(Rfast::Rnorm(K*N1,0,theta),nrow = N1)^{(mix.prop)}+matrix(Rfast::Rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  stat.val=HC.stat(x=x,k0=k0, k1=k1,thre=thre)
  log_g_norm <- rowsums(log(dnorm(x))) # numerator part
  
  log_g_mix <- rowsums(log(transpose(transpose(dnorm(x,0,theta))*grad)+transpose(transpose(dnorm(x))*(1-grad))))
  
  w=exp(log_g_norm-log_g_mix)
  return(list(stat.val=stat.val,weight=w,x=x))
}

mixed.par.update.grad=function(object,ro,K,idx=1) {
  
  stat.val=object[[1]]
  gamma=quantile(stat.val,1-ro)
  weight=object[[2]]
  x.mat=transpose(object[[3]])
  B=length(stat.val)
  grad=1/((1:K)^idx+1)
  
  opt=optimize(function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,function(x) sum(log(grad*dnorm(x,0,t)+(1-grad)*dnorm(x,0,1))))),
               interval = c(0,5),maximum = TRUE)
  #opt=optim(2,function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,function(x) sum(log(grad*dnorm(x,0,t)+(1-grad)*dnorm(x,0,1)))),method="L-BFGS-B"), 
  # control=list(fnscale=-1))
  par=opt$maximum
  #par=opt$par
  return(list(gamma=gamma,par=par))
}


CE.mixed.grad=function(K,ro,N=10^4,q,idx=1,k0=1,k1,thre=FALSE) {
  theta=5
  gamma=-Inf
  t=0
  while (gamma<q) {
    t=t+1
    obj=CE.HC.mixed.grad(K=K,theta=theta,k0=k0,k1=k1,N1=N,idx=idx,thre=thre)
    par=mixed.par.update.grad(object = obj,ro=ro,K=K,idx = idx)
    gamma=par$gamma
    theta=par$par
    message(paste0(gamma,", ",theta))
  }
  message("stop")
  obj.final=CE.HC.mixed.grad(K=K,theta=theta,k0=k0,k1=k1,N1=N,idx=idx,thre=thre)
  stat.val.final=obj.final[[1]]
  weight.final=obj.final[[2]]
  
  p.val=1/N*sum(weight.final[stat.val.final>=q])
  return(list(p.val=p.val,iter=t,theta=theta))
}

parsHC.CE.mixGrad=function(q,K=50,ro=0.01,N=10^4,idx,k0=1,k1,thre=FALSE){
  #temp.lib=mclapply(q.val.set, function(x) CE.mixed.grad(K=K,ro=ro,N=N,q=x, k0=k0, k1=k1,thre = thre),mc.cores = 30)
  temp.lib=lapply(idx, function(x) CE.mixed.grad(K=K,ro=ro,N=N,q=q, idx=x,k0=k0,k1=k1,thre=thre))
  temp.par=sapply(temp.lib,"[[",3)
  return(temp.par)
}

CE.HC.mixgrad.fix=function(q,K,theta,N=1000,idx,k0=1,k1,thre=FALSE) {
  obj.final=CE.HC.mixed.grad(K=K,theta=theta,k0=k0,k1=k1,N1=N,idx = idx,thre = thre)
  stat.val.final=obj.final[[1]]
  weight.final=obj.final[[2]]
  
  p.val=1/N*sum(weight.final[stat.val.final>=q])
  print(p.val)
  return(p.val)
}
libs50_HC.CE.fix=function(q,K,theta,N=10^4,B=50,idx,k0=0,k1,thre=FALSE) {
  temp.p=lapply(1:B, function(b) CE.HC.mixgrad.fix(q=q,K=K,theta=theta,N=N,idx=idx,k0=k0,k1=k1,thre=thre))
  #temp.p=sapply(temp.lib,"[[",1)
  return(temp.p)
}

q.vals50=exp(seq(log(10^3),log(10^7),length.out = 200))
qvals=exp(seq(log(30),17,length=100))
idx_thetas=lapply(c(10^4,10^5,10^6),function(N) parsHC.CE.mixGrad(q=q.vals50[100],K=100,ro=0.01,N=N,idx=c(0.3,0.5,1,2),k0=1,k1=100,thre=FALSE))
idx_thetas_THC=lapply(c(10^4,10^5,10^6),function(N) parsHC.CE.mixGrad(q=qvals[60],K=100,ro=0.01,N=N,idx=c(0.3,0.5,1,2),k0=1,k1=50,thre=FALSE))

idx_vec=c(0.3,0.5,1,2)
idxHC100N4_fix=lapply(1:4,function(i) libs50_HC.CE.fix(q=q.vals50[100],K=100,theta=idx_thetas[[1]][i],N=10^4,B=50,idx=idx_vec[i],k0=1,k1=100,thre=FALSE))
idxHC100N5_fix=lapply(1:4,function(i) libs50_HC.CE.fix(q=q.vals50[100],K=100,theta=idx_thetas[[2]][i],N=10^5,B=50,idx=idx_vec[i],k0=1,k1=100,thre=FALSE))
idxHC100N6_fix=lapply(1:4,function(i) libs50_HC.CE.fix(q=q.vals50[100],K=100,theta=idx_thetas[[3]][i],N=10^6,B=50,idx=idx_vec[i],k0=1,k1=100,thre=FALSE))
list.save(idxHC100N4_fix,"idxHC100N4_fix.RData")
list.save(idxHC100N5_fix,"idxHC100N5_fix.RData")
list.save(idxHC100N6_fix,"idxHC100N6_fix.RData")
saveRDS((1-phc.mod(q.vals50[100], M=diag(rep(1,100)), k0=1, k1=100)),"idxHC100analy.rds")

idxTHC100N4_fix=lapply(1:4,function(i) libs50_HC.CE.fix(q=qvals[60],K=100,theta=idx_thetas_THC[[1]][i],N=10^4,B=50,idx=idx_vec[i],k0=1,k1=50,thre=FALSE))
idxTHC100N5_fix=lapply(1:4,function(i) libs50_HC.CE.fix(q=qvals[60],K=100,theta=idx_thetas_THC[[2]][i],N=10^5,B=50,idx=idx_vec[i],k0=1,k1=50,thre=FALSE))
idxTHC100N6_fix=lapply(1:4,function(i) libs50_HC.CE.fix(q=qvals[60],K=100,theta=idx_thetas_THC[[3]][i],N=10^6,B=50,idx=idx_vec[i],k0=1,k1=50,thre=FALSE))
list.save(idxTHC100N4_fix,"idxTHC100N4_fix.RData")
list.save(idxTHC100N5_fix,"idxTHC100N5_fix.RData")
list.save(idxTHC100N6_fix,"idxTHC100N6_fix.RData")
saveRDS((1-phc.mod(q.vals50[100], M=diag(rep(1,100)), k0=1, k1=50)),"idxTHC100analy.rds")

CE.HC.mixed.Prop=function(K,theta,N1=1000,k0,k1,thre,prop=0.2) {
  if (prop==-Inf) {grad=rep(1,K)}
  else {grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}
  mix.prop=matrix(rbinom(K*N1,size=1,prob=grad),nrow=K)
  mix.prop=transpose(mix.prop)
  
  x=matrix(Rfast::Rnorm(K*N1,0,theta),nrow = N1)^{(mix.prop)}+matrix(Rfast::Rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  stat.val=HC.stat(x=x,k0=k0,k1=k1,thre=thre)
  log_g_norm <- rowsums(log(dnorm(x))) # numerator part
  
  log_g_mix <- rowsums(log(transpose(transpose(dnorm(x,0,theta))*grad)+transpose(transpose(dnorm(x))*(1-grad))))
  
  w=exp(log_g_norm-log_g_mix)
  return(list(stat.val=stat.val,weight=w,x=x))
}

mixed.par.update.Prop=function(object,ro,K,prop) {
  
  stat.val=object[[1]]
  gamma=quantile(stat.val,1-ro)
  weight=object[[2]]
  x.mat=transpose(object[[3]])
  B=length(stat.val)
  if (prop==-Inf) {grad=rep(1,K)}
  else {grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}
  opt=optimize(function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,function(x) sum(log(grad*dnorm(x,0,t)+(1-grad)*dnorm(x,0,1))))),
               interval = c(0,5),maximum = TRUE)
  par=opt$maximum
  return(list(gamma=gamma,par=par))
}


CE.mixed.Prop=function(K,ro,N=10^4,q,k0,k1,thre,prop=0.2,theta=1) {
  #theta=1
  gamma=-Inf
  t=0
  while (gamma<q) {
    t=t+1
    obj=CE.HC.mixed.Prop(K=K,theta=theta,N1=N, k0=k0,k1=k1,thre=thre,prop = prop)
    par=mixed.par.update.Prop(object = obj,ro=ro,K=K,prop = prop)
    gamma=par$gamma
    theta=par$par
    print(paste0(gamma,", ",theta))
  }
  print("stop")
  obj.final=CE.HC.mixed.Prop(K=K,theta=theta,N1=N,k0=k0,k1=k1,thre=thre,prop = prop)
  stat.val.final=obj.final[[1]]
  weight.final=obj.final[[2]]
  
  p.val=1/N*sum(weight.final[stat.val.final>=q])
  return(list(p.val=p.val,iter=t,theta=theta))
}

CE.MHC.mixed.fix=function(q,K,theta,N1=1000,k0,k1,thre,prop=0.2) {
  if (prop==-Inf) {grad=rep(1,K)}
  else {grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}
  mix.prop=matrix(rbinom(K*N1,size=1,prob=grad),nrow=K)
  mix.prop=transpose(mix.prop)
  
  x=matrix(Rfast::Rnorm(K*N1,0,theta),nrow = N1)^{(mix.prop)}+matrix(Rfast::Rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  stat.val=HC.stat(x=x,k0=k0,k1=k1,thre=thre)
  log_g_norm <- rowsums(log(dnorm(x))) # numerator part
  
  log_g_mix <- rowsums(log(transpose(transpose(dnorm(x,0,theta))*grad)+transpose(transpose(dnorm(x))*(1-grad))))
  
  w=exp(log_g_norm-log_g_mix)
  p.val=1/N1*sum(w[stat.val>=q])
  print(p.val)
  return(p.val)
}
libs50_MHC.CE.fix=function(q,K,theta,N=10^4,k0=1,k1, thre=FALSE,B=50,prop) {
  temp.p=lapply(1:B, function(b) CE.MHC.mixed.fix(q=q,K=K,theta=theta,N1=N,k0=k0,k1=k1,thre=thre,prop=prop))
  #temp.p=sapply(temp.lib,"[[",1)
  return(temp.p)
}


pars_HC.CE.mixProp=function(q,K=50,ro=0.01,N=10^4,k0=1,k1=25, thre=FALSE,prop=0.2,theta=1){
  #temp.lib=mclapply(q.val.set, function(x) CE.mixed.grad(K=K,ro=ro,N=N,q=x, k0=k0, k1=k1,thre = thre),mc.cores = 30)
  temp.lib=lapply(prop, function(x) CE.mixed.Prop(K=K,ro=ro,N=N,q=q, k0=k0, k1=k1,thre = thre,prop=x,theta = theta))
  temp.par=sapply(temp.lib,"[[",3)
  return(temp.par)
}

pars_MHC=lapply(c(10^4,10^5,10^6),function(N) pars_HC.CE.mixProp(q=exp(seq(log(3),log(13),length=100))[subInd][14],K=100,ro=0.01,N=N,k0=1,k1=100, thre=TRUE,prop=c(-Inf,0.01,0.1,0.2,0.5)))
pars_TMHC=lapply(c(10^4,10^5,10^6),function(N) pars_HC.CE.mixProp(q=exp(seq(log(3),log(13),length=100))[subInd][14],K=100,ro=0.01,N=N,k0=1,k1=50, thre=TRUE,prop=c(-Inf,0.01,0.1,0.2,0.5)))

prop_vec=c(-Inf,0.01,0.1,0.2,0.5)
propMHC100N4_fix=lapply(1:5,function(i) libs50_MHC.CE.fix(q=exp(seq(log(3),log(13),length=100))[subInd][14],theta=pars_MHC[[1]][i],K=100,N=10^4,k0=1,k1=100, thre=TRUE,B=50,prop=prop_vec[i]))
propMHC100N5_fix=lapply(1:5,function(i) libs50_MHC.CE.fix(q=exp(seq(log(3),log(13),length=100))[subInd][14],theta=pars_MHC[[2]][i],K=100,N=10^5,k0=1,k1=100, thre=TRUE,B=50,prop=prop_vec[i]))
propMHC100N6_fix=lapply(1:5,function(i) libs50_MHC.CE.fix(q=exp(seq(log(3),log(13),length=100))[subInd][14],theta=pars_MHC[[3]][i],K=100,N=10^6,k0=1,k1=100, thre=TRUE,B=50,prop=prop_vec[i]))

propTMHC100N4_fix=lapply(1:5,function(i) libs50_MHC.CE.fix(q=exp(seq(log(3),log(13),length=100))[subInd][14],theta=pars_TMHC[[1]][i],K=100,N=10^4,k0=1,k1=50, thre=TRUE,B=50,prop=prop_vec[i]))
propTMHC100N5_fix=lapply(1:5,function(i) libs50_MHC.CE.fix(q=exp(seq(log(3),log(13),length=100))[subInd][14],theta=pars_TMHC[[2]][i],K=100,N=10^5,k0=1,k1=50, thre=TRUE,B=50,prop=prop_vec[i]))
propTMHC100N6_fix=lapply(1:5,function(i) libs50_HC.CE.fix(q=exp(seq(log(3),log(13),length=100))[subInd][14],theta=pars_TMHC[[3]][i],K=100,N=10^6,k0=1,k1=50, thre=TRUE,B=50,prop=prop_vec[i]))

list.save(propMHC100N4_fix,"propMHC100N4_fix.RData")
list.save(propMHC100N5_fix,"propMHC100N5_fix.RData")
list.save(propMHC100N6_fix,"propMHC100N6_fix.RData")

list.save(propTMHC100N4_fix,"propTMHC100N4_fix.RData")
list.save(propTMHC100N5_fix,"propTMHC100N5_fix.RData")
list.save(propTMHC100N6_fix,"propTMHC100N6_fix.RData")



