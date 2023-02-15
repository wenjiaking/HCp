###############################################################################################
##################FULL MODIFICATIONS OF ZHEYANG WU'S METHOD (4 STEPS)######################
###############################################################################################
library(ggplot2)
library(cowplot)
library(rlist)
library(Rmpfr)
library(Rfast)
library(GHC)
#install.packages("SetTest")
library(SetTest)

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

PM_updateindep <-
  function(ind,t_mesh,cv,PM,tailprob){
    denom=sum(PM[1:(cv[ind-1]+1),ind-1])
    for(a_k in 0:cv[ind]){
      for(a_k1 in 0:cv[ind-1]){
        PM[a_k+1,ind] = PM[a_k+1,ind] + dbinom(a_k,a_k1,tailprob[ind]/tailprob[ind-1])*PM[a_k1+1,ind-1]
      }
    }
    PM[1:(cv[ind]+1),ind] = PM[1:(cv[ind]+1),ind]/denom  
    return(PM)
  }

XH=function(GHCstat,P) {
  t_mesh = rep(0,P)
  for(k in 1:P){
    t_mesh[k]=uniroot(function(x) GHCstat*sqrt(P*2*pnorm(x,lower.tail=F)*(1-2*pnorm(x,lower.tail=F)))+P*2*pnorm(x,lower.tail=F)-(P-k+1),interval=c(10^-10,50))$root
  }
  m=P;tailprob = 2*pnorm(t_mesh,lower.tail=F);cv=(P-1):0
  pval_prod=rep(0,m)
  PM = matrix(0,nrow=cv[1]+1,ncol=m) #P*P matrix
  PM[,1] = dbinom(0:cv[1],P,tailprob[1])
  pval_prod[1]=sum(PM[1:(cv[1]+1),1])
  for(j in 2:m){
    PM=PM_updateindep(j,t_mesh,cv,PM,tailprob)
    pval_prod[j]=sum(PM[1:(cv[j]+1),j])
  }
  return(1-prod(pval_prod))
}

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


CE.HC.mixed.grad=function(K,theta,N1=1000,k0,k1,thre,idx=1) {
  grad=1/((1:K)^idx+1)
  mix.prop=matrix(rbinom(K*N1,size=1,prob=grad),nrow=K)
  mix.prop=transpose(mix.prop)
  x=matrix(rnorm(K*N1,0,theta),nrow = N1)^{(mix.prop)}+matrix(rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  stat.val=HC.stat(x=x,k0=k0,k1=k1,thre=thre)
  log_g_norm <- rowsums(log(dnorm(x))) # numerator part
  
  log_g_mix <- rowsums(log(transpose(transpose(dnorm(x,0,theta))*grad)+transpose(transpose(dnorm(x))*(1-grad))))
  
  w=exp(log_g_norm-log_g_mix)
  #s=rowsums(x[,grad==1]^2)
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
  par=opt$maximum
  return(list(gamma=gamma,par=par))
}

CE.mixed.grad=function(K,ro,N0=10^4, N1=10^4,B=10^3,q,k0,k1,thre,idx=1,theta=1) {
  
  gamma=-Inf
  t=0
  while (gamma<q) {
    t=t+1
    obj=CE.HC.mixed.grad(K=K,theta=theta,N1=N0, k0=k0,k1=k1,thre=thre,idx=idx)
    par=mixed.par.update.grad(object = obj,ro=ro,K=K,idx = idx)
    gamma=par$gamma
    theta=par$par
    print(paste0(gamma,", ",theta))
  }
  print("stop")
  
  objs=lapply(1:B,function(b) CE.HC.mixed.grad(K=K,theta=theta,N1=N1,k0=k0,k1=k1,thre=thre,idx=idx))
  stat.val.final=as.vector(sapply(objs,"[[",1))
  weight.final=as.vector(sapply(objs,"[[",2))
  p.val=1/(N1*B)*sum(weight.final[stat.val.final>=q])
  return(p.val)
}

libs_HC.CE.mixGrad=function(q.val.set,K=50,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=50,thre=FALSE,idx=1,theta=1){
  temp.lib=lapply(q.val.set, function(x) CE.mixed.grad(K=K,ro=ro,N0=N0,N1=N1,B=B,q=x,k0=k0,k1=k1,thre=thre,idx=idx,theta=theta))
  temp.p=sapply(temp.lib,"[[",1)
  return(temp.p)
}



CE.HC.mixed.prop=function(K,theta,N1=1000,k0,k1,thre,prop=0.2) {
  if (prop==-Inf) {grad=rep(1,K)}
  else {grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}
  #else {grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}
  
  #grad=1/((1:K)^(0.1)+1)
  mix.prop=matrix(rbinom(K*N1,size=1,prob=grad),nrow=K)
  mix.prop=transpose(mix.prop)
  
  x=matrix(rnorm(K*N1,0,theta),nrow = N1)^{(mix.prop)}+matrix(rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  stat.val=HC.stat(x=x,k0=k0,k1=k1,thre=thre)
  log_g_norm <- rowsums(log(dnorm(x))) # numerator part
  
  log_g_mix <- rowsums(log(transpose(transpose(dnorm(x,0,theta))*grad)+transpose(transpose(dnorm(x))*(1-grad))))
  
  w=exp(log_g_norm-log_g_mix)
  #s=rowsums(x[,grad==1]^2)
  return(list(stat.val=stat.val,weight=w,x=x))
}

mixed.par.update.prop=function(object,ro,K,prop) {
  
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

CE.mixed.prop=function(K,ro,N0=10^4, N1=10^4,B=10^3,q,k0,k1,thre,prop=0.2,theta=1) {
  #theta=2
  gamma=-Inf
  t=0
  while (gamma<q) {
    t=t+1
    obj=CE.HC.mixed.prop(K=K,theta=theta,N1=N0, k0=k0,k1=k1,thre=thre,prop = prop)
    par=mixed.par.update.prop(object = obj,ro=ro,K=K,prop = prop)
    gamma=par$gamma
    theta=par$par
    print(paste0(gamma,", ",theta))
  }
  print("stop")
  
  objs=lapply(1:B,function(b) CE.HC.mixed.prop(K=K,theta=theta,N1=N1,k0=k0,k1=k1,thre=thre,prop = prop))
  stat.val.final=as.vector(sapply(objs,"[[",1))
  weight.final=as.vector(sapply(objs,"[[",2))
  p.val=1/(N1*B)*sum(weight.final[stat.val.final>=q])
  
  #p.val=1/(N1*B)*sum(weight.final[stat.val.final>=q])
  return(p.val)
}

libs_HC.CE.mixprop=function(q.val.set,K=50,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=50,thre=FALSE,prop=0.2,theta=1){
  #RNGkind("L'Ecuyer-CMRG")
  temp.lib=lapply(q.val.set, function(x) CE.mixed.prop(K=K,ro=ro,N0=N0,N1=N1,B=B,q=x,k0=k0,k1=k1,thre=thre,prop=prop,theta=theta))
  temp.p=sapply(temp.lib,"[[",1)
  return(temp.p)
}

########Li & Siegmund approximation######
Cx_HC=function(x,eps) {
  C=(2*x+eps^2-eps*sqrt(eps^2+4*(1-x)*x))/(2*(1+eps^2))
  deltaC=1/(1+eps^2)-eps*(1-2*x)/((1+eps^2)*sqrt(eps^2+4*(1-x)*x))
  return(list(C=C, deltaC=deltaC))
}
LiAppro_HC=function(K,b, k0, k1, thre=FALSE) {
  k_range=k0:k1
  k_Cs=Cx_HC(x=k_range/K,eps=b/sqrt(K))
  k_roots=k_Cs$C
  k_deltaC=k_Cs$deltaC
  k_coefs=1-((K-k_range+1)*k_deltaC)/(K*(1-k_roots))
  if (thre) {
    k_Bins=dbinom(k_range,K,sapply(k_roots,function(r) max(r,1/K)))-dbinom(k_range,K,1/K)*sapply(K*k_roots,function(r) max(r,1))
  }
  else {
    k_Bins=dbinom(k_range,K,k_roots)
  }
  return(sum(k_coefs*k_Bins))
}

MC.est=function(K,N1=10^4,B=10^3,k0,k1,thre,q) {
  stat.vals=lapply(1:B,function(b) {
    x=matrix(Rfast::Rnorm(K*N1,0,1),nrow = N1)
    stat.val=HC.stat(x=x,k0=k0,k1=k1,thre=thre)
    return(stat.val)
  })
  stat.vals=unlist(stat.vals)
  p=mean(stat.vals>=q)
  return(p)
}


HCp=function(HC_flibs,K,stats,cores=NULL) {
  #HC_flibs=list.load(flibsp_dir)
  f=HC_flibs[[K-1]]
  if (is.null(cores)) {
    ps=sapply(stats, function(x){
      return(exp(f(log(x))))
    } )
  }
  else {
    ps=mclapply(stats, function(i){
      return(exp(f(log(x))))
    } ,mc.cores=cores)
    ps=unlist(ps)
  }
  
  return(ps)
}

#####Prepare Data for Figure II (A) and (B)#########
q.vals50=exp(seq(log(10^3),log(10^7),length.out = 200))
ZY500.HCorig=sapply(q.vals50, function(x) 1-phc(x, M=diag(rep(1,500)), k0=1, k1=500))
saveRDS(ZY500.HCorig,"ZY500.HCorig.rds")

ZY500=sapply(q.vals50, function(x) 1-phc.mod(x, M=diag(rep(1,500)), k0=1, k1=500)) #Time difference of 48.01968 secs
saveRDS(ZY500,"ZY500.rds")

XH500=sapply(q.vals50,function(x) XH(GHCstat = x,P=500))
saveRDS(XH500,"XH500.rds")

libsHC500=lapply(1:15, function(i) libs_HC.CE.mixGrad(q.val.set=q.vals50,K=500,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=500, idx=1,thre=FALSE))
list.save(libsHC500,"libs500.RData")

ZY_modtrunc500=sapply(exp(seq(log(30),17,length=100)), function(x) 1-phc.mod(x, M=diag(rep(1,500)), k0=1, k1=250))
ZY_trunc500=sapply(exp(seq(log(30),17,length=100)), function(x) 1-phc(x, M=diag(rep(1,500)), k0=1, k1=250))
saveRDS(ZY_trunc500,"ZY_trunc500.rds")
saveRDS(ZY_modtrunc500,"ZY_modtrunc500.rds")
mixgradlibs_modHC500=lapply(1:15, function(i) libs_HC.CE.mixGrad(q.val.set=exp(seq(log(30),17,length=200)),K=500,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=250, idx=1,thre=FALSE))
list.save(mixgradlibs_modHC500,"mixgradlibs_modHC500.RData")

######Prepare data for Figure II (C) and (D)############3
subInd=seq(1,100,length.out = 15)
mixgrad5libs50_MHC100=mclapply(1:15, function(i) libs_HC.CE.mixprop(q.val.set=seq(log(3),log(13),length=100)[subInd],K=100,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=100, thre=TRUE), mc.cores=30)
saveRDS(do.call(cbind,mixgrad5libs50_MHC100),"mixgrad5libs50N7_MHC100.rds")

mixgrad5libs50_TMHC100=mclapply(1:15, function(i) libs_HC.CE.mixprop(q.val.set=seq(log(3),log(13),length=100)[subInd],K=100,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=50, thre=TRUE), mc.cores=30)
saveRDS(do.call(cbind,mixgrad5libs50_TMHC100),"mixgrad5libs50N7_TMHC100.rds")

ZY_MHC100=sapply(exp(seq(log(3),log(13),length=100)),function(x) 1-phc.mod(q=x,M=diag(1,100),k0=1,k1=100,MHC = TRUE))
saveRDS(ZY_MHC100,"ZY_MHC100.rds")

ZY_TMHC100=sapply(exp(seq(log(3),log(13),length=100)),function(x) 1-phc.mod(q=x,M=diag(1,100),k0=1,k1=50,MHC = TRUE))
saveRDS(ZY_TMHC100,"ZY_TMHC100.rds")

######Prepare data for Figure III ##############
mixgradSplib_MHC100=lapply(1:50,function(i) libs_HC.CE.mixprop(q.val.set=c(exp(seq(log(3.8),log(13),length.out = 150))[seq(80,150,length.out = 12)],exp(seq(log(13.5),log(14.5),length.out=3))),K=100,ro=0.01,N0=10^4, N1=10^4,B=10^3,k0=1,k1=100, thre=TRUE))
saveRDS(do.call(cbind,mixgradSplib_MHC100),"mixgradSplib_MHC100.rds")

mixgradSplib_TMHC100=lapply(1:50,function(i) libs_HC.CE.mixprop(q.val.set=c(exp(seq(log(3.8),log(13),length.out = 150))[seq(80,150,length.out = 12)],exp(seq(log(13.5),log(14.5),length.out=3))),K=100,ro=0.01,N0=10^4, N1=10^4,B=10^3,k0=1,k1=50, thre=TRUE))
saveRDS(do.call(cbind,mixgradSplib_TMHC100),"mixgradSplib_TMHC100.rds")

MCLpMHC100lib=lapply(1:50,function(i) sapply(c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5])), function(x) MC.est(100,N1=10^4,B=10^3,k0=1,k1=100,thre=TRUE,q=x)))
saveRDS(do.call(cbind,MCLpMHC100lib),"MCLpMHC100lib.rds")
MCLpTMHC100lib=mclapply(1:50, function(i) sapply(c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5])), function(x) MC.est(100,N1=10^4,B=10^3,k0=1,k1=50,thre=TRUE,q=x)),
                        mc.cores=30)
saveRDS(do.call(cbind, MCLpTMHC100lib),"MCLpTMHC100lib.rds")

LSapproSpMHC100=sapply(c(exp(seq(log(3.8),log(13),length.out = 150))[seq(80,150,length.out = 12)],exp(seq(log(13.5),log(14.5),length.out=3))),function(x) LiAppro_HC(K=100,b=x,k0=1,k1=100,thre = TRUE))
LSapproSpTMHC100=sapply(c(exp(seq(log(3.8),log(13),length.out = 150))[seq(80,150,length.out = 12)],exp(seq(log(13.5),log(14.5),length.out=3))),function(x) LiAppro_HC(K=100,b=x,k0=1,k1=50,thre = TRUE))
saveRDS(LSapproSpTMHC100,"LSapproSpTMHC100.rds")
saveRDS(LSapproSpMHC100,"LSapproSpMHC100.rds")

q_sm=exp(seq(log(0.5),log(3.7),length.out = 50)) 
ZY_SpMHC100=sapply(c(exp(seq(log(3.8),log(13),length.out = 150))[seq(80,150,length.out = 12)],exp(seq(log(13.5),log(14.5),length.out=3))),function(x) 1-phc.mod(x, M=diag(1,nrow = 100), k0=1, k1=100, LS = F, ZW = F, onesided=FALSE,MHC=TRUE))
ZY_LpMHC100=sapply(c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5])), function(x) 1-phc.mod(x, M=diag(1,nrow = 100), k0=1, k1=100, LS = F, ZW = F, onesided=FALSE,MHC=TRUE))
saveRDS(ZY_SpMHC100,"ZY_SpMHC100.rds")
saveRDS(ZY_LpMHC100,"ZY_LpMHC100.rds")

mod_SpTMHC100=sapply(c(exp(seq(log(3.8),log(13),length.out = 150))[seq(80,150,length.out = 12)],exp(seq(log(13.5),log(14.5),length.out=3))),function(x) 1-phc.mod(x, M=diag(1,nrow = 100), k0=1, k1=50, LS = F, ZW = F, onesided=FALSE,MHC=TRUE))
mod_LpTMHC100=sapply(c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5])), function(x) 1-phc.mod(x, M=diag(1,nrow = 100), k0=1, k1=50, LS = F, ZW = F, onesided=FALSE,MHC=TRUE))
saveRDS(mod_SpTMHC100,"mod_SpTMHC100.rds")
saveRDS(mod_LpTMHC100,"mod_LpTMHC100.rds")

######Prepare data for Figure IV ##############
set.seed(923)
Kset=c(30,100,500,1000,2000)
nset=c(1,10,100,1000,10000,100000)
q.vals.HC=c(exp(seq(log(1),log(9*10^2),length.out = 200)),exp(seq(log(10^3),log(10^7),length.out = 200)))
compHC.q=lapply(nset,function(x) sample(min(q.vals.HC):max(q.vals.HC),x))

set.seed(923)
Kset2=c(50,200,300,400,600,700,800,900,seq(1100,1900,by=100))
nset2=100
compHC.q2=sample(min(q.vals.HC):max(q.vals.HC),nset2)

libtimes=function(Kset,HC_flibs, compHC.q,cores=NULL) {
  res=lapply(compHC.q,function(x) lapply(Kset,function(k) {
    start.time=Sys.time()
    p=HCp(HC_flibs,k,x,cores=cores)
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  }))
  return(res)
}

libtimes2=function(Kset,HC_flibs, compHC.q,cores=NULL) {
  res=lapply(Kset, function(k) {
    start.time=Sys.time()
    p=HCp(HC_flibs,k,compHC.q,cores=cores)
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  })
  return(res)
}
HC_flibs=list.load("HC_flibs.RData")
HCres_libs=libtimes(Kset,HC_flibs,compHC.q,cores=NULL)
HCres2_libs=libtimes2(Kset2,HC_flibs,compHC.q2,cores=NULL)
list.save(HCres_libs,"HCres_libs.RData")
list.save(HCres2_libs,"HCres2_libs.RData")

HCIStimes=function(Kset,compHC.q) {
  res=lapply(compHC.q,function(x) lapply(Kset,function(k) {
    start.time=Sys.time() 
    p=lapply(x, function(q) CE.mixed.grad(K=k,ro=0.01,N=10^4,N1=10^4,B=1,q=q,idx=1,k0=1,k1=k,thre=FALSE,theta=1))
    t=Sys.time()-start.time
    return(list(time=t,ps=sapply(p,"[[",1)))
  }))
  return(res)
}
HCIStimes2=function(Kset,compHC.q) {
  res=lapply(Kset,function(k) {
    start.time=Sys.time()
    p=lapply(compHC.q, function(q) CE.mixed.grad(K=k,ro=0.01,N=10^4,N1=10^4,B=1,q=q,idx=1,k0=1,k1=k,thre=FALSE,theta=1))
    t=Sys.time()-start.time
    return(list(time=t,ps=sapply(p,"[[",1)))
  })
  return(res)
}

HCres_IS.sub=HCIStimes(Kset,compHC.q[1:5])
HCres_IS.lg=lapply(Kset,function(k) {
  start.time=Sys.time()
  p=lapply(compHC.q[[6]], function(q) CE.mixed.grad(K=k,ro=0.01,N=10^4,N1=10^4,B=1,q=q,idx=1,k0=1,k1=k,thre=FALSE,theta=1))
  t=Sys.time()-start.time
  return(list(time=t,ps=sapply(p,"[[",1))) #speed5 3291342
})
list.save(HCres_IS.sub,"HCres_IS.sub.RData")
list.save(HCres_IS.lg,"HCres_IS.lg.RData")
HCres2_IS=HCIStimes2(Kset2,compHC.q2) 
list.save(HCres2_IS,"HCres2_IS.RData")

HCmodtimes=function(Kset,compHC.q) {
  res=lapply(compHC.q,function(x) lapply(Kset,function(k) {
    start.time=Sys.time()
    p=lapply(x, function(q) 1-phc.mod(q, M=diag(rep(1,k)), k0=1, k1=k))
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  }))
  return(res)
}
HCmodtimes2=function(Kset,compHC.q) {
  res=lapply(Kset,function(k) {
    start.time=Sys.time()
    p=lapply(compHC.q, function(q) 1-phc.mod(q, M=diag(rep(1,k)), k0=1, k1=k))
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  })
  return(res)
}
HCres_mod=HCmodtimes(Kset,compHC.q)
HCres2_mod=HCmodtimes2(Kset2,compHC.q2)
list.save(HCres_mod,"HCres_mod.RData")
list.save(HCres2_mod,"HCres2_mod.RData")

HCLintimes=function(Kset,compHC.q) {
  res=lapply(compHC.q,function(x) lapply(Kset,function(k) {
    start.time=Sys.time()
    p=lapply(x, function(q) XH(GHCstat = q,P=k))
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  }))
  return(res)
}
HCLintimes2=function(Kset,compHC.q) {
  res=lapply(Kset, function(k) {
    start.time=Sys.time()
    p=lapply(compHC.q, function(q) XH(GHCstat = q,P=k))
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  })
  return(res)
}

HCres2_Lin=HCLintimes2(Kset2,compHC.q2) 
list.save(HCres2_Lin,"HCres2_Lin.RData")

HCZYtimes=function(Kset,compHC.q) {
  res=lapply(compHC.q,function(x) lapply(Kset,function(k) {
    start.time=Sys.time()
    p=lapply(x, function(q) 1-phc(q, M=diag(rep(1,k)), k0=1, k1=k))
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  }))
  return(res)
}
HCZYtimes2=function(Kset,compHC.q) {
  res=lapply(Kset,function(k) {
    start.time=Sys.time()
    p=lapply(compHC.q, function(q) 1-phc(q, M=diag(rep(1,k)), k0=1, k1=k))
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  })
  return(res)
}
HCres_ZY=HCZYtimes(Kset,compHC.q) 
HCres2_ZY=HCZYtimes2(Kset2,compHC.q2)
list.save(HCres_ZY,"HCres_ZY.RData")
list.save(HCres2_ZY,"HCres2_ZY.RData")

HCLStimes=function(Kset,compHC.q) {
  res=lapply(compHC.q,function(x) lapply(Kset,function(k) {
    start.time=Sys.time()
    p=lapply(x, function(q) LiAppro_HC(K=k,b=q,  k0=1, k1=k,thre = F))
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  }))
  return(res)
}
HCLStimes2=function(Kset,compHC.q) {
  res=lapply(Kset,function(k) {
    start.time=Sys.time()
    p=lapply(compHC.q, function(q) LiAppro_HC(K=k,b=q, k0=1, k1=k,thre = F))
    t=Sys.time()-start.time
    return(list(time=t,ps=unlist(p)))
  })
  return(res)
}
HCres_LS=HCLStimes(Kset,compHC.q) 
HCres2_LS=HCLStimes2(Kset2,compHC.q2) 
list.save(HCres_LS,"HCres_LS.RData")
list.save(HCres2_LS,"HCres2_LS.RData")

