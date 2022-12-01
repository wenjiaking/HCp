
library(Rmpfr)
########Prepare funs of MST for HC or THC##########
pphi.mod <- function(q, K,k0, k1, s=2, t=30, MHC=FALSE)
{
  qtemp = q
  M=diag(1,K)
  q = max(0.01, q)
  NODES=seq(-4,4,length.out=40)
  n = length(M[1,])
  nrep = n
  nthr = 5
  lower = q*0.9
  upper = q*1.1
  thr = seq(lower,upper,length.out=nthr)
  rho = sample(M[row(M)!=col(M)], nrep)
  rho = abs(rho)
  digi = rep(0.02, nrep)
  rho = round(rho/digi)*digi
  rho[rho==1] = 0.99
  tbl = table(rho)
  rho_unique = as.numeric(names(tbl))
  rho_freq = as.numeric(tbl)

  nrep_uniq = length(rho_unique)
  pvalue_cal = c()
  for(i in 1:nrep_uniq){
    freq = rho_freq[i]
    y = log(sapply(1:nthr, function(x)pphi.rho.mod(threshold=thr[x], n=n, rho=rho_unique[i], k0=k0, k1=k1, NODES=NODES, s=s, t=t,MHC = MHC)))
    #print(y)
    pvalue_cal = append(pvalue_cal,rep(y, freq)) #only depends on pphi.rho function
  }
  data = data.frame(y=pvalue_cal, x=rep(log(thr), nrep))

  fit = stats::loess(y~x, data = data)
  res = exp(stats::predict(fit, log(q)))
  if(qtemp>=0.01){
    return(res)
  }else{
    warning(paste0("Left-tail prob. < ",round(res,3), " or p-value > ", 1-round(res,3)))
    return(res)
  }
}

phi.f.inv.i.mod <- function(q, i, n, s){
  ep = 10^(-18)
  if(i!=n){
    if(q>0){
      if (s!=1&&s!=0){
        f <- function(x)(1-(i/n)^s*x^(1-s)-(1-i/n)^(s)*(1-x)^(1-s))*2*n/(s-s^2) - q^2
        if (f(ep)<0){
          CC=0
        }else{
          CC=stats::uniroot(f,c(ep,i/n),tol=1e-100)$root
        }
      }else if (s==1){
        f <- function(x)2*n*((i/n)*log(i/n/x)+(1-i/n)*log((1-i/n)/(1-x)))-q^2
        if (f(ep)<0){
          CC=0
        }else{
          CC=stats::uniroot(f,c(ep,i/n),tol=1e-100)$root
        }
      }else{
        f <- function(x)2*n*((x)*log(x/(i/n))+(1-x)*log((1-x)/(1-i/n)))-q^2
        if (f(ep)<0){
          CC=0
        }
        else{
          CC=stats::uniroot(f,c(ep,i/n),tol=1e-100)$root
        }
      }
    }else if(q<0){
      if (s!=1&&s!=0){
        f <- function(x)(1-(i/n)^s*x^(1-s)-(1-i/n)^(s)*(1-x)^(1-s))*2*n/(s-s^2) - q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=stats::uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
        }
      }else if (s==1){
        f <- function(x)2*n*((i/n)*log(i/n/x)+(1-i/n)*log((1-i/n)/(1-x)))-q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=stats::uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
        }
      }else{
        f <- function(x)2*n*((x)*log(x/(i/n))+(1-x)*log((1-x)/(1-i/n)))-q^2
        if (f(1-ep)<0){
          CC=1
        }else{
          CC=stats::uniroot(f,c(i/n, 1-ep),tol=1e-100)$root
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
        CC=stats::uniroot(f,c(ep,1),tol=1e-100)$root
      }else{
        f <- function(x)(1-x^(1-s))*2*n/(s-s^2) - q^2
        CC=stats::uniroot(f,c(ep,1),tol=1e-100)$root
      }
    }
  }
  return(CC)
}

pphi.rho.mod <- function(threshold, n, rho, k0, k1, NODES=seq(-4,4,length.out=40),
                         s=2, t=30, MHC=FALSE){
  t = k1-k0+1
  u = sapply(1:k1,function(x)phi.f.inv.i.mod(threshold, x, n, s))
  um = u[length(u)]
  if(rho!=0){
    SUM = sum(stats::dnorm(NODES))
    result = 0
    for(z in NODES){
      d = stats::pnorm((stats::qnorm(1-u/2)-sqrt(rho)*z)/sqrt(1-rho),lower.tail = F) +
        stats::pnorm((-stats::qnorm(1-u/2)-sqrt(rho)*z)/sqrt(1-rho))
      dm = stats::pnorm((stats::qnorm(1-um/2)-sqrt(rho)*z)/sqrt(1-rho),lower.tail = F) +
        stats::pnorm((-stats::qnorm(1-um/2)-sqrt(rho)*z)/sqrt(1-rho))
      result = result + stats::dnorm(z)*UnifCross_v1.mod(d,dm,t=t,n=n,k0=k0,k1=k1)/SUM
    }
    return(result)
  }else{
    if (MHC) return(UnifCross_v1_MHC(u,um,t=t,n=n,k0=k0,k1=k1))
    else return(UnifCross_v1.mod(u,um,t=t,n=n,k0=k0,k1=k1)) #for independent studies, only depends on UnifCross_v1 function
  }
}

UnifCross_v1.mod <- function(u, um, t, n, k0, k1){

  m = floor(k1)

  pp = rep(NA, t)
  pp[1] = stats::pbeta(um,m,n-m+1,lower.tail=F)
  pp[2:t] = (lfactorial(n) - lfactorial(n-(k0:(k0+t-2))))+log(stats::pbeta(um,m-(k0:(k0+t-2)),n-m+1,lower.tail=F))
  a = rep(1,t)
  a[2] = -stats::dpois(k0,u[k0])*exp(u[k0])
  if(t>1){
    for (i in 2:(t-1)){
      d = stats::dpois(c(1:(i-1), i+k0-1),u[i+k0-1])*exp(u[i+k0-1])

      a[i+1]=-a[i:1]%*%d
    }

    a_pos=which(a[2:t]>0)+1
    a_neg=which(a[2:t]<0)+1
    p_val=a[1]*pp[1]+sum(exp(pp[a_pos]+log(a[a_pos])))-sum(exp(pp[a_neg]+log(-a[a_neg])))
    return(p_val)

  }else{
    p = stats::pbeta(u[k0], k0, n-k0+1)
    return(1-drop(p))
  }
}


UnifCross_v1_MHC <- function(u, um, n, t,k0, k1){
  alpha=1/n #lower threshold
  m = floor(k1)
  u=sapply(u,function(x) max(x,alpha))
  um=max(um,alpha)

  pp = Rmpfr::mpfr(rep(NA, t),80)
  pp[1]=Rmpfr::mpfr(sum(stats::pbinom(0:(k0-1),size = n,prob=alpha)*stats::pbeta((um-alpha)/(1-alpha),m:(m-k0+1),n-m+1,lower.tail=F)),80)
  pp[2:t] = Rmpfr::mpfr(exp(lfactorial(n) - lfactorial(n-(k0:(k0+t-2))))*stats::pbeta(um,m-(k0:(k0+t-2)),n-m+1,lower.tail=F),80)
  a = Rmpfr::mpfr(rep(1,t),80)
  a[2] =Rmpfr::mpfr((alpha^k0-stats::pbinom(k0-1,k0,alpha/u[k0])*u[k0]^k0)/factorial(k0),80)
  if(t>1){
    for (i in 2:(t-1)){
      d = Rmpfr::mpfr(c((-stats::dpois(1:(i-1),u[i+k0-1])*exp(u[i+k0-1])),
                 (alpha^(i+k0-1)-stats::pbinom(k0-1,i+k0-1,alpha/u[i+k0-1])*u[i+k0-1]^(i+k0-1))/factorial(i+k0-1)),80)

      a[i+1]=a[i:1]%*%d
    }

    return(Rmpfr::asNumeric(drop(pp[1:t]%*%a[1:t])))
  }else{
    p = stats::pbeta(u[k0], k0, n-k0+1)
    return(1-drop(p))
  }
}


#####Prepare functions for IS#########
HC.stat=function(x, k0=1, k1,thre=FALSE) {
  t_mesh = Rfast::rowSort(abs(x),descending = TRUE)
  K=dim(t_mesh)[2]
  n=dim(t_mesh)[1]
  tailprob = 2*stats::pnorm(t_mesh,lower.tail=F)
  tailprob=tailprob[,k0:k1]
  if (thre) {
    tailprob[which(tailprob<1/K)]=k1/K
  }
  GHC.mesh=(matrix(rep(k0:k1,times=n),nrow = n,byrow = T)-K*tailprob)/sqrt(K*tailprob*(1-tailprob))
  GHCstat = Rfast::rowMaxs(GHC.mesh,value = T)
  return(stat=GHCstat)
}

CE.HC.mixed.prop=function(K,theta,N1=1000,k0,k1,thre,prop=0.2) {

  if (prop==-Inf) {grad=rep(1,K)}
  else {grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}
  mix.prop=matrix(stats::rbinom(K*N1,size=1,prob=grad),nrow=K)
  mix.prop=Rfast::transpose(mix.prop)

  x=matrix(Rfast::Rnorm(K*N1,0,theta),nrow = N1)^{(mix.prop)}+matrix(Rfast::Rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  stat.val=HC.stat(x=x,k0=k0,k1=k1,thre=thre)
  log_g_norm <- Rfast::rowsums(log(stats::dnorm(x))) # numerator part

  log_g_mix <- Rfast::rowsums(log(Rfast::transpose(Rfast::transpose(stats::dnorm(x,0,theta))*grad)+Rfast::transpose(Rfast::transpose(stats::dnorm(x))*(1-grad))))

  w=exp(log_g_norm-log_g_mix)
  #s=rowsums(x[,grad==1]^2)
  return(list(stat.val=stat.val,weight=w,x=x))
}

mixed.par.update.prop=function(object,ro,K,prop) {

  stat.val=object[[1]]
  gamma=stats::quantile(stat.val,1-ro)
  weight=object[[2]]
  x.mat=Rfast::transpose(object[[3]])
  B=length(stat.val)
  if (prop==-Inf) {grad=rep(1,K)}
  else {grad=c(rep(0,ceiling(K^(prop))),rep(1,floor(K-K^(prop))))}

  opt=stats::optimize(function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,function(x) sum(log(grad*stats::dnorm(x,0,t)+(1-grad)*stats::dnorm(x,0,1))))),
                      interval = c(0,5),maximum = TRUE)

  par=opt$maximum
  return(list(gamma=gamma,par=par))
}

CE.HC.mixed.grad=function(K,theta,k0=1,k1,N1=1000,idx=1,thre=FALSE) {
  grad=1/((1:K)^idx+1)
  mix.prop=matrix(stats::rbinom(K*N1,size=1,prob=grad),nrow=K)
  mix.prop=Rfast::transpose(mix.prop)

  x=matrix(Rfast::Rnorm(K*N1,0,theta),nrow = N1)^{(mix.prop)}+matrix(Rfast::Rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  stat.val=HC.stat(x=x,k0=k0, k1=k1,thre=thre)
  log_g_norm <- Rfast::rowsums(log(stats::dnorm(x))) # numerator part

  log_g_mix <- Rfast::rowsums(log(Rfast::transpose(Rfast::transpose(stats::dnorm(x,0,theta))*grad)+Rfast::transpose(Rfast::transpose(stats::dnorm(x))*(1-grad))))

  w=exp(log_g_norm-log_g_mix)
  return(list(stat.val=stat.val,weight=w,x=x))
}

mixed.par.update.grad=function(object,ro,K,idx=1) {

  stat.val=object[[1]]
  gamma=stats::quantile(stat.val,1-ro)
  weight=object[[2]]
  x.mat=Rfast::transpose(object[[3]])
  B=length(stat.val)
  grad=1/((1:K)^idx+1)

  opt=stats::optimize(function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,function(x) sum(log(grad*stats::dnorm(x,0,t)+(1-grad)*stats::dnorm(x,0,1))))),
                      interval = c(0,5),maximum = TRUE)
  #opt=optim(2,function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,function(x) sum(log(grad*dnorm(x,0,t)+(1-grad)*dnorm(x,0,1)))),method="L-BFGS-B"),
  # control=list(fnscale=-1))
  par=opt$maximum
  #par=opt$par
  return(list(gamma=gamma,par=par))
}



########Prepare funs of Li-Siegmund##########

Cx_HC=function(x,eps) {
  C=(2*x+eps^2-eps*sqrt(eps^2+4*(1-x)*x))/(2*(1+eps^2))
  deltaC=1/(1+eps^2)-eps*(1-2*x)/((1+eps^2)*sqrt(eps^2+4*(1-x)*x))
  return(list(C=C, deltaC=deltaC))
}

########Prepare funs of UFI##########


