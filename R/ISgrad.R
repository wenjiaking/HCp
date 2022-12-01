#' Calculate p-value of HC by cross-entropy based importance sampling method with mix Gussian gradient proposal family
#'
#' @param q  quantile or observed value of HC statistic, should be a scalar
#' @param K  dimension of HC, i.e. the total number of studies
#' @param ro  related to the heavy tailed proposal density, 0.01 by default
#' @param N0  sampling size for automatically selecting the proposal density
#' @param N1  sampling size of each epoch for estimating p-value, and the total size is N1*B
#' @param B  number of epoches
#' @param k0  search range starts from the k0th smallest p-value, default value is 1.
#' @param k1  search range ends at the k1th smallest p-value, default value is K.
#' @param thre  whether using thresholding supremum domain, Boolean variable equals False by default.
#' @param idx  hyper parameter of gradient degree to specify the Gaussian mixture family, 1 by default.
#' @param theta  initial parameter value of the proposal density family, 1 by default.
#'
#' @return The right-tail probability of the null distribution of HC statistic at the given quantile, i.e. p-value.
#' @export
#'
#' @examples
#' pset=runif(10)
#' hcstat <-HCstat(pset,k0=1)
#' CE.mixed.grad(q=hcstat, K=10)
CE.mixed.grad=function(q,K,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=NA,thre=FALSE,idx=1,theta=1) {
  k0=floor(k0)
  if(is.na(k1)){
    k1 = K
  }else{
    k1 = floor(k1)
  }
  gamma=-Inf
  t=0
  while (gamma<q) {
    t=t+1
    obj=CE.HC.mixed.grad(K=K,theta=theta,k0=k0,k1=k1,N1=N0,idx=idx,thre=thre)
    par=mixed.par.update.grad(object = obj,ro=ro,K=K,idx = idx)
    gamma=par$gamma
    theta=par$par
    #print(paste0(gamma,", ",theta))
  }
  #print("stop")

  objs=lapply(1:B,function(b) CE.HC.mixed.grad(K=K,theta=theta,k0=k0,k1=k1,N1=N1,idx=idx,thre=thre))
  stat.val.final=as.vector(sapply(objs,"[[",1))
  weight.final=as.vector(sapply(objs,"[[",2))
  p.val=1/(N1*B)*sum(weight.final[stat.val.final>=q])

  return(p.val)
}
