#' Barnett-Lin method of computing p-value for FHC
#'
#' @param q quantile or observed value of HC statistic, should be a scalar
#' @param K dimension of HC, i.e. the total number of studies
#'
#' @return The right-tail probability of the null distribution of HC statistic at the given quantile, i.e. p-value.
#' @export
#'
#' @examples
#' pset=runif(10)
#' hcstat <-HCstat(pset,k0=1)
#' XH(q=hcstat, K=10)
XH=function(q,K) {
  P=K
  GHCstat=q
  t_mesh = rep(0,P)
  for(k in 1:P){
    t_mesh[k]=stats::uniroot(function(x) GHCstat*sqrt(P*2*stats::pnorm(x,lower.tail=F)*(1-2*stats::pnorm(x,lower.tail=F)))+P*2*stats::pnorm(x,lower.tail=F)-(P-k+1),interval=c(10^-10,50))$root
  }
  m=P;tailprob = 2*stats::pnorm(t_mesh,lower.tail=F);cv=(P-1):0
  pval_prod=rep(0,m)
  PM = matrix(0,nrow=cv[1]+1,ncol=m) #P*P matrix
  PM[,1] = stats::dbinom(0:cv[1],P,tailprob[1])
  pval_prod[1]=sum(PM[1:(cv[1]+1),1])
  for(j in 2:m){
    PM=PM_updateindep(j,t_mesh,cv,PM,tailprob)
    pval_prod[j]=sum(PM[1:(cv[j]+1),j])
  }
  return(1-prod(pval_prod))
}
