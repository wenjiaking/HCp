#' Hybrid method of computing p-value for HC of four popular variations
#' @description This method is recommended for four popular HC variations of arbitrary dimension k: HC of full supremum domain (FHC), HC of half truncated supremum domain (THC), HC of thresholding supremum domain (MHC) or HC of thresholding and half truncated supremum domain (TMHC)
#' @param q quantile or observed value of HC statistic, should be a scalar
#' @param K dimension of HC, i.e. the total number of studies
#' @param flibs choose the spline function data set of one specific variation of HC form HC_flibs, THC_flibs, MHC_flibs or TMHC_flibs, and the default is HC_flibs.
#' @param N Samping size when sampling method is recommended, default value is 10^6.
#'
#' @return The right-tail probability of the null distribution of HC statistic at the given quantile, i.e. p-value.
#' @export
#'
#' @examples
#' pset=runif(2001)
#' hcstat <-HCstat(pset,k0=1,k1=NA, thre=FALSE)
#' hybridSpec(q=hcstat, K=2001, flibs=HC_flibs,N=10^6)
#'
hybridSpec=function(q,K,flibs=HC_flibs, N=10^6) {
  HCtype=deparse(substitute(flibs))
  if (HCtype %in% c("HC_flibs","THC_flibs","MHC_flibs","TMHC_flibs")) {
    if (HCtype=="HC_flibs") {
      if (K>1000) {
        return(mst(q=q, K=K, k0=1, k1=K, thre=FALSE))
      }
      else {
        f=flibs[[K-1]]
        return(exp(f(log(q))))
      }
    }
    if (HCtype=="THC_flibs") {
      if (K>1000) {
        return(mst(q=q, K=K, k0=1, k1=floor(K/2), thre=FALSE))
      }
      else {
        f=flibs[[K-1]]
        return(exp(f(log(q))))
      }
    }
    if (HCtype=="MHC_flibs") {
      if (K>1000 | K<20) {
        return(hybrid(q=q,K=K, k0=1, k1=K, thre=TRUE,N=10^6))
      }
      else {
        f=flibs[[K-19]]
        return(exp(f(log(q))))
      }
    }
    if (HCtype=="TMHC_flibs") {
      if (K>1000 | K<30) {
        return(hybrid(q=q,K=K, k0=1, k1=floor(K/2), thre=TRUE,N=10^6))
      }
      else {
        f=flibs[[K-29]]
        return(exp(f(log(q))))
      }
    }
  }
  else {
    stop("Please input the pre-calcuated data set of spline functions: flibs=HC_flibs, THC_flibs, MHC_flibs or TMHC_flibs.")
  }
}
