#' Hybrid method of computing p-value for HC of arbitrary variation
#' @description The function is recommended for HC of arbitrary variation except for FHC, THC, MHC, TMHC
#' @param q quantile or observed value of HC statistic, should be a scalar
#' @param K dimension of HC, i.e. the total number of studies
#' @param k0 search range starts from the k0th smallest p-value, default value is 1.
#' @param k1 search range ends at the k1th smallest p-value, default value is K.
#' @param thre whether using thresholding supremum domain, Boolean variable equals False by default.
#' @param N Samping size when sampling method is recommended, default value is 10^6.
#'
#' @return The right-tail probability of the null distribution of HC statistic at the given quantile, i.e. p-value.
#' @export
#'
#' @examples
#' pset=runif(200)
#' hcstat <-HCstat(pset,k0=1,k1=80, thre=FALSE)
#' hybrid(q=hcstat, K=200, k0=1,k1=80, thre=FALSE,N=10^6)

hybrid=function(q,K, k0=1, k1=NA, thre=FALSE,N=10^6) {
  if (!thre) {
    return(mst(q=q, K=K, k0=k0, k1=k1, thre=thre))
  }
  else {
    if (K<=100) {
      return(mst(q=q, K=K, k0=k0, k1=k1, thre=thre))
    }
    else {
      pest=ufi.p(flibs=MHC_flibs,K=K,q=q)
      if (pest<10^(-3)) {
        return(LiAppro_HC(q=q,K=K, k0=k0, k1=k1, thre=thre))
      }
      else {
        N1=ifelse(N%%(10^4)==N,N,10^4)
        B=ifelse((N/(10^4)) > 1,ceiling((N/(10^4))),1)
        return(MC.est(K=K,N1=N1,B=B,k0=k0,k1=k1,thre=thre,q=q))
      }
    }
  }
}
