#' Calculate p-value of HC by ultra-fast interpolation (UFI) method
#'
#' @description This method is recommended when K is no larger than 2000, for HC of full supremum domain (FHC), HC of half truncated supremum domain (THC), HC of thresholding supremum domain (MHC) or HC of thresholding and half truncated supremum domain (TMHC)
#' @param flibs choose the spline function data set of one specific variation of HC form HC_flibs, THC_flibs, MHC_flibs or TMHC_flibs, and the default is HC_flibs.
#' @param K a vector of HC dimensions i.e. the total number of studies
#' @param q a vector of quantiles or observed values of HC statistic
#'
#' @return a vector of p-values corresponding to the input HC statistics
#' @export
#'
#' @examples
#' q=sapply(10:15,function(i) {
#' pset=runif(i)
#' return(HCstat(pset,k0=1,thre=FALSE))
#' })
#' ufi.p(flibs=HC_flibs,K=10:15,q=q)
#'
ufi.p=function(flibs=HC_flibs,K,q) {
  HCtype=deparse(substitute(flibs))
  if (HCtype %in% c("HC_flibs","THC_flibs","MHC_flibs","TMHC_flibs")) {
    if (HCtype=="HC_flibs") {
      flist=flibs[K-1]
    }
    if (HCtype=="THC_flibs") {
      flist=flibs[K-1]
    }
    if (HCtype=="MHC_flibs") {
      flist=flibs[K-19]
    }
    if (HCtype=="TMHC_flibs") {
      flist=flibs[K-29]
    }
  }
  else {
    stop("Please input the pre-calcuated data set of spline functions: HC_flibs, THC_flibs, MHC_flibs or TMHC_flibs.")
  }

  ps=sapply(1:length(q), function(i){
    f=flist[[i]]
    return(exp(f(log(q[i]))))
  } )
  return(ps)
}
