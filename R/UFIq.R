#' Calculate threshold of HC by ultra-fast interpolation (UFI) method
#'
#' @description This method is recommended when K is no larger than 2000, for HC of full supremum domain (FHC), HC of half truncated supremum domain (THC), HC of thresholding supremum domain (MHC) or HC of thresholding and half truncated supremum domain (TMHC)
#' @param flibs_q choose the spline function data set of one specific variation of HC form HC_flibs_q, THC_flibs_q, MHC_flibs_q or TMHC_flibs_q, and the default is HC_flibs_q.
#' @param K a vector of HC dimensions i.e. the total number of studies
#' @param p a vector of signifiance level for each of the test
#' @return a vector of thresholds of HC tests corresponding to the input significance levels
#' @export
#'
#' @examples
#' ufi.q(flibs_q=HC_flibs_q,K=10:15,p=10^(seq(-3,-8)))
#'
ufi.q=function(flibs_q=HC_flibs_q,K,p) {
  HCtype=deparse(substitute(flibs_q))
  if (HCtype %in% c("HC_flibs_q","THC_flibs_q","MHC_flibs_q","TMHC_flibs_q")) {
    if (HCtype=="HC_flibs_q") {
      flist_q=flibs_q[K-1]
    }
    if (HCtype=="THC_flibs_q") {
      flist_q=flibs_q[K-1]
    }
    if (HCtype=="MHC_flibs_q") {
      flist_q=flibs_q[K-19]
    }
    if (HCtype=="TMHC_flibs_q") {
      flist_q=flibs_q[K-29]
    }
  }
  else {
    stop("Please input the pre-calcuated data set of spline functions: flibs_q=HC_flibs_q, THC_flibs_q, MHC_flibs_q or TMHC_flibs_q.")
  }

  qs=sapply(1:length(p), function(i){
    f=flist_q[[i]]
    return(exp(f(log(p[i]))))
  } )
  return(qs)
}
