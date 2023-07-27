
#' Calculate HC statitics given a set of individual p-values
#'
#' @param pset  a vector of p-values from individual studies
#' @param k0  search range starts from the k0th smallest p-value, default value is 1.
#' @param k1  search range ends at the k1th smallest p-value, default value is K.
#' @param thre  whether using thresholding supremum domain, Boolean variable equals False by default.
#'
#' @return HC statitic
#' @export
#'
#' @examples
#' pset=runif(10)
#' hcstat <-HCstat(pset,k0=1)
#'
HCstat <- function (pset, k0=1, k1=NA,thre=F) {
  K=length(pset)
  pset=sort(pset,decreasing = F)
  if(is.na(k1)){
    k1 = K
  }else{
    k1 = floor(k1)
  }
  tailprob=pset[k0:k1]
  if (thre) {
    tailprob[which(tailprob<1/K)]=k1/K
  }
  GHC.mesh=(k0:k1 - K*tailprob)/sqrt(K*tailprob*(1-tailprob))
  return(stat=max(GHC.mesh,na.rm=TRUE))
}
