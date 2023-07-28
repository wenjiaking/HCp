
#' Calculate the p-value of Higher Criticism statitic given an observed value
#'
#' @param q  quantile or observed value of HC statistic, should be a scalar
#' @param K  dimension of HC, i.e. the total number of studies
#' @param k0  search range starts from the k0th smallest p-value, default value is 1.
#' @param k1  search range ends at the k1th smallest p-value, default value is K.
#' @param thre  whether using thresholding supremum domain, Boolean variable equals False by default.
#'
#' @return The right-tail probability of the null distribution of HC statistic at the given quantile, i.e. p-value.
#' @export
#'
#' @examples
#' pset=runif(50)
#' hcstat <-HCstat(pset,k0=1)
#' mst(q=hcstat, K=length(pset), k0=1)

mst <- function(q, K, k0=1, k1=NA, thre=FALSE){
  k0=floor(k0)
  if(is.na(k1)){
    k1 = K
  }else{
    k1 = floor(k1)
  }
  if (thre & K>120) {warning("Only reliable for K<120 if thre=TRUE")}
  return(1-pphi.mod(q=q,K=K, k0=k0, k1=k1, s=2, MHC = thre))
}


