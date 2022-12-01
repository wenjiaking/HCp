#' Calculate p-value of HC by Li-Siegmund approximation
#'
#' @param q  quantile or observed value of HC statistic, should be a scalar
#' @param K  dimension of HC, i.e. the total number of studies
#' @param k0  search range starts from the k0th smallest p-value, default value is 1.
#' @param k1  search range ends at the k1th smallest p-value, default value is K.
#' @param thre  whether using thresholding supremum domain, Boolean variable equals TRUE by default.
#'
#' @return The right-tail probability of the null distribution of HC statistic at the given quantile, i.e. p-value.
#' @export
#'
#' @examples
#' pset=runif(10)
#' hcstat <-HCstat(pset,k0=1,thre=TRUE)
#' LiAppro_HC(q=hcstat, K=10)
LiAppro_HC=function(q,K, k0=1, k1=NA, thre=TRUE) {
  k0=floor(k0)
  if(is.na(k1)){
    k1 = K
  }else{
    k1 = floor(k1)
  }
  k_range=k0:k1
  k_Cs=Cx_HC(x=k_range/K,eps=q/sqrt(K))
  k_roots=k_Cs$C
  k_deltaC=k_Cs$deltaC
  k_coefs=1-((K-k_range+1)*k_deltaC)/(K*(1-k_roots))
  if (thre) {
    k_Bins=stats::dbinom(k_range,K,sapply(k_roots,function(r) max(r,1/K)))-stats::dbinom(k_range,K,1/K)*sapply(K*k_roots,function(r) max(r,1))
  }
  else {
    k_Bins=stats::dbinom(k_range,K,k_roots)
  }
  return(sum(k_coefs*k_Bins))
}
