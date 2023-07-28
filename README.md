# HCp R package

The R package "HCp" including functions of different methods for p-value computation of Higher Criticism test. The source code is under "master" branch.

## Install the pacakge

`library(devtools)` \
`install_github("wenjiaking/HCp@master")` \
`library(HCp)`

## Brief Description

The package provide functions for implementation of different methods for p-value computation of Higher Criticism test

* Barnett and Lin's method: `HCp::XH()`
 
Example: 

```
pset=runif(50)
hcstat <-HCstat(pset,k0=1,k1=length(pset),thre=F)
XH(q=hcstat,K=50)

```
* Modified SetTest method: `HCp::mst()` 
Example: 

```
mhcstat <-HCstat(pset,k0=1,k1=length(pset),thre=T)
mst(q=hcstat,K=50,k0=1,k1=50,thre=F)
mst(q=mhcstat,K=50,k0=1,k1=50,thre=T)

```
* Cross-entropy-based importance sampling method: `HCp::CE.mixed.grad()` and `HCp::CE.mixed.prop()` 
Example: 

```
CE.mixed.grad(q=hcstat,K=50,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=50,thre=F,idx=1,theta=1)
CE.mixed.prop(q=mhcstat,K=50,ro=0.01,N0=10^4, N1=10^4,B=1,k0=1,k1=50,thre=T,prop=0.2,theta=1)

```
* Li and Siegmund's method: `HCp::LiAppro_HC()` (Li-Siegmund approximation method is only reliable for small p-value computing (< 0.01).)
Example: 

```
pset=rbeta(100,1,2) #LiAppro_HC only reliable for small p-value computing
mhcstat <-HCstat(pset,k0=1,k1=length(pset),thre=T)
LiAppro_HC(q=mhcstat,K=length(pset))
mst(q=mhcstat,K=length(pset),k0=1,k1=length(pset),thre=T)

```
* Untra-fast interpolation method: `HCp::ufi.p()` to calcaulte p-value and `HCp::ufi.q()` to find quantile (This ultra-fast method is only for K no larger than 2000, otherwise please use function hybridSpec().)
Example: 

```
q=sapply(10:15,function(i) {
   pset=runif(i)
   return(HCstat(pset,k0=1,thre=FALSE))
   })
ufi.p(flibs=HC_flibs,K=10:15,q=q)
ufi.q(flibs_q=HC_flibs_q,K=10:15,p=10^(seq(-3,-8)))

```
* Hybrid computing strategy: `HCp::hybrid()` and `HCp::hybridSpec()`

`HCp::hybrid()`: Hybrid method of computing p-value for HC of arbitrary variation
`HCp::hybridSpec()`: Hybrid method of computing p-value for HC of four popular variations

Example: 

```
pset=runif(200)
hybrid(q=HCstat(pset,k0=1,k1=80, thre=FALSE), K=200, k0=1,k1=80, thre=FALSE,N=10^6)

pset=runif(2001)
hybridSpec(q=HCstat(pset,k0=1,k1=2001, thre=FALSE), K=2001, flibs=HC_flibs,N=10^6)

```
## Reference
* Barnett, R. Mukherjee, and X. Lin. The generalized higher criticism for testing SNP-Set effects in genetic association studies. Journal of the American Statistical Association, 112, 06 2016. doi: 10.1080/01621459.2016.1192039.
* J. Li and D. Siegmund. Higher criticism: p-values and criticism. The Annals of Statistics, 43(3): 1323 – 1350, 2015. doi: 10.1214/15-AOS1312.
* H. Zhang, J. Jin, and Z. Wu. Distributions and power of optimal signal-detection statistics in finite case. IEEE Transactions on Signal Processing, 68:1021–1033, 2020. doi: 10.1109/TSP. 2020.2967179.
