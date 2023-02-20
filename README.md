# HCp R package

The R package "HCp" including functions of different methods for p-value computation of Higher Criticism test.

## Install the pacakge

`library(devtools)` \
`install_github("wenjiaking/HCp@master")` \
`library(HCp)`

## Brief Description

The package provide functions for implementation of different methods for p-value computation of Higher Criticism test

* Barnett and Lin's method: `HCp::XH()`
* Li and Siegmund's method: `HCp::LiAppro_HC()` 
* Modified SetTest method: `HCp::mst()` 
* Cross-entropy-based importance sampling method: `HCp::CE.mixed.grad()` and `HCp::CE.mixed.prop()` 
* Untra-fast interpolation method: `HCp::ufi.p()` to calcaulte p-value and `HCp::ufi.q()` to find quantile 
* Hybrid computing strategy: `HCp::hybrid()` and `HCp::hybridSpec()`

## Reference
* Barnett, R. Mukherjee, and X. Lin. The generalized higher criticism for testing SNP-Set effects in genetic association studies. Journal of the American Statistical Association, 112, 06 2016. doi: 10.1080/01621459.2016.1192039.
* J. Li and D. Siegmund. Higher criticism: p-values and criticism. The Annals of Statistics, 43(3): 1323 – 1350, 2015. doi: 10.1214/15-AOS1312.
* H. Zhang, J. Jin, and Z. Wu. Distributions and power of optimal signal-detection statistics in finite case. IEEE Transactions on Signal Processing, 68:1021–1033, 2020. doi: 10.1109/TSP. 2020.2967179.
