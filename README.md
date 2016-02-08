# Source for "Hierarchical expectation propagation for Bayesian aggregation of average data"

Date: February 1st, 2016

This repository includes the R and Stan source files used to generate
the main outputs in our paper on Bayesian aggregation of average data.

The code is distributed under the BSD-3 license, see LICENSE.

Please cite our work when you use these codes or apply our approach:

_S. Weber, A. Gelman, B. Carpenter, D. Lee, M. Betancourt, A. Vehtari,
A. Racine (2016), "Hierarchical expectation propagation for Bayesian
aggregation of average data", arXiv:1602.02055_

To use the code, we recommend to start with one of the examples below
as template and modify it to your needs.

Updates will be posted to https://github.com/wds15/hep

# Contents

- `extrap_linear.R` Linear example
- `extrap_logistic.R` Logistic example
- `extrap_pkpd.R` Non-linear PK/PD example
- `extrap_ep_routine.R` Main EP routine
- `utils.R` Various utility functions

# News

2016-02-01: Initial release.

# Software used

R session info:

```
R version 3.1.3 (2015-03-09)
Platform: x86_64-unknown-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server release 6.5 (Santiago)

locale:
[1] C

attached base packages:
[1] methods   stats     graphics  grDevices utils     datasets  base     

other attached packages:
 [1] assertthat_0.1 abind_1.4-0    loo_0.1.2      arm_1.8-4      lme4_1.1-7    
 [6] Matrix_1.1-5   MASS_7.3-39    rstan_2.7.0    inline_0.3.13  Rcpp_0.11.6   

loaded via a namespace (and not attached):
 [1] BH_1.58.0-1         RcppEigen_0.3.2.3.0 StanHeaders_2.7.0  
 [4] coda_0.17-1         codetools_0.2-10    grid_3.1.3         
 [7] lattice_0.20-31     matrixStats_0.14.2  minqa_1.2.4        
[10] nlme_3.1-120        nloptr_1.0.4        parallel_3.1.3     
[13] splines_3.1.3       stats4_3.1.3        tools_3.1.3        
```
