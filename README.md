
# sgs <a><img src='man/figures/logo.avif' align="right" height="139" /></a>

[![CRAN status](https://www.r-pkg.org/badges/version/sgs)](https://CRAN.R-project.org/package=sgs)
[![CRAN downloads this month](https://cranlogs.r-pkg.org/badges/sgs)](https://cran.r-project.org/package=sgs)

Implementation of Sparse-group SLOPE (SGS), a sparse-group penalisation regression approach. SGS performs adaptive bi-level selection, controlling the FDR under orthogonal designs. Linear (Gaussian) and logistic (Binomial) regression are supported, both with dense and sparse matrix implementations. Cross-validation functionality is also supported. SGS is implemented using adaptive three operator splitting (ATOS) and the package also contains a general implementation of ATOS.

A detailed description of SGS can be found in [F. Feser, M. Evangelou (2023) "Sparse-group SLOPE: adaptive bi-level selection with FDR-control"](https://arxiv.org/abs/2305.09467).


## Installation

You can install the current stable release from [CRAN](https://cran.r-project.org/) with
``` r
install.packages("sgs")
```
Your R configuration must allow for a working Rcpp. To install a develop the development version from GitHub run
``` r
library(devtools)
install_github("ff1201/sgs")
```

## Example

The code for fitting a basic SGS model is:

``` r
library(sgs)

model = fit_sgs(X = X, y = y, groups = groups, vFDR=0.1, gFDR=0.1)
```

where `X` is the input matrix, `y` the response vector, `groups` a vector containing indices for the groups of the predictors, and `vFDR` and `gFDR` are the the target variable/group false discovery rates.

[A more extensive example can be found here](https://github.com/ff1201/sgs/blob/master/vignettes/reproducible_example.Rmd).
