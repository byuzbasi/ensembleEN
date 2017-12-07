
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/esmucler/ensembleEN.svg?branch=master)](https://travis-ci.org/esmucler/ensembleEN)

ensembleEN
==========

This package provides functions for computing ensembles of regularized linear regression estimators as defined in Christidis, Lakshmanan, Smucler and Zamar (2017).

------------------------------------------------------------------------

### Installation

You can install the **development** version from [GitHub](https://github.com/esmucler/ensembleEN)

``` r
library(devtools)
devtools::install_github("esmucler/ensembleEN")
```

### Usage

``` r
# A small example
library(MASS)
library(ensembleEN)
set.seed(1)
beta <- c(rep(5, 5), rep(0, 45))
Sigma <- matrix(0.5, 50, 50)
diag(Sigma) <- 1
x <- mvrnorm(50, mu = rep(0, 50), Sigma = Sigma)
y <- x %*% beta + rnorm(50)
fit <- cv.ensembleEN(x, y, num_groups=10) # Use 10 models
coefs <- predict(fit, type="coefficients")
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
