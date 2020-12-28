
<!-- README.md is generated from README.Rmd. Please edit that file -->

# acmvup

<!-- badges: start -->
[![License](https://badgen.net/github/license/PeppeSaccardi/acmvup)](LICENSE)
<!-- badges: end -->

The goal of acmvup is to provide some useful functions that are
implemented following the additive covariance modeling via unconstrained
parametrization idea. You will find more detail in my master thesis.

## Installation

You can install the released version of acmvup from
[Github](https://github.com/PeppeSaccardi/acmvup) with:

``` r
library("devtools")
install_github("PeppeSaccardi/acmvup")
```

## Example

This is a basic example which shows you how to compute the
log-likelihood function value, gradient and hessian

``` r
library(acmvup)
set.seed(1234)
# Let generate synthetic data, using p = 3
data <- matrix(rnorm(300), ncol=3)

# Then we need the lists containing the covariates for all the
# unconstrained parameters that we want to model
lista_d <- list()
lista_phi <- list()

lista_d[[1]] <- matrix(c(rep(1,100),rnorm(100)), byrow = FALSE, ncol=2)
lista_d[[2]] <- matrix(c(rep(1,100),rnorm(200)), byrow = FALSE, ncol=3)
lista_d[[3]] <- matrix(rep(1,100), byrow = FALSE, ncol=1)

lista_phi[[1]] <- matrix(c(rep(1,100),rnorm(200)),byrow = FALSE, ncol=3)
lista_phi[[2]] <- matrix(rep(1,100),ncol=1)
lista_phi[[3]] <- matrix(c(rep(1,100),rnorm(100)),byrow = FALSE, ncol=2)
```

Now we are able to use the functions implemented in the package to
compute the log-likelihood function value, gradient and hessian, for
instance in a suitable random point

``` r
par <- rnorm(12)

ll <- optimal_loglik(par,data,lista_phi,lista_d)
gll <- optimal_grad(par,data,lista_phi,lista_d)
hll <- optimal_hessian(par,data,lista_phi,lista_d)
```

Therefore for our example

``` r
print(paste0("The log-likelihood value is ",ll))
#> [1] "The log-likelihood value is -1082.81516785296"
print("The gradient vector is ")
#> [1] "The gradient vector is "
gll
#>  [1]  380.3925537   51.8455975  167.4419360 -215.2302577   71.1800255
#>  [6] -206.7230241  -17.9633814    0.5270398  725.7886307  503.6839641
#> [11] -266.7740084  262.2474045
print("The hessian matrix is ")
#> [1] "The hessian matrix is "
hll
#>              [,1]       [,2]        [,3]        [,4]        [,5]       [,6]
#>  [1,] -129.213953  -10.32009    3.795308    0.000000    0.000000    0.00000
#>  [2,]  -10.320090  -97.83978   28.219271    0.000000    0.000000    0.00000
#>  [3,]    3.795308   28.21927 -127.188731    0.000000    0.000000    0.00000
#>  [4,]    0.000000    0.00000    0.000000 -309.969137    9.849759  -13.36294
#>  [5,]    0.000000    0.00000    0.000000    9.849759 -320.007963  -51.85966
#>  [6,]    0.000000    0.00000    0.000000  -13.362941  -51.859657 -246.46343
#>  [7,]    0.000000    0.00000    0.000000    0.000000    0.000000    0.00000
#>  [8,]    0.000000    0.00000    0.000000    0.000000    0.000000    0.00000
#>  [9,] -380.392554  -51.84560 -167.441936    0.000000    0.000000    0.00000
#> [10,] -223.845450 -106.43793 -107.343399    0.000000    0.000000    0.00000
#> [11,]  112.712387   47.58699   16.062310    0.000000    0.000000    0.00000
#> [12,]    0.000000    0.00000    0.000000  215.230258  -71.180026  206.72302
#>             [,7]       [,8]      [,9]      [,10]      [,11]      [,12]
#>  [1,]   0.000000   0.000000 -380.3926  -223.8454  112.71239    0.00000
#>  [2,]   0.000000   0.000000  -51.8456  -106.4379   47.58699    0.00000
#>  [3,]   0.000000   0.000000 -167.4419  -107.3434   16.06231    0.00000
#>  [4,]   0.000000   0.000000    0.0000     0.0000    0.00000  215.23026
#>  [5,]   0.000000   0.000000    0.0000     0.0000    0.00000  -71.18003
#>  [6,]   0.000000   0.000000    0.0000     0.0000    0.00000  206.72302
#>  [7,] -32.036619  -0.121783    0.0000     0.0000    0.00000    0.00000
#>  [8,]  -0.121783 -36.842355    0.0000     0.0000    0.00000    0.00000
#>  [9,]   0.000000   0.000000 -775.7886  -502.5947  273.61786    0.00000
#> [10,]   0.000000   0.000000 -502.5947 -1196.5624  352.93666    0.00000
#> [11,]   0.000000   0.000000  273.6179   352.9367 -848.71568    0.00000
#> [12,]   0.000000   0.000000    0.0000     0.0000    0.00000 -312.24740
```
