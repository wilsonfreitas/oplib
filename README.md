
<!-- README.md is generated from README.Rmd. Please edit that file -->

# op

<!-- badges: start -->
<!-- badges: end -->

Option Volatility and Pricing Models.

## Installation

You can install the development version of op like so:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("wilsonfreitas/rb3")
```

## Example

Pricing options with Black-Scholes-Merton Models (BSM Model).

``` r
library(op)

bsmprice(c("call", "put"), c(50, 49, 48, 47), 45, 0.25, 0.13, 0.01, 0.2)
#> [1] 6.5031266 0.2854164 4.7405597 0.6079959
bsmprice("call", 50, 45, 0.25, 0.13, 0.01, seq(0.1, 0.5, 0.1))
#> [1] 6.316558 6.503127 7.012037 7.695351 8.464750
bsmprice("put", 20, 30, 0:4 / 2, 0.15, 0, 0.25)
#> [1] 10.000000  7.882056  6.272936  5.098353  4.203993
```
