
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MvNormTest

<!-- badges: start -->
<!-- badges: end -->

MvNormTest provides graphical techniques to find evidence of
non-normality of a multivariate dataset.

## Installation

You can install the development version of MvNormTest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("HuongTran53/MvNormTest")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MvNormTest)
set.seed(123)
x <- MASS::mvrnorm(1000, rep(0, 5), diag(5))
```

``` r
d3hCGF_plot(x); title("Using third derivatives of CGF")
#> [1] "accept"
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
d4hCGF_plot(x); title("Using fourth derivatives of CGF")
#> [1] "accept"
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
df <- Multi.to.Uni(x)
qqnorm(df$x.new, main = "Transfromation to nearly independent unvariate sample, Q-Q plot"); abline(0, 1)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r
# Maximum skewness under linear transformation
linear_transform(x, method = "skewness")$max_result
#> [1] 0.01160368
```
