# -- Filename: Multi_Uni.R -----#
#' @title Transformation to Independent Univariate Sample.
#' @name Independent_transformation
#' @description Leave-one-out method gives approximately independent sample of
#' standard multivariate normal distribution,
#' which then produces sample of standard univariate normal distribution.
#' @param x multivariate data matrix
#' @return Data frame containing univariate data and
#' the index from multivariate data.
#' @details
#' \deqn{S_{-k}^{-1/2} (X_k - \bar{X}_{-k}) , k = 1,... n}
#' are approximately independent sample of \eqn{N_p(0, I)},
#' which then gives \eqn{n \times p}
#' univariate sample of \eqn{N(0, 1)}.
# @note The transformation gives number of data points as a multiple of
# dimension \eqn{p}, density of data points is high and
# spline method may produce ill-posed matrices. Merging with \code{x.dist}
# parameter is recommended.
#' @export
#' @examples
#' set.seed(1)
#' x <- MASS::mvrnorm(100, mu = rep(0, 5), diag(5))
#' df <- Multi.to.Uni(x)
#' qqnorm(df$x.new); abline(0, 1)
Multi.to.Uni <- function(x){
  nx <- nrow(x)
  p <- ncol(x)
  xbar <- colMeans(x)
  S <- stats::cov(x)
  inv.S <- solve(S)
  temp <- ((nx - 1)**2)/nx
  ##############
  inv.Si <- function(i){
    #function calculate inverse S(-k)
    top <- inv.S %*% (x[i,] - xbar) %*% t((x[i, ] - xbar)) %*% inv.S
    bottom <- temp - t((x[i, ] - xbar)) %*% inv.S %*% (x[i, ] - xbar)
    inv.temp <- inv.S + top/sum(bottom)
    inv.temp <- (nx - 2) / (nx - 1) *inv.temp
    return(inv.temp)
  }
  x.new <- c()
  ##############
  for (i in 1: nx){
    xbari <- colMeans(x[-i, ])
    inv.S.not.i <- inv.Si(i)
    a <- chol(inv.S.not.i)
    xi <- a %*% (x[i,] - xbari)
    x.new <- c(x.new, xi)
  }
  ind <- kronecker(c(1:nx), rep(1, p))
  df <- data.frame(x.new, ind)
  df <- df[order(df$x.new),]
  return(newdata = df)
}





