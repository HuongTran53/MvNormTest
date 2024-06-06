# -- Filename: Multi_Uni.R -----#
#' @title Transformation to Independent Univariate Sample
#' @name Independent_transformation
#' @description Leave-one-out method gives approximately independent sample of
#' standard multivariate normal distribution,
#' which then produces sample of standard univariate normal distribution.
#' @param x multivariate data matrix
#' @return Data frame contains univariate data and
#' the index from multivariate data.
#' @details
#' Let \eqn{\bar{X}_{-k} } and \eqn{S_{-k}} are the sample mean sample variance
#' covariance matrix obtained by using all but \eqn{k^{th}} data point. Then
#' \eqn{S_{-k}^{-1/2} (X_k - \bar{X}_{-k}) , k = 1,... n} are approximately
#' independently distributed as \eqn{N_p(0, I)}. Thus all \eqn{n \times p}
#' entries in the data matrix so constructed can be treated as
#' univariate samples of size \eqn{n \time p} from \eqn{N(0, 1)}.
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





