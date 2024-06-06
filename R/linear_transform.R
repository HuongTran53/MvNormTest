# -- Filename: maximum_skewness.R --
#' @name skewness_kurtosis
#' @details Sample kurtosis is
#' \deqn{
#' \hat{\kappa}_4 =
#' \dfrac{1}{n-1} \sum_{i = 1}^n \left(\dfrac{X_i - \bar{X}}{S}\right)^4.
#' }
#' @return \code{kurtosis} returns sample kurtosis.
#' @export
#' @examples
#' set.seed(123)
#' y <- rnorm(100)
#' kurtosis(y)
kurtosis <- function(x){
  nx <- length(x)
  sdx <- stats::sd(x)
  y <- (x - mean(x))/sdx
  k <- sum(y^4)
  k <- 1/(nx-1) * k
  return(k - 3)
}
#' @rdname skewness_kurtosis
#' @title Sample skewness and Sample Kurtosis.
#' @param x univariate data sample
#' @details
#' Sample skewness is
#' \deqn{
#' \hat{\kappa}_3 =
#' \dfrac{1}{n-1} \sum_{i = 1}^n \left(\dfrac{X_i - \bar{X}}{S}\right)^3.
#' }
#' @return \code{skewness} returns sample skewness.
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' skewness(x)
skewness <- function(x) {
  nx <- length(x)
  sdx <- stats::sd(x)
  y <- (x - mean(x))/sdx
  k <- sum(y^3)
  k <- 1/(nx - 1) * k
  return(k)
}

###############################
######## max k^2 ##############
# -- Filename: maximum_skewness.R --
#' @title Best Linear Transformations
#' @name linear_transform
#' @description
#' The algorithm uses gradient descent algorithm to obtain the maximum of the
#' square of sample skewness, of the kurtosis or of their average under any
#' univariate linear transformation of the multivariate data.
# Departure from 0 of these values is an indication of non-normality.
#'
#' @param x multivariate data matrix.
#' @param l0 starting point for projection algorithm,
#' default is \code{rep(1, ncol(x))}.
#' @param method character strings,
#'  one of \code{c("skewness", "kurtosis", "both")}.
#' @param epsilon bounds on error of optimal solution, default is \code{1e-10}.
#' @param iter number of iteration of projection algorithm,
#' default is \code{5000}.
#' @param stepsize gradient descent stepsize, default is \code{.001}.
# @details
# The algorithm  looks for vector \eqn{l \in \mathbb{R}^p} that maximizes the
# univariate skewness, kurtosis, average of these both under any possible linear
# transformation.
#' @return
#' \itemize{
#' \item{\code{max_result}: The maximum value after linear transformation.}
#' \item{\code{x_uni}: Univariate data after transformation.}
#' \item{\code{vector_k}: Vector of the "best" linear transformation.}
#' \item{\code{error}: Error of projection algorithm.}
#' \item{\code{iteration}: Number of iteration. }
#' }
#' @seealso \code{\link[=skewness]{skewness()}},
#' \code{\link[=kurtosis]{kurtosis()}}
#' @export
#' @examples
#' set.seed(1)
#' x <- MASS::mvrnorm(100, mu = rep(0, 2), diag(2))
#' linear_transform(x, method = "skewness")$max_result
#' linear_transform(x, method = "kurtosis")$max_result
#' linear_transform(x, method = "both")$max_result

######## All method ##########
linear_transform <- function(x, l0  = rep(1, ncol(x)), method = "both",
                             epsilon = 1e-10, iter = 5e3, stepsize = 1e-3){
  nx <- nrow(x)
  p <- ncol(x)

  if (epsilon <= 0){
    stop("error should be positive number")
  }
  # if (sum(t(l0) %*% l0) != 1){
  #   stop(sQuote("l0"), "must has norm of 1")
  # }
  if (length(l0) != p){
    stop(sQuote("l0"), "wrong dimension")
  }
  if (iter < 0 || iter != round(iter) ){
    stop(sQuote("iter"), "must be positive integer")
  }
  if (stepsize < 0){
    stop(sQuote("stepsize"), "must be positive")
  }
  #### Data standardization ##################
  l0 <- l0 / sqrt(sum(t(l0) %*% l0))
  S <- stats::cov(x)
  m <- matrix(rep(colMeans(x), nx), byrow = T, ncol = p)
  chol.S <- chol(S)  # 50% faster
  A <-  backsolve(r = chol.S, x = diag(p))
  x <- (x - m) %*% A    # check var(x) is an identity matrix
  #### Choosing method #################
  if (method == "skewness"){
    fun <- function(y) skewness(y)^2
    gradG <- function(l, x){
      xl <- (x %*% l)
      v <- 3 * rowSums(sapply(1:nx, function(i) xl[i]^2 * x[i, ]))
      return(v)
    }
  } else{
    if (method == "kurtosis"){
      fun <- function(y) kurtosis(y)^2
      gradG <- function(l, x){
        xl <- (x %*% l)
        v <- -2 * (mean(xl^4) - 3) *4*rowMeans(
          sapply(1:nx, function(i) xl[i]^3 * x[i, ])
        )
        return(v)
      }
    } else {
      fun <- function(y) .5*(skewness(y)^2 + kurtosis(y)^2)
      gradG <- function(l, x){
        xl <- x %*% l
        v1 <- (mean(xl^4) - 3) * 4 * rowMeans(
          sapply(1:nx, function(i) xl[i]^3 * x[i, ])
        )
        v2 <- (mean(xl^3)) * 3 * rowMeans(
          sapply(1:nx, function(i) xl[i]^2 * x[i, ])
        )
        v <-  -(v1+v2)
        return(v)
      }
    }
  }

  lk <- l0; err <- 1; i <- 1;
  while (i <= iter && err >= epsilon ) {
    gradG.lk <- gradG(l = lk, x = x)
    gamma <- lk - stepsize * gradG.lk
    lkp1 <- gamma / sqrt(sum(gamma^2))
    # err <- sqrt(t(lk - lkp1) %*% (lk - lkp1))
    err <- sqrt(sum((lk - lkp1)^2))
    i <- i + 1
    lk <- lkp1
  }

  x.new <- x %*% lk
  re <- fun(x.new) / p #scale by p
  return(list(max_result = re, x_new = x.new, l = lk,
              error = err, iteration = i))


}

