# --- Filename: d_hCGF.R ---#
#' @title Calculation of derivatives of empirical
#' cumulant generating function (CGF).
#' @name  dCGF
#' @description
#' Get the third/fortth derivatives of sample CGF at a given point.
#' @details
#' Estimator of standardized cumulant function is
#' \deqn{ \log\hat{M}_X(t) = \log \left(\dfrac{1}{n}
#' \sum_{i = 1}^n \exp(t'S^{\frac{-1}{2}}(X_i - \bar{X})) \right)
#' }
#' and its \deqn{k^{th}} order derivatives is defined as
#' \deqn{
#' T_k(t) = \dfrac{\partial^k}{
#' \partial t_{j_1}t_{j_2} \dots t_{j_k}} \log(\hat{M}_X(t)), t \in \mathbb{R}^p
#'  }
#' where \eqn{t_{j_1}t_{j_2} \dots t_{j_k}} are the corresponding components
#' of vector \eqn{t \in \mathbb{R}^p}.
# The number of distinct third derivatives is:
# \deqn{
# l_T = p + 2 \times \begin{pmatrix}
# p\\2
# \end{pmatrix} + \begin{pmatrix}
# p \\ 3
# \end{pmatrix}
# }
#' @param myt,t numeric vector of length \code{p}.
#' @param x data matrix.
#' @return \code{d3hCGF} returns the sequence of third derivatives of
#' empirical CGF, ordered by index of \eqn{j_1 \leq j_2 \leq j_3 \leq p}.
#' @export
#' @examples
#' p <- 3
#' # Number of distinct derivatives
#' l_dhCGF(p)
#' set.seed(1)
#' x <- MASS::mvrnorm(100, rep(0, p), diag(p))
#' myt <- rep(.2, p)
#' d3hCGF(myt = myt, x = x)
#' d4hCGF(myt = myt, x = x)
d3hCGF <- function(myt, x){
  nx <- nrow(x)
  p <- ncol(x)
  #standardize data: not using Cholesky decomposition
  S <- stats::cov(x)
  xbar <- colMeans(x)
  m <- matrix(rep(xbar, nx), byrow = T, ncol = p)
  x <- x - m
  chol.S <- chol(S)  # 50% faster
  A <-  backsolve(r = chol.S, x = diag(p))
  x <- x %*% A    # check var(x) is an identity matrix
  #############################################
  pi <- c(exp(x %*% myt))
  pbar <- 1/sum(exp(x %*% myt))
  pi <- pbar*pi
  m1  <- colSums(x*pi)
  m2  <-  t(x*pi)%*%x
  t3 <-  lapply( 1:p, function(k){
    temp2 <- m1 %*% t(m2[k, ]);
    t(x*pi*x[,k])%*%x -  m1[k] *m2 -temp2 - t(temp2) +2* (m1[k]) * (m1) %*% t(m1)
    # t(x*pi*x[,k])%*%x -  m1[k] *m2 -temp2 - t(temp2) +  (m1[k]) * (m1) %*% t(m1)
  })
  v3 <- c()
  for (i in 1:p){
    t3i <- t3[[i]][i:p, i:p]
    t3i <- t3i[upper.tri(t3i, diag= T)]
    v3 <- c(v3, c(t3i))
  }
  return(v3)
}
#####################################################
#####################################################
#' @rdname dCGF
#' @return \code{d4hCGF} returns the sequence of fourth derivatives of empirical
#'  CGF ordered by index of \eqn{j_1 \leq j_2 \leq j_3 \leq j_4 \leq p}.
#' @export
d4hCGF <- function(myt, x){
  nx <- nrow(x)
  p <- ncol(x)
  #standardize data: not using Cholesky decomposition
  S <- stats::cov(x)
  xbar <- colMeans(x)
  m <- matrix(rep(xbar, nx), byrow = T, ncol = p)
  x <- x - m
  if (kappa(S) > 1e13){
    warning("Non-invertable covariance matrix")
  } else {
    chol.S <- chol(S)  # 50% faster
    A <-  backsolve(r = chol.S, x = diag(p))
    x <- x %*% A    # check var(x) is an identity matrix
  }
  #############################################
  pi <- c(exp(x %*% myt))
  pbar <- 1/sum(exp(x %*% myt))
  pi <- pbar*pi
  v4 <- c()
  for (j1 in (1:p)){
    for (j2 in (j1:p)){
      for (j3 in (j2:p)){
        for (j4 in (j3:p)){
          s0 <- sum(pi *x[, j1]*x[, j2]*x[,j3] *x[,j4])

          sj1 <- pi * x[,j1]; s.sj1 <- sum(sj1)
          sj2 <- pi * x[,j2]; s.sj2 <- sum(sj2)
          sj3 <- pi * x[,j3]; s.sj3 <- sum(sj3)
          sj4 <- pi * x[,j4]; s.sj4 <- sum(sj4)

          sj1j2 <- sj1*x[, j2]; s.sj1j2 <- sum(sj1j2)
          sj1j3 <- sj1*x[, j3]; s.sj1j3 <- sum(sj1j3)
          sj1j4 <- sj1*x[, j4]; s.sj1j4 <- sum(sj1j4)
          sj2j3 <- sj2*x[, j3]; s.sj2j3 <- sum(sj2j3)
          sj2j4 <- sj2*x[, j4]; s.sj2j4 <- sum(sj2j4)
          sj3j4 <- sj3*x[, j4]; s.sj3j4 <- sum(sj3j4)

          s.sj1j2j3 <- sum(pi *x[, j1]*x[, j2]*x[,j3])

          s_j4 <- s.sj1j2j3*s.sj4
          s_j3 <- sum(pi *x[, j1]*x[, j2]*x[,j4])*s.sj3
          s_j2 <- sum(pi *x[, j1]*x[, j3]*x[,j4])*s.sj2
          s_j1 <- sum(pi *x[, j2]*x[, j3]*x[,j4])*s.sj1


          tmp1 <- s.sj1j2 *s.sj3 * s.sj4 +
            s.sj1j3 * s.sj2 * s.sj4 +
            s.sj1j4 * s.sj3 * s.sj2 +
            s.sj2j3 * s.sj1 * s.sj4 +
            s.sj2j4 * s.sj3 * s.sj1 +
            s.sj3j4 * s.sj1 * s.sj2

          tmp2 <- s.sj1j2 *s.sj3j4 +
            s.sj1j3 * s.sj2j4 +
            s.sj2j3 * s.sj1j4

          tmp3 <- s.sj1*s.sj2*s.sj3*s.sj4

          s4 <- s0 - s_j4 - s_j3 - s_j2 - s_j1 + 2*tmp1 - tmp2 - 6*tmp3
          v4 <- c(v4, s4)
        }
      }
    }
  }
  return(v4)
}
###################################
#' @rdname dCGF
#' @param p Dimension.
#' @return \code{l_dhCGF} returns number of distinct third and
#' fourth derivatives.
#' @export
l_dhCGF <- function(p){
   l3 <- p + choose(p, 2)*2 + choose(p, 3)
   l4 <- p + 3 *choose(p, 2) + 3*choose(p, 3) + choose(p, 4)
   return(list("Third derivatives" = l3, "fourth derivatives" = l4))
}
############################
############################
#' @rdname dCGF
#' @return \code{dhCGF1D} returns third/fourth derivatives of univariate
#' empirical CGF, which are \code{d3hCGF} and \code{d4hCGF} when \eqn{p = 1}.
#' @export
#' @examples
#' #Univariate data
#' set.seed(1)
#' x <- rnorm(100)
#' t <- .3
#' dhCGF1D(t, x)
dhCGF1D <- function(t, x){
  z <- (x - mean(x))/stats::sd(x)
  nz <- length(z)
  pi <- c(exp(z*t))
  pbar <- 1/sum(exp(z * t))
  pi <- pbar*pi

  m1 <- sum(z*pi)
  m2 <- sum(z^2 * pi)
  m3 <- sum(z^3 *pi)
  m4 <- sum(z^4 * pi)

  t3 <- m3 - 3*m2*m1 + 2*m1^3
  t4 <- m4 - 4*m3 *m1 - 3*m2^2 + 12 *m2*m1^2 - 6 *m1^4

  return(list(t3 = t3, t4 = t4))
}



