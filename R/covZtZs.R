#' @name covZtZs
#' @title Covariance matrix of derivatives of sample moment generating function (MGF).
#'
#' @description
#' Stacking derivatives upto third/fourth order of sample MGF
#' together to obtain a vector, which (under normality assumption) approaches a normally distributed vector
#' with zero mean and a covariance matrix.
#' \code{covZtZs} calculates covariance between any two points \eqn{t} and \eqn{s} in \eqn{\mathbb{R}^p}.
#' @param t,s  a vector of length \eqn{p}.
# @param s a vector of length \eqn{p}.
#' @param pos.matrix matrix contains information of position of any derivatives. Default is \code{NULL},
#' where the function will call \code{\link[=mt3_pos]{mt3_pos()}} or \code{\link[=mt4_pos]{mt4_pos()}}.
#' @return \code{mt3_covZtZs} Covariance matrix relating to the use of third derivatives.
# @seealso \code{\link[=pos]{pos()}}.
#' @export
#' @examples
#' set.seed(1)
#' p <- 3
#' x <- MASS::mvrnorm(100, rep(0, p), diag(p))
#' t <- rep(0.2, p)
#' s <- rep(-.3, p)
#' # Using third derivatives
#' pos.matrix3 <- mt3_pos(p)
#' sZtZs3 <- mt3_covZtZs(t, s, pos.matrix = pos.matrix3)
#' dim(sZtZs3)
#' sZtZs3[1:5, 1:5]


mt3_covZtZs <- function(t, s, pos.matrix = NULL){
  p <- length(t)

  if (is.null(pos.matrix)){
    pos.matrix <- mt3_pos(p)
  }

  lT <- sum(unlist(lapply(2:(p +1), FUN = function(x) choose(x,k = 2))))
  lZ <- 1 + p + p*(p+1)/2 + lT

  ts <- t + s
  expt <- exp(sum(t*t)/2)
  exps <- exp(sum(s*s)/2)

  expts <- exp(sum(ts * ts)/2)

  sZtZs <- array(0, dim = c(lZ, lZ))
  #################
  for (i in 1:lZ){
    for (k in 1:lZ){
      j1 <- pos.matrix[i, 1]
      j2 <- pos.matrix[i, 2]
      j3 <- pos.matrix[i, 3]

      k1 <- pos.matrix[k, 1]
      k2 <- pos.matrix[k, 2]
      k3 <- pos.matrix[k, 3]

      tab1 <- as.data.frame(table(c(j1, j2, j3)))
      # muj <- prod(apply(tab1, MARGIN = 1, dMGF, t = t))
      muj <- prod(apply(tab1, MARGIN = 1, function(tab) dMGF(tab, t = t)))

      tab2 <-  as.data.frame(table(c(k1, k2, k3)))
      muk <- prod(apply(tab2, MARGIN = 1, dMGF, t = s))


      mytab <- as.data.frame(table(c(j1, j2, j3, k1, k2, k3)))
      sZtZs[i, k] <- prod(apply(mytab, MARGIN = 1, dMGF, t = ts))*expts - muj*expt * muk*exps
    }

  }

  return(sZtZs)
}


#####################################
#####################################
#' @rdname covZtZs
#' @return \code{mt4_covZtZs} Covariance matrix relating to the use of fourth derivatives,
#' which also contains information of the use of third derivatives \code{mt3_covZtZs}.
#' @export
#' @examples
#' # Using fourth derivatives
#' sZtZs4 <- mt4_covZtZs(t, s)
#' dim(sZtZs4)
#' sZtZs4[1:5, 1:5]
mt4_covZtZs <- function(t, s, pos.matrix = NULL){
  p <- length(t)

  if (is.null(pos.matrix)){
    pos.matrix <- mt4_pos(p)
  }

  # lT <- sum(unlist(lapply(2:(p +1), FUN = function(x) choose(x,k = 2))))
  # lZ <- 1 + p + p*(p+1)/2 + lT
  # lZ <- 5 # only univariate

  lT3 <- p + choose(p, 2)*2 + choose(p, 3)
  lZ3 <- 1 + p + p*(p+1)/2 + lT3
  lT4 <- p + 3 *choose(p, 2) + 3*choose(p, 3) + choose(p, 4)
  lZ4 <- (lT4 + lZ3)

  ts <- t + s
  expt <- exp(sum(t*t)/2)
  exps <- exp(sum(s*s)/2)

  expts <- exp(sum(ts * ts)/2)

  sZtZs <- array(0, dim = c(lZ4, lZ4))
  #################
  for (i in 1:lZ4){
    for (k in 1:lZ4){
      j1 <- pos.matrix[i, 1]
      j2 <- pos.matrix[i, 2]
      j3 <- pos.matrix[i, 3]
      j4 <- pos.matrix[i, 4]

      k1 <- pos.matrix[k, 1]
      k2 <- pos.matrix[k, 2]
      k3 <- pos.matrix[k, 3]
      k4 <- pos.matrix[k, 4]

      tab1 <- as.data.frame(table(c(j1, j2, j3, j4)))
      # muj <- prod(apply(tab1, MARGIN = 1, dMGF, t = t))
      muj <- prod(apply(tab1, MARGIN = 1, function(tab) dMGF(tab, t = t)))

      tab2 <-  as.data.frame(table(c(k1, k2, k3, k4)))
      muk <- prod(apply(tab2, MARGIN = 1, dMGF, t = s))


      mytab <- as.data.frame(table(c(j1, j2, j3, j4, k1, k2, k3, k4)))
      sZtZs[i, k] <- prod(apply(mytab, MARGIN = 1, dMGF, t = ts))*expts - muj*expt * muk*exps
    }

  }

  return(sZtZs)
}
