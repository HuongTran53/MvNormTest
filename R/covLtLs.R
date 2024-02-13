#' @name covLtLs
#' @title Linear combinations of distinct derivatives of empirical cumulant generating function (CGF).
#' @description
#' Linear combination of third/fourth derivatives of CGF results in an asymptotically
#' univariate Gaussian process with mean 0 and covariance between two points
#'  \eqn{t \in \mathbb{R}^p} and \eqn{s \in \mathbb{R}^p} is defined.
#' We are considering vector \eqn{t} and \eqn{s} as the form  \eqn{t = t^*1_p}
#' and \eqn{s = s^*1_p}.
#'
#'
#' @param l vector of linear combination with length is the number of distinct derivatives, see \code{\link[=l_dhCGF]{l_dhCGF()}}.
#' @param p dimension.
#' @param bigt array of value \eqn{t^*} and \eqn{s^*}.
#' @param sTtTs Covariance matrix of derivatives vector, see \code{\link[=covTsTt]{covTsTt()}}. Default is \code{NULL}, when the algorithm
#' will call \code{\link[=mt3_covTsTt]{mt3_covTsTt()}} or \code{\link[=mt4_covTsTt]{mt4_covTsTt()}}.
#' @param seed Random seed to get the estimated of the supremum of the univariate Gaussian process obtained from the linear combination.
#'
#' @return
#' \itemize{
#' \item{\code{sLtLs} covariance matrix of the linear combination of distinct derivatives, which is a zero-mean Gaussian process. }
#' \item{\code{m.supLt} Monte-Carlo estimates of supremum of this Gaussian process}
#' }
#' \code{mt3_covLtLs} returns values related to the use of third derivatives.
#' \code{mt4_covLtLs} returns values related to the use of fourth derivatives.
#'
#' @export
#' @examples
#' bigt <- seq(-1, 1, .4)
#' p <- 3
#' # Third derivatives
#' lT3 <- l_dhCGF(p)[[1]]
#' l3 <- rep(1/sqrt(lT3), lT3)
#' mt3_covLtLs(l = l3, p = p, bigt = bigt/sqrt(p), seed = 1)
#' #fourth derivatives
#' lT4 <- l_dhCGF(p)[[2]]
#' l4 <- rep(1/sqrt(lT4), lT4)
#' mt4_covLtLs(l = l4, p = p, bigt = bigt/sqrt(p), seed = 1)

mt3_covLtLs <- function(l, p, bigt = seq(-1, 1, 0.05)/sqrt(p),
                        sTtTs = NULL, seed = 1){

  # bigt = seq(-1, 1, by = 0.05)/sqrt(p)
  # bigt <- bigt/sqrt()

  if (is.null(l)){
    lT3 <- p + choose(p, 2)*2 + choose(p, 3)
    l <- rep(1/sqrt(lT3), lT3)
  }

  if (is.null(sTtTs)){
    sTtTs <- mt3_covTtTs(bigt, p = p)
  }

  numt <- length(bigt)

  numt <- nrow(sTtTs)
  sLtLs <- array(0, dim = c(numt, numt))

  for (i in 1:numt){
    for (j in i:numt){
      sLtLs[i, j] <- t(l) %*% sTtTs[[i, j]] %*% l
      sLtLs[j, i] <- sLtLs[i, j]
    }
  }

  til.Lt <- diag(1/sqrt(diag(sLtLs))) %*% sLtLs %*% diag(1/sqrt(diag(sLtLs)))


  set.seed(seed)
  sup.Lt <- replicate(8e3, expr= {
    Ltrep <- MASS::mvrnorm(1, mu = rep(0, numt), Sigma = til.Lt)
    max(abs(Ltrep))
  })
  set.seed(NULL)

  m.supLt <- mean(sup.Lt)
  return(list(sLtLs = sLtLs, m.supLt = m.supLt))

}

#############################################
#############################################
#' Title
#' @rdname covLtLs
#' @export
#'
mt4_covLtLs <- function(l, p, bigt = seq(-1, 1, 0.05)/sqrt(p),
                        sTtTs = NULL, seed = 1){

  if (is.null(l)){
    lT4 <- p + 3 *choose(p, 2) + 3*choose(p, 3) + choose(p, 4)
    l <- rep(1/sqrt(lT4), lT4)
  }

  if (is.null(sTtTs)){
    sTtTs <- mt4_covTtTs(bigt, p = p)
  }

  numt <- nrow(sTtTs)
  sLtLs <- array(0, dim = c(numt, numt))

  for (i in 1:numt){
    for (j in i:numt){
      sLtLs[i, j] <- t(l) %*% sTtTs[[i, j]] %*% l
      sLtLs[j, i] <- sLtLs[i, j]
    }
  }

  til.Lt <- diag(1/sqrt(diag(sLtLs))) %*% sLtLs %*% diag(1/sqrt(diag(sLtLs)))


  set.seed(seed)
  sup.Lt <- replicate(8e3, expr= {
    Ltrep <- MASS::mvrnorm(1, mu = rep(0, numt), Sigma = til.Lt)
    max(abs(Ltrep))
  })
  set.seed(NULL)

  m.supLt <- mean(sup.Lt)
  return(list(sLtLs = sLtLs, m.supLt = m.supLt))

}
