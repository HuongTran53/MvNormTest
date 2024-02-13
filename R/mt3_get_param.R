#-- Filename: mt3_get_params.R --#
#' @name get_params
#' @title Get parameters for plots of multivariate CGF to assess normality assumption.
#' @description
#' Obtain neccesaary paramaters to build a graphical test using third/fourth derivatives
#' of cumulant generating function. The set of availble parameters for test in dimensions \code{p = 1, 2,..., 10},
#' test points \eqn{t = t^*1_p} with value of \eqn{t^*} in
#' \code{bigt = seq(-1, 1, by = 0.05)/sqrt(p)} and \code{l} chosen as the average of distinct derivatives
#' is available, see
#' \code{\link[=data/mt3_lst_param.RDATA]{mt3_lst_param.RDATA}} and
#' \code{\link[=data/mt4_lst_param.RDATA]{mt4_lst_param.RDATA}}.
#'
#'
#' @param p Dimension.
#' @param l Linear transformation of vector of third/fourth distinct derivatives, default is their average.
#' @param bigt Array containing value of \eqn{t^*}.
#' @return
#' \itemize{
#' \item{\code{p}}{ Dimension.}
#' \item{\code{lT}}{ Number of distinct third/fourth order derivatives.}
#' \item{\code{sTtTs}}{ Two dimensional array, each element contains covariance
#'  matrix of vector of derivatives, the function called \code{\link[=mt3_covTtTs]{mt3_covTtTs()}}, \code{\link[=mt4_covTtTs]{mt4_covTtTs()}}}.
#' \item{\code{l.sTtTs}}{ Covariance matrix of linear combination of distinct derivatives, the function called \code{\link[=mt3_covLtLs]{mt3_covLtLs()}}, \code{\link[=mt4_covLtLs]{mt4_covLtLs()}}.}
#' \item{\code{m.supLT}}{ The Monte Carlo estimate of expected value supremum of the Gaussian process, see \code{\link[=mt4_covLtLs]{mt4_covLtLs()}}.
#' }
#' }
#' \code{mt3_get_param} returns necessary parameters for the 2D plot
#' relying on third derivatives.
#' \code{mt4_get_param} returns necessary parameters for the 2D plot
#' relying on fourth derivatives.
#'
#' @seealso \code{\link[=covZtZs]{covZtZs()}},
#' \code{\link[=covLtLs]{covLtLs()}}, \code{\link[=covTtTs]{covTtTs()}}
#' @export
#'
#' @examples
#' p <- 2
#' mt3 <- mt3_get_param(p, bigt = seq(-1, 1, .4)/sqrt(p))
#' names(mt3)
#' mt4 <- mt4_get_param(p, bigt = seq(-1, 1, .4)/sqrt(p))
#' names(mt4)
mt3_get_param <- function(p,
                          bigt = seq(-1, 1, by = .05)/sqrt(p),
                          l = NULL){
  print(paste("MT3 params for dimension: ", p))


  lT <- p + choose(p, 2)*2 + choose(p, 3)

  numt <- length(bigt)
  pos.matrix <- mt3_pos(p)
  l.sTtTs  <- mt3_covTtTs(bigt = bigt, p)

  ############################
  if (is.null(l)){
    l <- rep(1/sqrt(lT), lT)
    # l <- rep(1/, lT)
  }

  tmp <- mt3_covLtLs(l= l, p = p, sTtTs = l.sTtTs, bigt)

  sLtLs <- tmp$sLtLs
  m.supLt <- tmp$m.supLt
  ############################
  mt3 <- list(p = p, lT = lT, bigt = bigt, l = l,
              pos.matrix = pos.matrix,
              l.sTtTs = l.sTtTs,
              sLtLs = sLtLs,
              m.supLt = m.supLt
  )
  return(mt3)
}


################################
################################
#' @rdname get_params
#' @export
#'
mt4_get_param <- function(p, bigt = seq(-1, 1, by = .05)/sqrt(p),
                          l = NULL){
  print(paste("MT4 params for dimension: ", p))


  lT4 <- p + 3 *choose(p, 2) + 3*choose(p, 3) + choose(p, 4)

  numt <- length(bigt)
  pos.matrix <- mt4_pos(p)
  l.sTtTs  <- mt4_covTtTs(bigt = bigt, p)

  ############################
  if (is.null(l)){
    l <- rep(1/sqrt(lT4), lT4)
    # l <- rep(1/, lT)
  }

  tmp <- mt4_covLtLs(l= l, p = p, sTtTs = l.sTtTs)

  sLtLs <- tmp$sLtLs
  m.supLt <- tmp$m.supLt
  ############################
  mt4 <- list(p = p, lT = lT4, bigt = bigt, l = l,
              pos.matrix = pos.matrix,
              l.sTtTs = l.sTtTs,
              sLtLs = sLtLs,
              m.supLt = m.supLt
  )
  return(mt4)
}

