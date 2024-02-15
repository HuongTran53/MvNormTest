# -- Filename: hCGF_plot.R ---#
#' @title Graphical plots to assess multivariate normality assumption of data.
#' @name Multivariate_CGF_PLot
#' @description
#' Cumulant generating functions of normally distributed random variables has derivatives
#' of order higher than 3 are all 0. Hence, plots of empirical third/fourth order
#' derivatives with large value or high slope gives indication of non-normality.
#' \code{Multivariate_CGF_PLot} estimates and provides confidence region for average (or any linear combination)
#' of third/fourth derivatives of empirical cumulant function at the points \eqn{t = t^*1_p}.
#' Plots for \eqn{p = 2, 3, \dots, 10} will be faster to obtain, as confidence regions and other necessary parameters are available in
#'  \code{mt3_lst_param.rda} and \code{mt4_lst_param.rda}.
#'  Higher dimension requires expensive computational cost.
#'
#' @param x Data matrix of size \eqn{n \times p}
#' @param l Vector of linear combination, having length is the number of unique third/fourth
#' derivatives. Default is \code{NULL} where the algorithm will use the average of distinct derivatives.
#'
#' @return \code{d3hCGF_plot} returns plot relying in third derivatives.
#'
#' @seealso \code{\link[=dhCGF_plot1D]{dhCGF_plot1D()}}
#' @export
#'
#' @examples
#' # Normal distribution
#' set.seed(1234)
#' p <- 3
#' x <- MASS::mvrnorm(500, rep(0, p), diag(p))
#' d3hCGF_plot(x)
#' d4hCGF_plot(x)

d3hCGF_plot <- function(x, l = NULL, alpha = 0.05){

  bigt <- seq(-1, 1, by = 0.05)
  p <- ncol(x); n <- nrow(x)
  lT <- p + choose(p, 2)*2 + choose(p, 3)
  numt <- length(bigt)

  if (is.null(l)){
    l <- rep(1/sqrt(lT), lT)
  }

  if (!(p %in% 1:10)){
    lst <- mt3_get_param(p, l = l)
    warning("Heavy computation")
  } else {

    # if (p %in% 1:5){
    #   load("data/mt3_lst_param.rda")
    # } else {
    #   load("data/mt3_lst_param2.rda")
    # }

    lst <- mt3_lst_param[[p]]


    if (any(l != rep(1/sqrt(lT), lT))){
      tmp <- mt3_covLtLs(l= l, p = p, sTtTs = lst$l.sTtTs)
      lst$sLtLs <- tmp$sLtLs
      lst$m.supLt <- tmp$m.supLt
      lst$varLtLs <- diag(lst$sLtLs)
    }

  }

  v3 <- lapply(bigt/sqrt(p),function(u) d3hCGF(myt = rep(u, p), x = x))
  L3 <- unlist(lapply(v3, function(v) sqrt(n) * sum(l*v)))


  u <- sqrt(-2 *log(alpha))
  band1 <- (u + lst$m.supLt)*sqrt(lst$varLtLs)
  band2 <- -(u + lst$m.supLt)*sqrt(lst$varLtLs)
  ylim <- max(max(abs(band1)), abs(L3))

  plot(bigt, L3, ylim = c(-ylim, ylim), col = "blue", lty = 1, type = "l", lwd = 2,
       ylab = bquote(L[3]),
       xlab = "t")
  lines(bigt, band1, col = "orange", lty = 6, lwd = 2)
  lines(bigt, band2, col = "orange", lty = 6, lwd = 2)
  legend("top", c("Probability  bands", expression(L[3])),
         lty = c(6, 1), col = c("orange", "blue"),
         lwd = c(2, 2),
         merge = TRUE, y.intersp = 1.5, text.width = .6)
  title(main = bquote(MT[3] ~"plot"),
        adj = 0)

  til.L <- L3/sqrt(lst$varLtLs)
  deci <- ifelse(max(abs(til.L)) >= u + lst$m.supLt, "reject", "accept")

  return(deci)
}

######################################################
######################################################
#' Title
#' @rdname Multivariate_CGF_PLot
#'
#' @return \code{d4hCGF_plot} returns plot relying in forth derivatives.
#' @export
#'
#' @examples
d4hCGF_plot <- function(x, l = NULL, alpha = 0.05){

  bigt <- seq(-1, 1, by = 0.05)
  p <- ncol(x); n <- nrow(x)
  lT4 <- p + 3 *choose(p, 2) + 3*choose(p, 3) + choose(p, 4)
  numt <- length(bigt)

  if (is.null(l)){
    l <- rep(1/sqrt(lT4), lT4)
  }

  if (!(p %in% 1:10)){
    lst <- mt4_get_param(p, l = l)
    warning("Heavy computation")
  } else {

    # if (p %in% 1:6){
    #   load("data/mt4_lst_param.rda")
    # } else {
    #   if (p == 10){
    #     load("data/mt4_lst_param10.rda")
    #   } else {
    #     load("data/mt4_lst_param2.rda")
    #   }
    # }

    lst <- mt4_lst_param[[p]]

    if (any(l != rep(1/sqrt(lT4), lT4))){
      tmp <- mt4_covLtLs(l= l, p = p, sTtTs = lst$l.sTtTs)
      lst$sLtLs <- tmp$sLtLs
      lst$m.supLt <- tmp$m.supLt
      lst$varLtLs <- diag(lst$sLtLs)
    }

  }

  v4 <- lapply(bigt/sqrt(p),function(u) d4hCGF(myt = rep(u, p), x = x))
  L4 <- unlist(lapply(v4, function(v) sqrt(n) * sum(l*v)))

  u <- sqrt(-2 *log(alpha))
  band1 <- (u + lst$m.supLt)*sqrt(lst$varLtLs)
  band2 <- -(u + lst$m.supLt)*sqrt(lst$varLtLs)
  ylim <- max(max(abs(band1)), abs(L4))
  plot(bigt, L4, ylim = c(-ylim, ylim), col = "blue", lty = 1, type = "l", lwd = 2,
       ylab = bquote(L[4]),
       xlab = "t")
  lines(bigt, band1, col = "orange", lty = 6, lwd = 2)
  lines(bigt, band2, col = "orange", lty = 6, lwd = 2)
  legend("top", c("Probability  bands", expression(L[4])),
         lty = c(6, 1), col = c("orange", "blue"),
         lwd = c(2, 2),
         merge = TRUE, y.intersp = 1.5, text.width = .6)

  til.L <- L4/sqrt(lst$varLtLs)
  deci <- ifelse(max(abs(til.L)) >= u + lst$m.supLt, "reject", "accept")
  return(deci)
}


