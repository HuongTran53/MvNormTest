#' @title Graphical plots to assess multivariate univarite assumption of data.
#' @description
#' Plots the empirical third/fourth derivatives of cumulant generating function together with confidence probability region.
#' Indication of non-normality is either violation of probability bands or curves with high slope.
#'
#' @name Univariate_CGF_plot
#' @param x Univariate data
#' @param alpha Significant level (default is .05)
#' @param method string, \code{"T3"} used the third derivatives, and \code{"T4"} uses the fourth derivatives.
#'
#'
#' @return Plots
#' @export
#' @references
#' \insertRef{ref:ghoshuni}{MvNormTest}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' dhCGF_plot1D(x, method = "T3")
#' dhCGF_plot1D(x, method = "T4")
#'
dhCGF_plot1D <- function(x, alpha = 0.05, method){
  if (!(method %in% c("T3", "T4"))){
    stop("Method is either T3 or T4")
  }
  bigt = seq(-1, 1, by = 0.05)


  load("data/mt3_lst_param.RData")
  load("data/mt4_lst_param.RData")

  sT3 <- diag(mt3_lst_param$'1'$l.sTtTs)
  m.supLt3 <- mt3_lst_param$'1'$m.supLt

  sT4 <- diag(mt4_lst_param$'1'$l.sTtTs)
  m.supLt4 <- mt4_lst_param$'1'$m.supLt


  nx <- length(x)
  re <- lapply(bigt, dhCGF1D, x = x)
  T3 <- sqrt(nx)* unlist(lapply(re, function(i) i$t3))
  T4 <- sqrt(nx)* unlist(lapply(re, function(i) i$t4))

  u <- sqrt(-2 *log(alpha))
  bandt4 <- (u + m.supLt4)*sqrt(unlist(sT4))
  bandt3 <- (u + m.supLt3)*sqrt(unlist(sT3))

  ylimt3 <- max(40, abs(T3))
  ylimt4 <- max(95, abs(T4))

  if (method == "T3"){
    plot(bigt, T3, type = "l", lty = 1, lwd = 2, col = "blue",
         xlab = "t",
         ylab = bquote(T[3]),
         ylim = c(-ylimt3, ylimt3)
    )
    lines(bigt, -bandt3, lty = 6, lwd = 2, col = "darkorange")
    lines(bigt, bandt3, lty = 6, lwd = 2, col = "darkorange")
    legend("top", c("T3","Probability bands"),
           lwd = c(2, 2),
           lty = c(1, 6), col = c("blue", "darkorange"),
           merge = TRUE, y.intersp = 2,text.width = .6)
    title(main = bquote(T[3] ~"plot"),
          adj = 0)
  } else {
    plot(bigt, T4, type = "l", lty = 1, lwd = 2, col = "blue",
         xlab = "t", ylab = bquote(T[4]),
         ylim = c(-ylimt4, ylimt4)
         )
    lines(bigt, bandt4, lty = 6, lwd = 2, col = "darkorange")
    lines(bigt, -bandt4, lty = 6, lwd = 2, col = "darkorange")
    legend("top", c("T4","Probability bands"),
           lwd = c(2, 2),
           lty = c(1, 6), col = c("blue","orange"),
           merge = TRUE, y.intersp = 2, text.width = .6)
    title(main = bquote(T[4] ~"plot"),
          adj = 0)
  }
}

