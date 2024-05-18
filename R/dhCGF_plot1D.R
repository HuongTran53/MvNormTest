## -- Filename: dhCGF_plot1D.R
#' @title Graphical plots to assess multivariate univarite assumption of data.
#' @description
#' Plots the empirical third/fourth derivatives of cumulant generating function
#' together with confidence probability region.
#' Indication of non-normality is either violation of probability bands or
#'  curves with high slope.
#' @name Univariate_CGF_plot
#' @param x Univariate data
#' @param alpha Significant level (default is \eqn{.05})
#' @param method string, \code{"T3"} used the third derivatives,
#' and \code{"T4"} uses the fourth derivatives.
#' @return Plots
#' @export
#' @references
#' \insertRef{ref:ghoshuni}{PlotNormTest}
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

  # load("data/mt3_lst_param.RData")
  # load("data/mt4_lst_param.RData")

  sT3 <- c(
    92.42158217,	74.73877805,	60.7302487,	49.59500711,	40.7146832,
    33.61041073,	27.91023887,	23.3244267,	19.62665994,	16.63973104,
    14.22459407,	12.2719798,	10.69596014,	9.429002869,	8.418171,
    7.622206539,	7.00930269,	6.555417571,	6.243020088,	6.060187813,
    6
  )
  sT3 <- c(sT3, rev(sT3[-21]))
  # sT3 <- diag(mt3_lst_param$'1'$l.sTtTs)
  m.supLt3 <- 1.453986705

  sT4 <- c(568.1209021,	448.0515896,	355.0182299,	282.6795124,	226.2374173,
           182.0500981,	147.3446669,	120.0034577,	98.40450305,	81.30213209,
           67.73735645,	56.97044429,	48.4300794,	41.67496479,	36.3648053,
           32.23839726,	29.09714181,	26.79273777,	25.2181407,	24.30112719,
           24)

  sT4 <- c(sT4, rev(sT4[-21]))

  m.supLt4 <- 1.478985651

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
    graphics::lines(bigt, -bandt3, lty = 6, lwd = 2, col = "darkorange")
    graphics::lines(bigt, bandt3, lty = 6, lwd = 2, col = "darkorange")
    graphics::legend("top", c("T3","Probability bands"),
           lwd = c(2, 2),
           lty = c(1, 6), col = c("blue", "darkorange"),
           merge = TRUE, y.intersp = 2,text.width = .6)
    graphics::title(main = bquote(T[3] ~"plot"),
          adj = 0)
  } else {
    plot(bigt, T4, type = "l", lty = 1, lwd = 2, col = "blue",
         xlab = "t", ylab = bquote(T[4]),
         ylim = c(-ylimt4, ylimt4)
         )
    graphics::lines(bigt, bandt4, lty = 6, lwd = 2, col = "darkorange")
    graphics::lines(bigt, -bandt4, lty = 6, lwd = 2, col = "darkorange")
    graphics::legend("top", c("T4","Probability bands"),
           lwd = c(2, 2),
           lty = c(1, 6), col = c("blue","orange"),
           merge = TRUE, y.intersp = 2, text.width = .6)
    graphics::title(main = bquote(T[4] ~"plot"),
          adj = 0)
  }

}

