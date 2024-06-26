## -- Filename: hCGF_plot.R
#' @title Graphical plots to assess multivariate normality assumption.
#' @name Multivariate_CGF_PLot
#' @description Cumulant generating functions of normally distributed
#' random variables has derivatives of order higher than 3 are all 0.
#' Hence, plots of empirical third/fourth order derivatives with large value
#' or high slope gives indication of non-normality.
#' \code{Multivariate_CGF_PLot} estimates and provides confidence region for
#' average (or any linear combination) of third/fourth derivatives of empirical
#' cumulant function at the points \eqn{t = t^*1_p}. Plots for
#' \eqn{p = 2, 3, \dots, 10} will be faster to obtain, as confidence regions
#' and other necessary parameters are available in \code{mt3_lst_param.rda} and
#' \code{mt4_lst_param.rda}.
#' Higher dimension requires expensive computational cost.
#' @param x Data matrix of size \eqn{n \times p}
# @param l Vector of linear combination, having length is the number of
# unique third/fourth derivatives. Default is \code{NULL} where the algorithm
# will use the average of distinct derivatives.
#' @param alpha Significant level (default is \eqn{.05})
#' @return \code{d3hCGF_plot} returns plot relying in third derivatives.
#' @seealso \code{\link[=dhCGF_plot1D]{dhCGF_plot1D()}}
#' @export
#' @examples
#' set.seed(1234)
#' p <- 3
#' x <- MASS::mvrnorm(500, rep(0, p), diag(p))
#' d3hCGF_plot(x)
#' d4hCGF_plot(x)
d3hCGF_plot <- function(x,
                        # l = NULL,
                        alpha = 0.05){
  bigt <- seq(-1, 1, by = 0.05)
  p <- ncol(x); n <- nrow(x)
  lT <- p + choose(p, 2)*2 + choose(p, 3)
  numt <- length(bigt)
  l <- rep(1/sqrt(lT), lT)
  # This part allows other values of l
  # if (is.null(l)){
  #   l <- rep(1/sqrt(lT), lT)
  # }
  if (!(p %in% 1:10)){
    lst <- mt3_get_param(p, l = l)
    warning("Heavy computation")
  } else {
    Hvar <-
      matrix(
        nrow = 10, ncol = 21, byrow = T,
        data = c(
          92.42158217,	74.73877805,	60.7302487,	49.59500711,	40.7146832,
          33.61041073,	27.91023887,	23.3244267,	19.62665994,	16.63973104,
          14.22459407,	12.2719798,	10.69596014,	9.429002869,	8.418171,
          7.622206539,	7.00930269,	6.555417571,	6.243020088,	6.060187813,
          6,
          51.64735474,	42.06049707,	34.43383781,	28.34512273,	23.46757024,
          19.54747827,	16.3872706,	13.83262291,	11.76265797,	10.08245759,
          8.717328805,	7.608403151,	6.70925268,	5.983284647,	5.401735469,
          4.942128632,	4.58709455,	4.323475808,	4.141660713,	4.035103292,
          4,
          38.88149786,	31.69987493,	25.98318627,	21.41640716,	17.75565754,
          14.81153012,	12.43646553,	10.51516324,	8.957276844,	7.691833972,
          6.662962446,	5.826609195,	5.148016599,	4.599779034,	4.160346054,
          3.812871459,	3.544332294,	3.344860732,	3.207246338,	3.126577521,
          3.1,
          32.68733899,	26.64565406,	21.8369797,	17.99606145,	14.91756688,
          12.44204029,	10.4452675,	8.830198359,	7.520794105,	6.457327483,
          5.592783363,	4.890095636,	4.320022099,	3.859508021,	3.490425881,
          3.198606386,	2.973096793,	2.805598483,	2.690047987,	2.622315183,
          2.6,
          29.04678297,	23.66659881,	19.38571068,	15.96742076,	13.22853561,
          11.0268189,	9.251491708,	7.816018372,	6.652612096,	5.708038546,
          4.940402436,	4.316681083,	3.810827712,	3.402311133,	3.074991259,
          2.816254638,	2.616352869,	2.467900968,	2.365503721,	2.305486558,
          2.285714286,
          26.6578591,	21.70842882,	17.77156459,	14.62903647, 12.11197704,
          10.08929945,	8.458929531,	7.141148429,	6.0735229, 5.20703381,
          4.503111954,	3.931363246,	3.467819647,	3.093592715, 2.793836949,
          2.556952971,	2.373977783,	2.238122517,	2.144428152, 2.089517571,
          2.071428571,
          24.97305857,	20.3259158,	16.63063154,	13.68187755,	11.32079857,
          9.424102982,	7.895805196,	6.66095645,	5.6608713,	4.849482818,
          4.190552338,	3.655528423,	3.221900887,	2.871933871,	2.591690565,
          2.370283663,	2.199301874,	2.0723752,	1.984851205,	1.933561885,
          1.916666667,
          23.72267204,	19.29911891,	15.78260114,	12.97730357,	10.73175543,
          8.928416104,	7.475793238,	6.302459321,	5.352496366,	4.582018088,
          3.956507809,	3.448776016,	3.037390364,	2.705467401,	2.439742546,
          2.22985541,	2.067803051,	1.947525575,	1.864597561,	1.816005865,
          1.8,
          22.75872054,	18.50712799,	15.12814181,	12.43324677, 10.27663975,
          8.54519503,	7.150870282,	6.024942309,	5.113618989,	4.374689512,
          3.774960568,	3.288289533,	2.89407277,	2.576082348,	2.321570752,
          2.12058297,	1.965430277,	1.850291432,	1.770915745,	1.724409283,
          1.709090909,
          21.99337116,	17.87807412,	14.60812103,	12.00077286,	9.914711778,
          8.24030648,	6.892245795,	5.803947699,	4.923303532,	4.209429498,
          3.630180608,	3.160243679,	2.779671764,	2.472756467,	2.227160141,
          2.033249152,	1.883583916,	1.772532444,	1.69598262,	1.651135061,
          1.636363636
        )
      )

    Msup <- c(1.453986705, 1.445928528, 1.444276978, 1.44268609, 1.444282713,
              1.444054044, 1.446514826, 1.445421625, 1.44661582, 1.448280449
    )

    lst <- list(m.supLt = Msup[p],
                varLtLs = c(Hvar[p, ], Hvar[p, 20:1])
    )

    # lst <- mt3_lst_param[[p]]
    # This part costs memories, improve later.
    # if (any(l != rep(1/sqrt(lT), lT))){
    #   tmp <- mt3_covLtLs(l= l, p = p, sTtTs = lst$l.sTtTs)
    #   lst$sLtLs <- tmp$sLtLs
    #   lst$m.supLt <- tmp$m.supLt
    #   lst$varLtLs <- diag(lst$sLtLs)
    # }
  }
  v3 <- lapply(bigt/sqrt(p),function(u) d3hCGF(myt = rep(u, p), x = x))
  L3 <- unlist(lapply(v3, function(v) sqrt(n) * sum(l*v)))
  u <- sqrt(-2 *log(alpha))
  band1 <- (u + lst$m.supLt)*sqrt(lst$varLtLs)
  band2 <- -(u + lst$m.supLt)*sqrt(lst$varLtLs)
  ylim <- max(max(abs(band1)), abs(L3))
  plot(bigt, L3, ylim = c(-ylim, ylim), col = "blue",
       lty = 1, type = "l", lwd = 2, ylab = bquote(L[3]), xlab = "t")
  graphics::lines(bigt, band1, col = "orange", lty = 6, lwd = 2)
  graphics::lines(bigt, band2, col = "orange", lty = 6, lwd = 2)
  graphics::legend("top", c("Probability  bands", expression(L[3])),
                   lty = c(6, 1), col = c("orange", "blue"), lwd = c(2, 2),
                   merge = TRUE, y.intersp = 1.5, text.width = .6)
  graphics::title(main = bquote(MT[3] ~"plot"), adj = 0)
  til.L <- L3/sqrt(lst$varLtLs)
  deci <- ifelse(max(abs(til.L)) >= u + lst$m.supLt, "reject", "accept")
  return(deci)
}
######################################################
#' Title
#' @rdname Multivariate_CGF_PLot
#'
#' @return \code{d4hCGF_plot} returns plot relying in forth derivatives.
#' @export
d4hCGF_plot <- function(x,
                        # l = NULL,
                        alpha = 0.05){
  bigt <- seq(-1, 1, by = 0.05)
  p <- ncol(x); n <- nrow(x)
  lT4 <- p + 3 *choose(p, 2) + 3*choose(p, 3) + choose(p, 4)
  numt <- length(bigt)
  l <- rep(1/sqrt(lT4), lT4)
  # if (is.null(l)){
  #   l <- rep(1/sqrt(lT4), lT4)
  # }
  if (!(p %in% 1:10)){
    lst <- mt4_get_param(p, l = l)
    warning("Heavy computation")
  } else {
    Hvar <-
      matrix(nrow = 10, ncol = 21, byrow = T,
             data = c(
               568.1209021,	448.0515896,	355.0182299,	282.6795124,	226.2374173,
               182.0500981,	147.3446669,	120.0034577,	98.40450305,	81.30213209,
               67.73735645,	56.97044429,	48.4300794,	41.67496479,	36.3648053,
               32.23839726,	29.09714181,	26.79273777,	25.2181407,	24.30112719,
               24,
               280.078945,	220.4020183,	174.2856492,	138.5357285,	110.735366,
               89.05053643,	72.08563482,	58.77581847,	48.30654298,	40.05326137,
               33.53611834,	28.38583129,	24.31794483,	21.11337599,	18.60370546,
               16.66006845,	15.18479424,	14.10516504,	13.36883132,	12.94054861,
               12.8,
               192.017125,	150.4754251,	118.5211313,	93.85864813,	74.76032964,
               59.92208683,	48.35667067,	39.31457154,	32.22522209,	26.65317137,
               22.26533269,	18.80644783,	16.08066803,	13.93770524,	12.2624126,
               10.96695142,	9.984921882,	9.266999309,	8.777740015,	8.493314772,
               8.4,
               139.4637519,	109.4943184,	86.42252828,	68.59431806,	54.76680757,
               44.00306309,	35.5945326,	29.00370386,	23.82158724,	19.73609938,
               16.50848791,	13.95570699,	11.93721173,	10.34504666,	9.096400137,
               8.128014685,	7.392003839,	6.85274504,	6.484607316,	6.270339886,
               6.2,
               104.9909838,	82.87394461,	65.77466418,	52.50044615,	42.15382226,
               34.05709724,	27.69674951,	22.68233197,	18.71599557,	15.56982032,
               13.0689038,	11.0787099,	9.495581318,	8.239611218,	7.249281821,
               6.477434023,	5.888246794,	5.454990303,	5.15838038,	4.98541008,
               4.928571429,
               82.011931,	65.17873306,	52.08615944,	41.85866521,	33.83509445,
               27.51431271,	22.51507366,	18.54629239,	15.38495284,	12.85963481,
               10.83819206,	9.218508179,	7.921543093,	6.886091921,	6.064829952,
               5.421329472,	4.927816535,	4.563496973,	4.313326754,	4.167136587,
               4.119047619,
               66.62374719,	53.30475074,	42.87967482,	34.68359842,	28.21195052,
               23.08037733,	18.99500689,	15.73033808,	13.11274092,	11.00810133,
               9.31253868,	7.945410092,	6.8440239,	5.959636079,	5.254414896,
               4.699140766,	4.27146887,	3.954627172,	3.73645638,	3.608724244,
               3.566666667,
               56.33114824,	45.31277767,	36.64217968,	29.78935694,	24.35000213,
               20.01455816,	16.54517401,	13.75847828,	11.51265495,	9.69771301,
               8.228135973,	7.037312106,	6.073302562,	5.295619503,	4.67277072,
               4.180389742,	3.799816819,	3.51703094,	3.321859343,	3.207411075,
               3.16969697,
               49.50809201,	39.95575422,	32.4136072,	26.4335349,	21.67217263,
               17.86555133,	14.81010191,	12.34837317,	10.35825587,	8.744824767,
               7.434142751,	6.36854146,	5.503017571,	4.802476072,	4.239619901,
               3.793336062,	3.447466197,	3.189878105,	3.011776472,	2.907207735,
               2.872727273,
               45.07928983,	36.41542436,	29.56901672,	24.13680211,	19.80897198,
               16.34684669,	13.5660621,	11.32386571,	9.50948683,	8.036814083,
               6.838811252,	5.863247519,	5.069425651,	4.425671118, 3.907404411,
               3.495663057,	3.175973205,	2.937495809,	2.772391759, 2.675365267,
               2.643356643
             )
      )

    Msup <- c(1.478985651,
              1.454432793,
              1.448680764,
              1.449313815,
              1.453737199,
              1.4561804,
              1.463732497,
              1.470303197,
              1.476776159,
              1.484990068)

    lst <- list(m.supLt = Msup[p],
                varLtLs = c(Hvar[p, ], Hvar[p, 20:1])
    )

  }
  v4 <- lapply(bigt/sqrt(p),function(u) d4hCGF(myt = rep(u, p), x = x))
  L4 <- unlist(lapply(v4, function(v) sqrt(n) * sum(l*v)))
  u <- sqrt(-2 *log(alpha))
  band1 <- (u + lst$m.supLt)*sqrt(lst$varLtLs)
  band2 <- -(u + lst$m.supLt)*sqrt(lst$varLtLs)
  ylim <- max(max(abs(band1)), abs(L4))
  plot(bigt, L4, ylim = c(-ylim, ylim), col = "blue",
       lty = 1, type = "l", lwd = 2, ylab = bquote(L[4]), xlab = "t")
  graphics::lines(bigt, band1, col = "orange", lty = 6, lwd = 2)
  graphics::lines(bigt, band2, col = "orange", lty = 6, lwd = 2)
  graphics::legend("top", c("Probability  bands", expression(L[4])),
                   lty = c(6, 1), col = c("orange", "blue"), lwd = c(2, 2),
                   merge = TRUE, y.intersp = 1.5, text.width = .6)
  til.L <- L4/sqrt(lst$varLtLs)
  deci <- ifelse(max(abs(til.L)) >= u + lst$m.supLt, "reject", "accept")
  return(deci)
}


