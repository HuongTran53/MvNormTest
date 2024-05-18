## -- Filename: cox.R --#
#' @title Graphical plots to assess the univarite noramality assumption of data.
#' @name Univariate_Score_function
#' @description Score function of a univariate normal distribution is
#' a straight line. A non-linear graph of score function estimator shows
#' evidence of non-normality.
#' @param x univariate data.
#' @param P vector of weight.
#' @param lambda smoothing parameter, default is \eqn{0.5}.
#' @param x.dist the minimum distance between two data points in vector x.
#' @return \code{cox} returns the estimate of score function.
#' \itemize{
#' \item{\code{x}: The updated univariate data if merging happens.}
#' \item{\code{a}: Score value estimated at \code{x}.}
#' \item{\code{P}: Updated weight (if merging happens).}
#' \item{\code{slt}: Index of merged data point
#' (is \code{NULL} if \code{x.dist = NULL}).}
#' }
#' @details
#' To avoid the singularity of coefficient matrices in spline method, points
#' with distance less than \code{x.dist} are merged and weight of the
#' representative points is updated by the summation of weight of
#'discarded points.
#' @export
#' @import MatrixExtra
#' @import Matrix
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{ref:sp_estimate}{PlotNormTest}
#' @examples
#' set.seed(1)
#' x <- rnorm(100, 2, 4)
#' re <- cox(sort(x))
#' plot(re$x, re$a, xlab = "x", ylab = "Estimated Score",
#'  main = "Estimator of score function")
#' abline(0, 1)
#'
cox <- function(x, P = NULL, lambda = 0.5, x.dist = NULL) {
  nx <- length(x)
  if (is.unsorted(x) == T){
    stop("x should be sort")
  }
  if (lambda < 0) {
    stop("penalty parameters is non negative")
  }

  if (is.null(P)) {
    P <- rep(1/nx, nx)
  } else {
    if (length(P) < nx) {
      stop("Incorrect dimensions")
    }
  }
  x <- (x - mean(x))/stats::sd(x) #standardize data
  #merge all close data points and updata x
  if (!is.null(x.dist)){
    y <- x[-1] - x[-nx]
    dif <- which(y < x.dist)
    #slt: data point that we skipped
    slt <- split(dif, cumsum(c(1, diff(dif) != 1)))
    if(sum(slt[[1]]) != 0){
      for (i in length(slt):1){
        temp1 <- slt[[i]]
        P[max(temp1) + 1] <- P[max(temp1) + 1] + sum(P[temp1])
        P <- P[-temp1]
        x <- x[-temp1]
      }
    }
    slt <- unlist(slt)
  } else {
    slt = NULL
  }
  nx <- length(x)
  onethd <- 1/3
  twothd <- 2/3
  ift <- 0
  nxm1 <- nx - 1; nxm2 <- nx - 2; nxm3 <- nx - 3
  nxp1 <- nx + 1
  lambda <- lambda
  h <- x[-1] - x[-nx]
  # Matrix C:
  c <- Matrix::bandSparse(nx, nxm2, -1:0, list(twothd *h[-1],
                                               onethd *h[-nxm1]))
  c[nx, nxm2] <- -onethd * h[nxm1]
  # Store matrix T as w7:
  w7 <- Matrix::bandSparse(nxm2, nx, 1:2, list(2*h[-length(h)]/ 3,  h[-1]/3))
  w7[1, 1] <- -h[1]/3
  # Matrix Q:
  Q <- Matrix::bandSparse(nx, nxm2, -2: 0, list(1 / h[-1],
                                                - 1/h[-nxm1] - 1/h[-1],
                                                1 / h[-nxm1]))
  # Matrix P^-1Q is pq:
  pq <- Matrix::bandSparse(
    nx, nxm2, -2:0, list(1/(P[-c(1,2)] * h[-1]),
                         1/P[-c(1, nx)] * ( -1/h[-nxm1] - 1/h[-1]),
                         1/(P[-c(nxm1, nx)] * h[-nxm1])))
  # Matrix A is a:
  a <- Matrix::bandSparse(nx, nx, 0: 1, list(-1/c(h, -h[nxm1]) , 1 / h ))
  a[nx, nx -1] <- -1/h[nxm1]
  #Matrix P^-1A' is pa:
  pa <- Matrix::bandSparse(nx, nx, -1:0, list(1/(h * P[-1]),
                                              -1/(P * c(h, -h[nxm1]))))
  pa[nxm1, nx] <- -1/(P[nxm1] * h[nxm1])
  #Matrix R is r:
  r <- Matrix::bandSparse(nxm2, nxm2, (-1): 1,
                          list(h[-c(1, length(h))] / 3,
                               2 / 3 *(h[-1] + h[-length(h)]),
                               h[-c(1, length(h))] / 3))
  #Pseudo y:
  #Create C'P1 store in cp1
  cp <- onethd * (h[-nxm1] * P[-c(nxm1, nx)]) + twothd * (h[-1] * P[-c(1, nx)])
  cp[nxm2] <- cp[nxm2] - onethd * h[nxm1]*P[nx]
  #invert of R:
  invr <- solve(r)
  # Peuso y is in suy
  suy <- pa %*% P - pq%*%t(invr)%*%cp
  #compute R + 2lambda * Q'*P^-1 *Q
  temp <- r + 2*lambda* t(Q) %*% pq
  invtemp <- solve(temp)
  #compute Q'y:
  qy <- suy[-c(nxm1, nx)] * (1/h[-nxm1]) -
    suy[-c(1,nx)] *(1/h[-nxm1] + 1/h[-1]) + suy[-c(1,2)] * (1/h[-1])
  #compute a:
  a <- suy - 2 *lambda *pq %*% invtemp %*% qy
  return(list(x = x, a = a, P= P, slt = slt))
}
