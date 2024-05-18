# -- Filname: pos.R --#
#' @rdname CGF_transformations
#' @title Derivatives of empirical moment generating function (MGF).
# @description
# Given the position in the chain containg all derivatives of estimator
# \eqn{\hat{M}_X(t)}, get the information
# of which derivatives were taken.
#' @param j1 Index of the first variables
#' @param j2 Index of the first variables, should be at least \code{j1}
#' @param j3 Index of the first variables, should be at least \code{j2}
#' @param p Dimension
#' @return \code{mt3_rev_pos} returns the position of this particular derivative
#' in the chain of all derivatives, up to third order.
#' @details
#' The estimator of multivariate moment generating function is
#' \eqn{\hat{M}_X(t) = \dfrac{1}{n} \sum_{i = 1}^n \exp(t'X_i)}
#' The chain containing all derivatives up to the third order is
#' \deqn{
#' Z = \bigg(\hat{M}, \hat{M}^{001}, \dots \hat{M}^{00p},
#' \hat{M}^{011}, \hat{M}^{012}, \dots \hat{M}^{0pp}, \hat{M}^{111},
#' \hat{M}^{112}, \dots \hat{M}^{ppp}\bigg)'
#' }
#' and
#' \deqn{
#' \hat{M} = \hat{M}^{000}(t)= \hat{M}_X(t)
#' }
#' \deqn{
#' \hat{M}^{j_1j_2j_3}(t) =
#' \dfrac{\partial^k}{\partial t_{j_1} t_{j_2} t_{j_3}} \hat{M}(t)
#' }
#'  where \eqn{k} is the number of \eqn{j_1, j_2, j_3} different from 0.
#' Similar notation is applied when fourth derivatives is used.
#' @export
#' @examples
#' mt3_rev_pos(1, 2, 2, p = 3)
#' p <- 3
#' mt3_pos(p)
#' mt4_pos(p)
mt3_rev_pos <- function(j1, j2, j3, p){
  if (is.unsorted(c(j1, j2, j3))){
    stop("Index should be is descending order")
  }
  pos <- 0
  if (j1 >= 1){
    for (t in 0:(j1-1)){
      for (l in 0:(p-1-t)){
        pos <- pos + (p+1 - t-l)
      }
    }
  }
  pos <- pos + j1
  if (j1 <= (j2 - 1)){
    for (l in ((j1):(j2 - 1))){

      pos <- pos + p + 1 - l
    }
  }
  pos <- pos + j3 - (j2 -1)
  return(pos)
}
##################################
##################################
#' @rdname CGF_transformations
#' @description
#' Given dimension \eqn{p}, returns a dataframe containing the position of
#' all derivatives of
#' estimator of moment generating function \eqn{\hat{M}_X(t)},
#' upto third/fourth order.
# Reverse of \code{\link[=rev_pos]{rev_pos()}}.
#' @param p Dimension
#' @return \code{mt3_pos} an array contaning all position with respect
#' to index of \eqn{j_1, j_2, j_3}.
# @seealso \code{\link[=rev_pos]{rev_pos()}}
#' @export
mt3_pos <- function(p){
  j1j2j3 <- c()
  for (j1 in 0:p){
    for (j2 in j1:p){
      for (j3 in j2:p){
        j1j2j3 <- rbind(j1j2j3, c(j1, j2, j3))
      }
    }
  }
  return(j1j2j3)
}
##################################
##################################
#' @rdname CGF_transformations
#' @return \code{mt4_pos} an array contaning all position with respect to
#' the index of \eqn{j_1, j_2, j_3, j_4}.
#' @export
mt4_pos <- function(p){
  j1j2j3j4 <- c()
  for (j1 in 0:p){
    for (j2 in j1:p){
      for (j3 in j2:p){
        for (j4 in j3:p){
          j1j2j3j4  <- rbind(j1j2j3j4, c(j1, j2, j3, j4))
        }

      }
    }
  }
  return(j1j2j3j4)
}
