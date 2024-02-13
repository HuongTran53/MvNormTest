#' @name covTtTs
#' @title Covariance matrix of derivatives of sample cumulant generating function (CGF).
#'
#' @description Stacking third/fourth derivatives of sample CGF together to obtain a vector,
#' which (under normality assumption) approaches a normally distributed vector with zero mean and a covariance matrix.
#' More specifically, \code{covTsTs} computes covariance between any two
#' points as the form \eqn{t = t^*1_p} and \eqn{s = s^*1_p}.
#' @details Number of distinct third derivatives is
#' \deqn{
#' l_{T_3}= p + 2 \times \begin{pmatrix}
#' p\\2
#' \end{pmatrix} + \begin{pmatrix}
#' p \\ 3
#' \end{pmatrix}
#' }
#'
#' Number of distinct fourth derivatives is
#' \deqn{
#' l_{T_4} = p + 3 \times \begin{pmatrix}
#' p\\2
#' \end{pmatrix} + 3 \times \begin{pmatrix}
#' p \\ 3
#' \end{pmatrix} + \begin{pmatrix}
#' p \\ 4
#' \end{pmatrix}
#' }
#' For each pairs of \eqn{(t^*, s^*)}, \code{covTsTt} results a covariance matrix of size \eqn{l_{T_3} \times l_{T_3}} or \eqn{l_{T_4} \times l_{T_4}}.
#' @param bigt array contains value of \eqn{t^*}.
#' @param p dimension.
#' @param pos.matrix matrix containing information of position of any derivatives. Default is \code{NULL},
#' the function will call \code{\link[=mt3_pos]{mt3_pos()}} or \code{\link[=mt4_pos]{mt4_pos()}}.
#'
#' @return A 2 dimensional upper square triangular array, with size equals to length of \code{bigt},
#'  each position contains a covariance matrix of derivatives sequences between any two points \eqn{t = t^* 1_p} and
#'  \eqn{s = s^*1_p}.
#'  Diagonal of this array is the variance/covariance matrix of derivatives sequences, taking at points
#'  \eqn{t = t^*1_p}.
#'
#' \code{mt3_covTsTt} returns result relating the use of third derivatives.
#' @export
#' @examples
#' bigt <- seq(-1, 1, .4)
#' p <- 2
#' # Third derivatives
#' mt3_pos.matrix <- mt3_pos(p)
#' sTsTt3 <- mt3_covTtTs(bigt = bigt, p = p, pos.matrix = mt3_pos.matrix)
#' dim(sTsTt3)
#' sTsTt3[1:5, 1:5]
#' # Fourth derivatives
#' mt4_pos.matrix <- mt4_pos(p)
#' sTsTt4 <- mt4_covTtTs(bigt = bigt, p = p, pos.matrix = mt4_pos.matrix)
#' dim(sTsTt4)
#' sTsTt4[1:5, 1:5]

mt3_covTtTs <- function(bigt, p = 1, pos.matrix= NULL){
  if (is.null(pos.matrix)){
    pos.matrix<- mt3_pos(p)
  }

  numt <- length(bigt)

  re <- vector(mode = "list", length = numt * numt)
  dim(re) <- c(numt, numt)

  lst_matrixA <- lapply(bigt, function(t) mt3_matrix_A(rep(t, p)))


  for (i in 1:numt){
    myt <- rep(bigt[i], p)
    capAt <- lst_matrixA[[i]]
    for (j in i:numt){
      mys <- rep(bigt[j], p)
      # capAs <- mt3_matrix_A(t = mys)
      capAs <- lst_matrixA[[j]]
      sZtZs <- mt3_covZtZs(t = myt, s = mys, pos.matrix= pos.matrix)
      sTtTs <- capAt %*% sZtZs %*% t(capAs)
      re[[i, j]] <- sTtTs
    }
  }

  return(re)
}

#################################################
#################################################
#' @rdname covTtTs
#' @return \code{mt4_covTsTt} returns result relating the use of fourth derivatives.
#' @export
#'
mt4_covTtTs <- function(bigt, p = 1, pos.matrix= NULL){
  if (is.null(pos.matrix)){
    pos.matrix <- mt4_pos(p)
  }
  numt <- length(bigt)
  re <- vector(mode = "list", length = numt * numt)
  dim(re) <- c(numt, numt)

  lst_matrixA <- lapply(bigt, function(t) mt4_matrix_A(rep(t, p)))


  for (i in 1:numt){

    myt <- rep(bigt[i], p)
    capAt <- lst_matrixA[[i]]
    for (j in i:numt){
      mys <- rep(bigt[j], p)
      # capAs <- mt3_matrix_A(t = mys)
      capAs <- lst_matrixA[[j]]
      sZtZs <- mt4_covZtZs(t = myt, s = mys, pos.matrix= pos.matrix)
      sTtTs <- capAt %*% sZtZs %*% t(capAs)
      re[[i, j]] <- sTtTs
    }
  }

  return(re)
}

