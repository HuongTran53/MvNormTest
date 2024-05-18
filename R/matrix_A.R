## -- Filename: matrix_A.R
#' @name  matrix_A
#' @title From derivatives of MGF to derivatives of CGF.
#' @description
#' Taylor expansion implies that vectors of derivatives of
#' \eqn{\log(\hat{M}_X(t))} can be approximated
#' by a linear combination of vectors of derivatives of \eqn{\hat{M}_X(t)}.
#' \code{matrix_A} results the corresponding
#' linear combinations.
#' @param t vector of \eqn{\mathbb{R}^p}
#' @return \code{mt3_matrix_A} returns coefficient matrix relating to the use
#' of third derivatives.
#' @examples
#' p <- 3
#' t <- rep(.2, p)
#' A3 <- mt3_matrix_A(t)
#' dim(A3)
#' A3[1:5, 1:5]
#' A4 <- mt4_matrix_A(t)
#' dim(A4)
#' A4[1:5, 1:5]
#' @export
mt3_matrix_A <- function(t){
  p <- length(t)
  lT <- p + choose(p, 2)*2 + choose(p, 3)
  lZ <- 1 + p + p*(p+1)/2 + lT
  capA <- array(0, c(lT, lZ))
  i <- 1
  for (j1 in 1:p){
    vA001 <- mt3_rev_pos(0, 0, j1, p)
    vA011 <- mt3_rev_pos(0, j1, j1, p)
    vA111 <- mt3_rev_pos(j1, j1, j1, p)
    for (j2 in j1:p){
      vA002 <- mt3_rev_pos(0, 0, j2, p)
      vA022 <- mt3_rev_pos(0, j2, j2, p)
      vA222 <- mt3_rev_pos(j2, j2, j2, p)
      vA012 <- mt3_rev_pos(0, j1, j2, p)

      vA112 <- mt3_rev_pos(j1, j1, j2, p)
      vA122 <- mt3_rev_pos(j1, j2, j2, p)
      for (j3 in j2:p){
        vA003 <- mt3_rev_pos(0, 0, j3, p)
        vA033 <- mt3_rev_pos(0, j3, j3, p)
        vA333 <- mt3_rev_pos(j3, j3, j3, p)

        vA013 <- mt3_rev_pos(0, j1, j3, p)
        vA023 <- mt3_rev_pos(0, j2, j3, p)
        vA113 <- mt3_rev_pos(j1, j1, j3, p)
        vA123 <- mt3_rev_pos(j1, j2, j3, p)
        vA133 <- mt3_rev_pos(j1, j3, j3, p)
        vA223 <- mt3_rev_pos(j2, j2, j3, p)
        vA233 <- mt3_rev_pos(j2, j3, j3, p)

        if (j1 == j2){
          if (j2 == j3){
            capA[i, 1] <- 3*t[j1] - t[j1]^3
            capA[i, vA001] <- 3*(t[j1]^2 - 1)
            capA[i, vA011] <- -3*t[j1]
            capA[i, vA111] <- 1
          } else {
            capA[i, 1] <- -t[j3]*(t[j1]^2 - 1)
            capA[i, vA001] <- 2 * t[j1]*t[j3]
            capA[i, vA003] <- t[j1]^2 - 1
            capA[i, vA011] <- -t[j3]
            capA[i, vA013] <- -2*t[j1]
            capA[i, vA113] <- 1
          }
        } else{
          if (j2 == j3){
            capA[i, 1] <- -t[j1]*(t[j2]^2 - 1)
            capA[i, vA001] <- t[j2]^2 - 1
            capA[i, vA002] <- 2*t[j1]*t[j2]
            capA[i, vA012] <- -2*t[j2]
            capA[i, vA022] <- -t[j1]
            capA[i, vA122] <- 1
          } else {
            capA[i, 1] <- -t[j1]*t[j2]*t[j3]
            capA[i, vA001] <- t[j2]*t[j3]
            capA[i, vA002] <- t[j1]*t[j3]
            capA[i, vA003] <- t[j2]*t[j1]
            capA[i, vA012] <- -t[j3]
            capA[i, vA013] <- -t[j2]
            capA[i, vA023] <- -t[j1]
            capA[i, vA123] <- 1
          }
        }
        i <- i+1
      }
    }
  }
  capA <- exp(-sum(t*t)/2) * capA
  return(capA)
}
############################################
############################################
#' @rdname matrix_A
#' @return \code{mt4_matrix_A} returns coefficient matrix relating to the
#' use of fourth derivatives.
#' @export
mt4_matrix_A <- function(t){
  p <- length(t)
  lT3 <- p + choose(p, 2)*2 + choose(p, 3)
  lZ3 <- 1 + p + p*(p+1)/2 + lT3
  lT4 <- p + 3 *choose(p, 2) + 3*choose(p, 3) + choose(p, 4)
  lZ4 <- (lT4 + lZ3)
  capA <- array(0, c(lT4, lZ4))
  pos.matrix  <- mt4_pos(p)
  pos.matrixT4 <- matrix(pos.matrix[(lZ4 - lT4 +1):lZ4, ], ncol = 4)
  for (j1 in 1: p){
    indT4 <- which(apply(pos.matrixT4, 1,
                         function(r) all(r == rep(j1, 4)))
                   )
    vA0001 <- which(
      apply(pos.matrix, 1, function(r) all(r == c(0, 0, 0, j1)))
      )
    vA0011 <- which(
      apply(pos.matrix, 1, function(r) all(r == c(0, 0, j1, j1)))
                    )
    vA0111 <- which(
      apply(pos.matrix, 1, function(r) all(r == c(0, j1, j1, j1)))
      )
    vA1111 <- which(
      apply(pos.matrix, 1, function(r) all(r == c(j1, j1, j1, j1)))
      )
    t1 <- t[j1]
    capA[indT4, 1] <- t1^4 - 6*t1^2 + 3
    capA[indT4, vA0001] <- -4*t1^3 + 12*t1
    capA[indT4, vA0011] <- 6*t1^2 - 6
    capA[indT4, vA0111] <- -4*t1
    capA[indT4, vA1111] <- 1
    for (j2 in j1:p){
      #########################
      vA0002 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(0, 0, 0, j2)))
        )
      vA0012 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(0, 0, j1, j2)))
        )
      vA0022 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(0, 0, j2, j2)))
        )
      vA0112 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(0, j1, j1, j2)))
        )
      vA0122 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(0, j1, j2, j2)))
        )
      vA0222 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(0, j2, j2, j2)))
        )
      vA1112 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(j1, j1, j1, j2)))
        )
      vA1222 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(j1, j2, j2, j2)))
        )
      vA1122 <- which(
        apply(pos.matrix, 1, function(r) all(r == c(j1, j1, j2, j2)))
        )
      if (j2 > j1){
        #########################
        # when c(j1, j1, j1, j2)
        indT4 <- which(
          apply(pos.matrixT4, 1, function(r) all(r == c(j1, j1, j1, j2)))
          )
        capA[indT4, 1] <- (t[j1]^3 - 3*t[j1])*t[j2]
        capA[indT4, vA0001] <- (-3*t[j1]^2 + 6)*t[j2]
        capA[indT4, vA0002] <- -t[j1]^3 + 3*t[j1]
        capA[indT4, vA0011] <- 3*t[j1]*t[j2]
        capA[indT4, vA0012] <- 3*t[j1]^2 - 3
        capA[indT4, vA0111] <- -t[j2]
        capA[indT4, vA0112] <- -3*t[j1]
        capA[indT4, vA1112] <- 1
        #########################
        # when c(j1, j2, j2, j2)
        indT4 <- which(
          apply(pos.matrixT4, 1, function(r) all(r == c(j1, j2, j2, j2)))
          )
        capA[indT4, 1] <- t[j1] *(t[j2]^3 - 3*t[j2])
        capA[indT4, vA0001] <- -t[j2]^3 + 3*t[j2]
        capA[indT4, vA0002] <- t[j1]*(-3*t[j2]^2 + 6)
        capA[indT4, vA0012] <- 3*t[j2]^2 -3
        capA[indT4, vA0022] <- 3*t[j1]*t[j2]
        capA[indT4, vA0122] <- -3*t[j2]
        capA[indT4, vA0222] <- -t[j1]
        capA[indT4, vA1222] <- 1
        #########################
        # when c(j1, j1, j2, j2)
        indT4 <- which(
          apply(pos.matrixT4, 1, function(r) all(r == c(j1, j1, j2, j2)))
          )
        capA[indT4, 1] <- t[j1]^2*t[j2]^2 - t[j1]^2 - t[j2]^2 + 1 + 1
        capA[indT4, vA0001] <- 2*t[j1]*(1 - t[j2]^2)
        capA[indT4, vA0002] <- 2*t[j2]*(1 - t[j1]^2)
        capA[indT4, vA0011] <- t[j2]^2 -1
        capA[indT4, vA0012] <- 4*t[j1]*t[j2]
        capA[indT4, vA0022] <- t[j1]^2 -1
        capA[indT4, vA0112] <- -2 *t[j2]
        capA[indT4, vA0122] <- -2*t[j1]
        capA[indT4, vA1122] <- 1
      }

      for (j3 in j2:p){
        vA0003 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, 0, 0, j3)))
          )
        vA0013 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, 0, j1, j3)))
          )
        vA0023 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, 0, j2, j3)))
          )
        vA0033 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, 0, j3, j3)))
          )
        vA0113 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, j1, j1, j3)))
          )
        vA0123 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, j1, j2, j3)))
          )
        vA0133 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, j1, j3, j3)))
          )
        vA0223 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, j2, j2, j3)))
          )
        vA0233 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(0, j2, j3, j3)))
          )
        vA1123 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(j1, j1, j2, j3)))
          )
        vA1223 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(j1, j2, j2, j3)))
          )
        vA1233 <- which(
          apply(pos.matrix, 1, function(r) all(r == c(j1, j2, j3, j3)))
          )
        if ((j3 > j2)&(j2 > j1)){
          #########################
          # when c(j1, j1, j2, j3)
          indT4 <- which(
            apply(pos.matrixT4, 1, function(r) all(r == c(j1, j1, j2, j3)))
            )
          capA[indT4, 1] <- (t[j1]^2 - 1)*t[j2]*t[j3]
          capA[indT4, vA0001] <- -2 *t[j1] * t[j2]* t[j3]
          capA[indT4, vA0002] <- (1 - t[j2]^2)* t[j3]
          capA[indT4, vA0003] <- (1 - t[j1]^2)* t[j2]
          capA[indT4, vA0011] <- t[j2] * t[j3]
          capA[indT4, vA0012] <- 2 * t[j1] * t[j3]
          capA[indT4, vA0013] <- 2 * t[j1] * t[j2]
          capA[indT4, vA0023] <- t[j1]^2 - 1
          capA[indT4, vA0112] <- - t[j3]
          capA[indT4, vA0113] <- - t[j2]
          capA[indT4, vA0123] <- -2*t[j1]
          capA[indT4, vA1123] <- 1
          #########################
          # when c(j1, j2, j2, j3)
          indT4 <- which(
            apply(pos.matrixT4, 1, function(r) all(r == c(j1, j2, j2, j3)))
            )
          capA[indT4, 1] <- t[j1]*(t[j2]^2 - 1)* t[j3]
          capA[indT4, vA0001] <- (1 - t[j2]^2)*t[j3]
          capA[indT4, vA0002] <- -2*t[j1]* t[j2]*t[j3]
          capA[indT4, vA0003] <- t[j1]*(1 - t[j2]^2)
          capA[indT4, vA0012] <- 2*t[j2]*t[j3]
          capA[indT4, vA0013] <- t[j2]^2 - 1
          capA[indT4, vA0022] <- t[j1]*t[j3]
          capA[indT4, vA0023] <- 2*t[j1]* t[j2]
          capA[indT4, vA0122] <- -t[j3]
          capA[indT4, vA0123] <-  - 2*t[j2]
          capA[indT4, vA0223] <- -t[j1]
          capA[indT4, vA1223] <- 1
          #########################
          # when c(j1, j2, j3, j3)
          indT4 <- which(
            apply(pos.matrixT4, 1, function(r) all(r == c(j1, j2, j3, j3)))
            )
          capA[indT4, 1] <- t[j1]*t[j2]*(t[j3]^2 -1)
          capA[indT4, vA0001] <- t[j2]*(t[j3]^2 - 1)
          capA[indT4, vA0002] <- t[j1]*(t[j3]^2 - 1)
          capA[indT4, vA0003] <- -2 *t[j1]*t[j2]*t[j3]
          capA[indT4, vA0012] <- t[j3]^2 - 1
          capA[indT4, vA0013] <- 2 * t[j2]*t[j3]
          capA[indT4, vA0023] <- 2 * t[j1]* t[j3]
          capA[indT4, vA0033] <- t[j1]*t[j2]
          capA[indT4, vA0123] <- - 2*t[j3]
          capA[indT4, vA0133] <- -t[j2]
          capA[indT4, vA0233] <- - t[j1]
          capA[indT4, vA1233] <- 1
        }
        for (j4 in j3:p){
          vA0004 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(0, 0, 0, j4)))
            )
          vA0014 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(0, 0, j1, j4)))
            )
          vA0024 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(0, 0, j2, j4)))
            )
          vA0024 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(0, 0, j3, j4)))
            )
          vA0034 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(0, 0, j3, j4)))
            )
          vA0124 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(0, j1, j2, j4)))
            )
          vA0134 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(0, j1, j3, j4)))
            )
          vA0234 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(0, j2, j3, j4)))
            )
          vA1234 <- which(
            apply(pos.matrix, 1, function(r) all(r == c(j1, j2, j3, j4)))
            )
          # When j1 < j2 < j3 < j4
          if (all(!duplicated(c(j1, j2, j3, j4)))){
            indT4 <- which(
              apply(pos.matrixT4, 1, function(r) all(r == c(j1, j2, j3, j4)))
              )
            capA[indT4, 1] <- t[j1]*t[j2]*t[j3]*t[j4]
            capA[indT4, vA0001] <- -t[j2]*t[j3]*t[j4]
            capA[indT4, vA0002] <- -t[j1]*t[j3]*t[j4]
            capA[indT4, vA0003] <- -t[j1]*t[j2]*t[j4]
            capA[indT4, vA0004] <- -t[j1]*t[j2]*t[j3]

            capA[indT4, vA0012] <- t[j3]*t[j4]
            capA[indT4, vA0013] <- t[j2]*t[j4]
            capA[indT4, vA0014] <- t[j2]*t[j3]
            capA[indT4, vA0023] <- t[j1]*t[j4]
            capA[indT4, vA0024] <- t[j1]*t[j3]
            capA[indT4, vA0034] <- t[j1]*t[j2]

            capA[indT4, vA0123] <- -t[j4]
            capA[indT4, vA0124] <- -t[j3]
            capA[indT4, vA0134] <- -t[j2]
            capA[indT4, vA0234] <- -t[j1]

            capA[indT4, vA1234] <- 1

          }
        }
      }
    }
  }
  capA <- capA * exp(-sum(t*t)/2)
  return(capA)
}


