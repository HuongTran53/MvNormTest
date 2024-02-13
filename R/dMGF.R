# -- Filename: hCGF.R ---#
#' @title  Moment generating functions (MGF) of standard normal distribution.
#' @name dMGF
#' @description
#' Get the polynomial term in the expression of derivatives of  moment generating function of \eqn{N_p(0, I_p)}, with
#' respect to a given component and its exponent. Up to eighth order.
#'
#' @param tab dataframe whose whose columns are the components and the order of derivatives.
#' @param t vector in \eqn{\mathbb{R}^p}.
#' @param coef get only polynomial or whole expression by multiplying the polynomial term with the exponent term \eqn{\exp(.5 t't)}.
#'
#' @return Value of derivatives.
#'
#' @details
#' For a standard multivariate normal random variables \eqn{Y \sim N_p(0, I_p)}
#' \deqn{\mathbb{E}\left(Y_1^{k_1} ... Y_p^{k_p} \exp(t'X)\right) = \dfrac{\partial^{k_1}\dots
#'  \partial^{k_p}}{t_1^{k_1} \dots t_p^{k_p}} \exp(t't/2) =
#'  \mu^{(k_1)} (t_1) ... \mu^{(k_p)}(t_p) \exp(t't/2)}
#' For example,
#' \eqn{\mathbb{E}Y_2^4 \exp(t'Y) = \dfrac{\partial^4}{\partial t_2^4} \exp(t't/2)
#' = \mu^{(4)}(t_2) \exp(t't/2).}
#'
#'
#' @examples
#' #Calculation of above example
#' t <- rep(.2, 7)
#' tab <- data.frame(j = 2, exponent = 4)
#' dMGF(tab, t = t)
#' dMGF(tab, t = t, coef = F)
#'
#' @export
dMGF <- function(tab, t, coef = T){ #upto j.num == 8
  j <- as.numeric(tab[1])
  # j <- as.numeric(tab$Var1)
  j.num <- tab[2]
  # j.num <- as.numeric(tab$Freq)

  if (j == 0){
    mu <- 1
  } else {
    s <- t[j]
    if (j.num == 1){
      mu <- s
    } else if (j.num == 2) {
      mu <- (1 + s^2)
    } else if (j.num == 3){
      mu <- (3*s + s^3)
    } else if (j.num == 4){
      mu <- (3 + 6*s^2 + s^4)
    } else if (j.num == 5){
      mu <- (15*s + 10*s^3 + s^5)
    } else if (j.num == 6) {
      mu <- (15 + 45*s^2 + 15*s^4 + s^6)
    } else if (j.num == 7){
      mu <- 105*s + 105*s^3 + 21*s^5 + s^7
    } else {
      mu <- 105 + 420*s^2 + 210 * s^4 + 28*s^6 + s^8
    }
  }

  mu <- ifelse(coef, mu, mu * exp(sum(t*t)/2))

  return(mu)
}

