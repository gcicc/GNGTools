#' @title Probability for difference of x or more between two beta distributions
#'
#' @param x value
#' @param a1 alpha for first
#' @param b1 beta for first
#' @param a2 alpha for second
#' @param b2 beta for second
#'
#' @return returns the probability of a difference of x or more between two beta distributions.
#' This is a simplification of the original P2beta function simplified for vectorization.
#' Therre is a lot of space here for improved run times.
#' @export
#'
#' @author Randall Henner
p2beta_diff <- function( x, a1, b1, a2, b2 ){

  if (x < -1) {
    cdf <- 0
  }
  else if ((x > -1) & (x <= 0)) {
    cdf <- integrate(function(t) pbeta(x+t, a2, b2)*dbeta(t, a1, b1), -x, 1)$value
  }
  else if ((x > 0) & (x <= 1)) {
    cdf <- integrate(function(t) pbeta(x+t, a2, b2)*dbeta(t, a1, b1), 0, 1 - x)$value +
      integrate(function(t) dbeta(t, a1, b1), 1 - x, 1)$value
  }
  else if (x > 1) {
    cdf <- 1
  }

  return(cdf)
}

#' @title Probability for difference of x or more between two beta distributions
#' @param x value
#' @param a1 alpha for first
#' @param b1 beta for first
#' @param a2 alpha for second
#' @param b2 beta for second
#' @return a function is returned
p2beta_diff_Vector <- Vectorize(p2beta_diff)
