#' @title Density function for the Normal-Gamma distribution
#' @name normal.gamma
#' @param mu likelihood evaluated when mean takes on value mu
#' @param tau likelihood evaluated when precision takes on value tau
#' @param mu0 hyperparameter describing mu
#' @param n0 hyperparameter describing effective sample size associated with mu0
#' @param a0 hyperparameter describing shape parameter of precision parameter
#' @param b0 hyperparameter describing rate parameter of precision parameter
#' @param mu.0 hyperparameter describing mu
#' @param n.0 hyperparameter describing effective sample size associated with mu0
#' @param alpha.0 hyperparameter describing shape parameter of precision parameter
#' @param beta.0 hyperparameter describing rate parameter of precision parameter
#' @param n number of observations
#' @return Returns the value of the Normal-gamma density function at the point passed.
#' @examples
#' dnorgam(100, .25, 0, 10, 0, .25)
#' rnormgam()
#' @description Density function for the Normal-Gamma distribution
dnorgam <- function(mu, tau, mu0, n0, a0, b0) {
  Zng <- gamma(a0)/b0 ^ (a0) * sqrt(2 * pi / n0)
  den <- 1 / Zng * tau ^ (a0 - 1 / 2) *exp(-tau / 2 * (n0*(mu - mu0) ^ 2 + 2 * b0))
  return(den)
}

#' @rdname normal.gamma
rnormgam <- function(n=100000, mu.0 = 0, n.0 = 1, alpha.0 = .01, beta.0 = .01){
  tau <- rgamma(n=n, shape = alpha.0, rate = beta.0)
  x <- rnorm(n=n, mean=mu.0, sd=(1/n.0*tau)^.5)
  return(x)
}
