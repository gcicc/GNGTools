#' @name Reparameterized.normal
#' @title Reparameterized normal functions
#' @param mean mean
#' @param tau precision parameter
#' @param x likelihood function evaluates at point x
#' @param q quantile
#' @param p cumulative probability
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param log logical; if TRUE, probabilities p are given as log(p)
#' @param log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are P(X < x) otherwise, P(X>x).
#' @return Returns density values, cumulative probabilities, quantiles and random samples from a normal distribution.
#' @description Reparameterized normal functions
#' @examples
#' dnorm.rp(x=0)
#' pnorm.rp(q=.1)
#' qnorm.rp(p=.975)
#' rnorm.rp(10)

dnorm.rp <- function(x, mean = 0, tau = 1, log = FALSE){
  dnorm(x, mean = mean, sd = 1 / sqrt(tau), log = log)
}

#' @rdname Reparameterized.normal
pnorm.rp <- function(q, mean = 0, tau = 1, lower.tail = TRUE, log.p = FALSE){
  pnorm(q, mean = mean, sd = 1 / sqrt(tau), lower.tail = lower.tail, log.p = log.p)}


#' @rdname Reparameterized.normal
qnorm.rp <- function(p, mean = 0, tau = 1, lower.tail = TRUE, log.p = FALSE){
  qnorm(p, mean = mean, sd = 1 / sqrt(tau), lower.tail = lower.tail, log.p = log.p)}

#' @rdname Reparameterized.normal
rnorm.rp <- function(n, mean = 0, tau = 1){
  rnorm(n, mean = mean, sd = 1 / sqrt(tau))}


