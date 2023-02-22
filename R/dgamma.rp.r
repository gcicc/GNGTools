#' @name Reparameterized.gamma
#' @title Reparameterized Gamma distribution functions
#' @param mean mean
#' @param var variance
#' @param log as in stats::*gamma
#' @param x likelihood function evaluates at point x
#' @param q quantile
#' @param p cumulative probability
#' @param n sample size
#' @param lower.tail logical; if TRUE (default), probabilities are P(X < x), otherwise, P(X>x).
#' @param log.p logical; if TRUE, probabilities p are given as log(p)
#' @return Returns density values, cumulative probabilities, quantiles and random samples from a gamma distribution.
#' @examples
#' dgamma.rp(1, 1, 1)
#' pgamma.rp(1.96, 1, 1)
#' qgamma.rp(.975, 1, 1)
#' rgamma.rp(10, 1, 1)
#' @description Reparameterized Gamma distribution functions
dgamma.rp <- function(x, mean, var, log = FALSE){
  dgamma(x, shape = mean * mean / var, rate = mean / var)}

#' @rdname Reparameterized.gamma
pgamma.rp <- function(q, mean, var, lower.tail = TRUE, log.p = FALSE){
  pgamma(q, shape = mean * mean / var, rate = mean / var, lower.tail = lower.tail,
         log.p = log.p)}

#' @rdname Reparameterized.gamma
qgamma.rp <- function(p, mean, var, lower.tail = TRUE, log.p = FALSE){
  qgamma(p, shape = mean * mean / var, rate = mean / var, lower.tail = lower.tail,
         log.p = log.p)}

#' @rdname Reparameterized.gamma
rgamma.rp <- function(n, mean, var){
  rgamma(n, shape = mean * mean / var, rate = mean / var)}



