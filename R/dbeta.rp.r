#' @name Reparameterized.beta
#' @title Reparameterized Beta distribution functions
#' @param mean mean
#' @param effective.ss effective sample size
#' @param ncp non-centrality parameter
#' @param log as in stats::*beta
#' @param x likelihood function evaluates at point x
#' @param q quantile
#' @param p cumulative probability
#' @param lower.tail logical; if TRUE (default), probabilities are P(X < x), otherwise, P(X>x).
#' @param log.p logical; if TRUE, probabilities p are given as log(p)
#' @return Returns density values, cumulative probabilities, quantiles and random samples from a Beta distribution.
#' @examples
#' dbeta.rp(.5, .5, 1)
#' pbeta.rp(.5, .5, 1)
#' qbeta.rp(.975, .5, 1)
#' @description Reparameterized Beta distribution functions
dbeta.rp <- function(x, mean = .5, effective.ss = 1, ncp = 0, log = FALSE){
  dbeta(x = x,
        shape1 = mean * effective.ss,
        shape2 = effective.ss - mean * effective.ss,
        ncp = ncp, log = log)
}

#' @rdname Reparameterized.beta
pbeta.rp <- function(q, mean = .5, effective.ss = 1, ncp = 0, lower.tail = TRUE,
                     log.p = FALSE){
  pbeta(q = q,
        shape1 = mean * effective.ss,
        shape2 = effective.ss - mean * effective.ss,
        ncp = ncp, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname Reparameterized.beta
qbeta.rp <- function(p, mean = .5, effective.ss = 1, ncp = 0, lower.tail = TRUE,
                     log.p = FALSE){
  qbeta(p = p,
        shape1 = mean * effective.ss,
        shape2 = effective.ss - mean * effective.ss,
        ncp = ncp, lower.tail = lower.tail, log.p = log.p)
}


# dbeta.rp(.5, .5, 1)
# pbeta.rp(.5, .5, 1)
# qbeta.rp(.975, .5, 1)
