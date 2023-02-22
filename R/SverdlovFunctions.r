# SverdlovFunctions ------
## SUPPLEMENT for the paper "Exact Bayesian Inference Comparing Binomial Proportions,
##                           with Application to Proof of Concept Clinical Trials"
## R CODE
## Probability Density Function (PDF)
## Cumulative Distribution Function (CDF)
## Quantiles of the Distribution
## Random Generator
#  of the RELATION of two independent betas
#  RELATION means:
# -- ('DIFF') risk difference (X2 - X1)
# -- ('RR') relative risk (X2/X1)
# -- ('OR') odds ratio ((X2/(1 - X2))/(X1/(1 - X1)))

## Auxiliary functions used to calculate distribution functions and options

# 1) Appell's hypergeometric function:
#     appel.hypgeom(a, b1, b2, c, x, y) = integral from 0 to 1 of f(a, b1, b2, c, x, y) by dt
# where f(a, b1, b2, c, x, y) = 1/beta(a,c-a)*t^(a-1)*(1-t)^(c-a-1)*(1-t*x)^(-b1)*(1-t*y)^(-b2)

#' @name SverdlovFunction
#' @title Sverdlov's functions
#' @param a a
#' @param b b
#' @param a1 a1
#' @param b1 b1
#' @param a2 a2
#' @param b2 b2
#' @param c c
#' @param x x
#' @param y y
#' @param t t
#' @param fun fun
#' @param left left
#' @param right right
#' @param tol tol
#' @param relation relation
#' @param approach approach
#' @param method method
#' @param n n
#' @param left0 left0
#' @param right0 right0
#' @param alpha alpha
#' @description Sverdlov's f function
#' @return Various functions are offered.
#' @references Sverdlov O, Ryeznik Y, Wu S. Exact Bayesian Inference Comparing Binomial Proportions, With Application to Proof-of-Concept Clinical Trials. Therapeutic Innovation & Regulatory Science. 2015;49(1):163-174. doi:10.1177/2168479014547420
f <- function(a, b1, b2, c, x, y, t){
  return(gamma(c)/gamma(a)/gamma(c-a)*t^(a-1)*(1-t)^(c-a-1)*(1-t*x)^(-b1)*(1-t*y)^(-b2))
}

#' @rdname SverdlovFunction
appel.hypgeom <- function(a, b1, b2, c, x, y){
  res <- rep(0, length(x))
  for (r in 1:length(res)){
    res[r] <- integrate(function(t) f(a, b1, b2, c, x[r], y[r], t), 0, 1)$value
  }
  return(res)
}

#' @rdname SverdlovFunction
g <- function(a, b, c, x, t){
  return(1/beta(a, c-a)*t^(a-1)*(1-t)^(c-a-1)*(1-t*x)^(-b))
}

#' @rdname SverdlovFunction
gauss.hypgeom <- function(a, b, c, x){
  res <- rep(0, length(x))
  for (r in 1:length(res)){
    res[r] <- integrate(function(t) g(a, b, c, x[r], t), 0, 1)$value
  }
  return(res)
}

#' @rdname SverdlovFunction
root <- function(fun, left, right, tol){
  rt <- 0.5*(left + right)
  while (abs(fun(rt)) > tol){
    left <- rt*(fun(left)*fun(rt) > 0) + left*(fun(left)*fun(rt) < 0)
    right <- rt*(fun(right)*fun(rt) > 0) + right*(fun(right)*fun(rt) < 0)
    rt <- 0.5*(left + right)
  }
  return(rt)
}

#' @rdname SverdlovFunction
r2beta <- function(relation, n, a1, b1, a2, b2){

  X1 <- rbeta(n, shape1 = a1, shape2 = b1)
  X2 <- rbeta(n, shape1 = a2, shape2 = b2)
  # difference
  if (relation == 'DIFF'){
    rnd <- X2 - X1
  }
  # relative risk
  else if (relation == 'RR'){
    rnd <- X2/X1
  }
  # odds ratio
  else if (relation == 'OR'){
    rnd <- (X2/(1-X2))/(X1/(1-X1))
  }
  return(rnd)
}

#' @rdname SverdlovFunction
d2beta <- function(relation, x, a1, b1, a2, b2){
  pdf <- rep(0, length(x))
  # difference
  if (relation == 'DIFF'){
    x.neg <- x[x <= 0]
    if (length(x.neg) != 0){
      pdf[x <= 0] <- beta(a2,b1)/beta(a1,b1)/beta(a2,b2)*(-x.neg)^(b1+b2-1)*
        (1+x.neg)^(a2+b1-1)*
        appel.hypgeom(b1, a1+a2+b1+b2-2, 1-a1, a2+b1, 1+x.neg, 1-x.neg^2)
    }
    x.pos <- x[x > 0]
    if (length(x.pos) != 0){
      pdf[x > 0] <- beta(a1,b2)/beta(a1,b1)/beta(a2,b2)*(x.pos)^(b1+b2-1)*
        (1-x.pos)^(a1+b2-1)*
        appel.hypgeom(b2, a1+a2+b1+b2-2, 1-a2, a1+b2, 1-x.pos, 1-x.pos^2)
    }
  }
  # relative risk
  else if (relation == 'RR'){
    x.01 <- x[x <= 1]
    if (length(x.01) != 0){
      pdf[x <= 1] <- beta(a1+a2,b1)/beta(a1,b1)/beta(a2,b2)*(x.01)^(a2-1)*
        gauss.hypgeom(a1+a2, 1-b2, a1+a2+b1, x.01)
    }
    x.1 <- x[x > 1]
    if (length(x.1) != 0){
      pdf[x > 1] <- beta(a1+a2,b2)/beta(a1,b1)/beta(a2,b2)*(x.1)^(-a1-1)*
        gauss.hypgeom(a1+a2, 1-b1, a1+a2+b2, 1/x.1)
    }
  }
  # odds ration
  else if (relation == 'OR'){
    x.01 <- x[x <= 1]
    if (length(x.01) != 0){
      pdf[x <= 1] <- beta(a1+a2,b1+b2)/beta(a1,b1)/beta(a2,b2)*(x.01)^(a2-1)*
        gauss.hypgeom(a2+b2, a1+a2, a1+b1+a2+b2, 1-x.01)
    }
    x.1 <- x[x > 1]
    if (length(x.1) != 0){
      pdf[x > 1] <- beta(a1+a2,b1+b2)/beta(a1,b1)/beta(a2,b2)*(x.1)^(-b2-1)*
        gauss.hypgeom(a2+b2, b1+b2, a1+b1+a2+b2, 1-1/x.1)
    }
  }
  return(pdf)
}

#' @rdname SverdlovFunction
p2beta <- function(relation, approach, x, a1, b1, a2, b2, n = 1000000){
  cdf <- rep(0, length(x))
  for (r in 1:length(cdf)){
    if (approach == 'SIMULATION'){
      rnd <- r2beta(relation, n, a1, b1, a2, b2)
      cdf[r] <- sum(rnd < x[r])/length(rnd)
    }
    else if (approach == 'DIRECT'){
      # difference
      if (relation == 'DIFF'){
        if (x[r] < -1) {
          cdf[r] <- 0
        }
        else if ((x[r] > -1) & (x[r] <= 0)) {
          cdf[r] <- integrate(function(t) pbeta(x[r]+t, a2, b2)*dbeta(t, a1, b1), -x[r], 1)$value
        }
        else if ((x[r] > 0) & (x[r] <= 1)) {
          cdf[r] <- integrate(function(t) pbeta(x[r]+t, a2, b2)*dbeta(t, a1, b1), 0, 1 - x[r])$value +
            integrate(function(t) dbeta(t, a1, b1), 1 - x[r], 1)$value
        }
        else if (x[r] > 1) {
          cdf[r] <- 1
        }
      }
      # relative risk
      else if (relation == 'RR'){
        if (x[r] < 0) {
          cdf[r] <- 0
        }
        else if ((x[r] >= 0) & (x[r] <= 1)) {
          cdf[r] <- integrate(function(t) pbeta(x[r]*t, a2, b2)*dbeta(t, a1, b1), 0, 1)$value
        }
        else if (x[r] > 1) {
          cdf[r] <- integrate(function(t) pbeta(x[r]*t, a2, b2)*dbeta(t, a1, b1), 0, 1/x[r])$value +
            integrate(function(t) dbeta(t, a1, b1), 1/x[r], 1)$value
        }
      }
      # odds ratio
      else if (relation == 'OR'){
        if (x[r] < 0) {
          cdf[r] <- 0
        }
        else {
          cdf[r] <- integrate(function(t) pf(a1*b2/a2/b1*x[r]*t, 2*a2, 2*b2)*df(t, 2*a1, 2*b1), 0, Inf)$value
        }
      }
    }
  }
  return(cdf)
}

#' @rdname SverdlovFunction
q2beta <- function(relation, a1, b1, a2, b2, alpha, tol = 10^(-5)){
  quantile <- rep(0, length(alpha))
  for (r in 1:length(quantile)){
    if (relation == 'DIFF'){
      fun1 <- function(t, ...){p2beta(relation, approach = 'DIRECT', t, a1, b1, a2, b2) - alpha[r]}
      quantile[r] <- root(fun1, -0.9999, 0.9999, tol)
    }
    else {
      fun2 <- function(x, ...){p2beta(relation, approach = 'DIRECT', tan(pi*x/2), a1, b1, a2, b2) - alpha[r]}
      x <- root(fun2, -0.9999, 0.9999, tol)
      quantile[r] <- tan(pi*x/2)
    }

  }
  return(quantile)
}

## Credible interval based on Nelder-Mead algorithm or inverse CDF approach
#' @rdname SverdlovFunction
ci2beta <- function(relation, method, a1, b1, a2, b2, alpha, left0, right0){
  if (method == 'neldermead'){
    optim.fun <- function(x) {
      abs(p2beta(relation, approach = 'DIRECT', x[2], a1, b1, a2, b2) -
            p2beta(relation, approach = 'DIRECT', x[1], a1, b1, a2, b2) -
            1 + alpha) +
        abs(d2beta(relation, x[2], a1, b1, a2, b2) -
              d2beta(relation, x[1], a1, b1, a2, b2))
    }
    ci <- optim(c(left0, right0), optim.fun)$par
  }
  else if (method == 'inv.cdf'){
    ci <- q2beta(relation, a1, b1, a2, b2, c(alpha/2, 1 - alpha/2))
  }
  return(ci)
}
