#' @name Location.scale.t
#' @title Location-scale t distribution functions
#' @param mu mean
#' @param df degrees of freedom
#' @param sigma scale parameter
#' @param x likelihood function evaluates at point x
#' @param prob cumulative probability
#' @param n sample size
#' @return Returns density values, cumulative probabilities, quantiles and random samples from a Location-scale t distribution.
#' @author Greg Cicconetti
#' @examples
#' {
#' dt_ls(0, 100, 0, 1)
#' pt_ls(0, 100, 0, 1)
#' qt_ls(0.5, 100, 0, 1)
#' rt_ls(100, 100, 0, 1)
#' }
#' @description Location-scale t distribution functions
dt_ls <- function(x, df, mu, sigma) 1 / sigma * dt((x - mu)/sigma, df)

#' @rdname Location.scale.t
pt_ls <- function(x, df, mu, sigma) pt(q = (x - mu)/sigma, df)

#' @rdname Location.scale.t
qt_ls <- function(prob, df, mu, sigma) qt(p = prob, df)*sigma + mu

#' @rdname Location.scale.t
rt_ls <- function(n, df, mu, sigma) rt(n, df)*sigma + mu


