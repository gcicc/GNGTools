#' @title Get single sample normal-gamma treament OC data.frame
#'
#' @param mu.0.t prior mean
#' @param n.0.t prior effective sample size
#' @param alpha.0.t prior alpha parameter
#' @param beta.0.t prior beta parameter
#' @param s.t sample sd for treatment group
#' @param n.t sample size for treatment group
#' @param from.here treatment effect lower bound
#' @param to.here treatment effect upper bound
#' @param length.out number of points used to span (from.here, to.here)
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return Returns a data.frame ready to create a treatment effect OC curve
#' @export
#'
#' @examples
#' my.ss.ng.trt.oc.df <- get.ss.ng.trt.oc.df(mu.0.t = 0, n.0.t = 10, alpha.0.t = 0.25, beta.0.t = 1,
#' s.t = 2, n.t = 40, from.here = 0, to.here =4, length.out=20,
#' Delta.tv = 1.75, Delta.lrv = 1,
#' tau.tv = 0.1, tau.lrv = 0.8, tau.ng = 0.65)
#' head(my.ss.ng.trt.oc.df)
#' @author Greg Cicconetti
get.ss.ng.trt.oc.df <- function(mu.0.t = 0, n.0.t = 10, alpha.0.t = 0.25, beta.0.t = 1,
                                s.t = 2, n.t = 40, from.here = 0, to.here =4, length.out=20,
                                Delta.tv = 1.75, Delta.lrv = 1,
                                tau.tv = 0.1, tau.lrv = 0.8, tau.ng = 0.65)  {

  # Create a grid to pass to decision.ss.continuous
  # This grid really just runs over a sequence of values for underlying treatment effect

  results <- get.ss.ng.df(mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t,
                          beta.0.t = beta.0.t,
                          xbar.t = seq(from.here, to.here, length.out = length.out),
                          s.t = s.t, n.t = n.t,
                          Delta.tv = Delta.tv, Delta.lrv = Delta.lrv, tau.tv = tau.tv,
                          tau.lrv = tau.lrv, tau.ng = tau.ng)

  result.go <- results %>% dplyr::filter(result== "Go") %>% slice(1)
  # Determine max number of TRT responders for No-Go
  result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(n())


  my.df <- data.frame(xbar=seq(from.here, to.here, length.out = length.out)) %>%
    mutate(Go= 1 - pnorm(q=result.go$xbar, mean=xbar, sd=s.t/sqrt(n.t)),
           NoGo= pnorm(q = result.ng$xbar, mean=xbar, sd=s.t/sqrt(n.t))) %>%
    mutate(mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t,
           s.t = s.t, n.t = n.t,
           Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
           tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

  my.df
}
