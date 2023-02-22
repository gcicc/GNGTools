#' Get single sample normal-gamma study end GNG
#' @title Get single sample normal-gamma study-end GNG
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size for treatment group
#' @param xbar.t sample mean for treatment group
#' @param s.t sample sd for treatment group
#' @param n.t sample size for treatment group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return returns a list of data.frames holding what is needed from data for study-end Go/No-Go
#' @export
#'
#' @examples
#' my.ss.mg.studyend.GNG <- get.ss.ng.studyend.GNG(mu.0.t = 0, alpha.0.t=.25, beta.0.t = 1,
#' n.0.t = 10, xbar.t = 1.97, s.t = 2, n.t = 20, Delta.lrv = 1.25, Delta.tv = 1.75,
#' tau.tv=.1, tau.lrv=.8, tau.ng=.65)
#' my.ss.mg.studyend.GNG

get.ss.ng.studyend.GNG <- function(mu.0.t = 0, alpha.0.t=.25, beta.0.t = 1, n.0.t = 10,
                                   xbar.t = 1.97, s.t = 2, n.t = 20,
                                   Delta.lrv = 1.25, Delta.tv = 1.75,
                                   tau.tv=.1, tau.lrv=.8, tau.ng=.65){

  # Run this for fixed values of
  results <- get.ss.ng.df(mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t,
                          beta.0.t = beta.0.t,
                          xbar.t = seq(Delta.lrv*.75, Delta.tv*1.5, length.out=1000),
                          s.t = s.t, n.t = n.t, Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                          tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

  result.go <- results %>% dplyr::filter(result== "Go") %>% dplyr::slice(1)
  result.ng <- results %>% dplyr::filter(result== "No-Go") %>% dplyr::slice(dplyr::n())
  return(list(result.go = result.go, result.ng = result.ng))
}
