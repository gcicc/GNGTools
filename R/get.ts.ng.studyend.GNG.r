#' @title Get two-sample normal gamma study end GNG
#'
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size for control group
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size for treatment group
#' @param xbar.c sample mean for control group
#' @param s.c sample sd for control group
#' @param n.c sample size for control group
#' @param xbar.t sample mean for treatment group
#' @param s.t sample sd for treatment group
#' @param n.t sample size of treatment group
#' @param Delta.lrv Min TPP
#' @param Delta.tv Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param n.MC n for MC sampling
#'
#' @return A list is returned holding what is needed of data to achieve Go/No-Go
#' @export
#'
#' @examples
#' my.ts.ng.studyend.GNG <- get.ts.ng.studyend.GNG()
#' my.ts.ng.studyend.GNG
#' @author Greg Cicconetti
get.ts.ng.studyend.GNG <- function(mu.0.c = 0, alpha.0.c = .25, beta.0.c = 1, n.0.c = 1,
                                   mu.0.t = 0, alpha.0.t = .25, beta.0.t = 1, n.0.t = 1,
                                   xbar.c = 1.5, s.c = 4, n.c = 40,
                                   xbar.t = 26, s.t = 4, n.t = 40,
                                   Delta.lrv = 1, Delta.tv = 1.5,
                                   tau.ng = .65, tau.lrv = .8, tau.tv = .1,
                                   n.MC = 1000){

  # So xbar.t arguement is not really required!!
  my.means <- seq(from = -Delta.tv*4, to = Delta.tv*4, length.out=50)
  stage1 <- get.ts.ng.mc.df(mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c,
                            xbar.c = xbar.c, s.c = s.c, n.c = n.c, group.c = "Control",
                            mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t,
                            xbar.t = xbar.c + my.means, s.t = s.t, n.t = n.t,
                            group.t="treat",
                            Delta.tv = Delta.tv,Delta.lrv = Delta.lrv,tau.tv = tau.tv,
                            tau.lrv = tau.lrv,tau.ng = tau.ng, n.MC = n.MC)
  stage1.go <- stage1 %>% dplyr::filter(result== "Go") %>% dplyr::slice(1)
  # Determine max number of TRT responders for No-Go
  stage1.ng <- stage1 %>% dplyr::filter(result== "No-Go") %>% dplyr::slice(n())

  my.means <- c(seq(stage1.go$xbar.t - stage1.go$s.t/4, stage1.go$xbar.t, length.out=100),
                seq(stage1.ng$xbar.t , stage1.ng$xbar.t + stage1.ng$s.t/4, length.out=100))
  stage2 <- get.ts.ng.mc.df(mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c = alpha.0.c,
                            beta.0.c = beta.0.c,
                            xbar.c = xbar.c, s.c = s.c, n.c = n.c, group.c = "Control",
                            mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t,
                            beta.0.t = beta.0.t,
                            xbar.t = my.means, s.t = s.t, n.t = n.t, group.t="treat",
                            Delta.tv = Delta.tv,Delta.lrv = Delta.lrv,tau.tv = tau.tv,
                            tau.lrv = tau.lrv,tau.ng = tau.ng,
                            n.MC = n.MC)
  result.go <- stage2 %>% dplyr::filter(result== "Go") %>% dplyr::slice(1)
  result.ng <- stage2 %>% dplyr::filter(result== "No-Go") %>% dplyr::slice(dplyr::n())

  return(list(result.go = result.go, result.ng = result.ng, bind_rows(result.go, result.ng)))
}

