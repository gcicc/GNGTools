#' @title Get two-sample normal gamma MC sampling
#'
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size for control group
#' @param xbar.c sample mean for control group
#' @param s.c sample sd for control group
#' @param n.c sample size for control group
#' @param group.c group label for control group
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size for treatment group
#' @param xbar.t sample mean for treatment group
#' @param s.t sample sd for treatment group
#' @param n.t sample size for treatment group
#' @param group.t group label for treatment group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param seed random seed
#' @param n.MC n for MC sampling
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.ts.ng.mc <- get.ts.ng.mc()
#' my.ts.ng.mc
get.ts.ng.mc <- function(mu.0.c = 0, n.0.c = .0001, alpha.0.c=.25 , beta.0.c = 1 ,
                         xbar.c = 0, s.c = 2.3, n.c = 20000, group.c="Control",
                         mu.0.t = 1.5, n.0.t = .0001, alpha.0.t=.25,                 beta.0.t = 1 ,
                         xbar.t =1.5, s.t = 2.3, n.t = 20000, group.t="Treatment",
                         Delta.tv = 1.5, Delta.lrv = 1, tau.tv = 0,
                         tau.lrv = .8, tau.ng = .75,
                         seed = 1234, n.MC = 5000)
{
  # Compute individual posterior distribution parameters given priors and data
  CON.results <- get.ng.post(mu.0 = mu.0.c, n.0 = n.0.c, alpha.0 = alpha.0.c,
                             beta.0 = beta.0.c, xbar = xbar.c, s = s.c,
                             n = n.c, group = group.c)
  TRT.results <- get.ng.post(mu.0 = mu.0.t, n.0 = n.0.t, alpha.0 = alpha.0.t,
                             beta.0 = beta.0.t, xbar = xbar.t, s = s.t,
                             n = n.t, group = group.t)

  # if(is.null(seed) ==F) set.seed(seed = seed)

  # Estimating prior probabilities and posterior probabilities
  Z.0 <- rt_ls(n = n.MC, df = TRT.results$tdf.0, mu = TRT.results$mu.0,
               sigma = TRT.results$sigma.0) - rt_ls(n = n.MC, df = CON.results$tdf.0,
                                                    mu = CON.results$mu.0,
                                                    sigma = CON.results$sigma.0)
  # sampling from each posterior and taking differences
  Z.n <- rt_ls(n = n.MC, df = TRT.results$tdf.n, mu = TRT.results$mu.n,
               sigma = TRT.results$sigma.n) -
    rt_ls(n = n.MC, df = CON.results$tdf.n, mu = CON.results$mu.n,
          sigma = CON.results$sigma.n)
 my.df <- data.frame(P.R1 = mean(Z.n >= Delta.lrv),
                      P.R3 = mean(Z.n >= Delta.tv)) %>%
    mutate(result = ifelse(P.R1 > tau.lrv & P.R3 > tau.tv, "Go",
                           ifelse(P.R3 <= tau.tv & P.R1 <= tau.ng, "No-Go", "Consider")))

  return(my.df)
}

