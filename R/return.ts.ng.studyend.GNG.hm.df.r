#' @title Return two-sample normal-gamma study end GNG heatmap data.frame
#'
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size parameter for control group
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size parameter for treatment group
#' @param xbar.c.low lower bound for control group sample mean grid
#' @param xbar.c.high upper bound for control group sample mean grid
#' @param s.c control standard deviation assumed
#' @param n.c control sample size
#' @param s.t treatment standard deviation assumed
#' @param n.t treatment sample size
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param n.MC n for MC sampling
#' @param x_ng xbar needed for no-go; leave null for standard rule
#' @param x_go xbar needed for go; leave null for standard rule
#' @param npoints number of points to use in simulation
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.ng.studyend.GNG.hm.df <- return.ts.ng.studyend.GNG.hm.df()
#' head(my.ts.ng.studyend.GNG.hm.df)
#' }
#' @author Greg Cicconetti

return.ts.ng.studyend.GNG.hm.df <- function(mu.0.c = 0, alpha.0.c = .25, beta.0.c = 1, n.0.c = 1,
                                            mu.0.t = 0, alpha.0.t = .25, beta.0.t = 1, n.0.t = 1,
                                            xbar.c.low = -1, xbar.c.high = 1, s.c = 4, n.c = 40,
                                            s.t = 4, n.t = 40,
                                            npoints = 15,
                                            Delta.lrv = 1, Delta.tv = 1.5,
                                            tau.ng = .65, tau.lrv = .8, tau.tv = .1,
                                            n.MC = 1000,
                                            x_ng = NULL, x_go=NULL
                                            ){

  my.grid <- expand.grid(x.bar.c = seq(xbar.c.low, xbar.c.high, length.out = npoints))
  bind_rows(
    apply(X = matrix(1:nrow(my.grid)), MARGIN = 1, FUN=function(x) {

      results <- get.ts.ng.studyend.GNG(mu.0.c = mu.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c, n.0.c = n.0.c,
                                         mu.0.t = mu.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                                         xbar.c = my.grid$x.bar.c[x], s.c = s.c, n.c = n.c,
                                         xbar.t = my.grid$x.bar.c[x], s.t = s.t, n.t = n.t,
                                         Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                         tau.ng = tau.ng, tau.lrv = tau.lrv, tau.tv = tau.tv,
                                         n.MC = n.MC)
      # If we are in a case where either go or no/go is not possible, fill with missing

      results$result.go <- results$result.go %>% dplyr::rename(xbar.go=xbar.t, P.R1.Go = P.R1, P.R3.Go=P.R3, result.Go=result )
      results$result.ng <- results$result.ng %>% dplyr::rename(xbar.ng=xbar.t, P.R1.NG = P.R1, P.R3.NG=P.R3, result.NG=result )
      results$result.go %>% left_join(results$result.ng)
    } )  )

}

