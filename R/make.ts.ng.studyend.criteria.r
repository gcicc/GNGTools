#' @title Create lookup table for study-end
#'
#' @param mu.0.c prior parameter, mean, control
#' @param alpha.0.c prior parameter, precision, alpha, control
#' @param beta.0.c prior parameter, precision, beta, control
#' @param n.0.c prior parameter, effective sample size, control
#' @param mu.0.t prior parameter, mean, treatment
#' @param alpha.0.t prior parameter, precision, alpha, treatment
#' @param beta.0.t prior parameter, precision, beta, treatment
#' @param n.0.t prior parameter, effective sample size, treatment
#' @param xbar.c.LB lowerbound for xbar on control
#' @param xbar.c.UB upperbound for xbar on control
#' @param npoints number of points
#' @param s.c sample sd, control
#' @param n.c sample size, control
#' @param xbar.t sample mean - treatment
#' @param s.t sample sd, treatment
#' @param n.t sample size, treatment
#' @param Delta.lrv TTP info
#' @param Delta.tv TTP info
#' @param tau.ng thresholds
#' @param tau.lrv thresholds
#' @param tau.tv thresholds
#' @param n.MC MC size
#' @param go.parallel logical to engage parallel computing
#' @param cl cluster
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' make.ts.ng.studyend.criteria(go.parallel=FALSE)
#' }
make.ts.ng.studyend.criteria <- function(
  mu.0.c = 0, alpha.0.c = .25, beta.0.c = 1, n.0.c = 1,
  mu.0.t = 0, alpha.0.t = .25, beta.0.t = 1, n.0.t = 1,
  xbar.c.LB = 0, xbar.c.UB = 5,
  npoints=6, s.c = 7, n.c = 55,
  xbar.t = 0, s.t = 7, n.t = 55,
  Delta.lrv = 5.5, Delta.tv = 6.5,
  tau.ng = .65, tau.lrv = .8, tau.tv = .1,n.MC = 250, go.parallel=TRUE, cl=cl){


my.df <- data.frame(mu.0.c = mu.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c, n.0.c = n.0.c,
                    mu.0.t = mu.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                    xbar.c = seq(from = xbar.c.LB, to = xbar.c.UB, length.out=npoints), s.c = s.c, n.c = n.c,
                    xbar.t = xbar.t, s.t = s.t, n.t = n.t,
                    Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                    tau.ng = tau.ng, tau.lrv = tau.lrv, tau.tv = tau.tv,
                    n.MC = n.MC)


  if(go.parallel == FALSE){
    my.df <- data.frame(mu.0.c = mu.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c, n.0.c = n.0.c,
                        mu.0.t = mu.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                        xbar.c = seq(xbar.c.LB, xbar.c.UB, length.out=npoints), s.c = s.c, n.c = n.c,
                        xbar.t = xbar.t, s.t = s.t, n.t = n.t,
                        Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                        tau.ng = tau.ng, tau.lrv = tau.lrv, tau.tv = tau.tv,
                        n.MC = n.MC)
  return(
    bind_rows(
      apply(matrix(1:nrow(my.df)), 1, function(x) {
        get.ts.ng.studyend.GNG(mu.0.c = my.df$mu.0.c[x], alpha.0.c = my.df$alpha.0.c[x], beta.0.c = my.df$beta.0.c[x], n.0.c = my.df$n.0.c[x],
                               mu.0.t = my.df$mu.0.t[x], alpha.0.t = my.df$alpha.0.t[x], beta.0.t = my.df$beta.0.t[x], n.0.t = my.df$n.0.t[x],
                               xbar.c = my.df$xbar.c[x], s.c = my.df$s.c[x], n.c = my.df$n.c[x],
                               xbar.t = my.df$xbar.t[x], s.t = my.df$s.t[x], n.t = my.df$n.t[x],
                               Delta.lrv = my.df$Delta.lrv[x], Delta.tv = my.df$Delta.tv[x],
                               tau.ng = my.df$tau.ng[x], tau.lrv = my.df$tau.lrv[x], tau.tv = my.df$tau.tv[x],
                               n.MC = my.df$n.MC[x])[[3]]
      }))
  )
  }

  if(go.parallel == TRUE){
    my.df <- data.frame(mu.0.c = mu.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c, n.0.c = n.0.c,
                        mu.0.t = mu.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                        xbar.c = seq(xbar.c.LB, xbar.c.UB, length.out=npoints), s.c = s.c, n.c = n.c,
                        xbar.t = xbar.t, s.t = s.t, n.t = n.t,
                        Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                        tau.ng = tau.ng, tau.lrv = tau.lrv, tau.tv = tau.tv,
                        n.MC = n.MC)
    clusterEvalQ(cl = cl, expr = {requireNamespace("dplyr"); requireNamespace("tidyr")})
    clusterExport(cl = cl, varlist = c("get.ts.ng.studyend.GNG",  "get.ts.ng.mc.df","get.ts.ng.mc","get.ng.post","rt_ls", "my.df"), envir = environment())
    return(
      bind_rows(
        parApply(cl = cl, matrix(1:nrow(my.df)), 1, function(x) {
          return(
            get.ts.ng.studyend.GNG(mu.0.c = my.df$mu.0.c[x], alpha.0.c = my.df$alpha.0.c[x], beta.0.c = my.df$beta.0.c[x], n.0.c = my.df$n.0.c[x],
                                   mu.0.t = my.df$mu.0.t[x], alpha.0.t = my.df$alpha.0.t[x], beta.0.t = my.df$beta.0.t[x], n.0.t = my.df$n.0.t[x],
                                   xbar.c = my.df$xbar.c[x], s.c = my.df$s.c[x], n.c = my.df$n.c[x],
                                   xbar.t = my.df$xbar.t[x], s.t = my.df$s.t[x], n.t = my.df$n.t[x],
                                   Delta.lrv = my.df$Delta.lrv[x], Delta.tv = my.df$Delta.tv[x],
                                   tau.ng = my.df$tau.ng[x], tau.lrv = my.df$tau.lrv[x], tau.tv = my.df$tau.tv[x],
                                   n.MC = my.df$n.MC[x])[[3]]
          )
        })
        )
      )
    }
}
