#' @title Get two-sample normal-gamma treatment OC data.frame
#'
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size for control group
#' @param xbar.c sample mean for control group
#' @param s.c sample sd for control group
#' @param group.c label for control group
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size for treatment group
#' @param xbar.t sample mean for treatment group
#' @param s.t sample sd for treatment group
#' @param group.t label for treatment group
#' @param Delta.LB Lower bound for Delta
#' @param Delta.UB upper bound for delta
#' @param ARatio randomization ratio
#' @param N total sample size
#' @param Delta.tv Base TPP
#' @param Delta.lrv Min TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param npoints number of points
#' @param n.MC n for MC sampling
#' @param seed random seed
#' @param cl cluster
#' @param goparallel a logical to indicate if parallel computing is employed
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.ng.trt.oc.df <- get.ts.ng.trt.oc.df(goparallel=FALSE)
#' head(my.ts.ng.trt.oc.df)
#' }
#' @author Greg Cicconetti
get.ts.ng.trt.oc.df <- function(mu.0.c = 0, n.0.c = .0001, alpha.0.c=.25, beta.0.c = 1,
                                xbar.c = 0, s.c = 4,  group.c="Control",
                                mu.0.t = 0, n.0.t = .0001, alpha.0.t=.25, beta.0.t = 1,
                                xbar.t = 5, s.t = 4,  group.t="Treatment",
                                Delta.LB=0, Delta.UB=1.5, ARatio=1, N=50,
                                Delta.tv = 1.5, Delta.lrv = 1,
                                tau.tv = .1, tau.lrv = .65, tau.ng = 0,
                                npoints=2, n.MC = 500,
                                seed = 1234, cl=cl, goparallel=TRUE){
  message("\nstarting decision.ts.continuous\n")
  my.specs <- data.frame(treatment.effect=seq(Delta.LB, Delta.UB, length.out=npoints)) %>%
    mutate(mu.0.c = mu.0.c,
           n.0.c = n.0.c,
           alpha.0.c=alpha.0.c,
           beta.0.c = beta.0.c,
           xbar.c = xbar.c,
           s.c = s.c,
           n.c = N - floor(N*(ARatio/(1+ARatio))),
           mu.0.t = mu.0.t,
           n.0.t = n.0.t,
           alpha.0.t=alpha.0.t,
           beta.0.t = beta.0.t,
           s.t = s.t,
           n.t = floor(N*(ARatio/(1+ARatio))),
           Delta.tv = Delta.tv,
           Delta.lrv = Delta.lrv,
           tau.tv = tau.tv,
           tau.lrv = tau.lrv,
           tau.ng = tau.ng,
           n.MC = n.MC,
           npoints = npoints,
           seed = seed+1:n())

           # Fixed placebo data and all priors + observed standard deviation and sample size.

           my.func <- function(x){

             set.seed(my.specs$seed[x])
             # GC: 3/24/21: We were sampling the treatment effect... and the variance is not correct; it should be:
             #      *****************    sqrt((my.specs$s.t[x]/sqrt(my.specs$n.t[x]))^2 + (my.specs$s.c[x]/sqrt(my.specs$n.c[x])^2)       *****************
             #     **********************  In any case this is not the way to go because PBO is xbar is common across all sims


             # my.treatmeant.effect <- rnorm(n = my.specs$n.MC[x],
             #                               mean = my.specs$treatment.effect[x],
             #                               sd=my.specs$s.t[x]/sqrt(my.specs$n.t[x]))
             #
             #        *****************    The problem with this is xbar.c is fixed throughout - we want xbar.c to vary across simulated trials ***************

             #     for.return <- get.ts.ng.mc.df(mu.0.c = my.specs$mu.0.c[x], n.0.c = my.specs$n.0.c[x],
             #                                   alpha.0.c=my.specs$alpha.0.c[x], beta.0.c = my.specs$beta.0.c[x],
             #                                   xbar.c = my.specs$xbar.c[x],
             #                                   s.c = my.specs$s.c[x],
             #                                   n.c = my.specs$n.c[x], group.c="Control",
             #                                   mu.0.t = my.specs$mu.0.t[x], n.0.t = my.specs$n.0.t[x],
             #                                   alpha.0.t=my.specs$alpha.0.t[x], beta.0.t = my.specs$beta.0.t[x],
             #                                   xbar.t = my.specs$xbar.c[x] + my.treatmeant.effect ,
             #                                   s.t = my.specs$s.t[x], n.t = my.specs$n.t[x], group.t="Treatment",
             #                                   Delta.tv = my.specs$Delta.tv[x], Delta.lrv = my.specs$Delta.lrv[x],
             #                                   tau.tv = my.specs$tau.tv[x], tau.lrv = my.specs$tau.lrv[x],
             #                                   tau.ng = my.specs$tau.ng[x],
             #                                   n.MC = my.specs$n.MC[x]) %>%
             #       mutate(treatment.effect = my.specs$treatment.effect[x])

             # Instead we should be simulating PBO and control sample means
             # PBO population mean: Note xbar.c plays the role of the underlying treatment PBO mean
             my.PBO.means <- rnorm(n = my.specs$n.MC[x],
                                   mean = my.specs$xbar.c[x],
                                   sd=my.specs$s.c[x]/sqrt(my.specs$n.c[x]))

             # TRT population mean: Note xbar.c + my.specs$treatment.effect[x]; my.specs$xbar.c[x] is common throughout
             my.TRT.means <- rnorm(n = my.specs$n.MC[x],
                                   mean = my.specs$xbar.c[x] + my.specs$treatment.effect[x],
                                   sd=my.specs$s.t/sqrt(my.specs$n.t[x]))

             # GC 3/24/21: I added argument expand=F.
             for.return <- get.ts.ng.mc.df(mu.0.c = my.specs$mu.0.c[x], n.0.c = my.specs$n.0.c[x],
                                           alpha.0.c=my.specs$alpha.0.c[x], beta.0.c = my.specs$beta.0.c[x],
                                           xbar.c = my.PBO.means, s.c = my.specs$s.c[x],
                                           n.c = my.specs$n.c[x], group.c="Control",
                                           mu.0.t = my.specs$mu.0.t[x], n.0.t = my.specs$n.0.t[x],
                                           alpha.0.t=my.specs$alpha.0.t[x], beta.0.t = my.specs$beta.0.t[x],
                                           xbar.t = my.TRT.means,
                                           s.t = my.specs$s.t[x], n.t = my.specs$n.t[x], group.t="Treatment",
                                           Delta.tv = my.specs$Delta.tv[x], Delta.lrv = my.specs$Delta.lrv[x],
                                           tau.tv = my.specs$tau.tv[x], tau.lrv = my.specs$tau.lrv[x],
                                           tau.ng = my.specs$tau.ng[x],
                                           n.MC = my.specs$n.MC[x], expand=F) %>%
               mutate(treatment.effect = my.specs$treatment.effect[x])
             #message("\nget.NG.ts.df.mc finished\n")
             return(for.return)
           }

           # Add parallel parApply here
           if(goparallel==TRUE){
           parallel::clusterEvalQ(cl=cl, expr = {loadNamespace("tidyverse")})
           parallel::clusterExport(cl=cl, varlist=c( "get.ts.ng.mc.df", "my.func", "my.specs", "get.ts.ng.mc", "get.ng.post", "rt_ls"), envir = environment())
           for.plot <- bind_rows(parallel::parApply(cl=cl, X = matrix(1:nrow(my.specs)), MARGIN = 1, FUN = my.func)) %>%
             mutate(result = factor(result, c("Go", "Consider", "No-Go"))) %>%
             group_by(treatment.effect, result) %>%
             summarize(N = n()) %>%
             mutate(freq = N / sum(N)) %>% dplyr::select(result, treatment.effect, freq) %>%
             ungroup() %>%
             tidyr::complete(result, fill = list(N = 0, freq = 0))%>%
             mutate(mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c=alpha.0.c, beta.0.c = beta.0.c,
                    xbar.c = xbar.c, s.c = s.c,  group.c="Control",
                    mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t=alpha.0.t, beta.0.t = beta.0.t,
                    xbar.t = xbar.t, s.t = s.t,  group.t="Treatment",
                    Delta.LB=Delta.LB, Delta.UB=Delta.UB, ARatio=ARatio, N=N,
                    Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                    tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                    npoints=npoints, n.MC = n.MC,
                    seed = seed)
           } else {
                   for.plot <- bind_rows(apply( X = matrix(1:nrow(my.specs)), MARGIN = 1, FUN = my.func)) %>%
                           mutate(result = factor(result, c("Go", "Consider", "No-Go"))) %>%
                           group_by(treatment.effect, result) %>%
                           summarize(N = n()) %>%
                           mutate(freq = N / sum(N)) %>% dplyr::select(result, treatment.effect, freq) %>%
                           ungroup() %>%
                           complete(result, fill = list(N = 0, freq = 0))%>%
                           mutate(mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c=alpha.0.c, beta.0.c = beta.0.c,
                                  xbar.c = xbar.c, s.c = s.c,  group.c="Control",
                                  mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t=alpha.0.t, beta.0.t = beta.0.t,
                                  xbar.t = xbar.t, s.t = s.t,  group.t="Treatment",
                                  Delta.LB=Delta.LB, Delta.UB=Delta.UB, ARatio=ARatio, N=N,
                                  Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                  tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                                  npoints=npoints, n.MC = n.MC,
                                  seed = seed)
           }

           return(for.plot)
}


