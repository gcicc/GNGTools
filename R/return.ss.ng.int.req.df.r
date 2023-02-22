#' @title Get single sample normal-gamma interim treatment effect OC data.frame
#'
#' @param mu.0.t prior mean
#' @param n.0.t prior effective sample size
#' @param alpha.0.t prior alpha parameter
#' @param beta.0.t prior beta parameter
#' @param xbar.t sample mean at interim
#' @param s.t treatment standard deviation
#' @param interim.n.t interim sample sizes
#' @param final.n.t final sample size
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param xbar_ng Leave NULL to compute or provide user value
#' @param xbar_go Leave NULL to compute or provide user value
#' @param go.thresh interim go threshold
#' @param ng.thresh interim no-go threshold
#' @param include_nogo logical
#'
#' @return returns data.frame to assist with creating interim oc curve in single sample normal-gamma case
#' @export
#'
#' @examples
#' return.ss.ng.data.req.df()

return.ss.ng.data.req.df <-   function(mu.0.t = 0, n.0.t = .0001, alpha.0.t = 0.25, beta.0.t = 1,
                                 xbar.t=0.5, s.t = 2, interim.n.t =  10,
                                 final.n.t = 40,
                                 Delta.lrv = 1.25, Delta.tv = 1.75,
                                 tau.tv = .1, tau.lrv = .8, tau.ng = .65,
                                 xbar_ng = NULL, xbar_go=NULL,
                                 go.thresh=0.8, ng.thresh=0.8, include_nogo=TRUE){

        if(is.null(xbar_ng) | is.null(xbar_go)){
                results <- get.ss.ng.studyend.GNG(mu.0.t = mu.0.t, alpha.0.t=alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                                                  xbar.t = xbar.t, s.t = s.t, n.t = final.n.t,
                                                  Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                                  tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng)
                xbar_go <- ifelse(nrow(results$result.go) == 0, Inf, results$result.go$xbar)
                xbar_ng <- ifelse(nrow(results$result.ng) == 0, 0,  results$result.ng$xbar)
        }
        my.params <-   get.ng.post.df(mu.0 = mu.0.t, n.0 = n.0.t, alpha.0 = alpha.0.t, beta.0 = beta.0.t,
                                      xbar = xbar.t, s = s.t, n = interim.n.t, group="Treatment")  %>%
                rename(interim.n.t = n) %>%
                mutate(xbar_ng=xbar_ng, xbar_go=xbar_go,
                       complement.n.t = final.n.t - interim.n.t) %>%
                mutate(
                        t.go = xbar_go*(interim.n.t+complement.n.t)/complement.n.t - interim.n.t/complement.n.t*xbar.t,
                        t.ng = xbar_ng*(interim.n.t+complement.n.t)/complement.n.t - interim.n.t/complement.n.t*xbar.t,
                        Go = 1 - pt_ls(x = t.go, mu = mu.n, sigma = beta.n*(n.0.t+interim.n.t+complement.n.t)/((n.0.t+interim.n.t)*alpha.n), df=2*alpha.n),
                        NoGo =   pt_ls(x = t.ng, mu = mu.n, sigma = beta.n*(n.0.t+interim.n.t+complement.n.t)/((n.0.t+interim.n.t)*alpha.n), df=2*alpha.n),
                        decision=case_when(Go > go.thresh ~ "Go",
                                           NoGo > ng.thresh ~ "No-Go",
                                           TRUE ~ "Consider")
                )
        if (include_nogo==FALSE) my.params <- my.params %>% mutate(decision=case_when(Go > go.thresh ~ "Go",
                                                                                  TRUE ~ "Consider"))

        return(my.params)
}
