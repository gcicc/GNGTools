#' @title Return two-sample normal-gamma predictive probability
#'
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size parameter for control group
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size parameter for treatment group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param xbar.t treatment mean
#' @param s.t treatment sd
#' @param n.t treatment sample size
#' @param xbar.c control mean
#' @param s.c control sd
#' @param n.c control sample size
#' @param xbar_ng Leave NULL to determine what is required or supply a value
#' @param xbar_go Leave NULL to determine what is required or supply a value
#' @param go.thresh If the predictive probability that study will conclude as 'Go' is larger than this threshold: Declare 'Interim go'.
#' @param ng.thresh If the predictive probability that study will conclude as 'No-Go' is larger than this threshold: Declare 'Interim no-go'
#' @param n.MC Monte Carlo simulation size
#' @return A dataframe is returned
#' @export
#'
#' @examples
#' holdit <- return.ts.ng.int.req.df()
#' head(holdit)

# This version replaces expand grid with data.frame
# xbar.c should be viewed as study-end xbar
return.ts.ng.int.req.df <-   function(mu.0.t = 0, n.0.t = .0001, alpha.0.t = 0.25, beta.0.t = 1,
                                 mu.0.c = 0, n.0.c = .0001, alpha.0.c = 0.25, beta.0.c = 1,
                                 xbar.t=c(1.9, 2, 2.1, 2.05), s.t = c(2, 2.1, 1.9, 2.04), n.t =  c(10,20,30, 40),
                                 xbar.c=c(1, 1.1, 1.5, 1.25), s.c = c(1.9, 2, 2.5, 2.25), n.c =  c(10,20,30, 40),
                                 Delta.lrv = 1.25, Delta.tv = 1.75,
                                 tau.tv = .1, tau.lrv = .8, tau.ng = .65,
                                 xbar_ng = NULL, xbar_go=NULL,
                                 go.thresh=0.8, ng.thresh=0.8,
                                 n.MC = 1000){

  # Now we pass the study.end results here to learn what is needed at final GIVEN study-end CONTROL mean
  if(is.null(xbar_ng) | is.null(xbar_go)){
    # Single values should be passed
    results <- get.ts.ng.studyend.GNG(mu.0.c = mu.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c, n.0.c = n.0.c,
                                      mu.0.t = mu.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                                      xbar.c = xbar.c[length(xbar.c)], s.c = s.c[length(s.c)], n.c = n.c[length(n.c)],
                                      xbar.t = xbar.t, s.t = s.t[length(s.t)], n.t = n.t[length(n.t)],
                                      Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                      tau.ng = tau.ng, tau.lrv = tau.lrv, tau.tv = tau.tv,
                                      n.MC = 1000)
    # Given control, what is required from treatment...
    xbar_go <- ifelse(nrow(results$result.go) == 0, Inf, results$result.go$xbar.t)
    xbar_ng <- ifelse(nrow(results$result.ng) == 0, 0,  results$result.ng$xbar.t)
  }

  # These process the interim results
  my.params.t <- get.ng.post.df(mu.0 = mu.0.t, n.0 = n.0.t, alpha.0 = alpha.0.t, beta.0 = beta.0.t,
                              xbar = xbar.t, s = s.t, n = n.t) %>%
    dplyr::rename(mu.0.t = mu.0, n.0.t = n.0, alpha.0.t = alpha.0, beta.0.t = beta.0, xbar.t = xbar, s.t = s, n.t = n,
                  tdf.0.t = tdf.0, sigma.0.t= sigma.0, mu.n.t= mu.n, n.n.t= n.n, alpha.n.t= alpha.n, beta.n.t=   beta.n,
                  t.df.n.t = tdf.n, sigma.n.t = sigma.n) %>%
    mutate(complement.n.t = n.t[length(n.t)] - n.t)
  my.params.c <- get.ng.post.df(mu.0 = mu.0.c, n.0 = n.0.c, alpha.0 = alpha.0.c, beta.0 = beta.0.c,
                              xbar = xbar.c, s = s.c, n = n.c)%>%
    dplyr::rename(mu.0.c = mu.0, n.0.c = n.0, alpha.0.c = alpha.0, beta.0.c = beta.0, xbar.c = xbar, s.c = s, n.c = n,
                  tdf.0.c = tdf.0, sigma.0.c= sigma.0, mu.n.c= mu.n, n.n.c= n.n, alpha.n.c= alpha.n, beta.n.c=   beta.n,
                  t.df.n.c = tdf.n, sigma.n.c = sigma.n )%>%
    mutate(complement.n.c = n.c[length(n.c)] - n.c)
  for.return <- data.frame(cbind(my.params.c, my.params.t)) %>%
    mutate(xbar_ng=xbar_ng, xbar_go=xbar_go) %>%

    mutate(
      # What is needed in complement to mean study end criteria
      t.go = xbar_go*(max(n.t))/complement.n.t - n.t/complement.n.t*xbar.t,
      t.ng = xbar_ng*(max(n.t))/complement.n.t - n.t/complement.n.t*xbar.t,
      # This is the probability of doing better/worse than t.go, t.ng
      Go = 1 - pt_ls(x = t.go, mu = mu.n.t, sigma = beta.n.t*(n.0.t+max(n.t)+complement.n.t)/((n.0.t+max(n.t))*alpha.n.t), df=2*alpha.n.t),
      NoGo =   pt_ls(x = t.ng, mu = mu.n.t, sigma = beta.n.t*(n.0.t+max(n.t)+complement.n.t)/((n.0.t+max(n.t))*alpha.n.t), df=2*alpha.n.t),
      # Go = 1 - pnorm(q = t.go, mean = mu.n.t, sd = s.t/sqrt(complement.n.t)),
      # NoGo =   pnorm(q = t.ng, mean = mu.n.t, sd = s.t/sqrt(complement.n.t)),
      decision=case_when(Go > go.thresh & n.c < max(n.c) ~ "Go",
                         NoGo > ng.thresh & n.c < max(n.c) ~ "No-Go",
                         NoGo < ng.thresh & Go < go.thresh & n.c < max(n.c) ~ "Consider",
                         xbar.t > xbar_go & n.c == max(n.c) ~ "Go",
                         xbar.t < xbar_ng & n.c == max(n.c) ~ "No-Go",
                         TRUE ~ "Consider"))# %>%
        # dplyr::select( mu.0.c, n.0.c, alpha.0.c, beta.0.c,mu.0.t, n.0.t, alpha.0.t, beta.0.t,
        #                xbar.c,  s.c, n.c, xbar.t, s.t, n.t, xbar_ng,  xbar_go, Go, NoGo, decision)
  return(for.return)
}
