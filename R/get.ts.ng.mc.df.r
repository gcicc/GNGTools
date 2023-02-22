#' @title Get two-sample normal-gamma MC-based data.frame
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
#' @param n.MC number of MC sampling
#' @param seed random seed
#' @param expand logical; if true expand.grid is employed; else data.frame is employed. Former provides all combinations
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.ts.ng.mc.df <- get.ts.ng.mc.df(mu.0.c = 0, n.0.c = 10,
#' alpha.0.c=1, beta.0.c = 4,
#' xbar.c = seq(-3,3,length.out=20), s.c = 3, n.c = 25, group.c="Control",
#' mu.0.t = 0, n.0.t = 10, alpha.0.t=1, beta.0.t = 4,
#' xbar.t = seq(0, 6, length.out=20), s.t = 2, n.t = 25, group.t="Treatment",
#' Delta.tv = 1.75, Delta.lrv = 1.5,  tau.tv = .1, tau.lrv = .8, tau.ng = 0,
#' n.MC = 1000, seed=1234, expand=TRUE)
#' my.ts.ng.mc.df
#'
get.ts.ng.mc.df <- function(
  mu.0.c = 0, n.0.c = 10, alpha.0.c=.25 * 4, beta.0.c = 1 * 4,
  xbar.c = seq(-3,3,length.out=20), s.c = 3, n.c = 25, group.c="Control",
  mu.0.t = 0, n.0.t = 10, alpha.0.t=.25 * 4, beta.0.t = 1 * 4,
  xbar.t = seq(0, 6, length.out=20), s.t = 2, n.t = 25, group.t="Treatment",
  Delta.tv = 1.75, Delta.lrv = 1.5,  tau.tv = .1, tau.lrv = .8, tau.ng = 0,
  n.MC = 1000, seed=1234, expand=TRUE){

  # Create a simulation grid
  if(expand==TRUE) {
  my.grid <- expand.grid(
    mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c,
    xbar.c = xbar.c, s.c = s.c, n.c = n.c,
    mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t,
    xbar.t = xbar.t, s.t = s.t, n.t = n.t,
    Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
    tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
    n.MC = n.MC, seed = seed)
  } else {
    my.grid <- data.frame(
      mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c,
      xbar.c = xbar.c, s.c = s.c, n.c = n.c,
      mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t,
      xbar.t = xbar.t, s.t = s.t, n.t = n.t,
      Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
      tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
      n.MC = n.MC, seed = seed)
}
  # Create a function that makes a single evaluation
  my.function1 <- function(
    mu.0.c = my.grid$mu.0.c, n.0.c = my.grid$n.0.c, alpha.0.c = my.grid$alpha.0.c,
    beta.0.c = my.grid$beta.0.c,
    xbar.c = my.grid$xbar.c, s.c = my.grid$s.c, n.c = my.grid$n.c,
    mu.0.t = my.grid$mu.0.t, n.0.t = my.grid$n.0.t, alpha.0.t = my.grid$alpha.0.t,
    beta.0.t = my.grid$beta.0.t,
    xbar.t = my.grid$xbar.t, s.t = my.grid$s.t, n.t = my.grid$n.t,
    Delta.tv = my.grid$Delta.tv, Delta.lrv = my.grid$Delta.lrv,
    tau.tv = my.grid$tau.tv, tau.lrv = my.grid$tau.lrv, tau.ng = my.grid$tau.ng,
    seed = my.grid$seed, n.MC = my.grid$n.MC){

    return(get.ts.ng.mc(mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c = alpha.0.c,
                        beta.0.c = beta.0.c,
                        xbar.c = xbar.c, s.c = s.c, n.c = n.c, group.c="Control",
                        mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t,
                        beta.0.t = beta.0.t,
                        xbar.t = xbar.t, s.t = s.t, n.t = n.t, group.t="Treatment",
                        Delta.tv = Delta.tv, Delta.lrv = Delta.lrv, seed = seed,
                        # This next line was missing, so defaults were being used.
                        tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                        n.MC = n.MC))
  }

  # Apply the function to each row of the grid
  my.grid <- cbind(my.grid, list_to_dataframe(do.call(Map, c(f= my.function1,
                                                                    my.grid)))) %>%
    mutate(result = factor(result, c("Go", "Consider", "No-Go")))
  return(my.grid)
}


