#' @title Get single sample normal-gamma data.frame with decision output
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
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' my.ss.ng.df <- get.ss.ng.df()
#' head(my.ss.ng.df)
#' }
get.ss.ng.df <- function(mu.0.t = 0, n.0.t = 10, alpha.0.t = .25, beta.0.t = 1,
                         xbar.t = seq(-1,5,.25), s.t = 1, n.t = 50,
                         Delta.tv = 1.75, Delta.lrv = 1.5, tau.tv = 0.10,
                         tau.lrv = 0.80, tau.ng = 0.65)
{
  # Create a simulation grid
  my.grid <- expand.grid(
    mu.0 = mu.0.t, n.0 = n.0.t, alpha.0 = alpha.0.t, beta.0 = beta.0.t,
    xbar = xbar.t,  s = s.t, n = n.t,
    Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
    tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

  my.results <- bind_rows(apply(X = matrix(1:nrow(my.grid)), MARGIN = 1,
                                FUN = function(x){
                                  get.ng.post(mu.0 = my.grid$mu.0[x], n.0 = my.grid$n.0[x],
                                              alpha.0 = my.grid$alpha.0[x],
                                              beta.0 = my.grid$beta.0[x],
                                              xbar = my.grid$xbar[x], s = my.grid$s[x], n = my.grid$n[x],
                                              group = "Treatment")})) %>%
    mutate(Delta.tv =Delta.tv, Delta.lrv=Delta.lrv, tau.tv=tau.tv,
           tau.lrv=tau.lrv, tau.ng=tau.ng) %>%
    mutate(
      P.R1 = 1 - pt_ls(x = Delta.lrv, df = tdf.n, mu = mu.n, sigma = sigma.n),
      P.R3 = 1 - pt_ls(x = Delta.tv, df = tdf.n, mu = mu.n, sigma = sigma.n),
      result = ifelse(P.R1 > tau.lrv & P.R3 > tau.tv, "Go",
                      ifelse(P.R1 <= tau.ng  & P.R3 <= tau.tv, "No-Go", "Consider"))) %>%
    mutate(result = factor(result, c("Go", "Consider", "No-Go")))

  return(my.results)
}
