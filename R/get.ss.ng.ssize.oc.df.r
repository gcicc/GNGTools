#' @title Get single sample normal-gamma sample size Oc curve data.frame
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size for treatment group
#' @param s.t sample sd for treatment group
#' @param n.t sample size for treatment group
#' @param SS.OC.N.LB lower bound for OC curve
#' @param SS.OC.N.UB upper bound for OC curve
#' @param npoints number of points
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param Delta.user User's value for underlying treatment response
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' my.ss.ng.ssize.oc.df <- get.ss.ng.ssize.oc.df()
#' head(my.ss.ng.ssize.oc.df)
#' }
get.ss.ng.ssize.oc.df <- function(mu.0.t = 3, n.0.t = 10, alpha.0.t = 0.25,
                                  beta.0.t = 1,
                                  s.t = 5, n.t = 50, SS.OC.N.LB = floor(50*0.75),
                                  SS.OC.N.UB=floor(50*2), npoints=15,
                                  Delta.lrv=2.5, Delta.tv=4, Delta.user = 3,
                                  tau.tv=.1, tau.lrv=.8, tau.ng=.65){
  # Create a grid
  specs <- expand.grid(n.t=floor(seq(SS.OC.N.LB, SS.OC.N.UB, length.out=npoints)),
                       Deltas.from=c(0, Delta.lrv, Delta.tv, Delta.user)) %>%
    mutate(mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t,
           s.t = s.t, n.t = n.t,
           Delta.lrv=Delta.lrv, Delta.tv=Delta.tv, Delta.user = Delta.user,
           tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng)


  for.plot <- bind_rows(apply(X = matrix(1:nrow(specs)), MARGIN = 1,
                              FUN = function(x){
                                # For each row of specs, get all possible probs along a fine grid
                                results <- get.ss.ng.df(mu.0.t = specs$mu.0.t[x], n.0.t = specs$n.0.t[x],
                                                        alpha.0.t = specs$alpha.0.t[x], beta.0.t = specs$beta.0.t[x],
                                                        s.t=specs$s.t[x], n.t=specs$n.t[x],
                                                        xbar.t=seq(min(min(specs$Delta.lrv, specs$Delta.user, 0)
                                                                       -specs$s.t/sqrt(specs$n.t)*5),
                                                                   max(max(specs$Delta.tv, specs$Delta.user, 0)
                                                                       +specs$s.t/sqrt(specs$n.t)*5), length.out=75),
                                                        Delta.lrv=specs$Delta.lrv[x], Delta.tv = specs$Delta.tv[x],
                                                        tau.tv=specs$tau.tv[x], tau.lrv=specs$tau.lrv[x],
                                                        tau.ng=specs$tau.ng[x])
                                # Identify max criteria for Go and min criteria for no-go
                                result.go <- results %>% dplyr::filter(result== "Go") %>% slice(1)
                                # Determine max number of TRT responders for No-Go
                                result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(n())

                                refine.go <-  get.ss.ng.df(mu.0.t = specs$mu.0.t[x], n.0.t = specs$n.0.t[x],
                                                           alpha.0.t = specs$alpha.0.t[x], beta.0.t = specs$beta.0.t[x],
                                                           s.t=specs$s.t[x], n.t=specs$n.t[x],
                                                           xbar.t=seq(result.go$xbar - 0.25*result.go$s/sqrt(result.go$n),
                                                                      result.go$xbar + 0.25*result.go$s/sqrt(result.go$n),
                                                                      length.out=50),
                                                           Delta.lrv=specs$Delta.lrv[x], Delta.tv = specs$Delta.tv[x],
                                                           tau.tv=specs$tau.tv[x], tau.lrv=specs$tau.lrv[x],
                                                           tau.ng=specs$tau.ng[x])

                                refine.ng <-  get.ss.ng.df(mu.0.t = specs$mu.0.t[x], n.0.t = specs$n.0.t[x],
                                                           alpha.0.t = specs$alpha.0.t[x], beta.0.t = specs$beta.0.t[x],
                                                           s.t=specs$s.t[x], n.t=specs$n.t[x],
                                                           xbar.t=seq(result.ng$xbar - 0.25*result.ng$s/sqrt(result.ng$n),
                                                                      result.ng$xbar + 0.25*result.ng$s/sqrt(result.ng$n),
                                                                      length.out=50),
                                                           Delta.lrv=specs$Delta.lrv[x], Delta.tv = specs$Delta.tv[x],
                                                           tau.tv=specs$tau.tv[x], tau.lrv=specs$tau.lrv[x],
                                                           tau.ng=specs$tau.ng[x])
                                # Identify max criteria for Go and min criteria for no-go
                                result.go <- refine.go %>% dplyr::filter(result== "Go") %>% slice(1)
                                # Determine max number of TRT responders for No-Go
                                result.ng <- refine.ng %>% dplyr::filter(result== "No-Go") %>% slice(n())

                                my.df <- data.frame(Go = 1 - pnorm(q = result.go$xbar,
                                                                   mean = c(0, Delta.lrv, Delta.tv, Delta.user),
                                                                   sd=(result.go$s/sqrt(result.ng$n))),
                                                    NoGo = pnorm(q = result.ng$xbar,
                                                                 mean = c(0, Delta.lrv, Delta.tv, Delta.user),
                                                                 sd=result.ng$s/sqrt(result.ng$n)),
                                                    Delta=c(0, Delta.lrv, Delta.tv, Delta.user))  %>%
                                  mutate(n.t = specs$n.t[x], mu.0.t = specs$mu.0.t[x], n.0.t = specs$n.0.t[x],
                                         beta.0.t = specs$beta.0.t[x], s.t = specs$s.t[x],
                                         Delta.lrv=specs$Delta.lrv[x], Delta.tv=specs$Delta.tv[x],
                                         Delta.user = specs$Delta.user[x], tau.tv=specs$tau.tv[x],
                                         tau.lrv=specs$tau.lrv[x], tau.ng=specs$tau.ng[x])
                                my.df
                              }))

  for.plot$key <- factor(as.character(for.plot$Delta), as.character(for.plot$Delta)[1:4])
  levels(for.plot$key) <-  c(
    TeX("$\\Delta\\,$ = NULL = 0"),
    TeX(paste("$\\Delta\\,$ = Min TPP = $", Delta.lrv, "$ ")),
    TeX(paste("$\\Delta\\,$ = Base TPP = $", Delta.tv, "$ ")),
    TeX(paste("$\\Delta\\,$ = User defined = $", Delta.user, "$")))
  for.plot$key <- factor(for.plot$key, levels(for.plot$key)[order(c(0, Delta.lrv, Delta.tv,
                                                                    Delta.user))])

  for.plot
}
