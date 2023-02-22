#' @title Get time to event sample size OC data.frame
#'
#' @param m.con.prior number of prior events for control
#' @param m.trt.prior number of prior events for treamtent
#' @param HR.prior HR estimate
#' @param ARatio randomization ratio
#' @param m.obs observed number of events
#' @param m.lower lower bound on number of events
#' @param m.upper upper bound on number of events
#' @param HR.tv Base TPP for HR
#' @param HR.lrv Min TPP for HR
#' @param HR.user user's HR
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return A data frame is returned
#' @export
#'
#' @examples \donttest{
#' my.tte.ssize.oc.df <- get.tte.ssize.oc.df()
#' my.tte.ssize.oc.df
#' }
get.tte.ssize.oc.df <- function(m.con.prior=10, m.trt.prior=10, HR.prior=.75,
                                ARatio=1, m.obs=50,
                                m.lower=40, m.upper=120,
                                HR.lrv=0.75, HR.tv=.75, HR.user = 0.845, tau.tv=.1,
                                tau.lrv=.8, tau.ng=.65){
  # Create a grid
  specs <- expand.grid(m.obs=floor(seq(m.lower, m.upper, length.out=20)),
                       HRs.from=c(HR.lrv, HR.tv, HR.user, 1)) %>%
    mutate(m.con.prior=m.con.prior, m.trt.prior=m.trt.prior, HR.prior=HR.prior,
           ARatio=ARatio, HR.lrv=HR.lrv, HR.tv=HR.tv, HR.user = HR.user, tau.tv=tau.tv,
           tau.lrv=tau.lrv, tau.ng=tau.ng)

  # Apply this function to every row of the grid
  for.plot <- bind_rows(apply(X = matrix(1:nrow(specs)), MARGIN = 1,
                              FUN = function(x){
                                # For each row of specs, get all possible probs along a fine grid
                                results <- get.tte.df(m.con.prior = specs$m.con.prior[x], m.trt.prior = specs$m.trt.prior[x],
                                                      HR.prior = specs$HR.prior[x],
                                                      ARatio = specs$ARatio[x], m.obs=specs$m.obs[x],
                                                      HR.obs=seq(0.1, 2,.005), HR.tv=specs$HR.tv[x],
                                                      HR.lrv = specs$HR.lrv[x], tau.tv=specs$tau.tv[x],
                                                      tau.lrv=specs$tau.lrv[x], tau.ng=specs$tau.ng[x])
                                # Identify max criteria for Go and min criteria for no-go
                                result.go <- results %>% dplyr::filter(result== "Go") %>% slice(n())
                                # Determine max number of TRT responders for No-Go
                                result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(1)

                                my.df <- data.frame(Go = pnorm(q = log(result.go$HR.obs),
                                                               mean = log(c(HR.lrv, HR.tv, HR.user, 1)),
                                                               sd=sqrt(4/result.go$m.obs)),
                                                    NoGo = 1 - pnorm(q = log(result.ng$HR.obs),
                                                                     mean = log(c(HR.lrv, HR.tv, HR.user, 1)),
                                                                     sd=sqrt(4/result.go$m.obs)),
                                                    HR=c(HR.lrv, HR.tv, HR.user, 1))  %>%
                                  mutate(m.con.prior = specs$m.con.prior[x],
                                         m.trt.prior = specs$m.trt.prior[x], HR.prior = specs$HR.prior[x],
                                         ARatio = specs$ARatio[x], m.obs=specs$m.obs[x],
                                         HR.tv=specs$HR.tv[x], HR.lrv = specs$HR.lrv[x], HR.user=specs$HR.user[x],
                                         tau.tv=specs$tau.tv[x], tau.lrv=specs$tau.lrv[x], tau.ng=specs$tau.ng[x])
                                my.df
                              }))
  for.plot$key.label = c("HR = Max TPP = ", "HR = Base TPP = ", "HR = User's TPP = ", "HR = Null = ")
  for.plot$key = factor(paste0(for.plot$key.label, for.plot$HR),
                        paste0(for.plot$key.label, for.plot$HR)[1:4][order(c(HR.lrv, HR.tv, HR.user, 1))])

  for.plot
}
