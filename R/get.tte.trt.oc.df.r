#' Get time to event treatment oc data.frame
#' @title Get TTE treatment effect OC curve data.frame
#' @param m.con.prior number of prior events for control
#' @param m.trt.prior number of prior events for treamtent
#' @param HR.prior HR estimate
#' @param ARatio randomization ratio
#' @param m.obs observed events
#' @param HR.tv HR for Base TPP
#' @param HR.lrv HR for Min TPP
#' @param HR.lower Lower bound for OC curve
#' @param HR.upper Upper bound for OC curve
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.tte.trt.oc.df <- get.tte.trt.oc.df()
#' my.tte.trt.oc.df
#' @author Greg Cicconetti
get.tte.trt.oc.df <- function(m.con.prior = 10, m.trt.prior = 10, HR.prior = .8,
                              ARatio = .5, m.obs = 50,
                              HR.tv = .7, HR.lrv = .9,
                              HR.lower=0.3, HR.upper=2,
                              tau.tv = 0.1, tau.lrv = 0.8, tau.ng = 0.65){

  if(HR.tv < HR.lrv) {
    # Get results for fine grid of HR.obs values
    results <- get.tte.df(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior,
                          HR.prior = HR.prior,
                          HR.obs = seq(HR.lower, HR.upper, .005), m.obs = m.obs,
                          ARatio = ARatio,
                          HR.tv = HR.tv, HR.lrv = HR.lrv,
                          tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

    # Identify min/max needed for no-go, go
    result.go <- results %>% dplyr::filter(result== "Go") %>% slice(n())
    # Determine max number of TRT responders for No-Go
    result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(1)

    my.df <- results %>%
      mutate(Go=pnorm(q =log(result.go$HR.obs), mean=log(HR.obs), sd=sqrt(4/m.obs)),
             NoGo=1 - pnorm(q = log(result.ng$HR.obs), mean = log(HR.obs), sqrt(4/m.obs))) %>%
      mutate(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior, HR.prior = HR.prior,
             m.obs = m.obs, ARatio = ARatio,
             HR.tv = HR.tv, HR.lrv = HR.lrv,
             tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
             HR.upper=HR.upper, HR.lower=HR.lower)} else {
               results <- get.tte.df(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior,
                                     HR.prior = HR.prior,
                                     HR.obs = seq(HR.lower, HR.upper, .005), m.obs = m.obs,
                                     ARatio = ARatio,
                                     HR.tv = HR.tv, HR.lrv = HR.lrv,
                                     tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

               # Identify min/max needed for no-go, go
               result.go <- results %>% dplyr::filter(result== "Go") %>% slice(1)
               # Determine max number of TRT responders for No-Go
               result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(n())

               my.df <- results %>%
                 mutate(Go=1 - pnorm(q =log(result.go$HR.obs), mean=log(HR.obs), sd=sqrt(4/m.obs)),
                        NoGo= pnorm(q = log(result.ng$HR.obs), mean = log(HR.obs), sqrt(4/m.obs))) %>%
                 mutate(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior, HR.prior = HR.prior,
                        m.obs = m.obs, ARatio = ARatio,
                        HR.tv = HR.tv, HR.lrv = HR.lrv,
                        tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                        HR.upper=HR.upper, HR.lower=HR.lower)

             }


  my.df
}

# get.tte.trt.oc.df()
