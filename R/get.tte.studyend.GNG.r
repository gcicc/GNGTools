#' @title Get TTE studyend GNG decision
#' @param m.con.prior number of prior events for control
#' @param m.trt.prior number of prior events for treatment
#' @param HR.prior HR estimate
#' @param ARatio randomization ratio
#' @param HR.obs Observed HR
#' @param m.obs observed events
#' @param HR.tv Base TPP for HR
#' @param HR.lrv Min TPP for HR
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @return a list is returned holding what is needed from data to achieve Go/No-Go
#' @export
#'
#' @examples
#' my.tte.studyend.GNG <- get.tte.studyend.GNG(m.con.prior = 50,m.trt.prior = 50, HR.prior=1.2,
#' ARatio=1, HR.obs=1.3, m.obs = 200, HR.tv= 1.4, HR.lrv = 1.25, tau.tv=.1, tau.lrv=.8, tau.ng=.65)
#' my.tte.studyend.GNG

get.tte.studyend.GNG <- function(m.con.prior = 50,m.trt.prior = 50, HR.prior=1.2, ARatio=1, HR.obs=1.3, m.obs = 200,
                                 HR.tv= 1.4, HR.lrv = 1.25, tau.tv=.1, tau.lrv=.8, tau.ng=.65){

if(HR.tv <= HR.lrv){
  results <- get.tte.df(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior,
                        HR.prior = HR.prior, ARatio = ARatio, m.obs=m.obs,
                        HR.obs=seq(.0025, 5,.0025), HR.tv=HR.tv, HR.lrv = HR.lrv,
                        tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng)
  result.go <- results %>% dplyr::filter(result== "Go") %>% dplyr::slice(n())
  # Determine max number of TRT responders for No-Go
  result.ng <- results %>% dplyr::filter(result== "No-Go") %>% dplyr::slice(1)
} else {
  results <- get.tte.df(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior,
                        HR.prior = HR.prior, ARatio = ARatio, m.obs=m.obs,
                        HR.obs=seq(.0025, 5,.0025), HR.tv=HR.tv, HR.lrv = HR.lrv,
                        tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng)
  result.go <- results %>% dplyr::filter(result== "Go") %>% dplyr::slice(1)
  # Determine max number of TRT responders for No-Go
  result.ng <- results %>% dplyr::filter(result== "No-Go") %>% dplyr::slice(dplyr::n())
}

  return(list(result.go = result.go, result.ng = result.ng))
}


# get.tte.studyend.GNG()
