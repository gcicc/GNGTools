#' @title Get time to event data.frame
#'
#' @param m.con.prior number of prior events for control
#' @param m.trt.prior number of prior events for treamtent
#' @param HR.prior HR estimate
#' @param HR.obs Observed HR
#' @param m.obs Obeserved number of events
#' @param ARatio Randomization ratio
#' @param HR.tv Base TPP for HR
#' @param HR.lrv Min TPP for HR
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.tte.df <- get.tte.df()
#' head(my.tte.df)

get.tte.df <- function(m.con.prior = 50, m.trt.prior = 50, HR.prior = .845,
                       HR.obs = seq(0.3, 1, .01), m.obs = seq(10, 200, 5),
                       ARatio = 0.5,
                       HR.tv = .80, HR.lrv = 0.9,
                       tau.tv = .10, tau.lrv = .20, tau.ng = 0.35) {

  P = ARatio / (ARatio + 1)

  my.grid <- expand.grid(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior,
                         HR.prior = HR.prior, ARatio = ARatio, HR.obs = HR.obs, m.obs = m.obs)
  my.results <- bind_rows(do.call(Map, c(f = get.tte.post.param, my.grid)))
  my.grid <- data.frame(cbind(my.grid, my.results)) %>%
    dplyr::mutate(m.prior = m.con.prior + m.trt.prior,
           ARatio=ARatio,
           HR.tv = HR.tv,
           HR.lrv = HR.lrv,
           tau.tv = tau.tv,
           tau.lrv = tau.lrv,
           tau.ng = tau.ng,
           P.R1 = pnorm(log(HR.tv), post.mean, post.sd),
           P.R3 = pnorm(log(HR.lrv),post.mean, post.sd),
           result = factor(
             ifelse(P.R1 > tau.tv & P.R3 > tau.lrv, "Go",
                    ifelse(P.R1 < tau.tv & P.R3 < tau.ng, "No-Go", "Consider")),
             c("Go", "Consider", "No-Go")))

  return(my.grid)
}


