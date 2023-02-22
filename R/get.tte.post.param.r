#' Update parameters in time to event case
#' @title Get TTE posterior parameters
#' @param m.con.prior prior number of events for control group
#' @param m.trt.prior prior number of events for treatment group
#' @param HR.prior hazard ratio estimate
#' @param ARatio randomization ratio
#' @param HR.obs HR observed
#' @param m.obs number of events
#'
#' @return A data.frame is returned.
#' @export
#'
#' @examples
#' my.tte.post.param <- get.tte.post.param()
#' my.tte.post.param

get.tte.post.param <- function(m.con.prior=10, m.trt.prior=10, HR.prior=.7,
         ARatio=1, HR.obs=.8, m.obs=50){
  ARatio <- 1/(1+ARatio)
  data.frame(
    post.mean = (1 / m.con.prior + 1 / m.trt.prior) ^ (-1)/
      ((1 / m.con.prior + 1 / m.trt.prior) ^ (-1) +
         (1/(m.obs * ARatio*(1 - ARatio))) ^ (-1))*log(HR.prior) +
      (1/(m.obs * ARatio*(1 - ARatio))) ^ (-1)/((1 / m.con.prior + 1 / m.trt.prior) ^ (-1) +
                                                  (1/(m.obs * ARatio*(1 - ARatio))) ^ (-1))*log(HR.obs),
    post.sd = sqrt(((1 / m.con.prior + 1 / m.trt.prior)*
                      (1/(m.obs * ARatio*(1 - ARatio))))/((1 / m.con.prior + 1 / m.trt.prior)+
                                                            (1/(m.obs * ARatio*(1 - ARatio))))))
}


