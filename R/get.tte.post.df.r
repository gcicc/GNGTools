#' @title Get posterior parameters for two-sample time to event case using normal approximation to log(HR) - data.frame
#' @description Returns the parameters of the posterior distribution in the case of normal approximation to hazard ratio - this version does not expand grid
#' @param m.con.prior prior number of events on control
#' @param m.trt.prior prior number of events on treatment
#' @param HR.prior Prior Hazard ratio
#' @param ARatio Randomization ratio in current trial
#' @param HR.obs Observed HR in current trial
#' @param m.obs Total events in current trail
#' @return a data.frame is returned
#' @examples
#' my.tte.post.df <- get.tte.post.df()
#' head(my.tte.post.df)
#' @export
get.tte.post.df <- function(m.con.prior=10, m.trt.prior=10, HR.prior=.7,
                            ARatio=1, HR.obs=.8, m.obs=50){
        data.frame(m.con.prior=m.con.prior, m.trt.prior=m.trt.prior, HR.prior=HR.prior,
                   ARatio=ARatio, HR.obs=HR.obs, m.obs=m.obs) %>%
                mutate(
                        prior.mean = log(HR.prior),
                        prior.sd = sqrt(4/(m.con.prior + m.trt.prior)),
                        ARatio = 1/(1+ARatio),
                        post.mean = (1 / m.con.prior + 1 / m.trt.prior) ^ (-1)/((1 / m.con.prior + 1 / m.trt.prior) ^ (-1) +
                                                                                        (1/(m.obs * ARatio*(1 - ARatio))) ^ (-1))*log(HR.prior) +    (1/(m.obs * ARatio*(1 - ARatio))) ^ (-1)/((1 / m.con.prior + 1 / m.trt.prior) ^ (-1) +  (1/(m.obs * ARatio*(1 - ARatio))) ^ (-1))*log(HR.obs),
                        post.sd = sqrt(((1 / m.con.prior + 1 / m.trt.prior)*
                                                (1/(m.obs * ARatio*(1 - ARatio))))/((1 / m.con.prior + 1 / m.trt.prior)+
                                                                                            (1/(m.obs * ARatio*(1 - ARatio))))))
}
