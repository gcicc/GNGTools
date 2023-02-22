#' @title Get two-sample binary treatment OC data.frame
#'
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param dcurve.con Response rate for control group
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior beta parameter for treatment group
#' @param TE.OC.N total sample size for OC
#' @param Aratio Randomization ratio
#' @param TE.OC.Delta.LB Lower bound for OC curve
#' @param TE.OC.Delta.UB Upper bound for OC curve
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @references Sverdlov O, Ryeznik Y, Wu S. Exact Bayesian Inference Comparing Binomial Proportions, With Application to Proof-of-Concept Clinical Trials. Therapeutic Innovation & Regulatory Science. 2015;49(1):163-174. doi:10.1177/2168479014547420
#' @return A data.frame is returned that allows for creation of treatment effect operating characteristics curve
#' @export
#'
#' @examples \donttest{
#' my.ts.bin.trt.oc.df <- get.ts.bin.trt.oc.df()
#' head(my.ts.bin.trt.oc.df)
#' }

get.ts.bin.trt.oc.df  <- function(a.con = 1, b.con = 1,  dcurve.con = 0.12,
                                  a.trt = 1, b.trt = 1,
                                  TE.OC.N = 80, Aratio=2,
                                  TE.OC.Delta.LB = 0, TE.OC.Delta.UB = 1 - 0.12,
                                  Delta.tv =  .35, Delta.lrv = 0.2,
                                  tau.tv = .01, tau.lrv = 0.8, tau.ng = 0.65){

        # Determine number of control subjects
        n.con = round(TE.OC.N / (1 + Aratio))
        n.trt = TE.OC.N - n.con

        results <- get.ts.bin.dec.df(a.con = a.con, b.con = b.con, n.con = n.con,
                                          a.trt = a.trt, b.trt = b.trt, n.trt = n.trt,
                                          x.trt = 0:n.trt, x.con = 0:n.con,
                                          Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                                          tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)


        my.df <-  bind_rows(apply(X=matrix(seq(TE.OC.Delta.LB, TE.OC.Delta.UB, 0.01)),
                                  MARGIN = 1, FUN = function(x) {
                                          results %>% mutate(
                                                  prob = dbinom(x = x.con, size = n.con, prob = dcurve.con) *
                                                          dbinom(x = x.trt, size = n.trt, prob = dcurve.con + x)) %>%
                                                  group_by(result) %>%
                                                  summarize(prob=sum(prob)) %>% mutate(dcurve.con=dcurve.con, x = x)
                                  })) %>%       mutate(a.con = a.con, b.con = b.con,
                                                       a.trt = a.trt, b.trt = b.trt,
                                                       n.con = n.con,  n.trt = n.trt, Aratio = Aratio,
                                                       Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                                                       tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                                                       TE.OC.Delta.LB=TE.OC.Delta.LB, TE.OC.Delta.UB=TE.OC.Delta.UB)

        my.df
}

