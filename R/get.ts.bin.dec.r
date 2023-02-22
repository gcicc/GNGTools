
#' @title Get two sample binary decision
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param n.con observed sample size for control group
#' @param x.con observed number of responders for control group
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior beta parameter for treatment group
#' @param n.trt observed sample size for treatment group
#' @param x.trt observed number of responders for treatment group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @author Greg Cicconetti
#' @return returns a data.frame holding Posterior probabilties of interest and Go/No-Go result
#' @examples
#' get.ts.bin.dec()

get.ts.bin.dec = function(a.con = 1, b.con = 1, n.con = 40, x.con = 5,
                               a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 20,
                               Delta.tv = 0.25, Delta.lrv = 0.2,
                               tau.tv = 0.10, tau.lrv = 0.8, tau.ng = 0.65)
{
        # Reporting P(Delta > Min.tpp), P(Delta > Base.tpp)
        probs <- 1 - p2beta(relation="DIFF", approach="DIRECT", x=c(Delta.lrv, Delta.tv),
                            a1=a.con+x.con, b1=b.con+n.con-x.con,
                            a2=a.trt+x.trt, b2=b.trt+n.trt-x.trt)
        # NOW COMPARE HERE
        my.df <- data.frame(P.R1 = probs[1], P.R3 = probs[2]) %>%
                mutate(result = ifelse(P.R1 >= tau.lrv & P.R3 >= tau.tv, "Go",
                                       ifelse(P.R1 < tau.ng & P.R3 < tau.tv, "No-Go", "Consider")))
        return(my.df)
}



