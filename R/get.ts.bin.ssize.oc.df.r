#' @title Get two-sample binary sample size OC data.frame
#'
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior beta parameter for treatment group
#' @param dcurve.con underlying rate for control group
#' @param Aratio randomization ratio
#' @param SS.OC.N.LB lower bound for sample size
#' @param SS.OC.N.UB upper bound for sample size
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param SS.OC.Delta user's treatment effect
#' @param npoints number of points
#'
#' @return returns a data.frame ready to create two-sample binary sample size operating characteristics curve
#' @export
#'
#' @examples \donttest{
#' my.ts.bin.ssize.oc.df <- get.ts.bin.ssize.oc.df()
#' }
#' @author Greg Cicconetti
get.ts.bin.ssize.oc.df <- function(a.con = 1, b.con = 1,
                                   a.trt = 1, b.trt = 1,
                                   dcurve.con = .12,
                                   Aratio=2,
                                   SS.OC.N.LB = 40, SS.OC.N.UB = 160,
                                   Delta.lrv = .2, Delta.tv=.25, SS.OC.Delta = .25,
                                   tau.tv = 0.10, tau.lrv = .8, tau.ng = .65,
                                   npoints=3){

        # Run total sample size over all integers between bounds
        specs <- data.frame(TE.OC.N = seq(SS.OC.N.LB, SS.OC.N.UB, 1)) %>%
                mutate(n.con = round(TE.OC.N/(1+Aratio)),
                       n.trt = TE.OC.N - n.con,
                       dcurve.con = dcurve.con,
                       p.con1 = dcurve.con + Delta.lrv,
                       p.con2 = dcurve.con + Delta.tv,
                       p.con3 = dcurve.con + SS.OC.Delta,
                       a.con = a.con, b.con = b.con,
                       a.trt = a.trt, b.trt = b.trt,
                       Delta.lrv = Delta.lrv, Delta.tv=Delta.tv,
                       tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng, Aratio=Aratio) %>%
                dplyr::filter(n.trt == Aratio*n.con)
        # This function sets the observed placebo rate based on rounding the user's dcurve.con*n.con
        runit <- function(x){
                get.ts.bin.dec.df(a.con = specs$a.con[x], b.con = specs$b.con[x],
                                       a.trt = specs$a.trt[x], b.trt = specs$b.trt[x],
                                       n.con = specs$n.con[x], x.con = 0:specs$n.con[x],
                                       n.trt = specs$n.trt[x], x.trt = 0:specs$n.trt[x],
                                       Delta.tv = specs$Delta.tv[x], Delta.lrv = specs$Delta.lrv[x],
                                       tau.tv = specs$tau.tv[x], tau.lrv = specs$tau.lrv[x],
                                       tau.ng = specs$tau.ng[x]) %>%
                        mutate(result = factor(result, c("Go", "Consider", "No-Go")),
                               total.n=n.con+n.trt,
                               prob.null = dbinom(x=x.trt, size = n.trt, prob = specs$dcurve.con[x])*
                                       dbinom(x=x.con, size = n.con, prob = specs$dcurve.con[x]),
                               prob.lrv = dbinom(x=x.trt, size = n.trt, prob = specs$p.con1[x])*
                                       dbinom(x=x.con, size = n.con, prob = specs$dcurve.con[x]),
                               prob.tv = dbinom(x=x.trt, size = n.trt, prob = specs$p.con2[x])*
                                       dbinom(x=x.con, size = n.con, prob = specs$dcurve.con[x]),
                               prob.user = dbinom(x=x.trt, size = n.trt, prob = specs$p.con3[x])*
                                       dbinom(x=x.con, size = n.con, prob = specs$dcurve.con[x]),
                               Aratio = specs$Aratio[x]) %>%
                        group_by(result) %>%
                        summarize(prob.null = sum(prob.null),
                                  prob.lrv = sum(prob.lrv),
                                  prob.tv = sum(prob.tv),
                                  prob.user = sum(prob.user)) %>%
                        mutate(a.con = specs$a.con[x], b.con = specs$b.con[x],
                               a.trt = specs$a.trt[x], b.trt = specs$b.trt[x],
                               n.con = specs$n.con[x], n.trt = specs$n.trt[x],
                               n.total=n.trt + n.con, Delta.tv = specs$Delta.tv[x],
                               Delta.lrv = specs$Delta.lrv[x], tau.tv = specs$tau.tv[x],
                               tau.lrv = specs$tau.lrv[x], tau.ng = specs$tau.ng[x],
                               Aratio=specs$Aratio[x], dcurve.con=specs$dcurve.con[x]) }

        big.results <- bind_rows(apply(X = matrix(seq(1,nrow(specs), length.out = npoints)),
                                       MARGIN=1, FUN= function(x) runit(x) )) %>%
                gather(key=key, value=value, prob.null, prob.lrv, prob.tv, prob.user, factor_key = T)

        levels(big.results$key) <-c(TeX("$\\Delta\\,$ = NULL = 0%"),
                                    TeX(paste("$\\Delta\\,$ = Min TPP = ", Delta.lrv*100, "% ")),
                                    TeX(paste("$\\Delta\\,$ = Base TPP = ", Delta.tv*100, "% ")),
                                    TeX(paste("$\\Delta\\,$ = User defined = ", SS.OC.Delta*100,
                                              "%")))
        big.results$key <- factor(big.results$key,
                                  levels(big.results$key)[order(c(0, Delta.lrv, Delta.tv,
                                                                  SS.OC.Delta))])

        return(big.results)
}

