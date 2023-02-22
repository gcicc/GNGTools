#' @title Get two-sample binary interim decision
#'
#' @param runs the number of simulation runs
#' @param include_nogo logical
#' @param a.con alpha parameter for control
#' @param b.con beta parameter for control
#' @param a.trt alpha parameter for treatment
#' @param b.trt beta parameter for treatment
#' @param n.trt final sample size for treatment
#' @param n.con final sample size for control
#' @param n.int.c interim sample sizes for control
#' @param n.int.t interim sample sizes for treatment
#' @param Delta.lrv min TPP
#' @param Delta.tv Base TPP
#' @param tau.tv study-end threshold for Base TPP
#' @param tau.lrv study-end threshold for Min TPP
#' @param tau.ng study-end threshold for NG
#' @param go.thresh threshold for predictive probabilities at interim
#' @param ng.thresh threshold for predictive probabilities at interim
#'
#' @return a Data frame with the following columns
#' * **IntermR_C** the number of control responses
#' * **IntermR_T** the number of treatment responses
#' * **Interm** the index of which interim analysis is being assessed.
#' * **Go/No-Go/Consider** the proportion of simulations generating a final go/No-Go/Consider decision.
#' * **Decision** The interim decision based off of the simulated results.
#' @export
#' @examples
#' my.ts.bin.int.dec <- get.ts.bin.int.dec()
#' my.ts.bin.int.dec
get.ts.bin.int.dec <- function(
                a.con = 1, b.con = 1,
                a.trt = 1, b.trt = 1,
                n.trt = 40, n.con = 40,
                n.int.c = c(10, 20, 30),
                n.int.t = c(10, 20, 30),
                Delta.lrv = 0.15, Delta.tv = .30,
                tau.tv = 0.10, tau.lrv = 0.80, tau.ng = 0.65,
                go.thresh=.8, ng.thresh=.8,
                runs=500, include_nogo=TRUE){

        # simulate number of successes across the whole grid.
        DecisionTableSims <- expand.grid(
                IntermR_C = 0:max(n.int.c),
                IntermR_T = 0:max(n.int.t),
                Interim = 1:length(n.int.c),
                Run = 1:runs) %>%
                dplyr::filter(IntermR_C <= n.int.c[Interim],
                              IntermR_T <= n.int.t[Interim]) %>% # remove impossible outcomes for the lower interim analysis
                dplyr::mutate(
                        x.trt = rbinom(n(),n.trt-n.int.t[Interim],
                                            rbeta(n(),a.trt+IntermR_T,b.trt +n.int.t[Interim]-IntermR_T))+IntermR_T,
                        x.con = rbinom(n(),n.con-n.int.c[Interim],
                                             rbeta(n(),a.con+IntermR_C,b.con +n.int.c[Interim]-IntermR_C))+IntermR_C)

        # to save some CPU time, we only calculate the probability once for each outcome. then merge back in.
        # at least for the example used here this goes from 700k calculations to 1.6k.
        ProbFrame <-DecisionTableSims %>%
                dplyr::distinct(x.trt,x.con) %>%
                dplyr::mutate(
                        P.R3 =1- p2beta_diff_Vector(Delta.tv,
                                                       a.con + x.con,
                                                       b.con + n.con-x.con,
                                                       a.trt + x.trt,
                                                       b.trt + n.trt - x.trt),
                        P.R1 =1- p2beta_diff_Vector(Delta.lrv,
                                                        a.con + x.con,
                                                        b.con + n.con-x.con,
                                                        a.trt + x.trt,
                                                        b.trt + n.trt - x.trt),
                        Decision = case_when(
                                P.R1 >= tau.lrv & P.R3 >= tau.tv ~ 'Go',
                                P.R1 < tau.ng & P.R3 < tau.tv ~ 'NoGo',
                                TRUE ~'Consider'
                        ) %>% factor(levels = c('Go','Consider','NoGo'))
                )

        # Join the simulations with the probabilities and then count the number of each outcome.
        DecisionTable <- left_join(DecisionTableSims,ProbFrame, by =c('x.trt','x.con') ) %>%
                dplyr::group_by(IntermR_C, IntermR_T, Interim, Decision,.drop=F) %>%
                dplyr::summarise(Rate = n()/runs) %>%
                pivot_wider(names_from = Decision,values_from = Rate) %>%
                dplyr::mutate(
                        Decision = factor(case_when(
                                Go >= go.thresh ~ 'Go',
                                NoGo >= ng.thresh ~ 'NoGo',
                                TRUE ~ 'Consider'), levels = c('Go','Consider','NoGo')),
                        include_nogo=TRUE
                )


        if(include_nogo==FALSE){
                DecisionTable <- DecisionTable %>%
                        dplyr::mutate(
                                Decision = factor(case_when(
                                        Go >= go.thresh ~ 'Go',
                                        TRUE ~ 'Consider'), levels = c('Go','Consider')),
                                include_nogo=FALSE)
        }

        DecisionTable$n.int.c <- n.int.c[DecisionTable$Interim]
        DecisionTable$n.int.t <- n.int.t[DecisionTable$Interim]

        DecisionTable <-  DecisionTable %>% dplyr::mutate(
                a.con = a.con, b.con = b.con,
                a.trt = a.trt, b.trt = b.trt,
                n.trt = n.trt, n.con = n.con,
                Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                go.thresh=go.thresh, ng.thresh=ng.thresh,
                runs=runs, include_nogo=include_nogo
        )


        return(dplyr::ungroup(DecisionTable))
}
