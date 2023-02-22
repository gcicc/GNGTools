#' @title Get two sample binary interim OC
#'
#' @param DecisionTable results from get.ts.bin.int.dec
#' @param runs number of simulation runs
#' @param ControlRate assumed rate in the control group
#' @param TreatmentEffect vector of treatment effects relative to control group.
#' @param a.con alpha parameter for control
#' @param b.con beta parameter for control
#' @param a.trt alpha parameter for treatment
#' @param b.trt beta parameter for treatment
#' @param Delta.lrv min TPP
#' @param Delta.tv Base TPP
#' @param tau.tv study-end threshold for Base TPP
#' @param tau.lrv study-end threshold for Min TPP
#' @param tau.ng study-end threshold for NG
#' @param go.thresh threshold for predictive probabilities at interim
#' @param ng.thresh threshold for predictive probabilities at interim
#' @param n.trt final sample size for treatment
#' @param n.con final sample size for control
#' @param n.int.c interim sample sizes for control
#' @param n.int.t interim sample sizes for treatment
#' @param include_nogo logical
#' @return dataframe of simulation results with the following columns:
#' * **Effect** the difference between the control Rate and treatment rate. (TreatmentRate = ControlRate + Effect)
#' * **run** the grouping variable by run.
#' * **assessment** the assessment in order, the last assessment is the final analysis results.
#' * **TreatmentEffect** The true response rate for the treatment group. (ControlRate is the same for all runs, it is not in the output, only the function input. We may want to change this for completeness, or include it in Param)
#' * **TotalResponse** for control (_C) and treatment(_T) respectively the total number of responses observed at that particular assessment. This is the cumulative sum over the assessments.
#' * **Decision** The decision at the assessment given the total response up to that point.  Note that the decision is calculated differently depending the assessment, the last assessment uses final assessment calculation while the earlier assessments use the interim assessment calculations.
#' @export
#'
#' @examples
#' my.ts.bin.int.oc <-  get.ts.bin.int.oc()
#' head(my.ts.bin.int.oc)

get.ts.bin.int.oc <- function (
                a.con = 1, b.con = 1,
                a.trt = 1, b.trt = 1,
                Delta.tv = .3, Delta.lrv = .2,
                tau.tv = .1, tau.lrv = .8, tau.ng = .65,
                go.thresh = .8, ng.thresh = .8,
                n.con = 40, n.trt = 40,
                n.int.c = c(10, 20, 30),
                n.int.t = c(10, 20, 30),
                DecisionTable= NULL,
                runs=500,
                ControlRate=.2,
                TreatmentEffect=seq(0, 0.8,.1),
                include_nogo=TRUE){

        if(is.null(DecisionTable)) DecisionTable <- get.ts.bin.int.dec(
                a.con = a.con, b.con = b.con,
                a.trt = a.trt, b.trt = b.trt,
                Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                go.thresh = go.thresh, ng.thresh = ng.thresh,
                n.con = n.con, n.trt = n.trt,
                n.int.c = n.int.c,  n.int.t = n.int.t,
                runs=runs, include_nogo=include_nogo)

        #Set up the sample sizes for each of the assessments.
        N.c <- c(n.int.c, n.con)
        ndifc <- N.c- dplyr::lag(N.c,1,default = 0)
        N.t <- c(n.int.t,n.trt)
        ndift <- N.t- dplyr::lag(N.t,1,default = 0)

        #Run the Simulations
        OCsims <- crossing(
                Effect = TreatmentEffect,
                run = 1:runs,
                assessment = 1:(length(n.int.c)+1)) %>%
                dplyr::mutate(
                        ControlRate=ControlRate,
                        TreatmentRate = ControlRate + Effect,
                        x.con = rbinom(n(),ndifc[assessment],ControlRate),
                        x.trt = rbinom(n(),ndift[assessment],TreatmentRate)
                ) %>%
                dplyr::group_by(run,Effect) %>%
                dplyr::mutate(
                        cumsum.x.con = cumsum(x.con),
                        cumsum.x.trt = cumsum(x.trt)
                ) %>% ungroup() # %>% select(-x.con,-x.trt)

        # Using the same time saving cheat as above for the final analysis. Also saving time by repeateding use of the get.int.decision interim decision simulations.
        FinalDecsisions <- OCsims %>%
                dplyr::filter(assessment == length(N.c)) %>%
                dplyr::distinct(assessment,cumsum.x.con,cumsum.x.trt) %>%
                dplyr::mutate(
                        Prob.trtV =1- p2beta_diff_Vector(Delta.tv,
                                                         a.con + cumsum.x.con,
                                                         b.con + n.con-cumsum.x.con,
                                                         a.trt + cumsum.x.trt,
                                                         b.trt + n.trt - cumsum.x.trt),
                        Prob_LRV =1- p2beta_diff_Vector(Delta.lrv,
                                                        a.con + cumsum.x.con,
                                                        b.con + n.con-cumsum.x.con,
                                                        a.trt + cumsum.x.trt,
                                                        b.trt + n.trt - cumsum.x.trt),
                        Decision = case_when(
                                Prob_LRV >= tau.lrv & Prob.trtV >= tau.tv ~ 'Go',
                                Prob_LRV < tau.ng & Prob.trtV < tau.tv ~ 'NoGo',
                                TRUE ~'Consider'
                        ) %>% factor(levels = c('Go','Consider','NoGo'))
                ) %>% dplyr::select(-Prob.trtV,-Prob_LRV) %>%
                rbind(
                        DecisionTable %>%
                                select(assessment = Interim,cumsum.x.con = IntermR_C,cumsum.x.trt= IntermR_T,Decision)
                )

        left_join(OCsims, FinalDecsisions, by = c('assessment','cumsum.x.con','cumsum.x.trt')) %>%
                dplyr::mutate(a.con = a.con, b.con = b.con,
                       a.trt = a.trt, b.trt = b.trt,
                       Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                       tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                       go.thresh = go.thresh, ng.thresh = ng.thresh,
                       n.con = n.con, n.trt = n.trt,
                       runs=runs, n.int.c = list(n.int.c), n.int.t = list(n.int.t),
                       ControlRate=ControlRate) %>% rename(decision=Decision)
}

