#' Get single sample normal-gamma interim operating characteristic data.frame
#'
#' @param mu.0.t prior mean
#' @param n.0.t prior effective sample size
#' @param alpha.0.t prior alpha parameter
#' @param beta.0.t prior beta parameter
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param Delta.OC.LB OC Lower bound
#' @param Delta.OC.UB OC Uppder bound
#' @param npoints number of points to span
#' @param n.MC number of trials at each point
#' @param s.t standard deviation
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param interim.n.t interim sample sizes
#' @param final.n.t final sample size
#' @param go.thresh interim go threshold
#' @param ng.thresh interim no-go threshold
#' @param include_nogo logical
#'
#' @return A data.frame is returned.
#' @export
#'
#' @examples
#' holdit <- get.ss.ng.trt.int.oc.df(n.MC=100, npoints=3)
#' head(holdit)
get.ss.ng.trt.int.oc.df <- function(mu.0.t = 0,
                                    n.0.t = .0001,
                                    alpha.0.t = 0.25,
                                    beta.0.t = 1,
                                    Delta.lrv = 1.25,
                                    Delta.tv = 1.75,
                                    Delta.OC.LB= -5,
                                    Delta.OC.UB=5,
                                    npoints = 10,
                                    n.MC=1000,
                                    s.t = 2,
                                    tau.tv = 0.1,
                                    tau.lrv = 0.8,
                                    tau.ng = 0.65,
                                    interim.n.t = c(10, 20, 30),
                                    final.n.t = 40,
                                    go.thresh=0.8,
                                    ng.thresh=0.8,
                                    include_nogo=TRUE){

  # create a function to sample correlated means at interims.
  # It appears that brute force summary stats calculations beats group-sequential theory
  # We need to generate interim and final results for each trial - this builds summary stats


  # 1. Create grid of treatment effects
  my.deltas <- seq(Delta.OC.LB, Delta.OC.UB, length.out=npoints)

  # 2. function to get observed treatment effects at interim
  get.interim.results <- function(IAs = c(interim.n.t,final.n.t), delta=my.deltas[1], sd=s.t){
    t(apply(matrix(1:length(IAs)), MARGIN=1, function(x){
      temp <- rnorm(n = IAs[length(IAs)], mean=delta, sd=s.t)
      #data.frame(n.t=IAs[x], mean=mean(temp[1:IAs[x]]), sd=sd(temp[1:IAs[x]]))
      matrix(data = c(n.t=IAs[x], mean=mean(temp[1:IAs[x]]), sd=sd(temp[1:IAs[x]])), nrow=1, ncol=3)
    }))
  }

  # 3. Build trial data
  n.analyses = c(interim.n.t,final.n.t)
  my.trials <- bind_rows(apply(X = matrix(1:length(my.deltas)), MARGIN=1, FUN = function(y){
  data.frame(do.call(rbind,
          (lapply(X=1:n.MC, FUN = function(x) {
            get.interim.results(IAs = n.analyses, delta=my.deltas[y], sd=s.t)})))) %>% mutate(delta=my.deltas[y])})) %>%
    mutate(trial=rep(1:(nrow(.)/length(n.analyses)), each=length(n.analyses)))
  colnames(my.trials)[1:3] <- c("n.t", "mean", "sd")

  # 4. Subset simulation to interim
  my.trials.ints <- my.trials %>% dplyr::filter(n.t != final.n.t)

  # 5. Run this once to grab the no-go and go thresholds at study end to save recalculation
  run.once <- return.ss.ng.int.data.req(
    mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t,
    xbar.t=my.trials.ints$mean[1], s.t = my.trials.ints$sd[1], interim.n.t = my.trials.ints$n.t[1],
    final.n.t = final.n.t,
    Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
    tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = .65,
    xbar_ng = NULL, xbar_go=NULL,
    go.thresh=0.8, ng.thresh=0.8)
  # this is what is required at study-end for Go/No-go
  xbar.ng <- run.once$xbar_ng
  xbar.go <- run.once$xbar_go

  # 6. Process interim decisions
 process.sim <- return.ss.ng.data.req.df(mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t,
                    xbar.t=my.trials.ints$mean, s.t = my.trials.ints$sd, interim.n.t = my.trials.ints$n.t,
                    final.n.t = final.n.t,
                    Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                    tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = .65,
                    xbar_ng = xbar.ng[1], xbar_go=xbar.go[1],
                    go.thresh=0.8, ng.thresh=0.8, include_nogo=include_nogo) %>%
   dplyr::select(interim.n.t, Go, NoGo, decision) %>%
   mutate(trial = my.trials.ints$trial,
          delta=my.trials.ints$delta)

 # Here's study-end results per simulation
 study.end <- my.trials %>% dplyr::filter(n.t == final.n.t) %>%
         dplyr::mutate(decision=ifelse(mean < run.once$xbar_ng[1], "No-Go",
                                       ifelse(mean > run.once$xbar_go, "Go", "Consider"))) %>%
         group_by(delta) %>%
         dplyr::mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
         dplyr::count(decision, .drop=F) %>%
         dplyr::mutate(analysis="Study-end")

  # these report first choices
  anyGo.NoGo <- bind_rows(
  process.sim %>%
          group_by(trial) %>%
          dplyr::mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
          dplyr::mutate(dnum = as.numeric(decision)) %>%
          dplyr::filter(dnum != 3) %>% dplyr::select(trial, delta, decision) %>% slice(1),
  process.sim %>%
          group_by(trial) %>%
          dplyr::mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
          dplyr::mutate(dnum = as.numeric(decision)) %>%
          group_by(trial, delta) %>%
          dplyr::summarize(Consider = all(dnum==3)) %>%
          dplyr::filter(Consider == TRUE) %>%
          dplyr::rename(decision = Consider) %>%
          dplyr::mutate(decision="Consider")) %>%
          dplyr::rename(anyDecision=decision)

  # Broken down by trial with adjustments for previous interims.
  any.interim.summary <- anyGo.NoGo %>%
          group_by(delta) %>%
          dplyr::mutate(anyDecision = factor(anyDecision, c("Go", "No-Go", "Consider"))) %>%
          dplyr::count(anyDecision=anyDecision, .drop=FALSE) %>%
          dplyr::rename(decision=anyDecision) %>%
          dplyr::mutate(analysis="Any Interim")

  # Broken down by interim.n.t, with no adjustments
  interim.results.summary <- process.sim %>%
          group_by(delta, interim.n.t) %>%
          dplyr::mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
          dplyr::count(decision=decision, .drop=FALSE) %>%
          dplyr::mutate(analysis=as.character(interim.n.t)) %>%
          ungroup() %>% dplyr::select(-interim.n.t)

  study.end.results <- my.trials %>%
          dplyr::filter(n.t == final.n.t) %>%
          dplyr::mutate(decision=ifelse(mean < run.once$xbar_ng[1], "No-Go",
                                        ifelse(mean > run.once$xbar_go, "Go", "Consider"))) %>%
          group_by(delta) %>%
          dplyr::mutate(decision = factor(decision, c("Go", "No-Go", "Consider")))

 # any analysis
  any.analysis <- left_join(anyGo.NoGo, study.end.results) %>%
          dplyr::mutate(decision = case_when(
                  # Inteirm Go trumps final No-Go
                  anyDecision == "Go" & decision == "Go" ~ "Go",
                  anyDecision == "Go" & decision == "No-Go" ~ "Go",
                  anyDecision == "Go" & decision =="Consider" ~ "Go",
                  # interim No-go trumps final Go
                  anyDecision == "No-Go" & decision == "Go" ~ "No-Go",
                  anyDecision == "No-Go" & decision == "No-Go" ~ "No-Go",
                  anyDecision == "No-Go" & decision == "Consider" ~ "No-Go",
                  # Final Go/No-Go trumps interim consider
                  anyDecision == "Consider" & decision == "Go" ~ "Go",
                  anyDecision == "Consider" & decision =="No-Go" ~ "No-Go",
                  anyDecision == "Consider" & decision =="Consider" ~ "Consider"
                  )) %>% group_by(delta) %>% dplyr::count(decision, .drop=F) %>%
          dplyr::mutate(relFreq = n/sum(n),
                 analysis = "Any Analysis") %>% ungroup()

  for.return <- bind_rows(study.end, any.interim.summary, interim.results.summary, any.analysis) %>%
          group_by(delta, analysis) %>%
          dplyr::mutate(relFreq=n/sum(n))%>%
          dplyr::mutate(Delta.lrv = Delta.lrv,
                        Delta.tv = Delta.tv,
                        tau.tv = tau.tv,
                        tau.lrv = tau.lrv,
                        tau.ng = tau.ng)

  return(for.return)
  }

