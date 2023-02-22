#' @title Make time to event predictive probability OC curve
#'
#' @param m.con.prior prior number of events for control group
#' @param m.trt.prior prior number of events for treamtnet group
#' @param HR.prior hazard ratio estimate
#' @param ARatio Randomization ratio
#' @param interim.m interim number of events
#' @param final.m final number of events
#' @param HR.tv Base TPP for HR
#' @param HR.lrv Min TPP for HR
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param HR.lower lower bound for HR
#' @param HR.upper upper bound for HR
#' @param go.thresh predictive probability go threshold
#' @param ng.thresh predictive probability no-go threshold
#' @param npoints number of points for OC Curve
#' @param n.MC Monte Carlo sample size
#' @param include_nogo logical
#'
#' @return A data.frame is returned.
#' @export
#'
#' @examples \donttest{
#' my.tte.int.oc.df <- get.tte.int.oc.df()
#' my.tte.int.oc.df
#'}
get.tte.int.oc.df <- function(m.con.prior=50,
                              m.trt.prior=50,
                              HR.prior=.845,
                              ARatio=1,
                              HR.lower=.0025,
                              HR.upper=2,
                              npoints=20,
                              interim.m=c(428, 750, 1000),
                              final.m=1500,
                              HR.tv=0.75, HR.lrv = .9,
                              tau.tv=.1, tau.lrv=.8, tau.ng=.65,
                              n.MC=2000,
                              go.thresh= 0.8, ng.thresh =0.8, include_nogo=TRUE){

    # Run this once to grab the HR.ng and HR.go required for NOGO/GO at end of trial

  # 1. Create grid of HRs
  my.hrs <- seq(HR.lower, HR.upper, length.out=npoints)

  # 2. function to get observed HRs at interim
  get.interim.results <- function(n.MC = n.MC, HR = .1, IAs = c(interim.m, final.m)){
    # Univariate
    # log(HR)-hat ~ N(log(HR), 4/nEvents)
    # Group-sequential
    # Z_k ~ N(theta*sqrt(d/4), 1)
    # Cov(Z_m1, Z_m2) = sqrt( (d_m1/4)/(d_m2/4))

    # Distribution theory
    # Z_k ~ N(theta * sqrt(I_k))
    # Cov(Z_k1, Z_k2) = sqrt(I_k1/I_k2); k1 < k2
    # Z1, ..., Zk ~ MVN
    # Ik ~ d_k/4
    # check <- (rmvnorm(n = 10000, log(HR)*sqrt(bigK/4),my.matrix))

    # Determine number of IAs and final
    K <- length(IAs)

    # Build covariance matrix
    my.matrix <- matrix(nrow=K, ncol=K, 0)
    for(i in 1:K){
      for(j in 1:K){
        my.matrix[i,j] <- ifelse(j >= i, sqrt((IAs[i]/4)/(IAs[j]/4)), sqrt((IAs[j]/4)/(IAs[i]/4)))
      }
    }

    # Get standardized values
    check <- rmvnorm(n = n.MC, mean = log(HR)*sqrt(IAs/4), sigma = my.matrix)
    # colMeans(check)
    # exp(colMeans(check))
    # Transform
    check <- apply(matrix(1:length(IAs)), MARGIN=1, FUN = function(x) {
      check[,x]<- check[,x]/sqrt(IAs[x]/4)
    })

    check <- data.frame(check) %>% mutate(trial=1:nrow(.))
    return(check)
  }

  # 3. Build trial data
  my.trials <- bind_rows(apply(X = matrix(1:length(my.hrs)), MARGIN=1, FUN = function(x){
    get.interim.results(n.MC = n.MC, HR=my.hrs[x], IAs = c(interim.m, final.m)) %>% mutate(UHR=my.hrs[x])
    })) %>% pivot_longer(-c(trial, UHR)) %>% dplyr::rename(interim=name) %>% mutate(analysis=rep(c(interim.m, final.m), n.MC*length(my.hrs))) %>%
    dplyr::select(-interim) %>%
    dplyr::rename(HR=value) %>% mutate(HR=exp(HR))

  # 4. Subset simulation to interim
  my.trials.int <- my.trials %>% dplyr::filter(analysis != final.m)

  # 5. run once to get studyend GNG cutoffs
  run.once <- get.tte.studyend.GNG(m.con.prior = m.con.prior,
                                   m.trt.prior = m.trt.prior,
                                   HR.prior=HR.prior,
                                   ARatio=ARatio,
                                   HR.obs=1,
                                   m.obs = final.m,
                                   HR.tv= HR.tv,
                                   HR.lrv = HR.lrv,
                                   tau.tv=tau.tv,
                                   tau.lrv=tau.lrv,
                                   tau.ng=tau.ng)
  HR.ng <- run.once$result.ng$HR.obs
  HR.go <- run.once$result.go$HR.obs



  temp <- return.tte.int.data.req(m.con.prior=m.con.prior, m.trt.prior=m.trt.prior, HR.prior=HR.prior,
                           ARatio=ARatio,
                           interim.HR=my.trials.int$HR,
                           interim.m=my.trials.int$analysis,
                           final.m=final.m,
                           HR.tv=HR.tv, HR.lrv = HR.lrv,
                           tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng,
                           HR.ng = HR.ng, HR.go=HR.go,
                           go.thresh= go.thresh, ng.thresh =ng.thresh, include_nogo=include_nogo) %>%
    mutate(trial=my.trials.int$trial, UHR=my.trials.int$UHR)%>%
    mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
    group_by(trial) %>%
    mutate(dnum = as.numeric(decision))


  # 6. Process interim decisions
  process.sims <- return.tte.int.data.req(m.con.prior=m.con.prior, m.trt.prior=m.trt.prior, HR.prior=HR.prior,
                   ARatio=ARatio,
                   interim.HR=my.trials.int$HR,
                   interim.m=my.trials.int$analysis,
                   final.m=final.m,
                   HR.tv=HR.tv, HR.lrv = HR.lrv,
                   tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng,
                   HR.ng = HR.ng, HR.go=HR.go,
                   go.thresh= go.thresh, ng.thresh =ng.thresh, include_nogo=include_nogo) %>%
    mutate(trial=my.trials.int$trial, UHR=my.trials.int$UHR)%>%
    mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
    group_by(trial) %>%
    mutate(dnum = as.numeric(decision)) %>% arrange(trial, interim.m)

  # Get study.end

  study.end <- my.trials %>% dplyr::filter(analysis==final.m)

  # 7. get  results by interim
  interim.results.summary <- process.sims %>%
    group_by(UHR) %>%
    mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
    group_by(interim.m, UHR) %>%
    dplyr::count(decision, .drop=F) %>%
    mutate(relFreq = n/sum(n)) %>%
    mutate(analysis = as.character(interim.m)) %>%
    ungroup() %>% arrange(UHR, interim.m, decision) %>%
    dplyr::select(-interim.m)

  # Here we group trials by first occurance of Go, No-go, and streak of Consider
  any.interim.results <- bind_rows(
    process.sims %>% group_by(UHR, trial) %>% mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
      mutate(dnum = as.numeric(decision)) %>%
      dplyr::filter(dnum != 3) %>% dplyr::select(UHR, trial, decision) %>% slice(1),
    process.sims %>% group_by(UHR, trial) %>%
      mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
      mutate(dnum = as.numeric(decision)) %>%
      dplyr::filter(all(dnum == 3)) %>%
      dplyr::select(UHR, trial, decision) %>% slice(1)) %>% mutate(analysis="any interim") %>% rename(interims.d = decision) %>% dplyr::select(-analysis)

  study.end.results <- study.end %>%
    mutate(decision=ifelse(HR < HR.go, "Go", ifelse(HR > HR.ng, "No-Go", "Consider"))) %>%
    mutate(analysis="Study End") %>%
    mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
    group_by(UHR) %>% dplyr::select(-HR, -analysis)%>% rename( final.d = decision)

  any.analysis.summary <- left_join(any.interim.results, study.end.results) %>% mutate(decision = case_when(
    # Inteirm Go trumps final No-Go
    interims.d == "Go" & final.d == "Go" ~ "Go",
    interims.d == "Go" & final.d == "No-Go" ~ "Go",
    interims.d == "Go" & final.d =="Consider" ~ "Go",
    # interim No-go trumps final Go
    interims.d == "No-Go" & final.d == "Go" ~ "No-Go",
    interims.d == "No-Go" & final.d == "No-Go" ~ "No-Go",
    interims.d == "No-Go" & final.d == "Consider" ~ "No-Go",
    # Final Go/No-Go trumps interim consider
    interims.d == "Consider" & final.d == "Go" ~ "Go",
    interims.d == "Consider" & final.d =="No-Go" ~ "No-Go",
    interims.d == "Consider" & final.d =="Consider" ~ "Consider"
  )) %>% group_by(UHR) %>% dplyr::count(decision, .drop=F) %>%
    mutate(relFreq = n/sum(n),
           analysis = "any analysis") %>% ungroup()

  any.interim.summary <- bind_rows(
    process.sims %>% group_by(UHR, trial) %>% mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
      mutate(dnum = as.numeric(decision)) %>%
      dplyr::filter(dnum != 3) %>% dplyr::select(UHR, trial, decision) %>% slice(1),
    process.sims %>% group_by(UHR, trial) %>%
      mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
      mutate(dnum = as.numeric(decision)) %>%
      dplyr::filter(all(dnum == 3)) %>%
      dplyr::select(UHR, trial, decision) %>% slice(1)) %>%
    group_by(UHR) %>% dplyr::count(decision, .drop=F) %>%
    mutate(relFreq = n/sum(n),
           analysis = "any interim") %>% ungroup()

  # GC: study end results need to be merged (by trial and UHR) with interim results. Then for each trial/UHR we need to ask: 'Any Go' and 'Any No-Go'
  # temp <- study.end %>%
  #   mutate(decision=ifelse(HR < HR.go, "Go", ifelse(HR > HR.ng, "No-Go", "Consider"))) %>%
  #   mutate(analysis="Study End") %>%
  #   mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
  #   group_by(UHR) %>% arrange(trial, UHR)

  study.end.summary <- study.end %>%
    mutate(decision=ifelse(HR < HR.go, "Go", ifelse(HR > HR.ng, "No-Go", "Consider"))) %>%
    mutate(analysis="Study End") %>%
    mutate(decision = factor(decision, c("Go", "No-Go", "Consider"))) %>%
    group_by(UHR) %>%
    dplyr::count(decision, .drop=F) %>%
    mutate(relFreq = n/sum(n)) %>%
    mutate(analysis = "Study End")

  for.return <- bind_rows(study.end.summary, any.analysis.summary, any.interim.summary, interim.results.summary) %>% group_by(UHR, analysis) %>%
    mutate(m.con.prior = m.con.prior,
           m.trt.prior = m.trt.prior,
           HR.prior=HR.prior,
           ARatio=ARatio,
           HR.obs=1,
           m.obs = final.m,
           HR.tv= HR.tv,
           HR.lrv = HR.lrv,
           tau.tv=tau.tv,
           tau.lrv=tau.lrv,
           tau.ng=tau.ng)

  return(for.return)
}





