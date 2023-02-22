#' @title interim.calculator
#'
#' @param a.trt alpha parameter
#' @param b.trt beta parameter
#' @param Delta.lrv TPP info
#' @param Delta.tv TPP info
#' @param tau.tv thresholds
#' @param tau.lrv info
#' @param tau.ng info
#' @param interim.n.t interim sample size
#' @param final.n.t final sample size
#' @param p.success probability of success
#' @param responses number of responders
#' @param x_ng number needed for no-go
#' @param x_go number needed for go
#' @param go.thresh go threshold
#' @param ng.thresh no-go threshold
#'
#' @return a data.frame is returned
#' @export
#'
#' @examples
#' interim.calculator(a.trt = 1, b.trt = 1,Delta.lrv = .3, Delta.tv = .45,tau.tv = 0.10,
#' tau.lrv = .80, tau.ng = .65,interim.n.t = 15,final.n.t = 100, p.success = .4,
#' responses=10,x_ng = NULL, x_go=NULL,go.thresh=0.8, ng.thresh=0.8)
interim.calculator <- function(a.trt = 1, b.trt = 1,
                             Delta.lrv = .3, Delta.tv = .45,
                             tau.tv = 0.10, tau.lrv = .80, tau.ng = .65,
                             interim.n.t = 15,
                             final.n.t = 100, p.success = .4, responses=10,
                             x_ng = NULL, x_go=NULL,
                             go.thresh=0.8, ng.thresh=0.8){

  # a.trt = input$a_trt_calc;
  # b.trt = input$b_trt_calc;
  # Delta.lrv = input$Delta_lrv_calc/100;
  # Delta.tv = input$Delta_tv_calc/100;
  # tau.tv = input$tau_tv_calc/100;
  # tau.lrv = input$tau_lrv_calc/100;
  # tau.ng = input$tau_ng_calc/100;
  # interim.n.t = input$int_n_t_calc;
  # final.n.t = input$final_n_t_calc; p.success = input$p_success_calc/100; responses=input$response_calc;
  # x_ng = NULL; x_go=NULL;
  # go.thresh=input$go_thresh_calc/100; ng.thresh=input$ng_thresh_calc/100

  # If the responses required at study-end to declare Go/NOGo are not provided, then compute based on rule
  if(is.null(x_ng) | is.null(x_go)){
    results <- get.ss.bin.studyend.GNG(a.trt = a.trt, b.trt = b.trt,  n.trt = final.n.t,
                                       Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                       tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)
    x_go <- ifelse(nrow(results$result.go) == 0, Inf, results$result.go$x)
    x_ng <- ifelse(nrow(results$result.ng) == 0, 0,  results$result.ng$x)
  }

  # Create a grid
  my.df <- expand.grid(a0 = a.trt, b0 = b.trt, interim.n.t = interim.n.t, observed=0:final.n.t) %>%
    # Trim so it makes sense at interims
    dplyr::filter(observed <= interim.n.t) %>%
    arrange(interim.n.t) %>%
    # Here's the number of trials remaining after interim
    mutate(n.comp = final.n.t - interim.n.t) %>%
    # Update to beta parameters
    mutate(an = a0 + observed, bn = b0 + (interim.n.t - observed)) %>%
    # These are the probabilities of observing the interim data given
    # P.interim1 assigns Probability of interim results according to prior
    # P.interim2 assigns probability of interim results according to binomial distribution with users' p.success
    mutate(P.interim1 = dbbinom(x = observed, size=interim.n.t, alpha=a0, beta=b0),
           P.interim2 = dbinom(x=observed, size=interim.n.t, prob = p.success)) %>%
    mutate(need.for.go = pmax(0, x_go-observed),
           need.for.ng = pmax(0, x_ng - observed))

  # Probability of Go is predictive probability that complementary number of successes is between need.at.least and n.comp.
  # If we already have the required number P.Go = 1
  my.df$P.Go <- apply(matrix(1:nrow(my.df)), MARGIN = 1, FUN = function(x)
    # We go if we are larger than needed.for.go
    ifelse(my.df$need.for.go[x] == 0, 1, sum(dbbinom(x=my.df$need.for.go[x]:my.df$n.comp[x], size = my.df$n.comp[x], alpha = my.df$an[x], beta = my.df$bn[x]))))
  # Probability of NoGo is predictive probability that complementary number of successes between 0 and need.at.least.
  # If number needed for ng ==0, then no.go is not possible, set this to 0.
  my.df$P.NG <- apply(matrix(1:nrow(my.df)), MARGIN = 1, FUN = function(x)
    ifelse(my.df$need.for.ng[x] == 0, 0, sum(dbbinom(x=0:my.df$need.for.ng[x], size = my.df$n.comp[x], alpha = my.df$an[x], beta = my.df$bn[x]))))

  my.df$P.Go2 <- apply(matrix(1:nrow(my.df)), MARGIN = 1, FUN = function(x)
    # We go if we are larger than needed.for.go
    ifelse(my.df$need.for.go[x] == 0, 1, sum(dbinom(x=my.df$need.for.go[x]:my.df$n.comp[x], size = my.df$n.comp[x], prob = p.success))))
  # Probability of NoGo is predictive probability that complementary number of successes between 0 and need.at.least.
  # If number needed for ng ==0, then no.go is not possible, set this to 0.
  my.df$P.NG2 <- apply(matrix(1:nrow(my.df)), MARGIN = 1, FUN = function(x)
    ifelse(my.df$need.for.ng[x] == 0, 0, sum(dbinom(x=0:my.df$need.for.ng[x], size = my.df$n.comp[x],  prob = p.success))))

  my.df <- my.df %>% mutate(
    go.thresh=go.thresh,
    ng.thresh=ng.thresh,
    decision = case_when(P.Go > go.thresh ~ "Go",
                         P.NG > ng.thresh ~ "No-Go",
                         TRUE ~ "Consider"),
    study.end.Go = x_go, study.end.NG = x_ng)


  my.df <- my.df %>% dplyr::select(
    a0, b0,  # Prior parameters
    interim.n.t, observed, # Interim data in hand
    an, bn,  # Update to parameters
    P.interim1, P.interim2, # Probability of seeing interim data under prior and under user's prob of success
    study.end.Go, study.end.NG,   #  What's required of data at study end for Go/NG per rule parameters
    n.comp, # Number of trials left after interim
    need.for.go, need.for.ng, # min/max number needed in order to trigger Go/NG among n.comp
    go.thresh,ng.thresh, # thresholds
    P.Go, P.NG, # Predictive probability of Go/NG
    P.Go2, P.NG2, # Predictive probability of Go/NG
    decision # decision based on comparison of P.GO/P.NG and go.thresh/ng.thresh
  ) %>% dplyr::filter(observed ==responses)

  for.return <- list(
    'Study end requirements' = my.df %>% dplyr::select(study.end.Go, study.end.NG),
    'priors parameters' = my.df %>% dplyr::select(a0, b0),
    'interim data' = my.df %>% dplyr::select(interim.n.t, observed) %>% rename('interim size' = interim.n.t, responders = observed),
    'posterior parameters' = my.df %>% dplyr::select(an, bn),
    'Probability of interim data' = round(my.df %>% dplyr::select(P.interim1, P.interim2) %>% rename("Under user's prior" = P.interim1, "Under user's probability of success"=P.interim2),4),
    'Probability of Go' = round(my.df %>% dplyr::select(P.Go, P.Go2) %>% rename("Under user's prior" = P.Go, "Under user's probability of success" = P.Go2),4),
    'Probability of No-go' = round(my.df %>% dplyr::select(P.NG, P.NG2) %>% rename("Under user's prior" = P.NG, "Under user's probability of success" = P.NG2),4)
  )

  return(for.return)
  }
