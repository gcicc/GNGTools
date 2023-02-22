#' Return single sample binary predictive probability
#'
#' @param a.trt prior alpha parameter
#' @param b.trt prior beta parameter
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param interim.n.t number of trials at interim
#' @param final.n.t number of trials at final
#' @param x.go responses for go; leave null for standard rule
#' @param x.ng responses for go; leave null for standard rule
#' @param p.success probability of success
#' @param go.thresh go threshold for predictive probability at interim
#' @param ng.thresh no-go threshold for predictive probability at interim
#'
#' @return Return single sample binary predictive probability
#' @export
#'
#' @examples
#' holdit <- return.ss.bin.int.req()
#' head(holdit)
return.ss.bin.int.req <- function(a.trt = 1, b.trt = 1,
                             Delta.lrv = .2, Delta.tv = .35,
                             tau.tv = 0.10, tau.lrv = .80, tau.ng = .65,
                             interim.n.t = c(0:39), final.n.t = 40, p.success = .5,
                             x.ng = NULL, x.go=NULL,
                             go.thresh=0.8, ng.thresh=0.8){

  if(is.null(x.ng) | is.null(x.go)){
    results <- get.ss.bin.studyend.GNG(a.trt = a.trt, b.trt = b.trt,  n.trt = final.n.t, x.trt = floor(final.n.t/2),
                                       Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                       tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)
    x.go <- ifelse(nrow(results$result.go) == 0, Inf, results$result.go$x)
    x.ng <- ifelse(nrow(results$result.ng) == 0, 0,  results$result.ng$x)
  }


  my.df <- expand.grid(a0 = a.trt, b0 = b.trt, interim.n.t = interim.n.t, observed=0:final.n.t) %>%
    dplyr::filter(observed <= interim.n.t) %>%
    arrange(interim.n.t) %>%
    mutate(n.comp = final.n.t - interim.n.t) %>%
    mutate(an = a0 + observed, bn = b0 + (interim.n.t - observed)) %>%
    # P.interim1 assigns Probability of interim results according to prior
    # P.interim2 assigns probability of interim results according to binomial distribution with users' p.success
    mutate(P.interim1 = dbbinom(x = observed, size=interim.n.t, alpha=a0, beta=b0),
           P.interim2 = dbinom(x=observed, size=interim.n.t, prob = p.success)) %>%
    mutate(need.for.go = pmax(0, x.go-observed),
           need.for.ng = pmax(0, x.ng - observed))
  # Probability of Go is predictive probability that complementary number of successes is between need.at.least and n.comp.
  # If we already have the required number PGo = 1
    my.df$PGo <- apply(matrix(1:nrow(my.df)), MARGIN = 1, FUN = function(x)
          ifelse(my.df$need.for.go[x] == 0, 1, sum(dbbinom(x=my.df$need.for.go[x]:my.df$n.comp[x], size = my.df$n.comp[x], alpha = my.df$an[x], beta = my.df$bn[x]))))
  # Probability of NoGo is predictive probability that complementary number of successes between 0 and need.at.least.
    # If number needed for ng ==0, then no.go is not possible, set this to 0.
    my.df$PNG <- apply(matrix(1:nrow(my.df)), MARGIN = 1, FUN = function(x)
      ifelse(my.df$need.for.ng[x] == 0, 0, sum(dbbinom(x=0:my.df$need.for.ng[x], size = my.df$n.comp[x], alpha = my.df$an[x], beta = my.df$bn[x]))))
    my.df <- my.df %>% mutate(decision = case_when(PGo > go.thresh ~ "Go",
                                                   PNG > ng.thresh ~ "No-Go",
                                                   TRUE ~ "Consider"),
                              study.end.Go = x.go, study.end.NG = x.ng,
                              go.thresh=go.thresh, ng.thresh=ng.thresh)


 return(my.df)
}

