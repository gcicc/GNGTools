#' @title Get single sample binary data.frame
#' @description Get Go/No-go/Continue result
#' @param a.trt beta prior hyperparameter
#' @param b.trt beta prior hyperparameter
#' @param beta.mean mean of beta prior
#' @param eff.ss effective sample size of beta prior
#' @param x.trt sample responses
#' @param n.trt sample size
#' @param Delta.tv Base TPP
#' @param Delta.lrv Min TPP
#' @param tau.tv Base TPP threshold
#' @param tau.lrv Min TPP threshold
#' @param tau.ng No-Go threshold
#' @param rp logical for reparameterized beta
#'
#' @return returns a data.frame with GO/No-Go probabilities and decisions
#' @export
#'
#' @examples
#' my.ss.bin.df <- get.ss.bin.df()
#' head(my.ss.bin.df)
get.ss.bin.df <- function(
  a.trt = seq(0.5,1,0.5), b.trt = seq(0.5,1,0.5),
  beta.mean = seq(0.3,0.7,.01), eff.ss = 1:40,
  x.trt = 0:80, n.trt = c(40:80),
  Delta.tv = .4, Delta.lrv = .3,
  tau.tv = 0.10, tau.lrv = .80, tau.ng = .65,
  rp = FALSE)
{
  if(rp==FALSE) {
    my.grid <- expand.grid(n.trt = n.trt,
                           x.trt = x.trt,
                           a.prior = a.trt,
                           b.prior = b.trt,
                           Delta.tv = Delta.tv,
                           Delta.lrv = Delta.lrv,
                           tau.tv = tau.tv,
                           tau.lrv = tau.lrv,
                           tau.ng = tau.ng) %>%
      mutate(beta.mean = a.prior / (a.prior + b.prior),
             eff.ss = (a.prior + b.prior))
  } else {
    my.grid <- expand.grid(n.trt = n.trt,
                           x.trt = x.trt,
                           beta.mean = beta.mean,
                           eff.ss = eff.ss,
                           Delta.tv = Delta.tv,
                           Delta.lrv = Delta.lrv,
                           tau.tv = tau.tv,
                           tau.lrv = tau.lrv,
                           tau.ng = tau.ng) %>%
      mutate(a.prior= beta.mean * eff.ss,
             b.prior= eff.ss - beta.mean * eff.ss)
  }

  my.grid <- my.grid %>%
    mutate(a.post = a.prior + x.trt, b.post = b.prior + (n.trt - x.trt)) %>%  # Get posterior parameters
    mutate(
      P.R1 = 1 - pbeta(Delta.lrv, a.post, b.post),
      P.R3 = 1 - pbeta(Delta.tv, a.post, b.post),

      result = ifelse(P.R1 >= tau.lrv & P.R3 >= tau.tv, "Go",
                      ifelse(P.R1 < tau.ng  & P.R3 < tau.tv, "No-Go", "Consider")))

  my.grid$result <- factor(my.grid$result,
                           c("Go", "Consider", "No-Go",
                             "NA: Successes > Subjects"))
  my.grid$result[is.na(my.grid$result)==TRUE] <- "NA: Successes > Subjects"
  return(my.grid)
}


