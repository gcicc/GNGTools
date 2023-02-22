#' @title Get single sample binary treatment OC data.frame
#' @param a.trt prior alpha parameter
#' @param b.trt prior beta parameter
#' @param n.trt observed sample size
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return A data.frame is returned.
#' @export
#'
#' @examples
#' my.ss.bin.trt.oc.df <- get.ss.bin.trt.oc.df()
#' my.ss.bin.trt.oc.df
get.ss.bin.trt.oc.df <- function(a.trt = 1, b.trt = 1, n.trt = 40,
                                 Delta.tv = 0.35, Delta.lrv = 0.2,
                                 tau.tv = 0.1, tau.lrv = 0.8, tau.ng = 0.65) {
  my.df <- expand.grid(successes=0:n.trt, Delta=seq(0,1,.01)) %>%
    mutate(prob=dbinom(successes, n.trt, Delta)) %>%
    mutate(a.post = a.trt + successes, b.post = b.trt + n.trt - successes) %>%
    mutate(P.R1 = 1 - pbeta(Delta.lrv, a.post, b.post),
           P.R3 = 1 - pbeta(Delta.tv, a.post, b.post),
           result = ifelse(P.R1 >= tau.lrv & P.R3 >= tau.tv, "Go",
                           ifelse(P.R1 < tau.ng & P.R3 < tau.tv, "No-Go",
                                  "Consider"))) %>%
    mutate(result = factor(result, c("Go", "Consider", "No-Go"))) %>%
    group_by(Delta, result, .drop=FALSE) %>%
    summarize(p=sum(prob)) %>%
    mutate(a.trt = a.trt, b.trt = b.trt, n.trt = n.trt,
           Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
           tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

  return(my.df)
}


