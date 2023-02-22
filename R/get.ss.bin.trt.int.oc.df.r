#' @title Get single sample binary interim treatment OC data.frame
#' @param a.trt prior alpha parameter
#' @param b.trt prior beta parameter
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param interim.n.t number of trials at interim
#' @param final.n.t number of trials at final
#' @param x.ng responses needed for no-go; leave null for standard rule
#' @param x.go responses needed for go; leave null for standard rule
#' @param go.thresh go threshold for predictive probability
#' @param ng.thresh no-go threshold for predictive probability
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' my.ss.bin.trt.int.oc.df <- get.ss.bin.trt.int.oc.df()
#' my.ss.bin.trt.int.oc.df[[1]]
#' my.ss.bin.trt.int.oc.df[[2]]
#' }
get.ss.bin.trt.int.oc.df <- function(a.trt = 1, b.trt = 1,
                                     Delta.tv = 0.35, Delta.lrv = 0.2,
                                     tau.tv = 0.1, tau.lrv = 0.8, tau.ng = 0.65,
                                     interim.n.t = c(10), final.n.t = 100,
                                     x.ng = NULL, x.go=NULL,
                                     go.thresh=0.8, ng.thresh=0.8){

# This will call return.ss.bin.pp for a spectrum of treatment effects, group by outcome and sum up probabilities
my.df <- bind_rows(
  apply(X = matrix(seq(0,1,.01)), MARGIN = 1, function(x) {
    return.ss.bin.int.req(
      a.trt = a.trt, b.trt = b.trt,
      Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
      tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
      interim.n.t = interim.n.t, final.n.t = final.n.t,
      x.ng = x.ng, x.go=x.go,
      go.thresh=go.thresh, ng.thresh=ng.thresh, p.success = x) %>%
      group_by(interim.n.t, decision) %>%
      summarize(prob=sum(P.interim2)) %>%
      mutate(p.success=x)
  }
  ))

# This function returns the data.frame associated with study end treatment operating characteristic curve
my.df2 <- get.ss.bin.trt.oc.df(a.trt = a.trt, b.trt = b.trt, n.trt = final.n.t,
                               Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                               tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng) %>% rename(decision=result)

return(list(my.df, my.df2))
}



