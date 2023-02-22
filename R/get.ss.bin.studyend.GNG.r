#' @title Get single sample binary study end GNG criteria
#' @param a.trt prior alpha parameter
#' @param b.trt prior beta parameter
#' @param n.trt observed sample size
#' @param x.trt observed number of responders
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return Returns a list of data.frames holding what is needed from data to acheive Go/No-Go
#' @export
#'
#' @examples
#' get.ss.bin.studyend.GNG()
get.ss.bin.studyend.GNG <- function(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
                                    Delta.lrv = .2, Delta.tv = .35,
                                    tau.tv = 0.10, tau.lrv = .80, tau.ng = .65){
  results <- get.ss.bin.df(a.trt = a.trt, b.trt = b.trt, n.trt = n.trt, x.trt = 0:n.trt,
                           Delta.lrv = Delta.lrv, Delta.tv=Delta.tv, tau.lrv=tau.lrv,
                           tau.tv=tau.tv, tau.ng=tau.ng)

  result.go <- results %>% dplyr::filter(result== "Go") %>% dplyr::slice(1)
  result.ng <- results %>% dplyr::filter(result== "No-Go") %>% dplyr::slice(dplyr::n())
  return(list(result.go = result.go, result.ng = result.ng))
}
