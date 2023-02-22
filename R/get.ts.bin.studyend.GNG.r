#' @title Get Two-sample binary study end GNG
#'
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param n.con number of trials for control group
#' @param x.con number of responses for control group
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior beta parameter for treatment group
#' @param n.trt number of trials for treatment group
#' @param x.trt number of responses for treatment group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return returns a list holding data.frames alterting to what is needed from data to acheive go/no-go
#' @export
#'
#' @examples
#' get.ts.bin.studyend.GNG()
get.ts.bin.studyend.GNG <- function(a.con = 1, b.con = 1, n.con = 30, x.con = 18,
                                    a.trt = 1, b.trt = 1, n.trt = 30, x.trt = 14,
                                    Delta.lrv = 0.3, Delta.tv = .4,
                                    tau.tv = 0.10, tau.lrv = 0.80, tau.ng = 0.65){

        results <- get.ts.bin.dec.df(a.con = a.con, b.con = b.con, n.con = n.con, x.con = x.con,
                                          a.trt = a.trt, b.trt = b.trt, n.trt = n.trt, x.trt = 0:n.trt,
                                          Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                                          tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

        # Determine mi/max number of TRT responders for Go/No G0 ----
        result.go <- results %>% dplyr::filter(result== "Go") %>% dplyr::slice(1)
        result.ng <- results %>% dplyr::filter(result== "No-Go") %>% dplyr::slice(dplyr::n())

        if(nrow(result.go) == 0) result.go <- result.ng %>% mutate(x.con=x.con, x.trt = NA, P.R1=NA, P.R3=NA, result="Not possible")
        if(nrow(result.ng) == 0) result.ng <- result.go %>% mutate(x.con=x.con, x.trt = NA, P.R1=NA, P.R3=NA, result="Not possible")


        return(list(result.go = result.go, result.ng = result.ng))
}
