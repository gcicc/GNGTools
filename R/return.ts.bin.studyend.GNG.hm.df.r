#' @title Return two-sample binary study end GNG heatmap data.frame
#'
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param n.con sample size for control group
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior beta parameter for treatment group
#' @param n.trt sample size for control treatment group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param x_ng responses needed for no-go; leave null for standard rule
#' @param x_go responses needed for go; leave null for standard rule
#' @param go.thresh go threshold for predictive probability
#' @param ng.thresh no-go threshold for predictive probability
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.ts.bin.studyend.GNG.hm.df <- return.ts.bin.studyend.GNG.hm.df()
#' head(my.ts.bin.studyend.GNG.hm.df)
return.ts.bin.studyend.GNG.hm.df <- function(a.con = 1, b.con = 1, n.con = 30,
                                             a.trt = 1, b.trt = 1, n.trt = 30,
                                             Delta.lrv = 0.3, Delta.tv = .4,
                                             tau.tv = 0.10, tau.lrv = 0.80, tau.ng = 0.65,
                                             x_ng = NULL, x_go=NULL,
                                             go.thresh=0.8, ng.thresh=0.8){

        # For every possible final number of success on the control arm
        bind_rows(
                apply(X = matrix(0:n.con), MARGIN = 1, FUN=function(x) {

                        results <- get.ts.bin.studyend.GNG(a.con = a.con, b.con = b.con, n.con = n.con, x.con = x,
                                                           a.trt = a.trt, b.trt = b.trt, n.trt = n.trt, x.trt = floor(n.trt/2),
                                                           Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                                           tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)
                        # If we are in a case where either go or no/go is not possible, fill with missing

                        results$result.go <- results$result.go %>% rename(x.trt.go=x.trt, P.R1.Go = P.R1, P.R3.Go=P.R3, result.Go=result )
                        results$result.ng <- results$result.ng %>% rename(x.trt.ng=x.trt, P.R1.NG = P.R1, P.R3.NG=P.R3, result.NG=result )
                        results$result.go %>% left_join(results$result.ng) %>% dplyr::select(
                                n.con, n.trt,a.con, a.trt, b.con, b.trt, Delta.tv, Delta.lrv, tau.tv, tau.lrv, tau.ng,  x.con, x.trt.ng,  x.trt.go,   P.R1.Go, P.R3.Go, result.Go, P.R1.NG, P.R3.NG, result.NG
                        )
                } )  )

}
