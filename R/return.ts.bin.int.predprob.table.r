#' @title Return two-sample binary predictive probability table
#'
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior alpha parameter for treatment group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param n.int.con sample size for control group at interim
#' @param n.int.trt sample size for treatment group at interim
#' @param n.trt sample size for control group at interim
#' @param n.con sample size for treatment group at interim
#' @param p.con probability of success for control group
#' @param p.trt probability of success for treatment group
#' @param studyend keep null
#' @param goparallel Option to use parallel computing
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.bin.int.predprob.table <- return.ts.bin.int.predprob.table()
#' }
#' @author Greg Cicconetti
return.ts.bin.int.predprob.table <- function(a.con = .5, b.con = .5, a.trt = .5, b.trt = .5, # Priors
                                             Delta.lrv = .08, Delta.tv = .08,
                                             tau.tv = .7, tau.lrv = .7, tau.ng = 0, # Decision rule components
                                             n.int.con = 5, n.int.trt = 5, # Interim data
                                             n.trt = 10, n.con = 10, # Final sample sizes
                                             p.con=.2, p.trt=.2 + seq(0, .5,.1),
                                             studyend = NULL, goparallel=FALSE){

        if(is.null(studyend) == T) {
                # NOTE: If studyend gets used repeatedly, consider passing this as an argument since it takes a few seconds to compute
                studyend <- return.ts.bin.studyend.GNG.hm.df(a.con = a.con, b.con = b.con, n.con = n.con,
                                                             a.trt = a.trt, b.trt = b.trt, n.trt = n.trt,
                                                             Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                                             tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                                                             x_ng = x_ng, x_go=x_go,
                                                             go.thresh=go.thresh, ng.thresh=ng.thresh)
        }

        my.specs <- expand.grid(x.int.con=0:n.int.con, x.int.trt = 0:n.int.trt) %>%
                mutate(n.int.con = n.int.con, n.int.trt = n.int.trt, n.trt = n.trt, n.con = n.con)

        if(goparallel==F){
                holdit <- bind_rows(
                        apply(X = matrix(1:nrow(my.specs)), MARGIN = 1, FUN = function(x){
                                return.ts.bin.int.predprob(a.con = a.con[1],
                                                 b.con = b.con[1],
                                                 a.trt = a.trt[1],
                                                 b.trt = b.trt[1], # Priors
                                                 Delta.lrv = Delta.lrv[1],
                                                 Delta.tv = Delta.tv[1],
                                                 tau.tv = tau.tv[1],
                                                 tau.lrv = tau.lrv[1],
                                                 tau.ng = tau.ng[1], # Decision rule components
                                                 n.int.con = my.specs$n.int.con[x],
                                                 n.int.trt = my.specs$n.int.trt[x],
                                                 x.int.con=my.specs$x.int.con[x],
                                                 x.int.trt = my.specs$x.int.trt[x], # Interim data
                                                 n.final.trt = my.specs$n.trt[x],
                                                 n.final.con = my.specs$n.con[x])
                        }))
        } else {
                parallel::clusterEvalQ(cl = cl, expr = {require(tidyverse);
                        require(extraDistr)})
                parallel::clusterExport(cl = cl, varlist = c("my.specs","return.ts.bin.int.predprob", "studyend",
                                                             "a.con" , "b.con" , "a.trt" , "b.trt", "n.con","n.trt",
                                                             "Delta.lrv" , "Delta.tv" , "tau.tv" , "tau.lrv" ,
                                                             "tau.ng" ), envir=environment())
                holdit <- bind_rows(
                        parallel::parApply(cl=cl,X = matrix(1:nrow(my.specs)), MARGIN = 1, FUN = function(x){
                                return.ts.bin.int.predprob(a.con = a.con[1],
                                                 b.con = b.con[1],
                                                 a.trt = a.trt[1],
                                                 b.trt = b.trt[1], # Priors
                                                 Delta.lrv = Delta.lrv[1],
                                                 Delta.tv = Delta.tv[1],
                                                 tau.tv = tau.tv[1],
                                                 tau.lrv = tau.lrv[1],
                                                 tau.ng = tau.ng[1], # Decision rule components
                                                 n.int.con = my.specs$n.int.con[x],
                                                 n.int.trt = my.specs$n.int.trt[x],
                                                 x.int.con=my.specs$x.int.con[x],
                                                 x.int.trt = my.specs$x.int.trt[x], # Interim data
                                                 n.final.trt = my.specs$n.trt[x],
                                                 n.final.con = my.specs$n.con[x])
                        })
                )


        }

        return(holdit)

}
