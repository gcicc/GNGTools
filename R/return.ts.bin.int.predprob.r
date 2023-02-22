#' Return two-sample binary interim predictive probability
#'
#' @param a.con alpha parameter for control
#' @param b.con beta parameter for control
#' @param a.trt alpha parameter for trt
#' @param b.trt beta parameter for trt
#' @param Delta.lrv TPP info
#' @param Delta.tv TPP info
#' @param tau.tv Threshold
#' @param tau.lrv Threshold
#' @param tau.ng Threshold
#' @param n.int.con Sample size, control group at interim
#' @param n.int.trt Sample size, treatment group at interim
#' @param x.int.con Responders, control group at interim
#' @param x.int.trt Responders, treatment group at interim
#' @param n.final.trt Final sample size treatment
#' @param n.final.con Final sample size control
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' return.ts.bin.int.predprob()
#' }
return.ts.bin.int.predprob <- function(a.con = 1, b.con = 1,a.trt = 1, b.trt = 1,
                             Delta.lrv = 0.3, Delta.tv = .4, tau.tv = 0.10, tau.lrv = 0.80, tau.ng = 0.8,
                             n.int.con = 20, n.int.trt = 20, x.int.con=17, x.int.trt = 17,
                             n.final.trt = 40, n.final.con = 40){

        # The number of success on treatment arm needed for Go/No-Go depends on the number of success observed on control arm
        # NOTE: x.con and x.trt are not doing anything in return.ts.bin.studyend.GNG.hm.df - consider removing
        studyend <- return.ts.bin.studyend.GNG.hm.df(a.con = a.con, b.con = b.con, n.con = n.final.con,
                                                     a.trt = a.trt, b.trt = b.trt, n.trt = n.final.trt,
                                                     Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                                     tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,

                                                     x_ng = x_ng, x_go=x_go,
                                                     go.thresh=go.thresh, ng.thresh=ng.thresh)

        # Suppose we are at an interim with n.int.con = 20, n.int.trt = 20, x.int.con=5, x.int.trt = 10
        # There are (n.final.con - n.int.con) and (n.final.trt - n.int.trt) more bernoulli trials
        # possible number of final successes on CON: x.int.con + 0:(n.final.con - n.int.con)
        # possible number of final successes on TRT: x.int.trt + 0:(n.final.trt - n.int.trt)
        # Probability that we observed x.int.con among n.int.con: dbbinom(x=x.int.con, n = n.int.con, alpha=a.con, beta=beta.con)
        # Probability that we observed x.int.trt among n.int.trt: dbbinom(x=x.int.trt, n = n.int.trt, alpha=a.trt, beta=beta.trt)
        # Given interim data, prob that we observe k success among remaining CON trials:
        # dbbinom(x=k, n = (n.final.con - n.int.con), alpha=a.con+x.int.con, beta=beta.con + (n.int.con - x.int.con))
        # Given interim data, prob that we observe k success among remaining trt trials:
        # dbbinom(x=k, n = (n.final.trt - n.int.trt), alpha=a.trt+x.int.trt, beta=beta.trt + (n.int.trt - x.int.trt))

        # Suppose we are at an interim with n.int.con = 20, n.int.trt = 20, x.int.con=5, x.int.trt = 10
        # For each possible outcome on CON group we have the studyend Go and NOgo criteria
        # Suppose we have 5 successes at interim. We can end with 5, 6, .., 20 success on control


        # return all possible ways to complement the control arm in part 2 (x.c.comp) along with
        # the probability based on dbbinom update and provide studyend final count on control

        build <-  data.frame(x.int.con = x.int.con, n.int.con=n.int.con,
                             P.x.int.con = dbbinom(x=x.int.con, size = n.int.con, alpha=a.con, beta=b.con),
                             n.final.con=n.final.con,
                             x.c.comp = 0:(n.final.con - n.int.con))%>%
                mutate(P.x.c.comp = dbbinom(x=x.c.comp, size = (n.final.con - n.int.con), alpha=a.con+x.int.con, beta=b.con + (n.int.con - x.int.con)),
                       x.con=x.int.con+x.c.comp)

        # Now tack on what's required at study end for Go and No-go
        build <-  build %>% left_join(
                studyend %>% dplyr::filter(x.con %in% build$x.con ) %>% select(x.con, x.trt.go, x.trt.ng)
        ) %>%
                mutate(x.int.trt = x.int.trt,
                       n.int.trt = n.int.trt,
                       n.final.trt=n.final.trt,
                       P.int.t = dbbinom(x=x.int.trt, size = n.int.trt, alpha=a.trt, beta=b.trt))

        # Now for each row of build: compute the posterior predictive probability of reaching study-end Go criteria
        # If not possible based on rule, return 0.  If there are already enough, return 1. If we can't reach due to limited number of trials return 0,
        # else
        build$P.Go <- apply(X = matrix(1:nrow(build)), MARGIN = 1, FUN =
                                    function(x){
                                            case_when(is.na(build$x.trt.go[x]) == T ~ 0,         # not possible
                                                      is.na(build$x.trt.go[x]) == F &
                                                              (build$x.int.trt[x] > build$x.trt.go[x]) ~ 1,                    # Already have enought
                                                      is.na(build$x.trt.go[x]) == F &
                                                              (build$x.int.trt[x] + build$n.int.trt[x] < build$x.trt.go[x]) ~ 0,  # Can't get enough
                                                      TRUE ~
                                                              ifelse(is.na(build$x.trt.go[x]) == T, 0,
                                                                     sum(
                                                                             dbbinom(x= (build$x.trt.go[x]- build$x.int.trt[x]):(build$n.final.trt[x] - build$n.int.trt[x]),  # These will give a Go
                                                                                     size = build$n.final.trt[x] - build$n.int.trt[x], # These are trials left
                                                                                     alpha=a.trt+build$x.int.trt[x], beta=b.trt + (build$n.int.trt[x] - build$x.int.trt[x])))
                                                              ))
                                    })
        # Now for each row of build: compute the post pred prob of reaching study-end No-go
        build$P.NG <- apply(X = matrix(1:nrow(build)), MARGIN = 1, FUN =
                                    function(x){
                                            case_when(is.na(build$x.trt.ng[x]) == T ~ 0,         # not possible
                                                      is.na(build$x.trt.ng[x]) == F &
                                                              (build$x.int.trt[x] > build$x.trt.ng[x]) ~ 0,  # Already exceed what's needed for nogo
                                                      is.na(build$x.trt.go[x]) == F &
                                                              (build$x.int.trt[x] <= build$x.trt.ng[x]) ~  # We still have room
                                                              ifelse(is.na(build$x.trt.ng[x]) == T, 0,
                                                                     sum(
                                                                             dbbinom(x= 0:(build$x.trt.ng[x] - build$x.int.trt[x]),  # These will give a Go
                                                                                     size = build$n.final.trt[x] - build$n.int.trt[x], # These are trials left
                                                                                     alpha=a.trt+build$x.int.trt[x], beta=b.trt + (build$n.int.trt[x] - build$x.int.trt[x])))
                                                              ))
                                    })

        # Now the final posterior predictive prob of reaching Go/No Go is the product of   P(A AND B) = P(A)P(B|A)
        # a: Obtaining required specific complement on control arm: P.x.c.comp
        # b: the probability of reaching Go given that control arm
        build %>% group_by(x.int.con, x.int.trt) %>% summarize(P.Go = sum( P.x.c.comp * P.Go, na.rm=T),
                                                               P.NG = sum( P.x.c.comp * P.NG, na.rm=T)) %>%
                mutate(a.con = a.con, b.con = b.con,
                       a.trt = a.trt, b.trt = b.trt,
                       Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                       tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                       n.int.trt = n.int.trt, n.final.trt = n.final.trt, n.final.con = n.final.con) %>%
                select(a.con, b.con, a.trt, b.trt, Delta.lrv, Delta.tv, tau.tv, tau.lrv, tau.ng, n.int.trt,x.int.con, x.int.trt, n.final.trt, n.final.con, P.Go, P.NG)

}
