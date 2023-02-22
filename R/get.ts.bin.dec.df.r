#' @title Get Two-sample binary decision data.frame
#'
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param n.con sample size for control
#' @param x.con responders on control
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior beta parameter for treatment group
#' @param n.trt sample size for treatment
#' @param x.trt responders for treatment
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#'
#' @return a dataframe is returned
#' @export
#'
#' @examples
#' holdit <- get.ts.bin.dec.df()
#' head(holdit)

get.ts.bin.dec.df = function(a.con = 1, b.con = 1, n.con = 40, x.con = 0:40,
                                  a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 0:40,
                                  Delta.tv = .25, Delta.lrv = .2,
                                  tau.tv = 0.10, tau.lrv = 0.80, tau.ng = 0.65)
{

        # create simulation grid

        my.grid <- expand.grid(n.con = n.con, n.trt = n.trt,
                               x.con = x.con,
                               x.trt = x.trt,
                               a.con= a.con, a.trt = a.trt,
                               b.con = b.con, b.trt = b.trt,
                               Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                               tau.tv = tau.tv,tau.lrv = tau.lrv,tau.ng = tau.ng)

        # Define function to be run on each row of grid
        my.function1 <- function(n.con = my.grid$n.con[1], n.trt = my.grid$n.trt[1], x.con = my.grid$x.con[1], x.trt = my.grid$x1[1],
                                 a.con = my.grid$a.con[1], a.trt = my.grid$a.trt[1], b.con = my.grid$b.con[1], b.trt = my.grid$b1[1],
                                 Delta.tv = my.grid$Delta.tv[1], Delta.lrv = my.grid$Delta.lrv[1],
                                 tau.tv = my.grid$tau.tv[1], tau.lrv = my.grid$tau.lrv[1], tau.ng = my.grid$tau.ng[1]){

                return(get.ts.bin.dec(x.con = x.con, x.trt = x.trt, n.con = n.con, n.trt = n.trt, a.con = a.con,
                                           a.trt = a.trt, b.con = b.con, b.trt = b.trt,
                                           Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                                           tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng))
        }

        # Make call for each row of grid
        listOfDataFrames <- do.call(Map, c(f = my.function1, my.grid))
        # Collapse list of data.frames into a single data.frame
        df <- do.call("rbind", listOfDataFrames)
        my.grid <- cbind(my.grid, df)
        my.grid$result <- factor(my.grid$result, c("Go", "Consider", "No-Go"))
        return(my.grid)
}
