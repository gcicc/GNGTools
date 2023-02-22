#' @title Make single sample normal-gamma interim data requirement plot
#'
#' @param my.table output from make.ss.ng.pp.table
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' my.ss.ng.int.req.table <- make.ss.ng.int.req.table()
#' my.ss.ng.int.req.plot <- make.ss.ng.int.req.plot(my.table = my.ss.ng.int.req.table)
#' my.ss.ng.int.req.plot

make.ss.ng.int.req.plot <- function(my.table = make.ss.ng.int.req.table()){

  my.table %>% ggplot() +
    geom_line(aes(x=interim.n.t, y=xbar.go), color="green")+
    geom_point(aes(x=interim.n.t, y=xbar.go), color="green")+
    geom_line(aes(x=interim.n.t, y=xbar.ng), color="red") +
    geom_point(aes(x=interim.n.t, y=xbar.ng), color="red")+
    labs(x="Number of observations at interim",
         y="Observed sample mean",
         title="Plot of number of observations vs. observed sample mean",
         subtitle="Reported: Sample means required to meet/exceed posterior predictive probability thresholds for Go/No-Go",
         caption="Dotted horizontal lines are the study end Go/NoGo cutoffs") +
    geom_hline(yintercept=my.table$study.end.Go[1], color="green", linetype=2) +
    geom_hline(yintercept=my.table$study.end.NG[1], color="red", linetype=2)
}

