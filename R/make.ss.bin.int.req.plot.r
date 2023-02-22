#' Make single sample binary predictive probability plot
#' @title Make single sample binary predictive probability plot
#' @param my.table output from make.ss.bin.pp.table
#'
#' @return A ggplot object is returned.
#' @export
#'
#' @examples
#' my.ss.bin.int.data.req.plot<- make.ss.bin.int.req.plot()
#' my.ss.bin.int.data.req.plot

make.ss.bin.int.req.plot <- function(my.table = make.ss.bin.int.req.table()){

ggplot(data=my.table, aes(x=interim.n.t, y= prop.needed, color=decision)) + geom_point() + geom_line()+
  scale_color_manual(values=c("green", "red"))+
    labs(x="Number of observations at interim",
         y="Observed sample proportion needed to meet PPP Thresholds",
         title="Plot of number of observations vs. observed sample proportion",
         subtitle="Reported: Sample proportion required to meet/exceed Post Pred Prob thresholds for Go/No-Go",
         caption="Dotted horizontal lines are the study end Go/NoGo cutoffs")+
  geom_hline(yintercept = c(my.table$study.end.Go[1]/(my.table$interim.n.t[1] + my.table$n.comp[1])), linetype=2, color="green") +
  geom_hline(yintercept = c(my.table$study.end.NG[1]/(my.table$interim.n.t[1] + my.table$n.comp[1])), linetype=2, color="red")  +
  guides(color="none")
}


