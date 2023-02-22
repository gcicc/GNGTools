#' @title Make time to event predictive probability cutoff plot
#'
#' @param my.table output from make.tte.int.data.table
#' @param include_nogo logical
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' my.tte.int.data.table <- make.tte.int.data.table()
#' my.tte.int.data.plot <- make.tte.int.data.plot(my.table = my.tte.int.data.table)
#' my.tte.int.data.plot

make.tte.int.data.plot <- function(my.table = make.tte.int.data.table(), include_nogo=FALSE){
  my.table$ARatio <- 1/my.table$ARatio - 1
  caption.text = paste0(
                        "Prior Control Events = ", my.table$m.con.prior[1], "; ",
                        "Prior Treatment Events = ", my.table$m.trt.prior[1], "; ",
                        "Prior HR estimate = ",  my.table$HR.prior[1], "; ",
                        "Control:TRT Ratio = 1:", round(my.table$ARatio[1], 3), "\n",
                        "Dotted lines indicate Study End HR needed for (Go, No-Go): ",
                        "(",  my.table$study.end.Go[1], ", ",
                        "",  my.table$study.end.NG[1], ")"
  )
  # Reverting back to original scale (control:TRT = 1:TRT)


 my.plot<- my.table %>% ggplot() +
    geom_line(aes(x=interim.m, y=HR.Go), color="green")+
    geom_point(aes(x=interim.m, y=HR.Go), color="green")+
    # geom_line(aes(x=interim.m, y=HR.NG), color="red") +
    # geom_point(aes(x=interim.m, y=HR.NG), color="red")+
    labs(x="Number of events", y="Observed HR required to hit PPP Thresholds",
         title="Plot of number of events vs. observed hazard ratios",
         subtitle="Reported: Hazard ratios required to meet/exceed Post Pred Prob thresholds for Go/No-Go",
         caption=caption.text)+
    geom_hline(yintercept=my.table$study.end.Go[1], color="green", linetype=2)
   # geom_hline(yintercept=my.table$study.end.NG[1], color="red", linetype=2)

  if(include_nogo==TRUE)  my.plot <- my.plot +
    geom_line(aes(x=interim.m, y=HR.NG), color="red") +
    geom_point(aes(x=interim.m, y=HR.NG), color="red")+
    geom_hline(yintercept=my.table$study.end.NG[1], color="red", linetype=2)

 return(my.plot)

}
