#' @title Make single sample binary treatment oc curve version 2
#'
#' @param my.df output from get.ss.bin.trt.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param nlines.ria Control for text spacing
#' @param Delta_OC_LB Lower bound for OC curve
#' @param Delta_OC_UB Upper bound for OC curve
#' @param add.table provides extended output summaries

#' @return A ggplot object is returned.
#' @export
#'
#' @examples
#' make.ss.bin.trt.oc2()
make.ss.bin.trt.oc2 <- function(my.df = get.ss.bin.trt.oc.df(), tsize=4, nlines=25,
                                nlines.ria=20, Delta_OC_LB = 0, Delta_OC_UB = 1, add.table=TRUE){

  Delta.tv <- my.df$Delta.tv[1]
  Delta.lrv <- my.df$Delta.lrv[1]
  tau.tv <- my.df$tau.tv[1]
  tau.lrv <- my.df$tau.lrv[1]
  tau.ng <- my.df$tau.ng[1]
  my.df$result <- factor(my.df$result, c("Go", "Consider", "No-Go"))
  total.n <- my.df$n.trt[1]

  main.plot <- ggplot(data=my.df, aes(x=Delta, y=p, color=result))+
    geom_line(size=.75)
  dpb <- ggplot_build(main.plot)
  main.plot <- main.plot +
    scale_color_manual(values=c("green", "grey50", "red"))+
    scale_x_continuous(expand = c(0,0), breaks=seq(0,1,.1), limits=c(Delta_OC_LB, Delta_OC_UB)) +
    scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2),
                       limits=c(0,1), labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"),
          axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
    geom_vline(xintercept = c(Delta.tv, Delta.lrv), color="blue", linetype=2)+
    annotate("text", label = TeX(paste0("Min TPP = ", Delta.lrv*100,"%")),
             x = Delta.lrv, y = 0 + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria, size = tsize,
             colour = "black", hjust = 1)+
    annotate("text", label = TeX(paste0("Base TPP = ",  Delta.tv*100,"%")),
             x = Delta.tv, y = 0 + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria, size = tsize,
             colour = "black", hjust = 0)+
    labs(x = TeX("$\\Delta\\,$ = Underlying proportion of treatment responders"),
         y = "Probability",
         title = "Operating characteristics as a function of treatment effect",
         color = "Decision",
         subtitle = paste0("Total randomized to treatment: ", total.n)) +
    theme(legend.position = "bottom")

  table.plot2 <- ggplot()+
    annotate("text", label = paste0("Decision Criteria"),
             x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
    annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
             x = -1, y = 1-3/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta\\,$ >= Min TPP) > ", tau.lrv*100, "% &")),
             x = 0, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
    annotate("text", label =  TeX(paste0("P($\\Delta\\,$ >= Base TPP) > ", tau.tv*100,"%")),
             x = 1.5, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
    annotate("text", label = TeX(paste0("No-Go if:")), color="red",
             x = -1, y = 1-4/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta\\,$ >= Min TPP) <= ", tau.ng*100, "% &")),
             x = 0, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0) +
    annotate("text", label = TeX(paste0("P($\\Delta\\,$ >= Base TPP) <= ", tau.tv*100,"%")),
             x = 1.5, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0)+
    annotate("text", label = TeX(paste0("Consider:")), color="black",
             x = -1, y = 1-5/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = "Otherwise",
             x = 0, y = 1-5/nlines, size = tsize+1, colour = "black", hjust = 0)+
    scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
    scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
    labs(x="\n", y="")

  if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2))
  if(add.table==FALSE) return(main.plot)
  }


