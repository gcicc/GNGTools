#' @title Make single sample normal-gamma treatment OC curve
#'
#' @param my.df output from get.ss.ng.trt.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provides extended output summaries

#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' my.ss.ng.trt.oc.df <- get.ss.ng.trt.oc.df(
#' mu.0.t = 0, n.0.t = 10, alpha.0.t = 0.25, beta.0.t = 1,
#' s.t = 2, n.t = 40, from.here = 0, to.here =4, length.out=10,
#' Delta.tv = 1.75, Delta.lrv = 1,
#' tau.tv = 0.1, tau.lrv = 0.8, tau.ng = 0.65)
#' my.ss.ng.trt.oc1 <- make.ss.ng.trt.oc1(my.df =my.ss.ng.trt.oc.df, add.table=TRUE)
#' plot(my.ss.ng.trt.oc1[[1]])
make.ss.ng.trt.oc1 <- function(my.df, nlines=25, tsize=4, add.table=TRUE){


  my.plot <- my.df  %>%
    ggplot() + geom_line(aes(x = xbar, y = Go), color="green")+
    geom_line(aes(x = xbar, y = 1-NoGo), color="red")
  dpb <- ggplot_build(my.plot)

  main.plot <- my.plot+
    geom_ribbon(data = my.df, aes(x = xbar, ymin=0, ymax=Go), fill=alpha("lightgreen", .5))+
    geom_ribbon(data = my.df, aes(x = xbar, ymin=Go, ymax=1-NoGo), fill=alpha("grey", .5))+
    geom_ribbon(data = my.df, aes(x = xbar, ymin=1-NoGo, ymax=1), fill=alpha("red", .5))+
    geom_line(data = my.df, aes(x = xbar, y = Go), color="green", size=.75)+
    geom_line(data = my.df, aes(x = xbar, y = 1-NoGo), color="red", size=.75)+
    geom_vline(xintercept = c(my.df$Delta.tv[1], my.df$Delta.lrv[1]), color="blue", linetype=2)+
    annotate("text", label = paste0("Min TPP = ", my.df$Delta.lrv[1]), x = my.df$Delta.lrv[1],
             y = 0 + max(dpb$data[[2]]$y,dpb$data[[2]]$y)/nlines, size = tsize, colour = "black", hjust = 1)+
    annotate("text", label = paste("Base TPP = ", my.df$Delta.tv[1]), x = my.df$Delta.tv[1],
             y = 0 + max(dpb$data[[2]]$y,dpb$data[[2]]$y)/nlines, size = tsize,
             colour = "black", hjust = 0)+
    labs(x = "Underlying treatment effect",
         y="Probability")+
    scale_x_continuous(expand = c(0,0), breaks=pretty(my.df$xbar, 10))+
    scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2), limits=c(0,1), labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
    labs(title="Operating characteristics as a function of treatment effect")

  table.plot2 <- ggplot()+
    annotate("text", label = paste0("Decision Criteria"),
             x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
    annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
             x = -1, y = 1-3/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta\\,$ <= Min TPP) > ", my.df$tau.lrv[1]*100, "% &")),
             x = 0, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
    annotate("text", label =  TeX(paste0("P($\\Delta\\,$ <= Base TPP) > ", my.df$tau.tv[1]*100,"%")),
             x = 1.5, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
    annotate("text", label = TeX(paste0("No-Go if:")), color="red",
             x = -1, y = 1-4/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta\\,$ <= Min TPP) >= ", my.df$tau.ng[1]*100, "% &")),
             x = 0, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0) +
    annotate("text", label = TeX(paste0("P($\\Delta\\,$ <= Base TPP) >= ", my.df$tau.tv[1]*100,"%")),
             x = 1.5, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0)+
    annotate("text", label = TeX(paste0("Consider:")), color="black",
             x = -1, y = 1-5/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = "Otherwise",
             x = 0, y = 1-5/nlines, size = tsize+1, colour = "black", hjust = 0)+
    scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
    scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
    labs(x="\n", y="")
  # plot_grid(main.plot,
  #            table.plot2, ncol=1, align = "v", rel_heights=c(.78,.22), axis="l")

  if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2))
  if(add.table==FALSE) return(main.plot)
}
