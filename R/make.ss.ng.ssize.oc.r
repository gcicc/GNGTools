#' Make single sample normal gamma sample size OC
#'
#' @param for.plot output from get.ss.ng.ssize.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table use TRUE for verbose output
#'
#' @return a list is returned
#' @export
#' @examples \donttest{
#' my.ss.ng.ssize.oc.df <- get.ss.ng.ssize.oc.df()
#' my.ss.ng.ssize.oc <- make.ss.ng.ssize.oc(for.plot=my.ss.ng.ssize.oc.df, add.table=TRUE)
#' plot(my.ss.ng.ssize.oc)
#' }
make.ss.ng.ssize.oc <- function(for.plot=get.ss.ng.ssize.oc.df(), nlines=25, tsize=4, add.table=TRUE){

  main.plot <- ggplot() +
    geom_line(data=for.plot, aes(x=n.t, y=Go), color="green") +
    geom_line(data=for.plot, aes(x=n.t, y=1 - NoGo), color="red") +
    facet_wrap(~key, labeller="label_parsed")+
    geom_ribbon(data=for.plot, aes(x=n.t, ymin=0, ymax=Go), fill="green", alpha=.5)+
    geom_ribbon(data=for.plot, aes(x=n.t, ymin=Go, ymax=1-NoGo), fill="grey", alpha=.5)+
    geom_ribbon(data=for.plot, aes(x=n.t, ymin=1-NoGo, ymax=1), fill="red", alpha=.5)+
    labs(x="Total number of observed events",
         y="Decision Probabilities",
         title="Operating characteristics as a function of sample size",
         subtitle=TeX(paste0("Prior and nuisance parameter and standard deviation held fixed." )))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1), labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x = element_text(angle=45, hjust=1,vjust=1))


  table.plot2 <- ggplot()+
    annotate("text", label = paste0("Decision Criteria"),
             x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
    annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
             x = -1, y = 1-3/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta$ >= Min TPP) > ", for.plot$tau.lrv[1]*100, "% &")),
             x = 0, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
    annotate("text", label =  TeX(paste0("P($\\Delta$ >= Base TPP) > ", for.plot$tau.tv[1]*100,"%")),
             x = 1.5, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
    annotate("text", label = TeX(paste0("No-Go if:")), color="red",
             x = -1, y = 1-4/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta$ >= Min TPP) <= ", for.plot$tau.ng[1]*100, "% &")),
             x = 0, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0) +
    annotate("text", label = TeX(paste0("P($\\Delta$ >= Base TPP) <= ", for.plot$tau.tv[1]*100,"%")),
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

  if(add.table==TRUE) return(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)))
  if(add.table==FALSE) return(main.plot)
  }

