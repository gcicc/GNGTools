#' @title Make time-to-event sample size OC curve
#' @param for.plot output from get.tte.ssize.oc.df
#'
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provide GNG Rule tables and extensive output
#'
#' @export
#' @return A ggplot object is returned
#' @examples \donttest{
#' my.tte.ssize.oc.df <- get.tte.ssize.oc.df()
#' make.tte.ssize.oc(for.plot = my.tte.ssize.oc.df)
#' }
make.tte.ssize.oc <- function(for.plot=get.tte.ssize.oc.df(), tsize=4, nlines=25, add.table=TRUE){

  # Simply takes a dataframe from get.tte.ssize.oc.df and plots
  main.plot <- ggplot() + geom_line(data=for.plot, aes(x=m.obs, y=Go), color="green") +
    geom_line(data=for.plot, aes(x=m.obs, y=1 - NoGo), color="red") +
    facet_wrap(~key)+
    geom_ribbon(data=for.plot, aes(x=m.obs, ymin=0, ymax=Go), fill="green", alpha=.5)+
    geom_ribbon(data=for.plot, aes(x=m.obs, ymin=Go, ymax=1-NoGo), fill="grey", alpha=.5)+
    geom_ribbon(data=for.plot, aes(x=m.obs, ymin=1-NoGo, ymax=1), fill="red", alpha=.5)+
    labs(x="Total number of observed events", y="Decision Probabilities",
         title="Operating Characteristics as a function of number of events")+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1), labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x = element_text(angle=45, hjust=1,vjust=1))

  table.plot2 <- ggplot()+
    annotate("text", label = paste0("Decision Criteria"),
             x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
    annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
             x = -1, y = 1-2.5/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta$ > $HR_{Max}$) > ", for.plot$tau.lrv[1]*100, "% &")),
             x = 0, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
    annotate("text", label =  TeX(paste0("P($\\Delta$ > $HR_{Base}$) > ", for.plot$tau.tv[1]*100,"%")),
             x = 1.5, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
    annotate("text", label = TeX(paste0("No-Go if:")), color="red",
             x = -1, y = 1-3.25/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta$ > $HR_{Max}$) <= ", for.plot$tau.ng[1]*100, "% &")),
             x = 0, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0) +
    annotate("text", label =  TeX(paste0("P($\\Delta$ > $HR_{Base}$) <= ", for.plot$tau.tv[1]*100,"%")),
             x = 1.5, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0)+
    annotate("text", label = TeX(paste0("Consider:")), color="black",
             x = -1, y = 1-4./nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = "Otherwise",
             x = 0, y = 1-4./nlines, size = tsize+1, colour = "black", hjust = 0)+
    scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
    scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
    labs(x="", y="")

  # plot_grid(main.plot,
  #            table.plot2, ncol=1, align = "v", rel_heights=c(.78,.22), axis="l")

  if(add.table==TRUE) return(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)))
  if(add.table==FALSE) return(main.plot)
}

