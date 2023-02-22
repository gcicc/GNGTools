#' Make two-sample normal gamma sample size OC curve
#'
#' @param for.plot output from get.ts.ng.ssize.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provides extended output summaries
#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.ng.ssize.oc.df <- get.ts.ng.ssize.oc.df(goparallel=FALSE)
#' my.ts.ng.ssize.oc <- make.ts.ng.ssize.oc()
#' my.ts.ng.ssize.oc
#' }
make.ts.ng.ssize.oc <- function(for.plot=get.ts.ng.ssize.oc.df(), tsize=4, nlines=25, add.table=TRUE){

        # Simply takes a dataframe from return.ssize.df and plots
        main.plot <-  for.plot %>%
                ggplot(aes(x=n.total, y=value, color=result)) +
                geom_line(alpha=.25, size=.75) +
                facet_wrap(~treatment.effect, labeller="label_parsed")+
                scale_color_manual(values=c("green", "red"))+
                geom_ribbon(data=for.plot %>% dplyr::filter(result=="Go"), aes(x=n.total, ymin=0, ymax=value), fill=alpha("green",.5), color=alpha("green",.25))+
                geom_ribbon(data=for.plot %>% dplyr::filter(result=="Go"), aes(x=n.total, ymin=value, ymax=1), fill=alpha("grey",.5), color=alpha("grey",0))+
                geom_ribbon(data=for.plot %>% dplyr::filter(result=="NoGo"), aes(x=n.total, ymin=value, ymax=1), fill=alpha("red",.5), color=alpha("red", .25))+
                guides(color="none", fill="none")+
                labs(x="Total sample size", y="Decision Probabilities",
                     title=paste0("Operating Characteristics as a function of sample size"))+
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

        if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2))
        if(add.table==FALSE) return(main.plot)
}
