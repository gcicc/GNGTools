#' @title Make single sample normal-gamma treatment OC curve V2
#' @param my.df output from get.ss.ng.trt.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provides extended output summaries

#' @return A ggplot object is returned.
#' @export
#' @examples \donttest{
#' make.ss.ng.trt.oc2()
#' }
#' @author Greg Cicconetti
make.ss.ng.trt.oc2 <- function(my.df = get.ss.ng.trt.oc.df(), nlines=25, tsize=4, add.table=TRUE){
        # Build main plot -----
        main.plot <-  my.df %>% mutate(Consider = 1 - Go - NoGo) %>%
                gather(value= value, key=key, factor_key = T, c(Go, NoGo, Consider)) %>%
                mutate(key=factor(key, c("Go", "Consider", "NoGo"))) %>%
                ggplot(aes(x=xbar, y=value, color=key))+
                geom_line(size=.75)
        dpb <- ggplot_build(main.plot)

        main.plot <- main.plot+
                labs(x = "Underlying treatment effect",
                     y="Probability",
                     color="Decision",
                     title="Operating characteristics as a function of treatment effect")+
                theme(legend.position = "bottom")+
                scale_color_manual(values=c( "green", "grey20","red"))+
                scale_x_continuous(expand = c(0,0), breaks=pretty(my.df$xbar, 10))+
                scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2), limits=c(-0.0005,1), labels=scales::percent)+
                theme(panel.spacing.x = unit(6, "mm"), axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
                geom_vline(xintercept = c(my.df$Delta.tv[1], my.df$Delta.lrv[1]), color="blue", linetype=2)+
                annotate("text", label = paste0("Min TPP = ", my.df$Delta.lrv[1]), x = my.df$Delta.lrv[1],
                         y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize, colour = "black", hjust = 1)+
                annotate("text", label = paste("Base TPP = ", my.df$Delta.tv[1]), x = my.df$Delta.tv[1],
                         y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize,
                         colour = "black", hjust = 0)+
                labs(x = "Underlying treatment effect",
                     y="Probability")
        # Build table plot ----
        table.plot2 <- ggplot()+
                annotate("text", label = paste0("Decision Criteria"),
                         x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
                annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
                         x = -1, y = 1-3/nlines, size = tsize+1, hjust = 0)+
                annotate("text", label = TeX(paste0("P($\\Delta\\,$ >= Min TPP) > ", my.df$tau.lrv[1]*100, "% &")),
                         x = 0, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
                annotate("text", label =  TeX(paste0("P($\\Delta\\,$ >= Base TPP) > ", my.df$tau.tv[1]*100,"%")),
                         x = 1.5, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
                annotate("text", label = TeX(paste0("No-Go if:")), color="red",
                         x = -1, y = 1-4/nlines, size = tsize+1, hjust = 0)+
                annotate("text", label = TeX(paste0("P($\\Delta\\,$ >= Min TPP) <= ", my.df$tau.ng[1]*100, "% &")),
                         x = 0, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0) +
                annotate("text", label = TeX(paste0("P($\\Delta\\,$ >= Base TPP) <= ", my.df$tau.tv[1]*100,"%")),
                         x = 1.5, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0)+
                annotate("text", label = TeX(paste0("Consider:")), color="black",
                         x = -1, y = 1-5/nlines, size = tsize+1, hjust = 0)+
                annotate("text", label = "Otherwise",
                         x = 0, y = 1-5/nlines, size = tsize+1, colour = "black", hjust = 0)+
                scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
                scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
                labs(x="\n", y="")

        # return ----
        if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2))
        if(add.table==FALSE) return(main.plot)
}
