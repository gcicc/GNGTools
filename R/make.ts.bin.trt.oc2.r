#' @title Make two-sample binary treatment OC curve V2
#'
#' @param for.plot output from get.ts.bin.trt.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provides extended output summaries

#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.bin.trt.oc.df <- get.ts.bin.trt.oc.df()
#' my.ts.bin.trt.oc2 <- make.ts.bin.trt.oc2(for.plot=my.ts.bin.trt.oc.df, add.table=TRUE)
#' plot(my.ts.bin.trt.oc2[[1]])
#' my.ts.bin.trt.oc2[[2]]
#' my.ts.bin.trt.oc2[[3]]
#' }
make.ts.bin.trt.oc2 <- function(for.plot=get.ts.bin.trt.oc.df(), tsize=4, nlines=25, add.table=TRUE){
        # Pull items for annotation from data.frame
        Delta.tv <- for.plot$Delta.tv[1]
        Delta.lrv <- for.plot$Delta.lrv[1]
        tau.tv <- for.plot$tau.tv[1]
        tau.lrv <- for.plot$tau.lrv[1]
        tau.ng <- for.plot$tau.ng[1]
        dcurve.con <- for.plot$dcurve.con[1]
        Aratio <- for.plot$Aratio[1]
        total.n <- for.plot$n.con[1] + for.plot$n.trt[1]

        my.plot <- ggplot(data = for.plot, aes(x = x, y = prob, color = result)) +
                geom_line(size=.75)
        dpb <- ggplot_build(my.plot)
        main.plot <-  my.plot +
                labs(x = "Underlying treatment effect",
                     y="Probability", color="Decision",
                     title = paste0("Operating characteristics as a function of treatment effect when the control rate is: ", round(for.plot$dcurve.con[1]*100, 2), "%"),
                     subtitle = paste0("Total randomized: ", total.n, " with ", for.plot$n.con[1],
                                       " to Control and ", for.plot$n.trt[1], " to Treatment"))+
                scale_color_manual(values = c("Go" = "green", "No-Go" = "red",
                                              "Consider" = "darkgrey"))+
                geom_vline(xintercept = c(Delta.tv, Delta.lrv), linetype=2, color="blue")+
                annotate("text", label = TeX(paste0("Min TPP = ", Delta.lrv*100, "%")),
                         x = Delta.lrv, y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize,
                         colour = "black", hjust = 1)+
                annotate("text", label = TeX(paste("Base TPP = ", Delta.tv*100, "%")),
                         x = Delta.tv, y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize,
                         colour = "black", hjust = 0)+
                theme(legend.position = "bottom")+
                scale_x_continuous(expand = c(0,0),
                                   breaks=seq(for.plot$TE.OC.Delta.LB[1],
                                              for.plot$TE.OC.Delta.UB[1], length.out=5),
                                   limits=c(head(dpb$data[[1]]$x,1),tail(dpb$data[[1]]$x,1)),
                                   labels=scales::percent) +
                scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                                   labels=scales::percent)+
                theme(panel.spacing.x = unit(6, "mm"),
                      axis.text.x = element_text(angle=45, hjust=1,vjust=1))
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
