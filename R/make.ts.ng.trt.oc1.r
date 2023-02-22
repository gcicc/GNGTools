#' @title Make two-sample normal-gamma treatment OC curve
#'
#' @param for.plot output from get.ts.ng.trt.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provides extended output summaries

#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.ng.trt.oc.df <- get.ts.ng.trt.oc.df(goparallel=FALSE)
#' my.ts.ng.trt.oc1 <- make.ts.ng.trt.oc1(for.plot = my.ts.ng.trt.oc.df, add.table=TRUE)
#' plot(my.ts.ng.trt.oc1[[1]])
#' my.ts.ng.trt.oc1[[2]]
#' my.ts.ng.trt.oc1[[3]]
#' }
make.ts.ng.trt.oc1 <- function(for.plot=get.ts.ng.trt.oc.df(goparallel=FALSE), nlines=25, tsize=4, add.table=TRUE){

        # for.plot <- for.return

        for.plot1 <- for.plot %>%
                dplyr::filter(result!="Consider") %>%
                mutate(freq=ifelse(result=="Go", freq, 1 - freq))


        my.plot <- ggplot() +
                geom_line(data=for.plot1 %>% dplyr::filter(result=="Go"), aes(x=treatment.effect, y=freq), color="darkgreen")
        dpb <- ggplot_build(my.plot)

        # Simply takes a dataframe from return.ssize.df and plots
        main.plot <- my.plot +
                geom_ribbon(data=for.plot1 %>% dplyr::filter(result=="Go"), aes(x=treatment.effect, ymin=0, ymax=freq), fill="lightgreen", alpha=.5)+
                geom_ribbon(data=for.plot1 %>% dplyr::filter(result=="Go"), aes(x=treatment.effect, ymin=freq, ymax=1), fill="grey", alpha=.5)+
                geom_line(data=for.plot1 %>% dplyr::filter(result=="No-Go"), aes(x=treatment.effect, y=freq), color="red")+
                geom_ribbon(data=for.plot1 %>% dplyr::filter(result=="No-Go"), aes(x=treatment.effect, ymin=freq, ymax=1), fill="red", alpha=.5, linetype=2)+
                labs(x="Underlying Treatment effect",
                     y="Probability",
                     title=paste0("Operating characteristics as a function of treatment effect when the control mean [sd] is ", round(for.plot$xbar.c[1],2), " [",round(for.plot$s.c[1],2),"]"),
                     subtitle = paste0("Total randomized with ",
                                       for.plot$N[1],
                                       " with ",
                                       for.plot$N[1] - floor(for.plot$N[1]*(for.plot$ARatio[1]/(1+for.plot$ARatio[1]))),
                                       " to Control and ",
                                       floor(for.plot$N[1]*(for.plot$ARatio[1]/(1+for.plot$ARatio[1]))),
                                       " to Treatment")) +
                scale_x_continuous(expand=c(0,0))+
                scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2), limits=c(0,1), labels=scales::percent)+
                theme(panel.spacing.x = unit(6, "mm"), axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
                geom_vline(xintercept = c(for.plot$Delta.tv[1], for.plot$Delta.lrv[1]), linetype=2, color="blue")+
                annotate("text", label = paste0("Min TPP = ", for.plot$Delta.lrv[1]),
                         x = for.plot$Delta.lrv[1], y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize, colour = "black", hjust = 1, vjust=0)+
                annotate("text", label = paste("Base TPP = ", for.plot$Delta.tv[1]),
                         x = for.plot$Delta.tv[1], y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize, colour = "black", hjust = 0, vjust=0)

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

