#' @title Make time to event treatment OC curve v2
#'
#' @param plot.df output from get.tte.trt.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provides extended output summaries
#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' my.tte.trt.oc.df <- get.tte.trt.oc.df()
#' make.tte.trt.oc2(plot.df = my.tte.trt.oc.df, add.table=TRUE)
make.tte.trt.oc2 <- function(plot.df = get.tte.trt.oc.df(), nlines=25, tsize=4, add.table=TRUE){
  HR.tv = plot.df$HR.tv[1]
  HR.lrv = plot.df$HR.lrv[1]
  tau.tv = plot.df$tau.tv[1]
  tau.lrv = plot.df$tau.lrv[1]
  tau.ng = plot.df$tau.ng[1]
  m.obs = plot.df$m.obs[1]
  ARatio = plot.df$ARatio[1]
  HR.upper = plot.df$HR.upper[1]
  HR.lower = plot.df$HR.lower[1]

  if(HR.tv < HR.lrv) {
    main.plot <- plot.df %>% mutate(Consider = 1 - Go - NoGo) %>%
      gather(key = key, value=value, c(Go, NoGo, Consider), factor_key = TRUE) %>%
      dplyr::filter(key !="NoGo2") %>%
      ggplot(aes(x= HR.obs, y=value, color=key))+
      geom_line(size=.75)+
      scale_color_manual(values = c("Go" = "green", "NoGo" = "red",
                                    "Consider" = "darkgrey"))
    dpb <- ggplot_build(main.plot)
    main.plot <- main.plot +
      geom_vline(xintercept = c(HR.tv, HR.lrv), color="blue", linetype=2)+
      labs(x="Underlying Hazard Ratio", y="Probability",
           color="Decision",
           title="Operating characteristics as a function of treatment effect",
           subtitle=paste0("Total events: ", m.obs,
                           ". Randomization ratio (Control:Treatment): (1:",
                           ARatio, ")"))+
      theme(legend.position = "bottom")+
      annotate("text", label = TeX(paste0("$HR_{Max}$ = ",HR.lrv)), x = HR.lrv,
               y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size =
                 tsize, colour = "black", hjust = 0)+
      annotate("text", label = TeX(paste0("$\\HR_{Base}$ = ", HR.tv)), x = HR.tv,
               y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size = tsize,
               colour = "black", hjust = 1)+
      scale_x_continuous(expand = c(0,0), breaks=seq(0,2,.2),
                         limits=c(HR.lower, HR.upper))+
      scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2),
                         minor_breaks=seq(0,1,.1), limits=c(0,1), labels=scales::percent)+
      theme(panel.spacing.x = unit(6, "mm"), axis.text.x = element_text(angle=45, hjust=1,vjust=1))

    table.plot2 <- ggplot()+
      annotate("text", label = paste0("Decision Criteria"),
               x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
      annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
               x = -1, y = 1-2.5/nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = TeX(paste0("P($\\Delta$ > $HR_{Max}$) > ", tau.lrv*100, "% &")),
               x = 0, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
      annotate("text", label =  TeX(paste0("P($\\Delta$ > $HR_{Base}$) > ", tau.tv*100,"%")),
               x = 1.5, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
      annotate("text", label = TeX(paste0("No-Go if:")), color="red",
               x = -1, y = 1-3.25/nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = TeX(paste0("P($\\Delta$ > $HR_{Max}$) <= ", tau.ng*100, "% &")),
               x = 0, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0) +
      annotate("text", label =  TeX(paste0("P($\\Delta$ > $HR_{Base}$) <= ", tau.tv*100,"%")),
               x = 1.5, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0)+
      annotate("text", label = TeX(paste0("Consider:")), color="black",
               x = -1, y = 1-4./nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = "Otherwise",
               x = 0, y = 1-4./nlines, size = tsize+1, colour = "black", hjust = 0)+
      scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
      scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
      labs(x="", y="")

  } else {
    main.plot <- plot.df %>% mutate(Consider = 1 - Go - NoGo) %>%
      gather(key = key, value=value, c(Go, NoGo, Consider), factor_key = T) %>%
      dplyr::filter(key !="NoGo2") %>%
      ggplot(aes(x= HR.obs, y=value, color=key))+
      geom_line(size=.75)+
      scale_color_manual(values = c("Go" = "green", "NoGo" = "red",
                                    "Consider" = "darkgrey"))
    dpb <- ggplot_build(main.plot)
    main.plot <- main.plot +
      geom_vline(xintercept = c(HR.tv, HR.lrv), color="blue", linetype=2)+
      labs(x="Underlying Hazard Ratio", y="Probability",
           color="Decision",
           title="Operating characteristics as a function of treatment effect",
           subtitle=paste0("Total events: ", m.obs,
                           ". Randomization ratio (Control:Treatment): (1:",
                           ARatio, ")"))+
      theme(legend.position = "bottom")+
      annotate("text", label = TeX(paste0("$HR_{Max}$ = ",HR.lrv)), x = HR.lrv,
               y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size =
                 tsize, colour = "black", hjust = 1)+
      annotate("text", label = TeX(paste0("$\\HR_{Base}$ = ", HR.tv)), x = HR.tv,
               y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size = tsize,
               colour = "black", hjust = 0)+
      scale_x_continuous(expand = c(0,0), breaks=seq(0,2,.2),
                         limits=c(HR.lower, HR.upper))+
      scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2),
                         minor_breaks=seq(0,1,.1), limits=c(0,1), labels=scales::percent)+
      theme(panel.spacing.x = unit(6, "mm"), axis.text.x = element_text(angle=45, hjust=1,vjust=1))

    table.plot2 <- ggplot()+
      annotate("text", label = paste0("Decision Criteria"),
               x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
      annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
               x = -1, y = 1-2.5/nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = TeX(paste0("P($\\Delta$ > $HR_{Max}$) > ", tau.lrv*100, "%$\\;\\;\\;\\;\\,$ &")),
               x = 0, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
      annotate("text", label =  TeX(paste0("P($\\Delta$ > $HR_{Base}$) > $", tau.tv*100,"$%")),
               x = 1.5, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
      annotate("text", label = TeX(paste0("No-Go if:")), color="red",
               x = -1, y = 1-3.25/nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = TeX(paste0("P($\\Delta$ > $HR_{Max}$) <= ", tau.ng*100, "%$\\;\\;\\;\\,\\,$ &")),
               x = 0, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0) +
      annotate("text", label =  TeX(paste0("P($\\Delta$ > $HR_{Base}$) <= ", tau.tv*100,"$%")),
               x = 1.5, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0)+
      annotate("text", label = TeX(paste0("Consider:")), color="black",
               x = -1, y = 1-4./nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = "Otherwise",
               x = 0, y = 1-4./nlines, size = tsize+1, colour = "black", hjust = 0)+
      scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
      scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
      labs(x="", y="")

  }


  if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2))
  if(add.table==FALSE) return(main.plot)}
