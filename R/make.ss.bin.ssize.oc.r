#' @title Make single sample binary sample size OC curve
#'
#' @param for.plot ouput from get.ss.bin.ssize.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provides extended output summaries

#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' my.ss.bin.ssize.oc.df <- get.ss.bin.ssize.oc.df()
#' my.ss.bin.ssize.oc <- make.ss.bin.ssize.oc(for.plot= my.ss.bin.ssize.oc.df, add.table=TRUE)
#' plot(my.ss.bin.ssize.oc)
make.ss.bin.ssize.oc <- function(for.plot=get.ss.bin.ssize.oc.df(), nlines=25, tsize=4, add.table=TRUE){
  Delta.lrv <- for.plot$Delta.lrv[1]
  Delta.tv <- for.plot$Delta.tv[1]
  Delta.user <- for.plot$Delta.user[1]
  tau.tv <- for.plot$tau.tv[1]
  tau.lrv <- for.plot$tau.lrv[1]
  tau.ng <- for.plot$tau.ng[1]
  # Simply takes a dataframe from get.ss.bin.ssize.oc.df and plots
  main.plot <-  ggplot() +
    geom_line(data=for.plot %>% dplyr::filter(result=="Go"),
              aes(x=n.trt, y=value, group=result), color="green") +
    geom_line(data=for.plot %>% dplyr::filter(result=="No-Go"),
              aes(x=n.trt, y=value, group=result), color="red")+
    facet_wrap(~key, labeller="label_parsed")+
    geom_ribbon(data=for.plot %>% dplyr::filter(result=="Go"),
                aes(x=n.trt, ymin=0, ymax=value), fill="green", alpha=.5)+
    geom_ribbon(data=for.plot %>% dplyr::filter(result=="Go"),
                aes(x=n.trt, ymin=value, ymax=1), fill="grey", alpha=.5)+
    geom_ribbon(data=for.plot %>% dplyr::filter(result=="No-Go"),
                aes(x=n.trt, ymin=value, ymax=1), fill="red", alpha=.5)+
    labs(x="Sample size per arm", y="Decision Probabilities",
         title="Operating characteristics as a function of sample size")+
    scale_x_continuous(expand=c(0,0))+
    theme(panel.spacing = unit(2, "lines")) +
    scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2), limits=c(0,1.01),
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

  if(add.table==TRUE) return(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)))
  if(add.table==FALSE) return(main.plot)
  }


