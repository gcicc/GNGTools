#' Make two-sample binary sample size operating characteristics curve
#'
#' @param for.plot output from get.ts.bin.ssize.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param add.table provides extended output summaries
#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.bin.ssize.oc.df <- get.ts.bin.ssize.oc.df()
#' make.ts.bin.ssize.oc(for.plot = my.ts.bin.ssize.oc.df)
#' }
make.ts.bin.ssize.oc <- function(for.plot = get.ts.bin.ssize.oc.df(), tsize=4,
                                 nlines=25, add.table=TRUE){

  for.plot <- for.plot %>% mutate(value = ifelse(result=="No-Go", 1 - value, value))
  main.plot <-  for.plot %>%
    dplyr::filter(result!="Consider") %>%
    ggplot(aes(x=n.total, y=value, color=result)) +
    geom_line(alpha=.25, size=.75) +
    facet_wrap(~key, labeller="label_parsed")+
    scale_color_manual(values=c("green", "red"))+
    geom_ribbon(data=for.plot %>% dplyr::filter(result=="Go"),
                aes(x=n.total, ymin=0, ymax=value),
                fill=alpha("green",.5), color=alpha("green",.25))+
    geom_ribbon(data=for.plot %>% dplyr::filter(result=="Go"),
                aes(x=n.total, ymin=value, ymax=1),
                fill=alpha("grey",.5), color=alpha("grey",0))+
    geom_ribbon(data=for.plot %>% dplyr::filter(result=="No-Go"),
                aes(x=n.total, ymin=value, ymax=1), fill=alpha("red",.5),
                color=alpha("red", .25))+
    guides(color="none", fill="none")+
    labs(x="Total sample size", y="Decision Probabilities",
         title=paste0("Operating Characteristics as a function of sample size when control response rate is ",
                      round(for.plot$dcurve.con[1]*100,2), "%"))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                       labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
            element_text(angle=45, hjust=1,vjust=1))

  head(for.plot)
  tau.lrv <- for.plot$tau.lrv[1]
  tau.tv <- for.plot$tau.tv[1]
  tau.ng <- for.plot$tau.lrv[1]

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

  if(add.table==TRUE) return(  grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)))
  if(add.table==FALSE) return(main.plot)
}


