#' @title Make single sample binary treatment oc curve
#'
#' @param my.df output from get.ss.bin.trt.oc.df
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param nlines.ria Control for text spacing
#' @param Delta_OC_LB Lower bound for OC curve
#' @param Delta_OC_UB Upper bound for OC curve
#' @param add.table provides extended output summaries
#' @return A ggplot object is returned
#' @export
#'
#' @examples
#'  my.ss.bin.trt.oc.df <- get.ss.bin.trt.oc.df()
#'  make.ss.bin.trt.oc1(my.df = my.ss.bin.trt.oc.df)

make.ss.bin.trt.oc1 <- function(my.df = get.ss.bin.trt.oc.df(),
                                tsize=4, nlines=25, nlines.ria=20,
                                Delta_OC_LB = 0, Delta_OC_UB = .41, add.table=TRUE){
  Delta.tv <- my.df$Delta.tv[1]
  Delta.lrv <- my.df$Delta.lrv[1]
  tau.tv <- my.df$tau.tv[1]
  tau.lrv <- my.df$tau.lrv[1]
  tau.ng <- my.df$tau.ng[1]
  my.df2 <- subset(my.df, result=="Go") %>%
    dplyr::filter(Delta <= Delta_OC_UB & Delta >= Delta_OC_LB)
  my.df3 <- subset(my.df, result=="No-Go") %>%
    dplyr::filter(Delta <= Delta_OC_UB & Delta >= Delta_OC_LB)
  # In the case where No-Go is never attained...
  if(nrow(my.df3) ==0){
    my.df3 <- my.df2
    my.df3$result <- "No-Go"
    my.df3$p <- 0
  }
  my.df3$p <- 1 - my.df3$p

  total.n <- my.df$n.trt[1]

  my.plot <- ggplot()+geom_line(data=my.df2, aes(x=Delta, y=p))+
    geom_line(data=my.df3, aes(x=Delta, y=p))
  dpb <- ggplot_build(my.plot)
  main.plot <- my.plot +
    geom_area(data = dpb$data[[1]], aes(x = x, y = y), fill=alpha("lightgreen",0.5))+
    geom_ribbon(data = dpb$data[[1]], aes(x = x, ymin = y, ymax=1),
                fill=alpha("grey", .5))+
    geom_ribbon(data = dpb$data[[2]], aes(x = x, ymin = y, ymax = 1),
                fill=alpha("red", 0.5))+
    geom_line(data=my.df2, aes(x=Delta, y=p), size=.75, color="green")+
    geom_line(data=my.df3, aes(x=Delta, y=p), size=.75, color="red")+
    geom_vline(xintercept = c(Delta.tv, Delta.lrv), linetype=2, color="blue")+
    annotate("text", label = TeX(paste0("Min TPP = ", Delta.lrv*100,"%")),
             x = Delta.lrv, y = 0 + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria,
             size = tsize, colour = "black", hjust = 1)+
    annotate("text", label = TeX(paste0("Base TPP = ",  Delta.tv*100,"%")),
             x = Delta.tv, y = 0 + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria,
             size = tsize, colour = "black", hjust = 0)+
    labs(x = TeX("$\\Delta\\,$ = Underlying proportion of treatment responders"),
         y = "Probability",
         title = "Operating characteristics as a function of treatment effect",
         subtitle = paste0("Total randomized to treatment: ", total.n))+
    scale_x_continuous(expand = c(0, 0), limits=c(Delta_OC_LB,
                                                  min(Delta_OC_UB, max(dpb$data[[1]]$x))),
                       labels=scales::percent)+
    scale_y_continuous(expand = c(0,0), breaks=seq(0,1,.2), labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"),
          axis.text.x = element_text(angle=45, hjust=1,vjust=1))


  table.plot2 <-ggplot()+
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

  # plot_grid(main.plot,
  #            table.plot2, ncol=1, align = "v", rel_heights=c(.78,.22), axis="l")

  if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2))
  if(add.table==FALSE) return(main.plot)
  }



