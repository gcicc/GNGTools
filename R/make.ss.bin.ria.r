#' @title Make single sample binary rule in action plot
#'
#' @param a.trt prior alpha parameter
#' @param b.trt prior beta parameter
#' @param n.trt observed sample size
#' @param x.trt observed number of responders
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param seed random seed
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param nlines.ria number of lines
#' @param add.table provides extended output summaries
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' my.ss.bin.ria <- make.ss.bin.ria(x.trt=10, add.table=TRUE)
#' plot(my.ss.bin.ria[[1]])
#' my.ss.bin.ria[[2]]
#' my.ss.bin.ria[[3]]
#' my.ss.bin.ria[[4]]
make.ss.bin.ria <- function(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
                            Delta.lrv = .2, Delta.tv = .35,
                            tau.tv = 0.10, tau.lrv = .80, tau.ng = .65,
                            seed = 1234, nlines = 25, tsize = 4, nlines.ria=20, add.table=TRUE){

  results <- get.ss.bin.df(a.trt = a.trt, b.trt = b.trt, n.trt = n.trt, x.trt = 0:n.trt,
                           Delta.lrv = Delta.lrv, Delta.tv=Delta.tv, tau.lrv=tau.lrv,
                           tau.tv=tau.tv, tau.ng=tau.ng)

  result.go <- results %>% dplyr::filter(result== "Go") %>% slice(1)

  # Determine max number of TRT responders for No-Go
  result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(n())
  my.df <- gcurve(expr = dbeta(x,shape1 = a.trt + x.trt, shape2 = b.trt + n.trt - x.trt),
                  from = 0, to =1, category="Posterior") %>%
    mutate(alpha = a.trt + x.trt, beta = b.trt + n.trt - x.trt, group="Posterior")
  my.df$facet <- paste0("Prior: Beta(", a.trt, ", ", b.trt, ")")

  P.R1 = 1 - pbeta(Delta.lrv, a.trt + x.trt, b.trt + (n.trt - x.trt))
  P.R3 = 1 - pbeta(Delta.tv, a.trt + x.trt, b.trt + (n.trt - x.trt))
  result = ifelse(P.R1 >= tau.lrv & P.R3 >= tau.tv, "Go",
                  ifelse(P.R1 < tau.ng  & P.R3 < tau.tv, "No-Go", "Consider"))

  result.color <- ifelse(result=="Go", "darkgreen",
                         ifelse(result=="No-Go", "red", "black"))

  for.subtitle <- paste0("Treatment data: ", round(x.trt/n.trt,2)*100,
                         "% (",  x.trt, "/", n.trt, "). ",
                         "Given control data, responders needed for Go: ",
                         result.go$x.trt[1],
                         " (", round(result.go$x.trt[1]/n.trt,2)*100,
                         "%). Needed for No-Go: ",  result.ng$x.trt[1], " (",
                         round(result.ng$x.trt[1]/n.trt, 2)*100,"%)")

  # Annotation line 1: Decision ----
  for.decision <- paste0("Decision: ", result)

  # Annotation line 2: Decision interval ----
  if(P.R3 > tau.tv){
    for.decision.interval <- paste0("Decision interval: (",
                                    round(qbeta(p = 1 - tau.lrv, shape1 = a.trt + x.trt,
                                                shape2 = b.trt + n.trt - x.trt),3)*100,
                                    "%, ",
                                    round(qbeta(p = 1 - tau.tv, shape1 = a.trt + x.trt,
                                                shape2 = b.trt + n.trt - x.trt),3)*100,
                                    "%)")
  } else {
    for.decision.interval <- paste0("Decision interval: (",
                                    round(qbeta(p = 1 - tau.ng, shape1 = a.trt + x.trt,
                                                shape2 = b.trt + n.trt - x.trt),3)*100, "%, ",
                                    round(qbeta(p = 1 - tau.tv, shape1 = a.trt + x.trt,
                                                shape2 = b.trt + n.trt - x.trt),3)*100, "%)")
  }

  # Annotation line 3: P(Delta >= Min TPP)

  annotate.P1 <- ifelse(
    result=="Go",
    TeX(paste0("$P(\\Delta$ >= Min TPP) = ",
               round(1 - pbeta(q = Delta.lrv, shape1 = a.trt + x.trt,
                               shape2 = b.trt + n.trt - x.trt),3)*100,
               "% > ", tau.lrv*100, "%")),
    ifelse(result=="No-Go",
           TeX(paste0("$P(\\Delta$ >=  Min TPP) = ",
                      round(1 - pbeta(q = Delta.lrv, shape1 = a.trt + x.trt,
                                      shape2 = b.trt + n.trt - x.trt),3)*100,
                      "% <= ", tau.ng*100, "%")),
           TeX(paste0("$P(\\Delta$ >= Min TPP) = ",
                      round(1 - pbeta(q = Delta.lrv, shape1 = a.trt + x.trt,
                                      shape2 = b.trt + n.trt - x.trt),3)*100, "%"))))

  # Annotation line 4: P(Delta >= Base TPP)
  annotate.P2 <- ifelse(
    result=="Go",
    TeX(paste0("$P(\\Delta$ >= Base TPP) = $",
               round(1 - pbeta(q = Delta.tv, shape1 = a.trt + x.trt,
                               shape2 = b.trt + n.trt - x.trt),3)*100, "$% > ",
               tau.tv*100, "%")),
    ifelse(result=="No-Go", TeX(paste0("$P(\\Delta\\$ >= Base TPP) = ",
                                       round(1 - pbeta(q = Delta.tv, shape1 = a.trt +
                                                         x.trt, shape2 = b.trt + n.trt
                                                       - x.trt),3)*100, "% <=",
                                       tau.tv*100, "%")),
           TeX(paste0("$P(\\Delta$ >= Base TPP) = ",
                      round(1 - pbeta(q = Delta.tv, shape1 = a.trt + x.trt,
                                      shape2 = b.trt + n.trt - x.trt),3)*100, "%"))))

  # Initialize a ggplot
  dplot <- ggplot() + geom_line(data = my.df, aes(x = x,y = y)) + facet_wrap(~facet)
  # Access the ggplot to get goodies to help accomplish shading
  dpb <- ggplot_build(dplot)
  x1.1 <- min(which(dpb$data[[1]]$x >=Delta.lrv))
  x2.1 <- max(which(dpb$data[[1]]$x <=Delta.tv))+1
  x1.2 <- min(which(dpb$data[[1]]$x >=Delta.tv))
  x2.2 <- max(which(dpb$data[[1]]$x <=1))

  if(which(dpb$data[[1]]$y== max(dpb$data[[1]]$y)) > length(dpb$data[[1]]$x)/2){
    annotate.x = min(min(dpb$data[[1]]$x), Delta.lrv, 0)
    annotate.j = 0
  } else {annotate.x = max(max(dpb$data[[1]]$x), Delta.tv, 1)
  annotate.j = 1}

  go.segment <- data.frame(x = qbeta(p = 1 - tau.lrv, shape1 = a.trt + x.trt,
                                     shape2 = b.trt + n.trt - x.trt),
                           xend = qbeta(p = 1 - tau.tv, shape1 = a.trt + x.trt,
                                        shape2 = b.trt + n.trt - x.trt),
                           y= min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines * 3,
                           yend = min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines * 3,
                           group=my.df$group[1])
  nogo.segment <- data.frame(x = qbeta(p = 1 - tau.ng, shape1 = a.trt + x.trt,
                                       shape2 = b.trt + n.trt - x.trt),
                             xend = qbeta(p = 1 - tau.tv, shape1 = a.trt + x.trt,
                                          shape2 = b.trt + n.trt - x.trt),
                             y= min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines * 3,
                             yend = min(dpb$data[[1]]$y) +
                               max(dpb$data[[1]]$y)/nlines * 3,
                             group=my.df$group[1])

  # Introduce shading
  main.plot <- dplot +
    geom_area(data = data.frame(x = dpb$data[[1]]$x[1:x1.1],
                                y = dpb$data[[1]]$y[1:x1.1]),
              aes(x = x, y = y), fill=alpha("red", 0))+
    geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.1:x2.1],
                                y = dpb$data[[1]]$y[x1.1:x2.1]),
              aes(x = x, y = y), fill=alpha("grey80", 0))+
    geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.2:x2.2],
                                y = dpb$data[[1]]$y[x1.2:x2.2]),
              aes(x = x, y = y), fill=alpha("green", 0))+
    geom_line(data = my.df, aes(x = x,y = y))+
    scale_x_continuous(limits = c(0,1),
                       breaks = unique(c(
                         pretty(x = c(0, Delta.lrv-.06),
                                n = diff(c(0, Delta.lrv))/2 * 10)[
                                  pretty(x = c(0, Delta.lrv-.06),
                                         n = diff(c(0, Delta.lrv))/2 * 10) <
                                    Delta.lrv-.06],
                         c(Delta.lrv, Delta.tv),
                         pretty(x = c(Delta.tv, 1), n = diff(c(Delta.tv,1))/2 * 10)[
                           pretty(x = c(Delta.tv, 1), n = diff(c(Delta.tv,1))/2 * 10) >
                             Delta.tv+.06]
                       )),
                       labels = scales::percent)+
    scale_y_continuous(breaks=NULL, labels = NULL)

  # Add Annotations -----
  main.plot <- main.plot+
    labs(title = TeX("Posterior distribution for the proportion of treatment responders"),
         subtitle = for.subtitle,
         x = TeX("$\\Delta\\,$ = Proportion of treatment responders"), y = NULL,
         caption = NULL)+
    annotate("text", label = for.decision,
             x = annotate.x, y = max(dpb$data[[1]]$y)-max(dpb$data[[1]]$y)/nlines.ria * 0,
             size = tsize+1, colour = result.color, hjust = annotate.j)+
   # annotate("text", label = for.decision.interval,
   #          x = annotate.x, y = max(dpb$data[[1]]$y)-max(dpb$data[[1]]$y)/nlines.ria * 1,
   #          size = tsize, colour = result.color, hjust = annotate.j)+
    annotate("text", label = annotate.P1, color=result.color,
             x = annotate.x, y = max(dpb$data[[1]]$y)-max(dpb$data[[1]]$y)/nlines.ria * 2,
             size = tsize,  hjust = annotate.j)+
    annotate("text", label = annotate.P2, color=result.color,
             x = annotate.x, y = max(dpb$data[[1]]$y)-max(dpb$data[[1]]$y)/nlines.ria * 3,
             size = tsize,  hjust = annotate.j)

  # Add reference lines and Credible interval

  # if(result.color == "red"){
  #   main.plot <- main.plot +
  #     geom_segment(data=nogo.segment, aes(x=x, xend=xend, y=y, yend=yend, group=group),
  #                  arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)
  # } else {
  #   main.plot <- main.plot +
  #     geom_segment(data=go.segment, aes(x=x, xend=xend, y=y, yend=yend, group=group),
  #                  arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)}
  main.plot <- main.plot +
    geom_vline(xintercept = c(Delta.lrv, Delta.tv), linetype = 2,
               color = c("blue", "blue")) +
    annotate("label", label = TeX(paste0("Min TPP = $", Delta.lrv*100,"$%")),
             x = Delta.lrv, y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize,
             colour = "black", hjust = 1, label.size=NA, fill="grey92", alpha=.8)+
    annotate("text", label = TeX(paste0("Base TPP = $",  Delta.tv*100,"$%")),
             x = Delta.tv, y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize,
             colour = "black", hjust = 0)

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

  # plot_grid(main.plot,
  #            table.plot2, ncol=1, align = "v", rel_heights=c(.78,.22), axis="l")

  if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2,
                          data.frame(a.trt = a.trt, b.trt = b.trt,
                                     Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                     tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                                      Needed.for.NG=result.ng$x.trt[1], Needed.for.GO = result.go$x.trt[1])))
  if(add.table==FALSE) return(main.plot)
}


