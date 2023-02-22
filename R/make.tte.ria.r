#' @title Make time to event rule in action plot
#'
#' @param m.con.prior prior number of control events
#' @param m.trt.prior prior number of treatment events
#' @param HR.prior prior estimate for HR
#' @param ARatio Allocation ratio
#' @param HR.obs observed hazard ratio
#' @param m.obs observed number of events
#' @param HR.lrv TPP Lower Reference Value aka Max TPP (large HRs lead to No-Go)
#' @param HR.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param nlines.ria Control for text spacing
#' @param add.table provides extended output summaries
#' @return A ggplot object is returned.
#' @export
#'
#' @examples
#' my.tte.ria <- make.tte.ria(add.table=TRUE)
#' plot(my.tte.ria[[1]])
#' my.tte.ria[[2]]
#' my.tte.ria[[3]]
#' my.tte.ria[[4]]
make.tte.ria <- function(m.con.prior = 1,m.trt.prior = 1, HR.prior=1, ARatio=.5, HR.obs=.845, m.obs = 1500,
                         HR.tv= .8, HR.lrv = .9, tau.tv=.1, tau.lrv=.7, tau.ng=.7, tsize = 4, nlines = 25, nlines.ria=20, add.table=TRUE){
  # Run this for fixed values of

  results <- get.tte.studyend.GNG(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior, HR.prior=HR.prior,
                                  ARatio=ARatio, HR.obs=HR.obs, m.obs = m.obs,
                                  HR.tv= HR.tv, HR.lrv = HR.lrv, tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng)
  result.go <- results$result.go
  result.ng <- results$result.ng
  # Determine min number of TRT responders for Go


  post.params <- get.tte.post.param(m.con.prior = m.con.prior,m.trt.prior = m.trt.prior, HR.prior = HR.prior, ARatio = ARatio, HR.obs = HR.obs, m.obs = m.obs)

  my.df <- rbind(
    gcurve(expr = dnorm(x, mean = post.params[1,1],sd = post.params[1,2]), from = log(0.01), to = log(3),
           n = 1001, category=paste0("Posterior hazard ratio: ", round(exp(post.params[1]),2), " worth ", m.trt.prior+ m.con.prior + m.obs, " events")) %>%
      mutate(HR = post.params[1,1], events = m.con.prior + m.trt.prior + m.obs, group="Posterior"))

  if(HR.tv <= HR.lrv){
    P.R1 = pnorm(log(HR.tv), post.params[1,1], post.params[1,2])
    P.R3 = pnorm(log(HR.lrv),post.params[1,1], post.params[1,2])

  } else {
    P.R1 = 1 - pnorm(log(HR.tv), post.params[1,1], post.params[1,2])
    P.R3 = 1 - pnorm(log(HR.lrv),post.params[1,1], post.params[1,2])
  }

  # Annotation line 1: Decision ----

  result = ifelse(P.R1 > tau.tv  & P.R3 > tau.lrv, "Go",
                  ifelse(P.R1 < tau.tv  & P.R3 < tau.ng, "No-Go", "Consider"))
  for.decision <- paste0("Decision: ", result)
  result.color <- ifelse(result=="Go", "darkgreen", ifelse(result=="No-Go", "red", "black"))


  # Annotation line 2: Decision interval ----
  if(HR.tv <= HR.lrv){
    if(P.R1 > tau.tv){
      a <- exp(qnorm(p = tau.tv,mean = post.params[1,1], sd = post.params[1,2]))
      b <- exp(qnorm(p = tau.ng,mean = post.params[1,1], sd = post.params[1,2]))
      for.decision.interval <- paste0("Decision interval: (",
                                      round(a,3), ", ",
                                      round(b,3),")")
    } else {
      a <- exp(qnorm(p = tau.tv,mean = post.params[1,1], sd = post.params[1,2]))
      b <- exp(qnorm(p = tau.lrv,mean = post.params[1,1], sd = post.params[1,2]))
      for.decision.interval <- paste0("Decision interval: (",
                                      round(a,3), ", ",
                                      round(b,3),")")
    }
  } else {
    if(P.R1 > tau.tv){
      a <- exp(qnorm(p = 1-tau.ng,mean = post.params[1,1], sd = post.params[1,2]))
      b <- exp(qnorm(p = 1-tau.tv,mean = post.params[1,1], sd = post.params[1,2]))
      for.decision.interval <- paste0("Decision interval: (",
                                      round(a,3), ", ",
                                      round(b,3),")")
    } else {
      a <- exp(qnorm(p = 1- tau.lrv,mean = post.params[1,1], sd = post.params[1,2]))
      b <- exp(qnorm(p = 1 - tau.tv,mean = post.params[1,1], sd = post.params[1,2]))
      for.decision.interval <- paste0("Decision interval: (",
                                      round(a,3), ", ",
                                      round(b,3),")")
    }

  }
  for.subtitle <- paste0("Data ", "(Obs HR, Total events): (", round(HR.obs, 4), ", ", m.obs, "). Observed HR needed for Go: ",
                         round(result.go$HR.obs,4), ". Needed for No-Go: ", round(result.ng$HR.obs,4))

  # Annotation line 3: P(Delta < HR Max)

  if(HR.tv <= HR.lrv){
    annotate.P1 <- ifelse(result=="Go",
                          TeX(paste0("P($\\Delta$ <= $\\HR_{Max}$) = ", round(P.R3,2), " >", tau.lrv)),
                          ifelse(result=="No-Go",
                                 TeX(paste0("P($\\Delta$ <= $\\HR_{Max}$) = ", round(P.R3,2), " < ", tau.ng)),
                                 TeX(paste0("P($\\Delta$ <= $\\HR_{Max}$) = ", round(P.R3,2),""  ))))
    # Annotation line 4: P(Delta < HR Base)

    annotate.P2 <- ifelse(result=="Go",
                          TeX(paste0("P($\\Delta$ <= $\\HR_{Base}$) = ",  round(P.R1 ,2), " > ", tau.tv)),
                          ifelse(result=="No-Go", TeX(paste0("P($\\Delta$ <= $\\HR_{Base}$) = ",  round(P.R1,2), " < ", tau.tv)),
                                 TeX(paste0("P($\\Delta$ <= $\\HR_{Base}$) = ",  round(P.R1,2),""))))
  } else {
    annotate.P1 <- ifelse(result=="Go",
                          TeX(paste0("P($\\Delta$ >= $\\HR_{Max}$) = ", round(P.R3,2), " >", tau.lrv)),
                          ifelse(result=="No-Go",
                                 TeX(paste0("P($\\Delta$ >= $\\HR_{Max}$) = ", round(P.R3,2), " < ", tau.ng)),
                                 TeX(paste0("P($\\Delta$ >= $\\HR_{Max}$) = ", round(P.R3,2),""  ))))
    # Annotation line 4: P(Delta < HR Base)

    annotate.P2 <- ifelse(result=="Go",
                          TeX(paste0("P($\\Delta$ >= $\\HR_{Base}$) = ",  round(P.R1 ,2), " > ", tau.tv)),
                          ifelse(result=="No-Go", TeX(paste0("P($\\Delta$ >= $\\HR_{Base}$) = ",  round(P.R1,2), " < ", tau.tv)),
                                 TeX(paste0("P($\\Delta$ >= $\\HR_{Base}$) = ",  round(P.R1,2),""))))

  }
  # Initialize a ggplot

  dplot <- ggplot(data = my.df, aes(x = exp(x),y = y)) + geom_line() + facet_wrap(~category)

  # Access the ggplot to get goodies to help accomplish shading
  dpb <- ggplot_build(dplot)
  x1.1 <- max(which(dpb$data[[1]]$x <= HR.lrv))+1
  x2.1 <- max(which((dpb$data[[1]]$x) <= (HR.tv)))
  x1.2 <- max(which((dpb$data[[1]]$x) <= (HR.tv)))
  x2.2 <- min(which((dpb$data[[1]]$x) >= 0))

  if(which(dpb$data[[1]]$y== max(dpb$data[[1]]$y)) < length(dpb$data[[1]]$x)/2){
    annotate.x = min(min(dpb$data[[1]]$x), HR.tv, 0)
    annotate.j = 0
  } else {annotate.x = max(max(dpb$data[[1]]$x), HR.lrv, 3)
  annotate.j = 1}


  for.decision.interval.df <- data.frame(x = round(a,3),
                                         xend = round(b,3),
                                         y= min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3,
                                         yend = min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3,
                                         group=my.df$group[1])


  # Introduce shading
  main.plot <- dplot +
    geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.1:length(dpb$data[[1]]$x)],
                                y = dpb$data[[1]]$y[x1.1:length(dpb$data[[1]]$y)]),
              aes(x=x, y = y), fill=alpha("red", 0))+
    geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.1:x2.1],
                                y = dpb$data[[1]]$y[x1.1:x2.1]),
              aes(x=(x), y = y), fill=alpha("grey80", 0))+
    geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.2:x2.2],
                                y = dpb$data[[1]]$y[x1.2:x2.2]),
              aes(x=(x), y = y), fill=alpha("green", 0))+
    geom_line(data = my.df, aes(x = exp(x),y = y), size=.75)


  main.plot <- main.plot +
    labs(title="Posterior distribution for the hazard ratio",
         subtitle = for.subtitle,
         x="Hazard Ratio", y = NULL,
         caption = NULL)+
    annotate("text", label = for.decision,
             x = annotate.x, y = max(dpb$data[[1]]$y)-max(dpb$data[[1]]$y)/nlines.ria * 0, size = tsize+1, colour = result.color, hjust = annotate.j)+
    # annotate("text", label = for.decision.interval,
    #          x = annotate.x, y = max(dpb$data[[1]]$y)-max(dpb$data[[1]]$y)/nlines.ria * 1, size = tsize, colour = result.color, hjust = annotate.j)+
    annotate("text", label = annotate.P1, color=result.color,
             x = annotate.x, y = max(dpb$data[[1]]$y)-max(dpb$data[[1]]$y)/nlines.ria * 2, size = tsize,  hjust = annotate.j)+
    annotate("text", label = annotate.P2, color=result.color,
             x = annotate.x, y = max(dpb$data[[1]]$y)-max(dpb$data[[1]]$y)/nlines.ria * 3, size = tsize,  hjust = annotate.j) +
    scale_y_continuous(breaks=NULL, labels=NULL)+
    scale_x_continuous(limits = c(0,3.01),
                       breaks = pretty(x=c(0,3), n=10))
  # Add reference lines and Credible interval

  # if(result.color == "red"){
  #   main.plot <- main.plot +
  #     geom_segment(data=for.decision.interval.df, aes(x=x, xend=xend, y=y, yend=yend, group=group), arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)
  # } else {
  #   main.plot <- main.plot +
  #     geom_segment(data=for.decision.interval.df, aes(x=x, xend=xend, y=y, yend=yend, group=group),  arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)
  # }

  if(HR.tv < HR.lrv){
    main.plot <- main.plot +
      geom_vline(xintercept = c(HR.lrv, HR.tv), linetype = 2, color = c("blue", "blue"))+
      annotate("text", label =paste0("Max TPP = ", HR.lrv), x = HR.lrv, y = 0 + max(dpb$data[[1]]$y)/nlines.ria, size = tsize, colour = "black", hjust = 0)+
      annotate("text", label = paste0("Base TPP = ",  HR.tv), x = HR.tv, y = 0 + max(dpb$data[[1]]$y)/nlines.ria, size = tsize, colour = "black", hjust = 1)
  } else {
    main.plot <- main.plot +
      geom_vline(xintercept = c(HR.lrv, HR.tv), linetype = 2, color = c("blue", "blue"))+
      annotate("text", label =paste0("Max TPP = ", HR.lrv), x = HR.lrv, y = 0 + max(dpb$data[[1]]$y)/nlines.ria, size = tsize, colour = "black", hjust = 1)+
      annotate("text", label = paste0("Base TPP = ",  HR.tv), x = HR.tv, y = 0 + max(dpb$data[[1]]$y)/nlines.ria, size = tsize, colour = "black", hjust = 0)
  }
  if(HR.tv < HR.lrv){
    table.plot2 <- ggplot()+
      annotate("text", label = paste0("Decision Criteria"),
               x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
      annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
               x = -1, y = 1-2.5/nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = TeX(paste0("P($\\Delta < HR_{Max}$) > ", tau.lrv*100, "%\\;\\;\\;\\;\\, &")),
               x = 0, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
      annotate("text", label =  TeX(paste0("P($\\Delta < HR_{Base}$) > ", tau.tv*100,"%")),
               x = 1.5, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
      annotate("text", label = TeX(paste0("No-Go if:")), color="red",
               x = -1, y = 1-3.25/nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = TeX(paste0("P($\\Delta < HR_{Max}$) <= ", tau.ng*100, "%\\;\\;\\;\\,\\, &")),
               x = 0, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0) +
      annotate("text", label =  TeX(paste0("P($\\Delta < HR_{Base}$) <= ", tau.tv*100,"%")),
               x = 1.5, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0)+
      annotate("text", label = TeX(paste0("Consider:")), color="black",
               x = -1, y = 1-4./nlines, size = tsize+1, hjust = 0)+
      annotate("text", label = "Otherwise",
               x = 0, y = 1-4./nlines, size = tsize+1, colour = "black", hjust = 0)+
      scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
      scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
      labs(x="", y="")} else {
        table.plot2 <- ggplot()+
          annotate("text", label = paste0("Decision Criteria"),
                   x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
          annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
                   x = -1, y = 1-2.5/nlines, size = tsize+1, hjust = 0)+
          annotate("text", label = TeX(paste0("P($\\Delta > HR_{Max}$) > ", tau.lrv*100, "% &")),
                   x = 0, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
          annotate("text", label =  TeX(paste0("P($\\Delta > HR_{Base}$) > ", tau.tv*100,"%")),
                   x = 1.5, y = 1-2.5/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
          annotate("text", label = TeX(paste0("No-Go if:")), color="red",
                   x = -1, y = 1-3.25/nlines, size = tsize+1, hjust = 0)+
          annotate("text", label = TeX(paste0("P($\\Delta > HR_{Max}$) <= ", tau.ng*100, "% &")),
                   x = 0, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0) +
          annotate("text", label =  TeX(paste0("P($\\Delta > HR_{Base}$) <= ", tau.tv*100,"%")),
                   x = 1.5, y = 1-3.25/nlines, size = tsize+1, colour = "red", hjust = 0)+
          annotate("text", label = TeX(paste0("Consider:")), color="black",
                   x = -1, y = 1-4./nlines, size = tsize+1, hjust = 0)+
          annotate("text", label = "Otherwise",
                   x = 0, y = 1-4./nlines, size = tsize+1, colour = "black", hjust = 0)+
          scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
          scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
          labs(x="", y="")
      }



  # plot_grid(main.plot,
  #            table.plot2, ncol=1, align = "v", rel_heights=c(.78,.22), axis="l")

  if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2,
                               data.frame(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior, HR.prior=HR.prior, ARatio=ARatio, m.obs = m.obs,
                                          HR.tv= HR.tv, HR.lrv = HR.lrv, tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng,
                                          Needed.for.NG = result.ng$HR.obs[1], Needed.for.GO = result.go$HR.obs[1]
                                          )
                               ))
  if(add.table==FALSE) return(main.plot)
}
