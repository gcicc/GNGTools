#' @title Make two-sample normal-gamma rule in action plot
#'
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size parameter for control group
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size parameter for treatment group
#' @param xbar.c control mean
#' @param s.c control sd
#' @param n.c control sample size
#' @param xbar.t treatment mean
#' @param s.t treatment sd
#' @param n.t treatment sample size
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param seed random seed
#' @param n.MC n for MC sampling
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param nlines.ria Control for text spacing
#' @param add.table provides extended output summaries
#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' my.ts.ng.ria <- make.ts.ng.ria(add.table=TRUE)
#' plot(my.ts.ng.ria[[1]])
#' my.ts.ng.ria[[2]]
#' my.ts.ng.ria[[3]]
#' my.ts.ng.ria[[4]]

make.ts.ng.ria <- function(mu.0.c = 0, alpha.0.c = 2.5, beta.0.c = 10, n.0.c = 10,
                           mu.0.t = 0, alpha.0.t = 0.25, beta.0.t = 1,n.0.t = 1e-04,
                           xbar.c = 0.25, s.c = 1.5, n.c = 40,
                           xbar.t = 2.5,  s.t = 1.5, n.t = 40,
                           Delta.lrv = 1.5,  Delta.tv = 3,
                           tau.ng = 0.65,  tau.lrv = 0.8,  tau.tv = 0.2,
                           seed = 1234,
                           n.MC = 1000,
                           nlines = 25,
                           nlines.ria = 20,
                           tsize = 4,
                           add.table = TRUE){
  set.seed(seed=seed)
  # Smart two-stage search-----
  my.means <- seq(-Delta.tv*4, Delta.tv*4, length.out=50)
  stage1 <- get.ts.ng.mc.df(mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c = alpha.0.c,
                            beta.0.c = beta.0.c,
                            xbar.c = xbar.c, s.c = s.c, n.c = n.c, group.c = "Control",
                            mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t,
                            beta.0.t = beta.0.t,
                            xbar.t = xbar.c + my.means, s.t = s.t, n.t = n.t,
                            group.t="treat",
                            Delta.tv = Delta.tv,Delta.lrv = Delta.lrv,tau.tv = tau.tv,
                            tau.lrv = tau.lrv,tau.ng = tau.ng, n.MC = n.MC)

  stage1.go <- stage1 %>% dplyr::filter(result== "Go") %>% slice(1)
  # Determine max number of TRT responders for No-Go
  stage1.ng <- stage1 %>% dplyr::filter(result== "No-Go") %>% slice(n())

  my.means <- c()
  if (nrow(stage1.go) > 0 ) my.means <- c(my.means, seq(stage1.go$xbar.t - stage1.go$s.t/4, stage1.go$xbar.t, length.out=100))
  if (nrow(stage1.ng) > 0 ) my.means <- c(my.means, seq(stage1.ng$xbar.t , stage1.ng$xbar.t + stage1.ng$s.t/4, length.out=100))

  stage2 <- get.ts.ng.mc.df(mu.0.c = mu.0.c, n.0.c = n.0.c, alpha.0.c = alpha.0.c,
                            beta.0.c = beta.0.c,
                            xbar.c = xbar.c, s.c = s.c, n.c = n.c, group.c = "Control",
                            mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t,
                            beta.0.t = beta.0.t,
                            xbar.t = my.means, s.t = s.t, n.t = n.t, group.t="treat",
                            Delta.tv = Delta.tv,Delta.lrv = Delta.lrv,tau.tv = tau.tv,
                            tau.lrv = tau.lrv,tau.ng = tau.ng,
                            n.MC = n.MC)
  result.go <- stage2 %>% dplyr::filter(result== "Go") %>% slice(1)
  result.ng <- stage2 %>% dplyr::filter(result== "No-Go") %>% slice(n())

  # Get post parameters
  pp.c <- get.ng.post(mu.0 = mu.0.c, n.0 = n.0.c, alpha.0 = alpha.0.c, beta.0 = beta.0.c,
                      xbar = xbar.c, s = s.c, n = n.c, group="Control")
  pp.t <- get.ng.post(mu.0 = mu.0.t, n.0 = n.0.t, alpha.0 = alpha.0.t, beta.0 = beta.0.t,
                      xbar = xbar.t, s = s.t, n = n.t, group="Treatment")

  # Get samples from marginal t-distribution
  samples <- data.frame(samples =
                          rt_ls(n = 100000, df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n))-
                          rt_ls(n = 100000, df = 2 * pp.c$alpha.n, mu = pp.c$mu.n,
                                sigma = pp.c$beta.n/(pp.c$alpha.n * pp.c$n.n))
  )

  # # Get samples from marginal t-distribution- all similar
  # samples1 <- data.frame(samples =
  #                         rt_ls(n = 100000, df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
  #                               sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n))-
  #                         rt_ls(n = 100000, df = 2 * pp.c$alpha.n, mu = pp.c$mu.n,
  #                               sigma = pp.c$beta.n/(pp.c$alpha.n * pp.c$n.n))
  # )
  #
  # # Get samples from normal distribution with marginal-t based parameters
  # samples2 <- data.frame(samples =
  #                         rnorm(n = 100000, mean = pp.t$mu.n,
  #                               sd = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n))-
  #                         rnorm(n = 100000, mean = pp.c$mu.n,
  #                               sd = pp.c$beta.n/(pp.c$alpha.n * pp.c$n.n))
  # )
  #
  # # get samples from normal distribution with marginal t-based mean and sampling distribution variability
  # samples3 <- data.frame(samples = rnorm(n = 100000, mean = pp.t$mu.n,
  #                                        sd = s.t/sqrt(n.t))-
  #                          rnorm(n = 100000, mean = pp.c$mu.n,
  #                                sd = s.c/sqrt(n.c)))
  # bind_rows(samples1 %>% mutate(V=1), samples2 %>% mutate(V=2), samples3 %>% mutate(V=3)) %>% ggplot(aes(x=samples, fill=factor(V)))+geom_density(alpha=.25)

  samples$group <- paste0(
    paste0("Control Prior: NG(",  round(pp.c$mu.0, 2), ", ", round(pp.c$n.0, 2), ", ",
           round(pp.c$alpha.0,2), ", ", round(pp.c$beta.0,2), ")"),
    paste0("; Treatment Prior: NG(",  round(pp.t$mu.0, 2), ", ", round(pp.t$n.0, 2), ", ",
           round(pp.t$alpha.0,2), ", ", round(pp.t$beta.0,2), ")")
  )

  P.R1 = mean(samples$samples > Delta.lrv)
  P.R3 = mean(samples$samples > Delta.tv)

  # comput result -----

  result = ifelse(P.R1 > tau.lrv & P.R3 > tau.tv, "Go",
                  ifelse(P.R1 <= tau.ng & P.R3 <= tau.tv, "No-Go",
                         "Consider"))
  result.color <- ifelse(result=="Go", "darkgreen", ifelse(result=="No-Go",
                                                           "red", "black"))

  # Subtitle ----

  for.subtitle <- TeX(paste0("Control, Treatment Data (n, $\\bar{x}_{T}$, s): (",
                             round(n.c,2), ", ", round(xbar.c,2), ", ", round(s.c,2),"), ",
                             "(", round(n.t,2), ", ", round(xbar.t,2), ", ", round(s.t,2),
                             "); $\\bar{x}_{T}$ needed for Go: ",round(result.go$xbar.t[1],2),
                             "; for No-Go: ", round(result.ng$xbar.t[1],2) , ""))

  # Annotation line 1: Decision ----
  for.decision <- paste0("Decision: ", result)

  # Annotation line 2: Decision interval----

  for.decision.interval <- ifelse(mean(samples$samples > Delta.tv) > tau.tv,
                                  paste0("Decision interval: (",
                                         round(quantile(samples$samples,
                                                        1 - tau.lrv), 2), ", ",
                                         round(quantile(samples$samples,
                                                        1 - tau.tv), 2), ")"),
                                  paste0("Decision interval: (",
                                         round(quantile(samples$samples,
                                                        1 - tau.ng), 2), ", ",
                                         round(quantile(samples$samples,
                                                        1 - tau.tv), 2), ")"))

  # Annotation line 3: P(Delta >= Min TPP)

  annotate.P1 <- ifelse(result=="Go",
                        TeX(paste0("P($\\Delta$ >= Min TPP) = ",
                                   round(mean(samples$samples > Delta.lrv),4)*100,
                                   "% > ", tau.lrv*100,"%")),
                        ifelse(result=="No-Go",
                               TeX(paste0("P($\\Delta$ >= Min TPP) = ",
                                          round(mean(samples$samples > Delta.lrv),4)*100 ,
                                          "% <= ", tau.ng*100,"%")),
                               TeX(paste0("$P(\\Delta$ >= Min TPP) = ",
                                          round(mean(samples$samples > Delta.lrv),4)*100,
                                          "%" ))))

  # Annotation line 3: P(Delta >= Base TPP)

  annotate.P2 <- ifelse(result=="Go",
                        TeX(paste0("P($\\Delta$ >= Base TPP) = ",
                                   round(mean(samples$samples > Delta.tv),4)*100 ,
                                   "% >=", tau.tv*100, "%")),
                        ifelse(result=="No-Go",
                               TeX(paste0("P($\\Delta$ >= Base TPP) = ",
                                          round(mean(samples$samples > Delta.tv),
                                                4)*100
                                          , "% <", tau.tv*100, "%")),
                               TeX(paste0("P($\\Delta$ >= Base TPP) = ",
                                          round(mean(samples$samples > Delta.tv),
                                                4)*100, "%"  ))))

  dplot <- ggplot() + geom_density(data = samples, aes(x = samples))+
    facet_wrap(~group)
  # Access the ggplot to get goodies to help accomplish shading
  dpb <- ggplot_build(dplot)

  x1.1 <- min(which(dpb$data[[1]]$x >=Delta.lrv))
  x2.1 <- max(which(dpb$data[[1]]$x <=Delta.tv))+1
  x1.2 <- pmin(min(which(dpb$data[[1]]$x >=Delta.tv)),
               max(which(dpb$data[[1]]$x <=Inf)))
  x2.2 <- max(which(dpb$data[[1]]$x <=Inf))

  go.segment <- data.frame(x = quantile(samples$samples, 1 - tau.lrv),
                           xend =quantile(samples$samples, 1 - tau.tv),
                           y= min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3,
                           yend = min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3) %>%
    mutate(group= paste0(
      paste0("Control Prior: NG(",  round(pp.c$mu.0, 2), ", ", round(pp.c$n.0, 2), ", ",
             round(pp.c$alpha.0,2), ", ", round(pp.c$beta.0,2), ")"),
      paste0("; Treatment Prior: NG(",  round(pp.t$mu.0, 2), ", ", round(pp.t$n.0, 2),
             ", ", round(pp.t$alpha.0,2), ", ", round(pp.t$beta.0,2), ")")
    ))

  nogo.segment <- data.frame(x = quantile(samples$samples, 1 - tau.ng),
                             xend =quantile(samples$samples, 1 - tau.tv),
                             y= min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3,
                             yend = min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3) %>%
    mutate(group= paste0(
      paste0("Control Prior: NG(",  round(pp.c$mu.0, 2), ", ", round(pp.c$n.0, 2), ", ",
             round(pp.c$alpha.0,2), ", ", round(pp.c$beta.0,2), ")"),
      paste0("; Treatment Prior: NG(",  round(pp.t$mu.0, 2), ", ", round(pp.t$n.0, 2),
             ", ", round(pp.t$alpha.0,2), ", ", round(pp.t$beta.0,2), ")")
    ))

  # Custom graphic parameters
  x.limits <- c(min(min(dpb$data[[1]]$x), Delta.lrv-s.c/5), max(max(dpb$data[[1]]$x),
                                                          Delta.tv+s.c/5))
  ticks <- pretty(x = c(x.limits[1], x.limits[2]), n = 15)

  # When the mound is on the right... we want to display annotation on the left
  if(dpb$data[[1]]$x[which(dpb$data[[1]]$y== max(dpb$data[[1]]$y))] > mean(x.limits)){
    annotate.x = min(min(dpb$data[[1]]$x), Delta.lrv-s.c/5)
    annotate.j = 0
  } else {annotate.x = max(max(dpb$data[[1]]$x), Delta.tv+s.t/5)
  annotate.j = 1}

  # Introduce shading ----
  main.plot <- dplot +
    # NOTE: GREG SHUT THIS OPTION OFF 3/27/20 BECAUSE LARGE XBAR_TRT VALUES WOULD CAUSE GRAPHIC TO CRASH
    # ALSO ISSUE IS RELATED TO SHADING WHICH HAD BEEN SHUT OFF WITH CODE IN PLACE BY SETTING FILL ALPHA TO 0.
    # geom_area(data = data.frame(x = dpb$data[[1]]$x[1:x1.1],
    #                             y = dpb$data[[1]]$y[1:x1.1]),
    #           aes(x = x, y = y), fill=alpha("red", 0))+
    # geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.1:x2.1],
    #                             y = dpb$data[[1]]$y[x1.1:x2.1]),
    #           aes(x = x, y = y), fill=alpha("grey80", 0))+
    # geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.2:x2.2],
    #                             y = dpb$data[[1]]$y[x1.2:x2.2]),
    #           aes(x = x, y = y), fill=alpha("green", 0))+
  scale_x_continuous(limits = x.limits, breaks=pretty(x=x.limits, n=10))+
    scale_y_continuous(breaks=NULL)

  # Add Annotations -----
  main.plot <- main.plot +
    labs(title = TeX("Posterior Density of Treatment Effect"),
         x = TeX("$\\Delta$ = Mean treatment difference (Treatment - Control)"),
         y = NULL,
         subtitle=for.subtitle)+
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
  # if(mean(samples$samples > Delta.tv) > tau.tv) {
  #   main.plot <- main.plot +
  #     geom_segment(data=go.segment, aes(x=x, xend=xend, y=y, yend=yend),
  #                  arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)
  # } else {
  #   main.plot <- main.plot +
  #     geom_segment(data=nogo.segment, aes(x=x, xend=xend, y=y, yend=yend, group=group),
  #                  arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)
  #
  # }
  main.plot <- main.plot +
    geom_vline(xintercept = c(Delta.lrv, Delta.tv), linetype = 2, color = c("blue", "blue")) +
    annotate("text", label = TeX(paste0("", Delta.lrv," = Min TPP")), x = Delta.lrv, y = 0 + max(dpb$data[[1]]$y)/nlines.ria, size = tsize, colour = "black", hjust = 1)+
    annotate("text", label = TeX(paste0("Base TPP = ",  Delta.tv,"")), x = Delta.tv, y = 0 + max(dpb$data[[1]]$y)/nlines.ria, size = tsize, colour = "black", hjust = 0)


  table.plot2 <- ggplot()+
    annotate("text", label = paste0("Decision Criteria"),
             x = -1, y = .95, size = tsize+2, colour = "black", hjust=0)+
    annotate("text", label = TeX(paste0("Go if:")), color="darkgreen",
             x = -1, y = 1-3/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta$ >= Min TPP) > ", tau.lrv*100, "% &")),
             x = 0, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0) +
    annotate("text", label =  TeX(paste0("P($\\Delta$ >= Base TPP) > ", tau.tv*100,"%")),
             x = 1.5, y = 1-3/nlines, size = tsize+1, colour = "darkgreen", hjust = 0)+
    annotate("text", label = TeX(paste0("No-Go if:")), color="red",
             x = -1, y = 1-4/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = TeX(paste0("P($\\Delta$ >= Min TPP) <= ", tau.ng*100, "% &")),
             x = 0, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0) +
    annotate("text", label = TeX(paste0("P($\\Delta$ >= Base TPP) <= ", tau.tv*100,"%")),
             x = 1.5, y = 1-4/nlines, size = tsize+1, colour = "red", hjust = 0)+
    annotate("text", label = TeX(paste0("Consider:")), color="black",
             x = -1, y = 1-5/nlines, size = tsize+1, hjust = 0)+
    annotate("text", label = "Otherwise",
             x = 0, y = 1-5/nlines, size = tsize+1, colour = "black", hjust = 0)+
    scale_y_continuous(limits=c(0.75,.975), expand=c(0,0), labels=NULL, breaks=NULL)+
    scale_x_continuous(limits=c(-1,3.5), expand=c(0,0), labels=NULL, breaks=NULL)+
    labs(x="\n", y="")


  if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2,
                               data.frame(mu.0.c = mu.0.c, alpha.0.c = alpha.0.c, beta.0.c = beta.0.c, n.0.c = n.0.c,
                                          mu.0.t = mu.0.t, alpha.0.t = alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                                          xbar.c = xbar.c, s.c = s.c, n.c = n.c,
                                          Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                          tau.ng = tau.ng, tau.lrv = tau.lrv, tau.tv = tau.tv,
                                          seed = seed, n.MC = n.MC,
                                          Needed.for.NG = result.ng$xbar.t[1],
                                          Needed.for.GO = result.go$xbar.t[1] )))
  if(add.table==FALSE) return(main.plot)
}

