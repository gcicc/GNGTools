#' @title Make single sample normal-gamma rule in action plot
#' @param mu.0.t prior mean
#' @param n.0.t prior effective sample size
#' @param alpha.0.t prior alpha parameter
#' @param beta.0.t prior beta parameter
#' @param xbar.t observed sample mean
#' @param s.t observed sample standard deviation
#' @param n.t sample size
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param nlines.ria Control for text spacing
#' @param add.table provides extended output summaries

#' @return A ggplot object is returned
#' @export
#' @examples \donttest{
#' my.ss.ng.ria <- make.ss.ng.ria(add.table=TRUE)
#' plot(my.ss.ng.ria[[1]])
#' my.ss.ng.ria[[2]]
#' my.ss.ng.ria[[3]]
#' my.ss.ng.ria[[4]]
#' }
make.ss.ng.ria <- function(mu.0.t = 0, alpha.0.t=.25, beta.0.t = 1, n.0.t = 1,
                           xbar.t = -.05, s.t = 3, n.t = 10,
                           Delta.lrv = 0, Delta.tv = 1,
                           tau.tv=.1, tau.lrv=.8, tau.ng=.65,
                           tsize=4, nlines=25, nlines.ria=20, add.table=TRUE){

        results <- get.ss.ng.studyend.GNG(mu.0.t = mu.0.t, n.0.t = n.0.t, alpha.0.t = alpha.0.t,
                                          beta.0.t = beta.0.t,
                                          xbar.t = seq(Delta.lrv*.75, Delta.tv*1.5, length.out=50),
                                          s.t = s.t, n.t = n.t, Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                                          tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)
        result.go <- results$result.go
        result.ng <- results$result.ng

        # Get post parameters
        pp.t <- get.ng.post(mu.0 = mu.0.t, n.0 = n.0.t, alpha.0 = alpha.0.t,
                            beta.0 = beta.0.t,
                            xbar = xbar.t, s = s.t, n = n.t, group = "Treatment")
        my.df <- gcurve(expr = dt_ls(x = x,df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                     sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),
                        from = pp.t$mu.n - 5 * s.t / sqrt(pp.t$n.n), to = pp.t$mu.n +
                                5 * s.t / sqrt(pp.t$n.n),n = 1001)

        # Compute result
        P.R1 = 1 - pt_ls(x = Delta.lrv, df = pp.t$tdf.n, mu = pp.t$mu.n, sigma = pp.t$sigma.n)
        P.R3 = 1 - pt_ls(Delta.tv, df = pp.t$tdf.n, mu = pp.t$mu.n, sigma = pp.t$sigma.n)
        result = ifelse(P.R1 > tau.lrv & P.R3 > tau.tv, "Go",
                        ifelse(P.R1 < tau.ng  & P.R3 < tau.tv, "No-Go", "Consider"))
        result.color <- ifelse(result=="Go", "darkgreen", ifelse(result=="No-Go", "red", "black"))

        for.subtitle <- TeX(paste0("Treatment data ", "(n, $\\bar{x}$, s): (",
                                   round(n.t,2), ", $", round(xbar.t,2), "$, ",
                                   round(s.t,2),"). Sample mean needed for Go: $",
                                   round(result.go$xbar,2), "$. Needed for No-Go: $",
                                   round(result.ng$xbar,2), "$."))

        my.df$group = paste0("Prior: NG(", round(mu.0.t,2), ", ", round(n.0.t,2),
                             ", ", round(alpha.0.t,2), ", ", round(beta.0.t,2), "); ")

        # Annotation line 1: Decision ----
        for.decision <- paste0("Decision: ", result)

        # Annotation line 2: Decision interval ----
        if(P.R3 > tau.tv){
                for.decision.interval <- paste0(
                        "Decision interval: (",
                        round(qt_ls(prob = 1 - tau.lrv,
                                    df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                    sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),3), ", ",
                        round(qt_ls(prob = 1 - tau.tv,
                                    df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                    sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),3), ")")
        } else {
                for.decision.interval <- paste0(
                        "Decision interval: (",
                        round(qt_ls(prob =  1 - tau.ng,
                                    df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                    sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),3)*100, ", ",
                        round(qt_ls(prob = 1 - tau.tv, df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                    sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),3)*100, ")")
        }

        # Annotation line 3: P(Delta >= Min TPP) ----
        annotate.P1 <- ifelse(
                result=="Go",
                TeX(paste0("P($\\Delta\\,$ >= Min TPP) = ",
                           round(1 - pt_ls(x = Delta.lrv,
                                           df = 2 * pp.t$alpha.n,
                                           mu = pp.t$mu.n,
                                           sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),
                                 4)*100, "% >", tau.lrv*100, "%")),
                ifelse(result=="No-Go",
                       TeX(paste0("P($\\Delta\\,$ >= Min TPP) = ",
                                  round(1 - pt_ls(x = Delta.lrv, df = 2 * pp.t$alpha.n,
                                                  mu = pp.t$mu.n,
                                                  sigma = pp.t$beta.n/(pp.t$alpha.n* pp.t$n.n)),4)*100,
                                  "% < ", tau.ng*100, "%")),
                       TeX(paste0("P($\\Delta\\,$ >= Min TPP) = ",
                                  round(1 - pt_ls(x = Delta.lrv, df = 2 * pp.t$alpha.n,
                                                  mu = pp.t$mu.n,
                                                  sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),
                                        4)*100, "%"))))

        annotate.P2 <- ifelse(
                result=="Go",
                TeX(paste0("P($\\Delta\\,$ >= Base TPP) = ",
                           round(1 - pt_ls(x = Delta.tv,
                                           df = 2 * pp.t$alpha.n,
                                           mu = pp.t$mu.n,
                                           sigma = pp.t$beta.n/(pp.t$alpha.n *pp.t$n.n)),4)*100,
                           " > ", tau.tv*100, "%")),
                ifelse(
                        result=="No-Go",
                        TeX(paste0("P($\\Delta\\,$ > Base TPP) = ",
                                   round(1 - pt_ls(x = Delta.tv,
                                                   df = 2 * pp.t$alpha.n,
                                                   mu = pp.t$mu.n,
                                                   sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),4)*100,
                                   " < ", tau.tv*100,"%")),
                        TeX(paste0("P($\\Delta\\,$ >= Base TPP) = ",
                                   round(1 - pt_ls(x = Delta.tv,
                                                   df = 2 * pp.t$alpha.n,
                                                   mu = pp.t$mu.n,
                                                   sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),
                                         4)*100, "%"))))

        # creat initial build ----
        dplot <- ggplot() + geom_line(data = my.df, aes(x = x, y = y)) +
                facet_wrap(~group)
        # Access the ggplot to get goodies to help accomplish shading
        dpb <- ggplot_build(dplot)
        x1.1 <- min(which(dpb$data[[1]]$x >=Delta.lrv))
        x2.1 <- max(which(dpb$data[[1]]$x <=Delta.tv))+1
        x1.2 <- pmin(min(which(dpb$data[[1]]$x >=Delta.tv)),
                     max(which(dpb$data[[1]]$x <=Inf)))
        x2.2 <- max(which(dpb$data[[1]]$x <=Inf))

        # Custom graphic parameters ----
        x.limits <- c(min(min(dpb$data[[1]]$x), Delta.lrv),
                      max(max(dpb$data[[1]]$x), Delta.tv))
        ticks <- pretty(x = c(x.limits[1], x.limits[2]), n = 15)
        # When the mound is on the right... we want to display annotation on the left ----
        if(which(dpb$data[[1]]$y== max(dpb$data[[1]]$y)) > length(dpb$data[[1]]$x)/2){
                annotate.x = min(min(dpb$data[[1]]$x), Delta.lrv)
                annotate.j = 0
        } else {annotate.x = max(max(dpb$data[[1]]$x), Delta.tv)
        annotate.j = 1}

        # Note these go and no-go segments differ on p = 1 - tau.lrv vs. p = 1 - tau.ng ----
        go.segment <- data.frame(x = qt_ls(prob = 1 - tau.lrv,
                                           df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                           sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),
                                 xend = qt_ls(prob = 1 - tau.tv,
                                              df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                              sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),
                                 y= min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3,
                                 yend = min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3,
                                 group=my.df$group[1])
        nogo.segment <- data.frame(x = qt_ls(prob = 1 - tau.ng,
                                             df = 2 * pp.t$alpha.n,
                                             mu = pp.t$mu.n,
                                             sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),
                                   xend = qt_ls(prob = 1 - tau.lrv,
                                                df = 2 * pp.t$alpha.n, mu = pp.t$mu.n,
                                                sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n)),
                                   y= min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3,
                                   yend = min(dpb$data[[1]]$y) + max(dpb$data[[1]]$y)/nlines.ria * 3,
                                   group=my.df$group[1])
        # Introduce shading ----
        main.plot <- dplot +
                geom_area(data = data.frame(x = dpb$data[[1]]$x[1:x1.1],
                                            y = dpb$data[[1]]$y[1:x1.1]),
                          aes(x = x, y = y), fill=alpha("red",0))+
                geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.1:x2.1],
                                            y = dpb$data[[1]]$y[x1.1:x2.1]),
                          aes(x = x, y = y), fill=alpha("grey80",0))+
                geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.2:x2.2],
                                            y = dpb$data[[1]]$y[x1.2:x2.2]),
                          aes(x = x, y = y), fill=alpha("green",0))+
                geom_vline(xintercept = c(Delta.lrv, Delta.tv), linetype = 2, color = c("blue", "blue"))+
                scale_y_continuous(breaks=NULL) +
                scale_x_continuous(limits = x.limits, breaks=pretty(x = x.limits, n=10))

        # Add Annotations -----
        main.plot <- main.plot +
                labs(title = "Posterior distribution for the treatment effect",
                     x = "Mean treatment response", y = NULL,
                     subtitle= for.subtitle,
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

        # Add reference lines and Credible interval----
        # if(result.color == "red"){
        #         main.plot <- main.plot +
        #                 geom_segment(data=nogo.segment, aes(x=x, xend=xend, y=y, yend=yend, group=group),
        #                              arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)
        # } else {
        #         main.plot <- main.plot +
        #                 geom_segment(data=go.segment, aes(x=x, xend=xend, y=y, yend=yend, group=group),
        #                              arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)
        # }

        main.plot <- main.plot +
                geom_vline(xintercept = c(Delta.lrv, Delta.tv), linetype = 2, color = c("blue", "blue")) +
                annotate("text", label = paste0("Min TPP = ", Delta.lrv), x = Delta.lrv,
                         y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize, colour = "black", hjust = 1)+
                annotate("text", label = paste0("Base TPP = ",  Delta.tv), x = Delta.tv,
                         y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize, colour = "black", hjust = 0)
        # Create table plot ----
        # In the shiny app we replicated the table with HTML - but retain this here
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

        # return ----
        if(add.table==TRUE) return(list(grid.arrange(main.plot, table.plot2, nrow=2, heights=c(.78,.22)), main.plot, table.plot2,
                                     data.frame(mu.0.t = mu.0.t, alpha.0.t=alpha.0.t, beta.0.t = beta.0.t, n.0.t = n.0.t,
                                                s.t = s.t, n.t = n.t,
                                                Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                                tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng,
                                                Needed.for.NG = result.ng$xbar[1], Needed.for.GO = result.go$xbar)))
        if(add.table==FALSE) return(main.plot)
}
