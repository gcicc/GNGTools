#' @title Make two-sample binary rule in action plot
#'
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param n.con sample size for control group
#' @param x.con number of responders for control group
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior beta parameter for treatment group
#' @param n.trt sample size for control treatment group
#' @param x.trt number of responders for treatment group
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#' @param nlines.ria Control for text spacing
#' @param add.table provides extended output summaries
#' @return a ggplot object is returned
#' @export
#'
#' @examples
#' my.ts.bin.ria <- make.ts.bin.ria(add.table=TRUE)
#' plot(my.ts.bin.ria[[1]])
#' my.ts.bin.ria[[2]]
#' my.ts.bin.ria[[3]]
#' my.ts.bin.ria[[4]]

make.ts.bin.ria <- function(a.con = 1, b.con = 1, n.con = 40, x.con = 9,
                            a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 17,
                            Delta.lrv = 0.15, Delta.tv = .30,
                            tau.tv = 0.10, tau.lrv = 0.80, tau.ng = 0.65,
                            nlines.ria = 20, tsize = 4, nlines = 25, add.table=TRUE)
{
        # Run this for fixed value of CON
        results <- get.ts.bin.dec.df(a.con = a.con, b.con = b.con, n.con = n.con, x.con = x.con,
                                          a.trt = a.trt, b.trt = b.trt, n.trt = n.trt, x.trt = 0:n.trt,
                                          Delta.tv = Delta.tv, Delta.lrv = Delta.lrv,
                                          tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

        # Determine mi/max number of TRT responders for Go/No G0 ----
        result.go <- results %>% dplyr::filter(result== "Go") %>% slice(1)
        result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(n())

        probs <- 1 - p2beta(relation="DIFF", approach="DIRECT", x=c(Delta.lrv, Delta.tv),
                            a1=a.con+x.con, b1=b.con+n.con-x.con,
                            a2=a.trt+x.trt, b2=b.trt+n.trt-x.trt, n = 1000000)

        # comput result -----
        result <- ifelse(probs[1] >= tau.lrv & probs[2] >= tau.tv, "Go",
                         ifelse(probs[1] < tau.ng & probs[2] < tau.tv, "No-Go", "Consider"))
        result.color <- ifelse(result=="Go", "darkgreen", ifelse(result=="No-Go", "red", "black"))

        # Subtitle ----
        for.subtitle <- paste0("Control data: ", round(x.con/n.con*100,2), "% (",  x.con, "/", n.con, "). ",
                               "Treatment data: ",round(x.trt/n.trt*100,2), "% (",  x.trt, "/", n.trt, "). ",
                               "Observed difference: ", round(x.trt/n.trt*100,2) - round(x.con/n.con*100,2),
                               "%\nGiven control data, treatment responders needed for Go: ", result.go$x.trt[1], " (",round(result.go$x.trt[1]/n.trt*100,2), "%).",
                               "Needed for No-Go ", result.ng$x.trt[1]," (", round(result.ng$x.trt/n.trt*100,2), "%)")

        # Annotation line 1: Decision ----
        for.decision <- paste0("Decision: ", result)


        # Annotation line 3: P(Delta >= Min TPP)
        annotate.P1 <- ifelse(result=="Go",
                              TeX(paste0("$P(\\Delta\\,$ >= Min TPP) = $",  round((probs[1])*100,1), "$% > ", (tau.lrv)*100,"%" )),
                              ifelse(result=="No-Go",
                                     TeX(paste0("$P(\\Delta\\,$ >= Min TPP) = ",  round((probs[1])*100,1), "% <= ", (tau.ng)*100,"%")),
                                     TeX(paste0("$P(\\Delta\\,$ >= Min TPP) = ",  round((probs[1])*100,1), "%"))))
        # Annotation line 4: P(Delta >= Base TPP)

        annotate.P2 <- ifelse(result=="Go",
                              TeX(paste0("$P(\\Delta\\,$ >= Base TPP) = ",  round(probs[2]*100,1), "% >", tau.tv*100,"%")),
                              ifelse(result=="No-Go", TeX(paste0("$P(\\Delta\\,$ >= Base TPP) = ",  round(probs[2]*100, 1), "% <=", tau.tv*100,"%")),
                                     TeX(paste0("$P(\\Delta\\,$ >= Base TPP) = ",  round(probs[2]*100, 1),"%"))))


        my.df <- try(data.frame(x = c(seq(-1, -0.001, by = 0.001), seq(0.001, 1, by = 0.001))) %>%
                             mutate(y = (d2beta(relation='DIFF', x=x, a1=a.con+x.con, b1=b.con+n.con-x.con, a2=a.trt+x.trt, b2=b.trt+n.trt-x.trt)))
        )

        if(as.character(class(my.df)) != "try-error"){
                my.df$group <- c(paste0("Control Prior: Beta(", a.con, ", ", b.con, "). Treatment Prior: Beta(", a.trt, ", ", b.trt, ")"))

                # Compute probabilities
                # Probs are reported as exceeding
                probs <- 1 - p2beta(relation="DIFF", approach="DIRECT", x=c(Delta.lrv, Delta.tv),
                                    a1=a.con+x.con, b1=b.con+n.con-x.con,
                                    a2=a.trt+x.trt, b2=b.trt+n.trt-x.trt, n = 100000)

                # comput result -----
                result <- ifelse(probs[1] >= tau.lrv & probs[2] >= tau.tv, "Go",
                                 ifelse(probs[1] < tau.ng & probs[2] < tau.tv, "No-Go", "Consider"))
                result.color <- ifelse(result=="Go", "darkgreen", ifelse(result=="No-Go", "red", "black"))

                # Subtitle ----
                for.subtitle <- paste0("Control data: ", round(x.con/n.con*100,2), "% (",  x.con, "/", n.con, "). ",
                                       "Treatment data: ",round(x.trt/n.trt*100,2), "% (",  x.trt, "/", n.trt, "). ",
                                       "Observed difference: ", round(x.trt/n.trt*100,2) - round(x.con/n.con*100,2),
                                       "%\nGiven control data, treatment responders needed for Go: ", result.go$x.trt[1], " (",round(result.go$x.trt[1]/n.trt*100,2), "%). ",
                                       "Needed for No-Go ", result.ng$x.trt[1]," (", round(result.ng$x.trt/n.trt*100,2), "%)")

                # Annotation line 1: Decision ----
                for.decision <- paste0("Decision: ", result)


                # Annotation line 3: P(Delta >= Min TPP)
                annotate.P1 <- ifelse(result=="Go",
                                      TeX(paste0("$P(\\Delta\\,$ >= Min TPP) = ",  round((probs[1])*100,1), "% > ", (tau.lrv)*100,"%" )),
                                      ifelse(result=="No-Go",
                                             TeX(paste0("$P(\\Delta\\,$ >= Min TPP) = ",  round((probs[1])*100,1), "% <= ", (tau.ng)*100,"%")),
                                             TeX(paste0("$P(\\Delta\\,$ >= Min TPP) = ",  round((probs[1])*100,1), "%"))))
                # Annotation line 4: P(Delta >= Base TPP)

                annotate.P2 <- ifelse(result=="Go",
                                      TeX(paste0("$P(\\Delta\\,$ >= Base TPP) = ",  round(probs[2]*100,1), "% >", tau.tv*100,"%")),
                                      ifelse(result=="No-Go", TeX(paste0("$P(\\Delta$ >= Base TPP) = ",  round(probs[2]*100, 1), "% <=", tau.tv*100,"%")),
                                             TeX(paste0("$P(\\Delta\\,$ >= Base TPP) = ",  round(probs[2]*100, 1),"%"))))

                # Initialize a ggplot
                dplot <- ggplot() + geom_line(data = my.df, aes(x = x, y=y)) + facet_wrap(~group)
                # Access the ggplot to get goodies to help accomplish shading
                dpb <- ggplot_build(dplot)
                x1.1 <- min(which(dpb$data[[1]]$x >=Delta.lrv))
                x2.1 <- max(which(dpb$data[[1]]$x <=Delta.tv))+1
                x1.2 <- min(which(dpb$data[[1]]$x >=Delta.tv))
                x2.2 <- max(which(dpb$data[[1]]$x <=1))


                beta.diff.quants <- q2beta(relation="DIFF", a1=a.con+x.con, b1=b.con+n.con-x.con,
                                           a2=a.trt+x.trt, b2=b.trt+n.trt-x.trt, alpha=c(1-tau.ng, 1-tau.lrv, 1 - tau.tv), tol = 10^(-5))

                # if(beta.diff.quants[3]>=Delta.tv){
                #         # Annotation line 2: Decision interval ----
                #         for.decision.interval <- paste0("Decision interval: (", round(beta.diff.quants[2]*100,1), "%, ", round(beta.diff.quants[3]*100,1), "%)")
                #         report.segment <- data.frame(x = beta.diff.quants[2],
                #                                      xend = beta.diff.quants[3],
                #                                      y= min(dpb$data[[1]]$y, na.rm=T) + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3,
                #                                      yend = min(dpb$data[[1]]$y, na.rm=T) + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3,
                #                                      group=my.df$group[1])
                # } else {  for.decision.interval <- paste0("Decision interval: (", round(beta.diff.quants[1]*100,1), "%, ", round(beta.diff.quants[3]*100,1), "%)")
                # report.segment <- data.frame(x = beta.diff.quants[1],
                #                              xend = beta.diff.quants[3],
                #                              y= min(dpb$data[[1]]$y, na.rm=T) + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3,
                #                              yend = min(dpb$data[[1]]$y, na.rm=T) + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3,
                #                              group=my.df$group[1])
                # }

                if(which(dpb$data[[1]]$y== max(dpb$data[[1]]$y)) > length(dpb$data[[1]]$x)/2){
                        annotate.x = min(min(dpb$data[[1]]$x), Delta.lrv,-1)
                        annotate.j = 0
                } else {annotate.x = max(max(dpb$data[[1]]$x), Delta.tv, 1)
                annotate.j = 1}

                # Introduce shading ----
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
                        geom_line(data = my.df, aes(x = x, y=y))+
                        scale_x_continuous(limits = c(-1,1),
                                           breaks = pretty(x=c(-1,1), n=10),
                                           labels = scales::percent)+
                        scale_y_continuous(breaks=NULL, labels=NULL)

                # Add Annotations -----
                main.plot <- main.plot +
                        labs(title = TeX("Posterior distribution for the treatment effect"),
                             subtitle = for.subtitle,
                             x = TeX("$\\Delta\\,$ = Treatment difference (Treatment - Control)"), y = NULL,
                             caption = NULL)+
                        annotate("text", label = for.decision,
                                 x = annotate.x, y = max(dpb$data[[1]]$y, na.rm=T)-max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 0, size = tsize+1, colour = result.color, hjust = annotate.j)+
                        # annotate("text", label = for.decision.interval,
                        #          x = annotate.x, y = max(dpb$data[[1]]$y, na.rm=T)-max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 1, size = tsize, colour = result.color, hjust = annotate.j)+
                        annotate("text", label = annotate.P1, color=result.color,
                                 x = annotate.x, y = max(dpb$data[[1]]$y, na.rm=T)-max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 2, size = tsize,  hjust = annotate.j)+
                        annotate("text", label = annotate.P2, color=result.color,
                                 x = annotate.x, y = max(dpb$data[[1]]$y, na.rm=T)-max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3, size = tsize,  hjust = annotate.j)

                # Add reference lines and Credible interval
                # main.plot <- main.plot +
                #         geom_segment(data=report.segment, aes(x=x, xend=xend, y=y, yend=yend, group=group), arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)

                main.plot <- main.plot +
                        geom_vline(xintercept = c(Delta.lrv, Delta.tv), linetype = 2, color = c("blue", "blue")) +
                        annotate("text", label = TeX(paste0("Min TPP = ", Delta.lrv*100,"%")), x = Delta.lrv, y = 0 + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria, size = tsize, colour = "black", hjust = 1)+
                        annotate("text", label = TeX(paste0("Base TPP = ",  Delta.tv*100,"%")), x = Delta.tv, y = 0 + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria, size = tsize, colour = "black", hjust = 0)
        } else {
                X1 <- rbeta(n=10^6, shape1 = a.con+x.con, shape2 = b.con+n.con-x.con)
                X2 <- rbeta(n=10^6, shape1 = a.trt+x.trt, shape2 = b.trt+n.trt-x.trt)
                my.df <- data.frame(x=X2 - X1, group=c(paste0("Control Prior: Beta(", a.con, ", ", b.con, "). Treatment Prior: Beta(", a.trt, ", ", b.trt, ")")))

                dplot <- ggplot() + geom_density(data = my.df, aes(x = x)) + facet_wrap(~group)
                # Access the ggplot to get goodies to help accomplish shading
                dpb <- ggplot_build(dplot)
                x1.1 <- ifelse(length(which(dpb$data[[1]]$x >=Delta.lrv)) == 0, length(dpb$data[[1]]$x >=Delta.lrv),  min(which(dpb$data[[1]]$x >=Delta.lrv)))
                x2.1 <- max(which(dpb$data[[1]]$x <=Delta.tv))
                x1.2 <- ifelse(length(which(dpb$data[[1]]$x >=Delta.tv)) == 0, length(dpb$data[[1]]$x >=Delta.tv),  min(which(dpb$data[[1]]$x >=Delta.tv)))
                x2.2 <- max(which(dpb$data[[1]]$x <=1))

                beta.diff.quants <- quantile(x = X2 - X1, probs = c(1-tau.ng, 1-tau.lrv, 1 - tau.tv))
                # beta.diff.quants <- q2beta(relation="DIFF", a1=a.con+x.con, b1=b.con+n.con-x.con,
                #                            a2=a.trt+x.trt, b2=b.trt+n.trt-x.trt, alpha=c(1-tau.ng, 1-tau.lrv, 1 - tau.tv), tol = 10^(-5))

                if(beta.diff.quants[3]>=Delta.tv){
                        # Annotation line 2: Decision interval ----
                        for.decision.interval <- paste0("Decision interval: (", round(beta.diff.quants[2]*100,1), "%, ", round(beta.diff.quants[3]*100,1), "%)")
                        report.segment <- data.frame(x = beta.diff.quants[2],
                                                     xend = beta.diff.quants[3],
                                                     y= min(dpb$data[[1]]$y, na.rm=T) + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3,
                                                     yend = min(dpb$data[[1]]$y, na.rm=T) + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3,
                                                     group=my.df$group[1])
                } else {  for.decision.interval <- paste0("Decision interval: (", round(beta.diff.quants[1]*100,1), "%, ", round(beta.diff.quants[3]*100,1), "%)")
                report.segment <- data.frame(x = beta.diff.quants[1],
                                             xend = beta.diff.quants[3],
                                             y= min(dpb$data[[1]]$y, na.rm=T) + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3,
                                             yend = min(dpb$data[[1]]$y, na.rm=T) + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3,
                                             group=my.df$group[1])
                }

                if(which(dpb$data[[1]]$y== max(dpb$data[[1]]$y)) > length(dpb$data[[1]]$x)/2){
                        annotate.x = min(min(dpb$data[[1]]$x), Delta.lrv,-1)
                        annotate.j = 0
                } else {annotate.x = max(max(dpb$data[[1]]$x), Delta.tv, 1)
                annotate.j = 1}

                # Introduce shading ----
                main.plot <- dplot +
                        geom_area(data = data.frame(x = dpb$data[[1]]$x[1:x1.1],
                                                    y = dpb$data[[1]]$y[1:x1.1]),
                                  aes(x = x, y = y), fill=alpha("red", 0.0))+
                        geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.1:x2.1],
                                                    y = dpb$data[[1]]$y[x1.1:x2.1]),
                                  aes(x = x, y = y), fill=alpha("grey80", 0.0))+
                        geom_area(data = data.frame(x = dpb$data[[1]]$x[x1.2:x2.2],
                                                    y = dpb$data[[1]]$y[x1.2:x2.2]),
                                  aes(x = x, y = y), fill=alpha("green", 0.0))+
                        geom_density(data = my.df, aes(x = x))+
                        scale_x_continuous(limits = c(-1,1),
                                           breaks = pretty(x=c(-1,1), n=10),
                                           labels = scales::percent)+
                        scale_y_continuous(breaks=NULL, labels=NULL)

                # Add Annotations -----
                main.plot <- main.plot +
                        labs(title = TeX("Posterior distribution for the treatment effect"),
                             subtitle = for.subtitle,
                             x = TeX("$\\Delta\\,$ = Treatment difference (Treatment - Control)"), y = NULL,
                             caption = NULL)+
                        annotate("text", label = for.decision,
                                 x = annotate.x, y = max(dpb$data[[1]]$y, na.rm=T)-max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 0, size = tsize+1, colour = result.color, hjust = annotate.j)+
                        # annotate("text", label = for.decision.interval,
                        #          x = annotate.x, y = max(dpb$data[[1]]$y, na.rm=T)-max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 1, size = tsize, colour = result.color, hjust = annotate.j)+
                        annotate("text", label = annotate.P1, color=result.color,
                                 x = annotate.x, y = max(dpb$data[[1]]$y, na.rm=T)-max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 2, size = tsize,  hjust = annotate.j)+
                        annotate("text", label = annotate.P2, color=result.color,
                                 x = annotate.x, y = max(dpb$data[[1]]$y, na.rm=T)-max(dpb$data[[1]]$y, na.rm=T)/nlines.ria * 3, size = tsize,  hjust = annotate.j)

                # Add reference lines and Credible interval
                # main.plot <- main.plot +
                #         geom_segment(data=report.segment, aes(x=x, xend=xend, y=y, yend=yend, group=group), arrow = arrow(ends="both", angle = 90), color=result.color, size=.75)

                main.plot <- main.plot +
                        geom_vline(xintercept = c(Delta.lrv, Delta.tv), linetype = 2, color = c("blue", "blue")) +
                        annotate("text", label = TeX(paste0("Min TPP = ", Delta.lrv*100,"%")), x = Delta.lrv, y = 0 + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria, size = tsize, colour = "black", hjust = 1)+
                        annotate("text", label = TeX(paste0("Base TPP = ",  Delta.tv*100,"%")), x = Delta.tv, y = 0 + max(dpb$data[[1]]$y, na.rm=T)/nlines.ria, size = tsize, colour = "black", hjust = 0)
        }

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
                                     data.frame(a.con = a.con, b.con = b.con, a.trt = a.trt, b.trt = b.trt,
                                                Delta.lrv = Delta.lrv, Delta.tv = Delta.tv,
                                                tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng,
                                                Needed.for.NG=result.ng$x.trt[1], Needed.for.GO = result.go$x.trt[1])))
        if(add.table==FALSE) return(main.plot)
}
