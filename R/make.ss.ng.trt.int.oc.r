#' @title Make single sample normal-gamma interim treatment oc curve
#'
#' @param my.df output from get.ss.ng.trt.int.oc.df
#' @param include_nogo logical
#' @return A ggplot object is returned.
#' @export
#'
#' @examples \donttest{
#' my.ss.ng.trt.int.oc.df <- get.ss.ng.trt.int.oc.df(npoints = 20, n.MC = 1000,
#' include_nogo = FALSE)
#' my.ss.ng.trt.int.oc <- make.ss.ng.trt.int.oc(my.df = my.ss.ng.trt.int.oc.df,
#' include_nogo=FALSE)
#' my.ss.ng.trt.int.oc[[1]]
#' my.ss.ng.trt.int.oc[[2]]
#' my.ss.ng.trt.int.oc[[3]]
#' }
make.ss.ng.trt.int.oc <- function(my.df, include_nogo=FALSE){

        my.df$analysis <- factor(my.df$analysis)
        levels(my.df$analysis)[(length(levels(my.df$analysis))-2):length(levels(my.df$analysis))]
        tlevels1 <- as.character(sort(as.numeric(levels(my.df$analysis)[1:(length(levels(my.df$analysis))-3)])))
        tlevels2 <- c("Any Interim", "Study-end", "Any Analysis")
        my.df$analysis <- factor(my.df$analysis, c(tlevels1, tlevels2))

        my.df.wide <- my.df  %>% pivot_wider(names_from = decision, values_from = relFreq) %>% rename(NoGo = 'No-Go')


        my.df.wide$analysis <- factor(my.df.wide$analysis)
        t.levels2 <- c("Any Interim", "Study-end", "Any Analysis")
        t.levels1 <- as.character(sort(as.numeric(levels(my.df.wide$analysis)[1:(length(levels(my.df.wide$analysis))-3)])))
        my.df.wide$analysis <- factor(my.df.wide$analysis, c(t.levels1, t.levels2))


        if (include_nogo==T){
                p1 <- ggplot(my.df, aes(x=delta, y=relFreq, color=factor(decision))) +
                        scale_color_manual(values = c("grey", "lightgreen", "red"))+
                        geom_line(size=1) + facet_wrap(~analysis) + labs(color="Decision", y="Probability", title="Treatment effect operating characteristics by analysis")


                main.plot <- ggplot() +
                        geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE), aes(x=delta, y= Go), color="lightgreen", size=1) +
                        geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end",is.na(NoGo) == FALSE), aes(x=delta, y= 1 - NoGo), color="red", size=1)+
                        geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(Go) == FALSE), aes(x=delta, y= Go), color="lightgreen", linetype=2, size=1) +
                        geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(NoGo) == FALSE), aes(x=delta, y=1 - NoGo), color="red", linetype=2, size=1)
                dpb <- ggplot_build(main.plot)

                p2 <- main.plot + geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE), aes(x=delta, y=Go),
                                            color="lightgreen") +
                        geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE), aes(x=delta, y=1 - NoGo), color="red") +
                        geom_ribbon(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE), aes(x=delta, ymin=0, ymax=Go),
                                    fill="lightgreen", alpha=.5)+
                        # geom_ribbon(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE & is.na(NoGo) == FALSE), aes(x=delta, ymin=Go, ymax=1-NoGo),
                        #             fill="grey", alpha=.5)+
                        geom_ribbon(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(NoGo) == FALSE), aes(x=delta, ymin=1-NoGo, ymax=1),
                                    fill="red", alpha=.5)+
                        geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(Go) == FALSE), aes(x=delta, y= Go), color="lightgreen", linetype=2, size=1) +
                        geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(NoGo) == FALSE), aes(x=delta, y=1 - NoGo), color="red", linetype=2, size=1)+
                        # labs(
                        #      subtitle=paste0("Total events: ", study.end$m.obs[1],
                        #                      ". Randomization ratio (Control:Treatment): (1:", study.end$ARatio[1], ")"))+
                        # annotate("text", label = TeX(paste0("$\\HR_{Max}$ =", study.end$HR.lrv[1])), x = study.end$HR.lrv[1],
                        #          y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size =
                        #            tsize, colour = "black", hjust = 0)+
                        # annotate("text", label = TeX(paste("$\\HR_{Base}$ =", study.end$HR.tv[1])), x = study.end$HR.tv[1],
                        #          y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size = tsize,
                        #          colour = "black", hjust = 1)+
                        labs(title="Operating characteristics as a function of treatment effect",
                             x = "Underlying treatment effect",
                             y="Probability")+
                        geom_vline(xintercept=(c(my.df.wide$Delta.tv[1], my.df.wide$Delta.lrv[1])), color="blue",linetype=2)+
                        scale_x_continuous(expand = c(0,0)
                                           #, breaks=seq(0,2,.2), limits=c(HR.lower, HR.upper)
                        )+
                        scale_y_continuous(expand=c(0,0),
                                           #, breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                                           labels=scales::percent)+
                        theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
                                      element_text(angle=45, hjust=1,vjust=1))



                p3<- ggplot() + geom_line(data = my.df.wide %>% dplyr::filter(is.na(Go) == FALSE), aes(x=delta, y=Go), color="darkgreen")+
                        geom_ribbon(data = my.df.wide %>% dplyr::filter(is.na(Go) == FALSE), aes(x=delta, ymin=0, ymax=Go),  fill="lightgreen", alpha=.5)+
                        geom_line(data = my.df.wide %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=delta, y=1-NoGo), color="red")+
                        geom_ribbon(data = my.df.wide %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=delta, ymin=1-NoGo, ymax=1),  fill="red", alpha=.5)+
                        facet_wrap(~analysis) +
                        labs(x="Underlying treatment effect", y="Probability", color="Decision", title="Treatment Effect Operating Curves by Analysis")+
                        geom_vline(xintercept=(c(my.df.wide$Delta.tv[1], my.df.wide$Delta.lrv[1])), color="blue",linetype=2)

        } else {
                {
                        my.df <- bind_rows(
                                my.df %>% dplyr::filter(analysis %in% c("Study-end", "Any Analysis")),
                                my.df %>% dplyr::filter(!(analysis %in% c("Study-end", "Any Analysis"))) %>% dplyr::filter(decision !="No-Go"))

                        p1 <- ggplot(my.df , aes(x=delta, y=relFreq, color=factor(decision))) +
                                scale_color_manual(values = c("grey", "lightgreen", "red"))+
                                geom_line(size=1) + facet_wrap(~analysis) + labs(color="Decision", y="Probability", title="Treatment effect operating characteristics by analysis")+
                                theme(legend.position = "bottom")

                        my.df.wide <- my.df  %>% pivot_wider(names_from = decision, values_from = relFreq) %>% rename(NoGo = 'No-Go')

                        main.plot <- ggplot() +
                                geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE), aes(x=delta, y= Go), color="lightgreen", size=1) +
                                geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end",is.na(NoGo) == FALSE), aes(x=delta, y= 1 - NoGo), color="red", size=1)+
                                geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(Go) == FALSE), aes(x=delta, y= Go), color="darkgreen", linetype=2, size=1)

                        dpb <- ggplot_build(main.plot)

                        p2 <- main.plot + geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE), aes(x=delta, y=Go),
                                                    color="lightgreen") +
                                geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE), aes(x=delta, y=1 - NoGo), color="red") +
                                geom_ribbon(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == FALSE), aes(x=delta, ymin=0, ymax=Go),
                                            fill="lightgreen", alpha=.5)+
                                geom_ribbon(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(NoGo) == FALSE), aes(x=delta, ymin=1-NoGo, ymax=1),
                                            fill="red", alpha=.5)+
                                geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(Go) == FALSE), aes(x=delta, y= Go), color="darkgreen", linetype=2, size=1) +
                                geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(NoGo) == FALSE), aes(x=delta, y=1 - NoGo), color="red", linetype=2, size=1)+
                                labs(title="Operating characteristics as a function of treatment effect",
                                     x = "Underlying treatment effect",
                                     y="Probability")+
                                geom_vline(xintercept=(c(my.df.wide$Delta.tv[1], my.df.wide$Delta.lrv[1])), color="blue",linetype=2)+
                                scale_x_continuous(expand = c(0,0))+
                                scale_y_continuous(expand=c(0,0), labels=scales::percent)+
                                theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
                                              element_text(angle=45, hjust=1,vjust=1))



                        p3<- ggplot() + geom_line(data = my.df.wide %>% dplyr::filter(is.na(Go) == FALSE), aes(x=delta, y=Go), color="darkgreen")+
                                geom_ribbon(data = my.df.wide %>% dplyr::filter(is.na(Go) == FALSE), aes(x=delta, ymin=0, ymax=Go),  fill="lightgreen", alpha=.5)+
                                geom_line(data = my.df.wide %>% dplyr::filter(is.na(NoGo) == FALSE) %>%
                                                  mutate(NoGo=ifelse(NoGo==0 & !(analysis %in% c("Study-end", "Any Analysis")), NA, NoGo)
                                                  ), aes(x=delta, y=1-NoGo), color="red")+
                                geom_ribbon(data = my.df.wide %>% dplyr::filter(is.na(NoGo) == FALSE) %>%
                                                    mutate(NoGo=ifelse(NoGo==0 & !(analysis %in% c("Study-end", "Any Analysis")), NA, NoGo)
                                                    ), aes(x=delta, ymin=1-NoGo, ymax=1),  fill="red", alpha=.5)+
                                facet_wrap(~analysis) +
                                labs(x="Underlying treatment effect", y="Probability", color="Decision", title="Treatment Effect Operating Curves by Analysis")+
                                geom_vline(xintercept=(c(my.df.wide$Delta.tv[1], my.df.wide$Delta.lrv[1])), color="blue",linetype=2)


                }

        }



        list(p1,p2, p3)




}
