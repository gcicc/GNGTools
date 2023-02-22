#' @title Make two-sample binary interim oc plot
#'
#' @param ts.bin.int.oc Results from get.ts.bin.int.oc
#' @param include_nogo logical
#' @param nlines Control for text spacing
#' @param tsize Control for text size
#'
#' @return A ggplot showing the probability of each out come given a particular difference from control
#' @export
#'
make.ts.bin.int.oc.plot <- function(ts.bin.int.oc=get.ts.bin.int.oc(
        a.con = 1, b.con = 1,
        a.trt = 1, b.trt = 1,
        Delta.tv = .3, Delta.lrv = .2,
        tau.tv = .1, tau.lrv = .8, tau.ng = .65,
        go.thresh = .8, ng.thresh = .8,
        n.con = 40, n.trt = 40,
        n.int.c = c(10, 20, 30),
        n.int.t = c(10, 20, 30),
        DecisionTable= NULL,
        runs=5000,
        ControlRate=.2,
        TreatmentEffect=seq(0, 0.8,.1),
        include_nogo=TRUE),
        nlines=25, tsize=4, include_nogo=TRUE){

        Final <- max(ts.bin.int.oc$assessment)
        runs <- max(ts.bin.int.oc$run)

        # final results
        study.end <- ts.bin.int.oc %>% dplyr::filter(assessment==Final) %>%
                dplyr::group_by(Effect) %>%
                dplyr::count(decision, .drop=FALSE) %>%
                dplyr::mutate(relFreq = n/runs, analysis="Study End")%>%
                dplyr::ungroup()

        any.interim <- ts.bin.int.oc %>%
                dplyr::filter(assessment!=Final) %>%
                dplyr::group_by(Effect, run)  %>%
                dplyr::summarise(decision = dplyr::case_when(
                        all(decision == 'Consider') ~ 'Consider',
                        T ~ as.character(decision[ decision != 'Consider'][1])))%>%
                dplyr::count(decision, .drop=FALSE) %>%
                dplyr::mutate(relFreq = n/runs, analysis = "Any Interim")%>%
                dplyr::ungroup()

        interims <- ts.bin.int.oc %>% dplyr::filter(assessment!=Final) %>%
                dplyr::group_by(assessment, Effect, n.int.c, n.int.t) %>%
                dplyr::count(decision, .drop=FALSE) %>% dplyr::mutate(relFreq = n/runs)%>%
                dplyr::ungroup() %>%
                dplyr::mutate(analysis=as.character(n.int.c[[1]][assessment]),
                       n.c=n.int.c[[1]][assessment],
                       n.t=n.int.t[[1]][assessment]) %>% dplyr::select(-assessment)

        any.analysis <- ts.bin.int.oc %>%
                dplyr::group_by(Effect,run) %>%
                dplyr::summarise(interims.d = dplyr::case_when(
                        all(decision == 'Consider') ~ 'Consider',
                        TRUE ~ as.character(decision[ decision != 'Consider'][1])),
                        final.d = decision[assessment == Final]) %>% mutate(decision = case_when(
                                # Inteirm Go trumps final No-Go
                                interims.d == "Go" & final.d == "Go" ~ "Go",
                                interims.d == "Go" & final.d == "NoGo" ~ "Go",
                                interims.d == "Go" & final.d =="Consider" ~ "Go",
                                # interim No-go trumps final Go
                                interims.d == "NoGo" & final.d == "Go" ~ "NoGo",
                                interims.d == "NoGo" & final.d == "NoGo" ~ "NoGo",
                                interims.d == "NoGo" & final.d == "Consider" ~ "NoGo",
                                # Final Go/No-Go trumps interim consider
                                interims.d == "Consider" & final.d == "Go" ~ "Go",
                                interims.d == "Consider" & final.d =="NoGo" ~ "NoGo",
                                interims.d == "Consider" & final.d =="Consider" ~ "Consider"
                        )) %>%
                dplyr::group_by(Effect) %>%
                dplyr::count(decision, .drop=FALSE) %>%
                dplyr::mutate(relFreq = n/sum(n),
                       analysis = "Any Analysis") %>%
                dplyr::ungroup()

        all <- bind_rows(interims %>% dplyr::select(Effect, decision ,    n ,relFreq ,analysis),
                         any.interim,
                         study.end,
                         any.analysis)  %>%
                pivot_wider(names_from = decision, values_from = relFreq)%>% dplyr::rename(NoGo = 'NoGo')
        all$analysis <- factor(all$analysis)
        t.levels2 <- c("Any Interim", "Study End", "Any Analysis")
        t.levels1 <- as.character(sort(as.numeric(levels(all$analysis)[1:(length(levels(all$analysis))-3)])))
        all$analysis <- factor(all$analysis, c(t.levels1, t.levels2))

        all.long <- bind_rows(interims, any.interim, study.end, any.analysis)
        all.long$analysis <- factor(all.long$analysis, c(t.levels1, t.levels2))

        if(include_nogo == FALSE) all.long <- rbind(all.long%>% dplyr::filter(analysis %in% c("Study End")),
                                                all.long%>% dplyr::filter(!(analysis %in% c("Study End"))) %>%
                                                        dplyr::filter(decision!="No-Go"))




        p1 <- ggplot(all.long, aes(x=Effect, y=relFreq, color=factor(decision))) +
                scale_color_manual(values = c("grey", "lightgreen", "red", "black"))+
                geom_line(size=1) + facet_wrap(~analysis) +
                labs(color="Decision", y="Probability",
                     x = "Underlying treatment effect",
                     title="Treatment effect operating characteristics by analysis")+
                scale_y_continuous(expand=c(0,0), labels=scales::percent, limits=c(0,1))


        all.collpased <- all %>% dplyr::filter(analysis == "Study End") %>% group_by(analysis, Effect) %>%
                summarize(Go = Go[1], Consider=Consider[2], NoGo=NoGo[3])

        main.plot <- ggplot() +
                geom_line(data=all %>% dplyr::filter(analysis == "Study End", is.na(Go)== FALSE), aes(x=Effect, y= Go), color="lightgreen", size=1) +
                geom_line(data=all %>% dplyr::filter(analysis == "Study End",is.na(NoGo)== FALSE), aes(x=Effect, y= 1 - NoGo), color="red", size=1)+
                geom_line(data=all %>% dplyr::filter(analysis == "Any Analysis",is.na(Go)== FALSE), aes(x=Effect, y= Go), color="lightgreen", linetype=2, size=1)
        if(include_nogo==TRUE) main.plot <- main.plot + geom_line(data=all %>% dplyr::filter(analysis == "Any Analysis",is.na(NoGo)== FALSE),
                                                                  aes(x=Effect, y=1 - NoGo), color="red", linetype=2, size=1)


        dpb <- ggplot_build(main.plot)

        p2 <- main.plot + geom_line(data=all %>% dplyr::filter(analysis == "Study End", is.na(Go)== FALSE), aes(x=Effect, y=Go),
                                    color="lightgreen") +
                geom_line(data=all %>% dplyr::filter(analysis == "Study End", is.na(Go)== FALSE), aes(x=Effect, y=1 - NoGo), color="red") +
                geom_ribbon(data=all %>% dplyr::filter(analysis == "Study End", is.na(Go)== FALSE), aes(x=Effect, ymin=0, ymax=Go),
                            fill="lightgreen", alpha=.5)+
                geom_ribbon(data=all %>% dplyr::filter(analysis == "Study End", is.na(NoGo)== FALSE), aes(x=Effect, ymin=1-NoGo, ymax=1),
                            fill="red", alpha=.5)+
                geom_line(data=all %>% dplyr::filter(analysis == "Any Analysis",is.na(Go)== FALSE), aes(x=Effect, y= Go), color="lightgreen", linetype=2, size=1) +

                labs(title="Operating characteristics as a function of treatment effect",
                     subtitle="Solid [dotted] lines provide operating characteristics for study-end [any analysis].",
                     x = "Underlying treatment effect",
                     y="Probability")+
                geom_vline(xintercept=(c(ts.bin.int.oc$Delta.tv[1], ts.bin.int.oc$Delta.lrv[1])), color="blue",linetype=2)+
                annotate("text", label = TeX(paste0("Min TPP = ", ts.bin.int.oc$Delta.lrv[1]*100, "%")),
                         x = ts.bin.int.oc$Delta.lrv[1], y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize,
                         colour = "black", hjust = 1)+
                annotate("text", label = TeX(paste("Base TPP = ", ts.bin.int.oc$Delta.tv[1]*100, "%")),
                         x = ts.bin.int.oc$Delta.tv[1], y = 0 + max(dpb$data[[1]]$y)/nlines, size = tsize,
                         colour = "black", hjust = 0)+
                scale_x_continuous(expand = c(0,0)) +
                scale_y_continuous(expand=c(0,0), labels=scales::percent, limits=c(0,1))+
                theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
                              element_text(angle=45, hjust=1,vjust=1))
        if(include_nogo==TRUE) p2 <- p2+   geom_line(data=all %>% dplyr::filter(analysis == "Any Analysis",is.na(NoGo)== FALSE),
                                                     aes(x=Effect, y=1 - NoGo), color="red", linetype=2, size=1)


        p3<- ggplot() + geom_line(data = all %>% dplyr::filter(is.na(Go)== FALSE), aes(x=Effect, y=Go), color="darkgreen")+
                geom_ribbon(data = all %>% dplyr::filter(is.na(Go)== FALSE), aes(x=Effect, ymin=0, ymax=Go),  fill="lightgreen", alpha=.5)+
                facet_wrap(~analysis) +
                labs(x = "Underlying treatment effect",
                     y="Probability", color="Decision", title="Treatment Effect Operating Curves by Analysis")+
                geom_vline(xintercept=(c(ts.bin.int.oc$Delta.lrv[1], ts.bin.int.oc$Delta.tv[1])), color="blue",linetype=2) +
                scale_y_continuous(expand=c(0,0), labels=scales::percent, limits=c(0,1))
        if(include_nogo==FALSE) p3 <- p3 + geom_ribbon(data = all %>% dplyr::filter(is.na(NoGo)== FALSE & analysis=="Study End"), aes(x=Effect, ymin=1-NoGo, ymax=1),
                                                   fill="red", alpha=.5)

        if(include_nogo==TRUE) p3 <- p3+ geom_ribbon(data = all %>% dplyr::filter(is.na(NoGo)== FALSE), aes(x=Effect, ymin=1-NoGo, ymax=1),
                                                  fill="red", alpha=.5)


        list(p1,p2, p3)

}


