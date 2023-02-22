#' @title Make TTE interim treatment OC curve
#'
#' @param my.df output from get.tte.int.OC.df
#' @param nlines number of lines
#' @param tsize text size
#' @param HR.lower lower bound for HR
#' @param HR.upper upper bound for HR
#' @param include_nogo logical
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' make.tte.int.oc()
#' }
make.tte.int.oc <- function(my.df = get.tte.int.oc.df(npoints = 20, include_nogo=TRUE),
                            nlines=25, tsize=6, HR.lower=.0025, HR.upper=2, include_nogo=TRUE){

  my.df$analysis <- factor(my.df$analysis)
  levels(my.df$analysis)
  temp.levels.1 <- levels(my.df$analysis)[1:(length(levels(my.df$analysis))-3)]
  temp.levels.1 <- temp.levels.1[order((as.numeric(temp.levels.1)))]
  temp.levels.2 <- c("any interim", "Study End", "any analysis")
  my.df$analysis <- factor(my.df$analysis, c(temp.levels.1, temp.levels.2))
  levels(my.df$analysis)[(length(levels(my.df$analysis))-2):length(levels(my.df$analysis))] <- c("Any Interim", "Study End", "Any Analysis")



  if(include_nogo==FALSE) {
    my.df <- bind_rows(my.df %>% dplyr::filter(analysis%in% c("Study End", "Any Analysis")),
                       my.df %>% dplyr::filter(!(analysis%in% c("Study End", "Any Analysis"))) %>% dplyr::filter(decision!="No-Go")
                       )

    p3 <- ggplot(data=my.df, aes(x=UHR, y=relFreq, color=decision)) + geom_line(size=.75) + facet_wrap(~analysis) +
      scale_color_manual(values = c("grey", "lightgreen", "red"))+
      labs(x="Underlying Hazard Ratio", y="Probability", color="Decision", title="Treatment Effect Operating Curves by analysis")

    study.end <- my.df %>% dplyr::filter(analysis == "Study End") %>%
      dplyr::filter(decision != "Consider") %>%
      pivot_wider(names_from=decision, values_from=relFreq)
    colnames(study.end)[16] <- "NoGo"

    any.interim <- my.df %>% dplyr::filter(analysis == "Any Interim") %>%
      dplyr::filter(decision != "Consider") %>%
      pivot_wider(names_from=decision, values_from=relFreq)
    any.interim$NoGo <- 1


    any.analysis <- my.df %>% dplyr::filter(analysis == "Any Analysis") %>%
      dplyr::filter(decision != "Consider") %>%
      pivot_wider(names_from=decision, values_from=relFreq)
    colnames(any.analysis)[16] <- "NoGo"

    main.plot <- ggplot() +
      geom_line(data=study.end %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y= Go), color="lightgreen", size=1) +
      geom_line(data=study.end %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, y= 1 - NoGo), color="red", size=1)+
      geom_line(data=any.analysis %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y= Go), color="lightgreen", linetype=2, size=1) +
      geom_line(data=any.analysis %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, y=1 - NoGo), color="red", linetype=2, size=1)
    dpb <- ggplot_build(main.plot)


    for.plot <- left_join(study.end %>% dplyr::filter(is.na(Go)==T) %>% dplyr::select(-c(Go, n)),
    study.end %>% dplyr::filter(is.na(Go)==F)%>% dplyr::select(-c(NoGo,n)))

    p1 <- main.plot + geom_line(data=study.end %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y=Go),
                                color="lightgreen") +
      geom_line(data=for.plot %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y=1 - NoGo), color="red") +
      geom_ribbon(data=for.plot %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, ymin=0, ymax=Go),
                  fill="lightgreen", alpha=.5)+
      geom_ribbon(data=for.plot, aes(x=UHR, ymin=Go, ymax=1-NoGo),
                  fill="grey", alpha=.5)+
      geom_ribbon(data=study.end %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, ymin=1-NoGo, ymax=1),
                  fill="red", alpha=.5)+
      labs(title="Operating characteristics as a function of treatment effect",
           subtitle=paste0("Total events: ", study.end$m.obs[1],
                           ". Randomization ratio (Control:Treatment): (1:", study.end$ARatio[1], ")"))+
      annotate("text", label = TeX(paste0("$\\HR_{Max}$ =", study.end$HR.lrv[1])), x = study.end$HR.lrv[1],
               y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size =
                 tsize, colour = "black", hjust = 0)+
      annotate("text", label = TeX(paste("$\\HR_{Base}$ =", study.end$HR.tv[1])), x = study.end$HR.tv[1],
               y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size = tsize,
               colour = "black", hjust = 1)+
      labs(x = "Underlying Hazard Ratio",
           y="Probability")+
      geom_vline(xintercept=(c(study.end$HR.tv[1], study.end$HR.lrv[1])), color="blue",linetype=2)+
      scale_x_continuous(expand = c(0,0), breaks=seq(0,2,.2),
                         limits=c(HR.lower, HR.upper))+
      scale_y_continuous(expand=c(0,0), breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                         labels=scales::percent)+
      theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
              element_text(angle=45, hjust=1,vjust=1))

    working <- my.df %>% pivot_wider(names_from = decision,  values_from=relFreq) %>% rename(NoGo = 'No-Go')

    p2<- ggplot() + geom_line(data = working %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y=Go), color="darkgreen")+
      facet_wrap(~analysis) +
      geom_ribbon(data = working %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, ymin=0, ymax=Go),  fill="lightgreen", alpha=.5)+
      geom_line(data = working %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, y=1-NoGo), color="red")+
      geom_ribbon(data = working %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, ymin=1-NoGo, ymax=1),  fill="red", alpha=.5)+

      labs(x="Underlying Hazard Ratio", y="Probability", color="Decision", title="Treatment Effect Operating Curves by Analysis")+
      geom_vline(xintercept=(c(study.end$HR.tv[1], study.end$HR.lrv[1])), color="blue",linetype=2)
    }


  if(include_nogo==T){
  p1 <- ggplot(data=my.df, aes(x=UHR, y=relFreq, color=decision)) + geom_line(size=.75) + facet_wrap(~analysis) +
    scale_color_manual(values = c("grey", "lightgreen", "red"))+
    labs(x="Underlying Hazard Ratio", y="Probability", color="Decision", title="Treatment Effect Operating Curves by analysis")

  study.end <- my.df %>% dplyr::filter(analysis == "Study End") %>%
    dplyr::filter(decision != "Consider") %>%
    pivot_wider(names_from=decision, values_from=relFreq)
  colnames(study.end)[16] <- "NoGo"

  any.interim <- my.df %>% dplyr::filter(analysis == "Any Interim") %>%
    dplyr::filter(decision != "Consider") %>%
    pivot_wider(names_from=decision, values_from=relFreq)
  colnames(any.interim)[16] <- "NoGo"

  any.analysis <- my.df %>% dplyr::filter(analysis == "Any Analysis") %>%
    dplyr::filter(decision != "Consider") %>%
    pivot_wider(names_from=decision, values_from=relFreq)
  colnames(any.analysis)[16] <- "NoGo"

  main.plot <- ggplot() +
    geom_line(data=study.end %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y= Go), color="lightgreen", size=1) +
    geom_line(data=study.end %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, y= 1 - NoGo), color="red", size=1)+
    geom_line(data=any.analysis %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y= Go), color="lightgreen", linetype=2, size=1) +
    geom_line(data=any.analysis %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, y=1 - NoGo), color="red", linetype=2, size=1)
  dpb <- ggplot_build(main.plot)

  p2 <- main.plot + geom_line(data=study.end %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y=Go),
                              color="lightgreen") +
    geom_line(data=study.end %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y=1 - NoGo), color="red") +
    geom_ribbon(data=study.end %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, ymin=0, ymax=Go),
                fill="lightgreen", alpha=.5)+
    geom_ribbon(data=study.end, aes(x=UHR, ymin=Go, ymax=1-NoGo),
                fill="grey", alpha=.5)+
    geom_ribbon(data=study.end %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, ymin=1-NoGo, ymax=1),
                fill="red", alpha=.5)+
    labs(title="Operating characteristics as a function of treatment effect",
         subtitle=paste0("Total events: ", study.end$m.obs[1],
                         ". Randomization ratio (Control:Treatment): (1:", study.end$ARatio[1], ")"))+
    annotate("text", label = TeX(paste0("$\\HR_{Max}$ =", study.end$HR.lrv[1])), x = study.end$HR.lrv[1],
             y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size =
               tsize, colour = "black", hjust = 0)+
    annotate("text", label = TeX(paste("$\\HR_{Base}$ =", study.end$HR.tv[1])), x = study.end$HR.tv[1],
             y = 0 + 2*max(dpb$data[[1]]$y)/nlines, size = tsize,
             colour = "black", hjust = 1)+
    labs(x = "Underlying Hazard Ratio",
         y="Probability")+
    geom_vline(xintercept=(c(study.end$HR.tv[1], study.end$HR.lrv[1])), color="blue",linetype=2)+
    scale_x_continuous(expand = c(0,0), breaks=seq(0,2,.2),
                       limits=c(HR.lower, HR.upper))+
    scale_y_continuous(expand=c(0,0), breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                       labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
            element_text(angle=45, hjust=1,vjust=1))

  working <- my.df %>% pivot_wider(names_from = decision,  values_from=relFreq) %>% rename(NoGo = 'No-Go')


  p3<- ggplot() + geom_line(data = working %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, y=Go), color="darkgreen")+
    facet_wrap(~analysis) +
    geom_ribbon(data = working %>% dplyr::filter(is.na(Go) == FALSE), aes(x=UHR, ymin=0, ymax=Go),  fill="lightgreen", alpha=.5)+
    geom_line(data = working %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, y=1-NoGo), color="red")+
    geom_ribbon(data = working %>% dplyr::filter(is.na(NoGo) == FALSE), aes(x=UHR, ymin=1-NoGo, ymax=1),  fill="red", alpha=.5)+

    labs(x="Underlying Hazard Ratio", y="Probability", color="Decision", title="Treatment Effect Operating Curves by Analysis")+
    geom_vline(xintercept=(c(study.end$HR.tv[1], study.end$HR.lrv[1])), color="blue",linetype=2)

  }

 list(p1,p2, p3)
  }
