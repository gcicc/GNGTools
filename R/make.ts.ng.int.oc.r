#' Make two-sample normal-gamma interim OC plot
#'
#' @param for.plot call to get.ts.ng.trt.int.oc.df
#' @param nlines number of lines
#' @param tsize text size
#' @param include_nogo logical
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.ng.trt.int.oc.df <- get.ts.ng.trt.int.oc.df(npointsLookup = 2, npoints=3, n.MC.lookup=5,
#' n.MC=5, go.parallel=FALSE)
#' make.ts.ng.int.oc(for.plot = my.ts.ng.trt.int.oc.df)
#' }
make.ts.ng.int.oc <- function(for.plot = get.ts.ng.trt.int.oc.df(go.parallel = FALSE,
                                                                 include_nogo=TRUE),
                              nlines=25, tsize=6, include_nogo=TRUE){

  my.df.wide <- for.plot  %>% pivot_wider(names_from = decision, values_from = relFreq) %>% rename(NoGo = 'No-Go')
  my.df.wide$analysis <- factor(my.df.wide$analysis)
  t.levels2 <- c("Any Interim", "Study-end", "Any Analysis")
  t.levels1 <- as.character(sort(as.numeric(levels(my.df.wide$analysis)[1:(length(levels(my.df.wide$analysis))-3)])))
  my.df.wide$analysis <- factor(my.df.wide$analysis, c(t.levels1, t.levels2))
  my.df.wide <- my.df.wide %>% dplyr::filter(analysis != t.levels1[length(t.levels1)])


  for.plot$analysis <- factor(for.plot$analysis, c(t.levels1, t.levels2))
  # Supposed to grab user's choice
  if(length(unique(my.df.wide$mu.c))==3) my.df.wide2 <- my.df.wide %>% ungroup() %>% dplyr::filter(mu.c!=max(mu.c) & mu.c !=min(mu.c)) else my.df.wide2 <- my.df.wide


  focused.plot <- ggplot()+
    geom_line(data=my.df.wide2 %>% dplyr::filter(analysis == "Study-end", is.na(Go) == F), aes(x=mu.t-mu.c, y= Go), color="lightgreen", size=1) +
    geom_line(data=my.df.wide2 %>% dplyr::filter(analysis == "Study-end",is.na(NoGo) == F), aes(x=mu.t-mu.c, y= 1 - NoGo), color="red", size=1)+
    geom_line(data=my.df.wide2 %>% dplyr::filter(analysis == "Any Analysis",is.na(Go) == F), aes(x=mu.t-mu.c, y= Go), color="lightgreen", linetype=2, size=1) +
    geom_line(data=my.df.wide2 %>% dplyr::filter(analysis == "Any Analysis",is.na(NoGo) == F), aes(x=mu.t-mu.c, y=1 - NoGo), color="red", linetype=2, size=1) +
    facet_grid(~mu.c)
  dpb <- ggplot_build(focused.plot)

 p1 <-  focused.plot +
    labs(title="Treatment Effect Operating Characteristics",
         subtitle="Study-end (solid lines), Any analysis (dotted lines)",
         x = "Underlying treatment effect",
         y="Probability")+
    geom_vline(xintercept=(c(my.df.wide$Delta.tv[1], for.plot$Delta.lrv [1])), color="blue",linetype=2)+
    geom_ribbon(data=my.df.wide2 %>% dplyr::filter(analysis == "Study-end", is.na(Go) == F), aes(x=mu.t-mu.c, ymin= 0, ymax=Go), fill="lightgreen", alpha=.5)+
    geom_ribbon(data=my.df.wide2 %>% dplyr::filter(analysis == "Study-end", is.na(NoGo) == F), aes(x=mu.t-mu.c, ymin= 1-NoGo, ymax=1), fill="red", alpha=.5)+
    annotate("text", label = paste0("Min TPP = ", for.plot$Delta.lrv[1]), x = for.plot$Delta.lrv[1],
             y = 0 + max(dpb$data[[2]]$y,dpb$data[[2]]$y)/nlines, size = tsize, colour = "black", hjust = 1)+
    annotate("text", label = paste("Base TPP = ", for.plot$Delta.tv[1]), x = for.plot$Delta.tv[1],
             y = 0 + max(dpb$data[[2]]$y,dpb$data[[2]]$y)/nlines, size = tsize,
             colour = "black", hjust = 0)+
    scale_x_continuous(expand = c(0,0)
                       #, breaks=seq(0,2,.2), limits=c(HR.lower, HR.upper)
    )+
    scale_y_continuous(expand=c(0,0),
                       #, breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                       labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
            element_text(angle=45, hjust=1,vjust=1))


  main.plot <- ggplot() +
    geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == F), aes(x=mu.t-mu.c, y= Go), color="lightgreen", size=1) +
    geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end",is.na(NoGo) == F), aes(x=mu.t-mu.c, y= 1 - NoGo), color="red", size=1)+
    geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(Go) == F), aes(x=mu.t-mu.c, y= Go), color="lightgreen", linetype=2, size=1) +
    geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(NoGo) == F), aes(x=mu.t-mu.c, y=1 - NoGo), color="red", linetype=2, size=1) +
    facet_grid(~mu.c, scales="free_x")
  dpb <- ggplot_build(main.plot)


  p2 <- main.plot + geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == F), aes(x=mu.t-mu.c, y=Go),
                              color="lightgreen") +
    geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == F), aes(x=mu.t-mu.c, y=1 - NoGo), color="red") +
    geom_ribbon(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(Go) == F), aes(x=mu.t-mu.c, ymin=0, ymax=Go),
                fill="lightgreen", alpha=.5)+
    geom_ribbon(data=my.df.wide %>% dplyr::filter(analysis == "Study-end", is.na(NoGo) == F), aes(x=mu.t-mu.c, ymin=1-NoGo, ymax=1),
                fill="red", alpha=.5)+
    geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(Go) == F), aes(x=mu.t-mu.c, y= Go), color="lightgreen", linetype=2, size=1) +
    geom_line(data=my.df.wide %>% dplyr::filter(analysis == "Any Analysis",is.na(NoGo) == F), aes(x=mu.t-mu.c, y=1 - NoGo), color="red", linetype=2, size=1)+
    labs(title="Treatment Effect Operating Characteristics",
         subtitle="Study-end (solid lines), Any analysis (dotted lines)",
         x = "Underlying treatment effect",
         y="Probability")+
    geom_vline(xintercept=(c(my.df.wide$Delta.tv[1], my.df.wide$Delta.lrv [1])), color="blue",linetype=2)+
    scale_x_continuous(expand = c(0,0)
                       #, breaks=seq(0,2,.2), limits=c(HR.lower, HR.upper)
    )+
    scale_y_continuous(expand=c(0,0),
                       #, breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                       labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
            element_text(angle=45, hjust=1,vjust=1))


  p4<- ggplot() +
    geom_line(data = my.df.wide %>% dplyr::filter(is.na(Go) == F), aes(x=mu.t-mu.c, y=Go), color="darkgreen")+
    geom_ribbon(data = my.df.wide %>% dplyr::filter(is.na(Go) == F), aes(x=mu.t-mu.c, ymin=0, ymax=Go),  fill="lightgreen", alpha=.5)+
    geom_line(data = my.df.wide %>% dplyr::filter(is.na(NoGo) == F), aes(x=mu.t-mu.c, y=1-NoGo), color="red")+
    geom_ribbon(data = my.df.wide %>% dplyr::filter(is.na(NoGo) == F), aes(x=mu.t-mu.c, ymin=1-NoGo, ymax=1),  fill="red", alpha=.5)+
    labs(x="Underlying Treatment Effect", y="Probability", color="Decision", title="Treatment Effect Operating Characteristics Curves", subtitle="By analysis and assumed control mean")+
    facet_grid(~mu.c~ analysis)+
    geom_vline(xintercept=(c(my.df.wide$Delta.tv[1], my.df.wide$Delta.lrv [1])), color="blue",linetype=2) +
    scale_x_continuous(expand = c(0,0) #, breaks=seq(0,2,.2), limits=c(HR.lower, HR.upper)
    )+

    scale_y_continuous(expand=c(0,0),
                       #, breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                       labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
            element_text(angle=45, hjust=1,vjust=1))


  p3 <- for.plot %>% dplyr::filter(analysis !=  t.levels1[length(t.levels1)]) %>% ggplot(aes(x=mu.t-mu.c, y=relFreq, color=factor(decision))) +
    scale_color_manual(values = c("grey", "lightgreen", "red"))+
    geom_line(size=1) + facet_grid(mu.c~analysis) + labs(color="Decision", x="Underlying Treatment Effect", y="Probability",
                                                         title="Treatment Effect Operating Characteristics Curves", subtitle="By analysis and assumed control mean")+
    scale_y_continuous(
                       #, breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                       labels=scales::percent)+theme(legend.position = "bottom")

  if(include_nogo==F) {
    my.df.wide <- for.plot  %>% pivot_wider(names_from = decision, values_from = relFreq) %>% rename(NoGo = 'No-Go')
    my.df.wide$analysis <- factor(my.df.wide$analysis)
    t.levels2 <- c("Any Interim", "Study-end", "Any Analysis")
    t.levels1 <- as.character(sort(as.numeric(levels(my.df.wide$analysis)[1:(length(levels(my.df.wide$analysis))-3)])))
    my.df.wide$analysis <- factor(my.df.wide$analysis, c(t.levels1, t.levels2))
    my.df.wide <- my.df.wide %>% dplyr::filter(analysis != t.levels1[length(t.levels1)])


    for.plot$analysis <- factor(for.plot$analysis, c(t.levels1, t.levels2))
    # Supposed to grab user's choice
    if(length(unique(my.df.wide$mu.c))==3) my.df.wide2 <- my.df.wide %>% ungroup() %>% dplyr::filter(mu.c!=max(mu.c) & mu.c !=min(mu.c)) else my.df.wide2 <- my.df.wide



    for.plot <- bind_rows(
      for.plot %>% dplyr::filter(analysis %in% c("Study-end",    "Any Analysis")),
      for.plot %>% dplyr::filter(!(analysis %in% c("Study-end",    "Any Analysis"))) %>% dplyr::filter(decision != "No-Go")
    )
  p3 <- ggplot(for.plot %>% dplyr::filter(analysis !=  t.levels1[length(t.levels1)]), aes(x=mu.t-mu.c, y=relFreq, color=factor(decision))) +
    scale_color_manual(values = c("grey", "lightgreen", "red"))+
    geom_line(size=1) + facet_grid(mu.c~analysis) + labs(color="Decision", x="Underlying Treatment Effect", y="Probability",
                                                         title="Treatment Effect Operating Characteristics Curves", subtitle="By analysis and assumed control mean")+
    scale_y_continuous(
      #, breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
      labels=scales::percent)+theme(legend.position = "bottom")

  p4<- ggplot() +
    geom_line(data = my.df.wide %>% dplyr::filter(is.na(Go) == F), aes(x=mu.t-mu.c, y=Go), color="darkgreen")+
    facet_grid(~mu.c~ analysis)+
    geom_ribbon(data = my.df.wide %>% dplyr::filter(is.na(Go) == F), aes(x=mu.t-mu.c, ymin=0, ymax=Go),  fill="lightgreen", alpha=.5)+
    geom_line(data = my.df.wide %>% dplyr::filter(is.na(NoGo) == F), aes(x=mu.t-mu.c, y=1-NoGo), color="red")+
    geom_ribbon(data = my.df.wide %>% dplyr::filter(is.na(NoGo) == F), aes(x=mu.t-mu.c, ymin=1-NoGo, ymax=1),  fill="red", alpha=.5)+
    labs(x="Underlying Treatment Effect", y="Probability", color="Decision", title="Treatment Effect Operating Characteristics Curves", subtitle="By analysis and assumed control mean")+

    geom_vline(xintercept=(c(my.df.wide$Delta.tv[1], my.df.wide$Delta.lrv [1])), color="blue",linetype=2) +
    scale_x_continuous(expand = c(0,0) #, breaks=seq(0,2,.2), limits=c(HR.lower, HR.upper)
    )+

    scale_y_continuous(expand=c(0,0),
                       #, breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1),
                       labels=scales::percent)+
    theme(panel.spacing.x = unit(6, "mm"), axis.text.x =
            element_text(angle=45, hjust=1,vjust=1))


}

    return(list(p1,p2, p3, p4))

}
