#' @title Make single sample binary interim OC Split plot
#'
#' @param ss.bin.int.GNG call to get.ss.bin.interim.GNG
#' @param ss.bin.int.oc call to get.ss.bin.interim.oc
#' @param ss.bin.studyend.GNG call to get.ss.bin.studyend.GNG
#' @param ss.bin.int.df call to get.ss.bin.int.df
#' @param lower lower bound
#' @param upper upper bound
#' @param step stepsize
#' @param include_nogo logical
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' my.ss.bin.studyend.GNG = get.ss.bin.studyend.GNG(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
#'                                       Delta.lrv = .2, Delta.tv = .35,
#'                                       tau.tv = 0.10, tau.lrv = .80, tau.ng = .65)
#' my.ss.bin.int.df <- get.ss.bin.int.df(ss.bin.studyend.GNG = my.ss.bin.studyend.GNG,
#'                                       goThreshold = .8, nogoThreshold = .8, include_nogo =TRUE)
#' my.ss.bin.int.GNG <- get.ss.bin.int.GNG(ss.bin.int.df = my.ss.bin.int.df,
#'                                     Interims = 20,
#'                                     ss.bin.studyend.GNG = my.ss.bin.studyend.GNG)
#'
#'  my.ss.bin.int.oc <- get.ss.bin.int.oc(
#'  ss.bin.int.df = my.ss.bin.int.df,
#'  ss.bin.int.GNG=my.ss.bin.int.GNG)
#'  my.ss.bin.int.oc.plot <- make.ss.bin.int.oc.plotSplit(
#'                  ss.bin.studyend.GNG = my.ss.bin.studyend.GNG,
#'                  ss.bin.int.oc = my.ss.bin.int.oc,
#'                  ss.bin.int.GNG=my.ss.bin.int.GNG,
#'                  ss.bin.int.df = my.ss.bin.int.df,
#'                  include_nogo =TRUE, lower=0, upper=1)
#' }


make.ss.bin.int.oc.plotSplit <- function(ss.bin.studyend.GNG,
                                         ss.bin.int.oc,
                                         ss.bin.int.GNG,
                                         ss.bin.int.df,
                                         lower = 0,upper = .75,step = .025,
                                         include_nogo=FALSE){

  if(include_nogo==TRUE){
  lastAssessment = nrow(ss.bin.int.GNG)

  Singles <- expand.grid(targetRate = seq(from= lower,to = upper,by = step), Analysis = 1:lastAssessment) %>%
    dplyr::mutate(Go = 1-pbinom(ss.bin.int.GNG$MinGo[Analysis]-1,ss.bin.int.GNG$n[Analysis],targetRate),
                  `No Go` = pbinom(ss.bin.int.GNG$MaxNoGo[Analysis],ss.bin.int.GNG$n[Analysis],targetRate),
                  Consider = 1-Go-`No Go`) %>%
          tidyr::pivot_longer(c(-targetRate,-Analysis),'Decision',values_to = 'Probability') %>%
          dplyr::mutate(Analysis = paste0(ifelse(Analysis ==lastAssessment,'Final','Interim'),', N = ',ss.bin.int.GNG$n[Analysis])) %>%
          dplyr::group_by(targetRate,Analysis) %>%
          dplyr::mutate(
                  upper = case_when(
                          Decision == 'Go' ~ Probability,
                          Decision == 'No Go' ~ 1,
                          Decision == 'Consider' ~ Probability[Decision == 'Go']+Probability),
                  lower = case_when(
                          Decision == 'Go' ~ 0,
                          Decision == 'No Go' ~ 1-Probability,
                          Decision == 'Consider' ~ Probability[Decision == 'Go']),
                  line = case_when(
                          Decision == 'Go' ~ Probability,
                          Decision == 'No Go' ~ 1-Probability,
                          Decision == 'Consider' ~ Probability)) %>%
          select(-Probability)

    Plt.frm <- purrr::map_df(
            seq(from= lower,to = upper,by = step),function(x) ocTable_multi(DesignTable = ss.bin.int.GNG,TargetRate = x)) %>%
            dplyr::filter(Analysis %in% paste('Analysis',1:(lastAssessment-1))) %>%
            dplyr::select(targetRate, `First Go`,`First No Go`, Consider = `Grey Area`) %>%
            dplyr::group_by(targetRate) %>%
            dplyr::summarise(across(c(`First Go`,`First No Go`, Consider),sum)) %>%
            dplyr::mutate(
                    Go_upper = `First Go`,
                    Go_lower = 0,
                    Go_line = `First Go`,
                    Consider_upper = 1-`First No Go`,
                    Consider_lower = `First Go`,
                    `No Go_upper`  = 1,
                    `No Go_lower` = 1-`First No Go`,
                    `No Go_line` = `No Go_lower`,
                    Consider_line = `No Go_line`-`Go_line`,
                    Analysis = 'Any Interim'
                    )%>%
            dplyr::select(Analysis,targetRate,ends_with(c('upper','lower','line'))) %>%
            tidyr::pivot_longer(ends_with(c('upper','lower','line')),names_to = c("Decision", ".value"),names_pattern = "(.*)_(.*)") %>%
            bind_rows(
                    ss.bin.int.oc %>%
                            dplyr::filter(Analysis == 'Overall') %>%
                            dplyr::mutate(Analysis = 'Any Analysis') %>%
                            dplyr::select(-Type),
                    Singles
                    ) %>%
            dplyr::mutate(
                    Decision = factor(Decision,levels = c('Go','Consider','No Go')),
                    Analysis = factor(Analysis, levels = c(paste('Interim, N =',ss.bin.int.GNG$n[-lastAssessment]),paste('Final, N =',ss.bin.int.GNG$n[lastAssessment]),'Any Interim','Any Analysis'))
                    )


  txt <- tibble(
          targetRate = c(ss.bin.studyend.GNG[[1]]$Delta.lrv,ss.bin.studyend.GNG[[1]]$Delta.tv),
          label = paste(c('Min TPP =','Target TPP ='),targetRate),
          align = c(1.1,-.1),
          Rate = .05
          )

    RibbonPlot <- ggplot(Plt.frm,aes(x=targetRate))+
            geom_ribbon(aes(ymin = lower,ymax=upper,fill=Decision),alpha=.6)+
            facet_wrap(~Analysis)+
            geom_vline(xintercept = txt$targetRate,linetype = 'dashed', color="blue")+
            geom_text(data = txt,color='black', aes(label = label,hjust = align,y=Rate), size=6)+
            scale_x_continuous(expand= c(0,0), limits=c(lower,upper))+
            scale_y_continuous(expand= c(0,0),labels = scales::percent)+
            scale_color_manual(values = c('green','grey','red'),drop = F)+
            scale_fill_manual(values = c('green','grey','red'),drop = F)+
            scale_linetype_manual(values = c(3,1))+
            theme(legend.position = 'bottom')+
            theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
                  panel.spacing.x = unit(2.25,"line"))+
            labs(title="Operating Characteristics as a function of treatment effect by analysis",
                 x="Response Rate",
                 y="Probability",
                 color="Conclusion",
                 linetype="Analysis")

    Plt.frm <- Plt.frm %>%
            mutate(line = ifelse(Decision == 'No Go',1-line,line))

    LinePlot <- ggplot(Plt.frm,aes(x = targetRate,y=line)) +
            geom_line(aes(color = Decision), size=.75)+
            facet_wrap(~Analysis)+
            geom_vline(xintercept = txt$targetRate,linetype = 'dashed', color="blue")+
            geom_text(data = txt,color='black', aes(label = label,hjust = align,y=Rate), size=6)+
            scale_x_continuous(expand= c(0,0), limits=c(lower,upper))+
            scale_y_continuous(expand= c(0,0),labels = scales::percent)+
            scale_color_manual(values = c('green','grey','red'),drop = F)+
            theme(legend.position = 'bottom')+
            theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
                  panel.spacing.x = unit(2.25,"line"))+
            labs(title="Operating Characteristics as a function of treatment effect by analysis",
                 x="Response Rate",
                 y="Probability",
                 color="Conclusion")
    } else {
            lastAssessment = nrow(ss.bin.int.GNG)

            Singles <- crossing(targetRate = seq(from= lower,to = upper,by = step), Analysis = 1:lastAssessment) %>%
                    dplyr::mutate(Go = 1-pbinom(ss.bin.int.GNG$MinGo[Analysis]-1,ss.bin.int.GNG$n[Analysis],targetRate),
                                  `No Go` = pbinom(ss.bin.int.GNG$MaxNoGo[Analysis],ss.bin.int.GNG$n[Analysis],targetRate),
                                  Consider = 1-Go-`No Go`) %>%
                    pivot_longer(c(-targetRate,-Analysis),'Decision',values_to = 'Probability') %>%
                    dplyr::mutate(Analysis = paste0(ifelse(Analysis ==lastAssessment,'Final','Interim'),', N = ',ss.bin.int.GNG$n[Analysis])) %>%
                    dplyr::group_by(targetRate,Analysis) %>%
                    dplyr::mutate(
                            upper = dplyr::case_when(
                                    Decision == 'Go' ~ Probability,
                                    Decision == 'No Go' ~ 1,
                                    Decision == 'Consider' ~ Probability[Decision == 'Go']+Probability),
                            lower = case_when(
                                    Decision == 'Go' ~ 0,
                                    Decision == 'No Go' ~ 1-Probability,
                                    Decision == 'Consider' ~ Probability[Decision == 'Go']),
                            line = case_when(
                                    Decision == 'Go' ~ Probability,
                                    Decision == 'No Go' ~ 1-Probability,
                                    Decision == 'Consider' ~ Probability)) %>%
                    dplyr::select(-Probability)

                Plt.frm <- purrr::map_df(
                        seq(from= lower,to = upper,by = step),function(x) ocTable_multi(DesignTable = ss.bin.int.GNG,x)) %>%
                        dplyr::filter(Analysis %in% paste('Analysis',1:(lastAssessment-1))) %>%
                        dplyr::select(targetRate, `First Go`,`First No Go`, Consider = `Grey Area`) %>%
                        dplyr::group_by(targetRate) %>%
                        dplyr::summarise(across(c(`First Go`,`First No Go`, Consider),sum)) %>%
                        dplyr::mutate(
                                Go_upper = `First Go`,
                                Go_lower = 0,
                                Go_line = `First Go`,
                                Consider_upper = 1-`First No Go`,
                                Consider_lower = `First Go`,
                                `No Go_upper`  = 1,
                                `No Go_lower` = 1-`First No Go`,
                                `No Go_line` = `No Go_lower`,
                                Consider_line = `No Go_line`-`Go_line`,
                                Analysis = 'Any Interim'
                                )%>%
                        dplyr::select(Analysis, targetRate, ends_with(c('upper','lower','line'))) %>%
                        tidyr::pivot_longer(ends_with(c('upper','lower','line')),names_to = c("Decision", ".value"),names_pattern = "(.*)_(.*)") %>%
                        dplyr::bind_rows(
        ss.bin.int.oc %>%
                dplyr::filter(Analysis == 'Overall') %>%
                dplyr::mutate(Analysis = 'Any Analysis') %>%
                dplyr::select(-Type),
        Singles) %>%
                dplyr::mutate(
                        Decision = factor(Decision,levels = c('Go','Consider','No Go')),
                        Analysis = factor(Analysis, levels = c(paste('Interim, N =',ss.bin.int.GNG$n[-lastAssessment]),paste('Final, N =',ss.bin.int.GNG$n[lastAssessment]),'Any Interim','Any Analysis'))
                        )

                txt <- tibble(
                        targetRate = c(ss.bin.studyend.GNG[[1]]$Delta.lrv,ss.bin.studyend.GNG[[1]]$Delta.tv),
                        label = paste(c('Min TPP =','Target TPP ='),targetRate),
                        align = c(1.1,-.1),
                        Rate = .05
                        )

                Plt.frm <- bind_rows(
                        Plt.frm %>% dplyr::filter(Analysis %in% levels(Plt.frm$Analysis)[c(length(levels(Plt.frm$Analysis)), length(levels(Plt.frm$Analysis))-2)]),
                        Plt.frm %>% dplyr::filter(!(Analysis %in% levels(Plt.frm$Analysis)[c(length(levels(Plt.frm$Analysis)), length(levels(Plt.frm$Analysis))-2)])) %>% dplyr::filter(Decision!="No Go")
                        )

                RibbonPlot <- ggplot(Plt.frm,aes(x=targetRate))+
                        geom_ribbon(aes(ymin = lower,ymax=upper,fill=Decision),alpha=.6)+
                        facet_wrap(~Analysis)+
                        geom_vline(xintercept = txt$targetRate,linetype = 'dashed', color="blue")+
                        geom_text(data = txt,color='black', aes(label = label,hjust = align,y=Rate), size=6)+
                        scale_x_continuous(expand= c(0,0), limits=c(lower,upper))+
                        scale_y_continuous(expand= c(0,0),labels = scales::percent)+
                        scale_color_manual(values = c('green','grey','red'),drop = F)+
                        scale_fill_manual(values = c('green','grey','red'),drop = F)+
                        scale_linetype_manual(values = c(3,1))+
                        theme(legend.position = 'bottom')+
                        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
                              panel.spacing.x = unit(2.25,"line"))+
                        labs(title="Operating Characteristics as a function of treatment effect by analysis",
                             x="Response Rate",
                             y="Probability",
                             color="Conclusion",
                             linetype="Analysis")

                Plt.frm <- Plt.frm %>%
                        dplyr::mutate(line = ifelse(Decision == 'No Go',1-line,line))

                LinePlot <- ggplot(Plt.frm,aes(x = targetRate,y=line)) +
                        geom_line(aes(color = Decision), size=.75)+
                        facet_wrap(~Analysis)+
                        geom_vline(xintercept = txt$targetRate,linetype = 'dashed', color="blue")+
                        geom_text(data = txt,color='black', aes(label = label,hjust = align,y=Rate), size=6)+
                        scale_x_continuous(expand= c(0,0), limits=c(lower,upper))+
                        scale_y_continuous(expand= c(0,0),labels = scales::percent)+
                        scale_color_manual(values = c('green','grey','red'),drop = F)+
                        theme(legend.position = 'bottom')+
                        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
                              panel.spacing.x = unit(2.25,"line"))+
                        labs(title="Operating Characteristics as a function of treatment effect by analysis",
                             x="Response Rate",
                             y="Probability",
                             color="Conclusion")
    }

          list(RibbonPlot = RibbonPlot,LinePlot = LinePlot)
}
