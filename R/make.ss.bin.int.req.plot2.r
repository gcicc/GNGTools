#' @title Make single sample binary interim Group-sequential style plot
#' @param InterimDF call to get.ss.bin.interim.df
#' @param Delta.lrv min TPP
#' @param Delta.tv Base TPP
#' @param tsize text size - default is for app
#'
#' @return A ggplot object is returned
#' @export

make.ss.bin.int.req.plot2 <- function(InterimDF,
                                      Delta.lrv = .2, Delta.tv = .3, tsize=5){

  Plt.frm <- InterimDF %>%
    group_by(n) %>%
    summarise(
      minGo = suppressWarnings(min(x[Decision =='Go'])),
      MaxNoGo = suppressWarnings(max(x[Decision =='No Go'])),
      Go = minGo/n,
      `No Go` = MaxNoGo/n
      ) %>%
    select(n,Go,`No Go`) %>%
    pivot_longer(-n, names_to = 'Decision',values_to = 'Observed Proportion') %>%
    mutate(Decision = factor(Decision, levels = c('Go','No Go'))) %>%
    filter(!is.infinite(`Observed Proportion`))

  txt.tbl <-tibble(
    y = c(Delta.lrv,Delta.tv),
    label = c('Min TPP','Base TPP')
  )

  levels(Plt.frm$Decision) <- c("Accelerate", "Do not Accelerate")
  ggplot(Plt.frm,aes(x=n,y=`Observed Proportion`)) +
    geom_hline(data=txt.tbl,aes(yintercept = y),linetype = 'dashed',alpha = .8)+
    annotate("text", label = TeX(paste0("Min TPP = ", Delta.lrv*100,"%")),
             x = 0, y = Delta.lrv-.03, size = tsize, colour = "black", hjust = 0)+
    annotate("text", label = TeX(paste0("Base TPP = ", Delta.tv*100,"%")),
             x = 0, y = Delta.tv+.03, size = tsize, colour = "black", hjust = 0)+
    geom_line(aes(color=Decision))+
    geom_point(aes(color=Decision))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    scale_color_manual(values = c('green','red'))+
    scale_y_continuous(labels = scales::percent)+
    xlab('Number of Samples') +
    theme(legend.position = 'bottom')+
    labs(title='Plot of number of samples vs Observed Proportion of Success',
         x="Number of subjects",
         y="Observed proportion of responders",
         color="Decision",
         subtitle="Reported: Sample proportions required to meet/exceed posterior predictive probability thresholds for Go/No-Go")

}
