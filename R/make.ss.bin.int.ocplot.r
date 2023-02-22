#' make.ss.bin.interim.ocplot - creates OC curve with Final and Any Interim
#' @title Make single sample binary interim OC curve
#' @param ss.bin.int.GNG call to get.ss.bin.interim.GNG
#' @param ss.bin.int.oc call to get.ss.bin.interim.oc
#' @param ss.bin.studyend.GNG call to get.ss.bin.studyend.GNG
#' @param ss.bin.int.df call to get.ss.bin.int.df
#' @param include_nogo logical
#' @param lower lower bound
#' @param upper upper bound
#' @param goThreshold predictive probability threshold for interim
#' @param nogoThreshold predictive probability threshold for interim
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
#' my.ss.bin.int.oc <- get.ss.bin.int.oc(ss.bin.int.df = my.ss.bin.int.df,
#' ss.bin.int.GNG= my.ss.bin.int.GNG)
#'  my.ss.bin.int.ocplot <- make.ss.bin.int.ocplot(
#'                  ss.bin.studyend.GNG = my.ss.bin.studyend.GNG,
#'                  ss.bin.int.oc = my.ss.bin.int.oc,
#'                  ss.bin.int.GNG=my.ss.bin.int.GNG,
#'                  ss.bin.int.df = my.ss.bin.int.df,
#'                  goThreshold = .8,
#'                  nogoThreshold = .8,
#'                  include_nogo =TRUE, lower=0, upper=1)
#' my.ss.bin.int.ocplot
#' }
make.ss.bin.int.ocplot <- function(
                ss.bin.studyend.GNG,
                ss.bin.int.oc,
                ss.bin.int.GNG,
                ss.bin.int.df,
                goThreshold = .8,
                nogoThreshold = .8,
                include_nogo =TRUE, lower=0, upper=1){

  Plt.frm <- ss.bin.int.oc %>% dplyr::filter(!grepl(pattern = 'Analysis',x = Type))

  txt <- tibble(
    targetRate = c(ss.bin.studyend.GNG[[1]]$Delta.lrv,ss.bin.studyend.GNG[[1]]$Delta.tv),
    label = paste(c('Min TPP =','Target TPP ='),targetRate),
    align = c(1,0),
    Rate = .05
  )
   if(include_nogo==TRUE){
  for.return <- ggplot(Plt.frm,aes(x=targetRate))+
    geom_ribbon(data = dplyr::filter(Plt.frm,Type == 'Final'),aes(ymin = lower,ymax=upper,fill=Decision), alpha=.6, show.legend = F)+
    geom_line(aes(y = line, color= Decision,linetype=Type),size = 1.1,data = dplyr::filter(Plt.frm,Decision!="Consider"))+
    geom_vline(xintercept = txt$targetRate,linetype = 'dashed', color="blue")+
    geom_text(data = txt,color='black', aes(label = label,hjust = align,y=Rate), size=6)+
    scale_x_continuous(expand= c(0,0), limits=c(lower,upper), labels = scales::percent)+
    scale_y_continuous(expand= c(0,0),labels = scales::percent)+
    scale_color_manual(values = c('green','grey','red'),drop = F)+
    scale_fill_manual(values = c('green','grey','red'))+
    scale_linetype_manual(values = c(1,3))+
    theme(legend.position = 'bottom')+
    labs(title="Operating Characteristics as a function of treatment effect", x="Response Rate", y="Probability", color="Conclusion", linetype="Analysis")
  } else if(include_nogo==FALSE) {
    Plt.frm <- dplyr::filter(ss.bin.int.oc,!grepl('Analysis',Type))

    txt <- tibble(
      targetRate = c(ss.bin.studyend.GNG[[1]]$Delta.lrv,ss.bin.studyend.GNG[[1]]$Delta.tv),
      label = paste(c('Min TPP =','Target TPP ='),targetRate),
      align = c(1,0),
      Rate = .05
    )
  for.return <- ggplot(Plt.frm,aes(x=targetRate))+
    geom_ribbon(data = dplyr::filter(Plt.frm,Type == 'Final'),aes(ymin = lower,ymax=upper,fill=Decision), alpha=.6, show.legend = F)+
    geom_line(aes(y = line, color= Decision,linetype=Type),size = 1.1,data = dplyr::filter(Plt.frm,Decision=='Go' | (Decision=="No Go" & Type=="Final")))+
    geom_vline(xintercept = txt$targetRate,linetype = 'dashed', color="blue")+
    geom_text(data = txt,color='black', aes(label = label,hjust = align,y=Rate), size=6)+
    scale_x_continuous(expand= c(0,0), limits=c(lower,upper), labels = scales::percent)+
    scale_y_continuous(expand= c(0,0),labels = scales::percent)+
    scale_color_manual(values = c('green','grey','red'),drop = F)+
    scale_fill_manual(values = c('green','grey','red'))+
    scale_linetype_manual(values = c(1,3))+
    theme(legend.position = 'bottom')+
    labs(title="Operating Characteristics as a function of treatment effect", x="Response Rate", y="Probability", color="Conclusion", linetype="Analysis")

  }
  return(for.return)
}
