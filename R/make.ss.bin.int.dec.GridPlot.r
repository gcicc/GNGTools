#' @title Make single sample binary interim decision grid
#' @param InterimDF call to get.ss.bin.interim.df
#' @param lower.bound.go lower bound for Go
#' @param lower.bound.ng Lower (upper!) bound for no-go
#' @param goThreshold  predictive probability threshold
#' @param nogoThreshold predictive probability threshold
#' @param include_nogo logical
#' @param add.black logical
#' @return A ggplot object is returned.
#' @export
#'
make.ss.bin.int.dec.GridPlot <- function(InterimDF,
                                         goThreshold = .8,
                                         nogoThreshold = .8,
                                         include_nogo = TRUE,
                                         lower.bound.go = 5,
                                         lower.bound.ng=10,
                                         add.black=FALSE){

  if(add.black==TRUE){
  InterimDF <- InterimDF %>% mutate(
    Decision= as.character(Decision),
    Decision=ifelse(n<lower.bound.go & Decision == "Go", "More Data Required", Decision),
    Decision=ifelse(n<lower.bound.ng & Decision == "No Go", "More Data Required", Decision),
    Decision=factor(Decision, c("Go","Consider", "No Go",  "More Data Required")))
  }

  # Fix Updated 7/13/21
  if(InterimDF$threshold.NoGo[1] >= 1) levels(InterimDF$Decision) <- c("Go", "Consider", "Consider", "More Data Required")

  levels(InterimDF$Decision) <- c("Go", "Consider", "No-Go", "More Data Required")
  ggplot(InterimDF,aes(x=x,y=n,fill=Decision))+
    geom_tile(color = 'grey60', alpha=.75)+
    scale_y_reverse(expand = expansion(0,0),breaks = function(x){(x[1]+.5):(x[2]+.5)})+
    scale_fill_manual(values = c('green','grey','red', "black"))+
    scale_x_continuous(expand = expansion(0,0),breaks = function(x){(x[1]+.5):(x[2]+.5)})+
    theme(legend.position = 'bottom', panel.grid = element_blank())+
    labs(title="Number of responders required for Go/No-Go for each trial outcome",
         subtitle="Trial results are categorized by posterior predictive probability thresholds",
         y='Number of Subjects available',
         x='Number of responders',
         fill="Decision")
}
