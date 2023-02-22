#' @title Make two-sample binary study end GNG heatmap
#'
#' @param my.df output from return.ts.bin.studyend.GNG.hm.df
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' make.ts.bin.studyend.GNG.hm()
#' }
make.ts.bin.studyend.GNG.hm <- function(my.df = return.ts.bin.studyend.GNG.hm.df()){
  my.df %>% ggplot() +
    geom_stepribbon(aes(x=x.con, ymin=0, ymax=x.trt.ng), fill="red")+
    geom_stepribbon(aes(x=x.con, ymin=x.trt.go, ymax=n.trt), fill="green")+
    geom_hline(aes(yintercept=x.con), color="grey80")+
    geom_vline(aes(xintercept=x.con), color="grey80")+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0)) +
    labs(x="Number of success on control",
         y="Number of successes required on treatment arm",
         title="Number of successes required for study-end Go and No-Go given number of success on control successes.")

}



