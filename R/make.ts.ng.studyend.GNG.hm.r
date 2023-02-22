#' @title Make two-sample normal-gamma studyend GNG Heatmap
#'
#' @param my.df output from return.ts.ng.studyend.GNG.hm.df
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' make.ts.ng.studyend.GNG.hm()
#' }
#' @author Greg Cicconetti
make.ts.ng.studyend.GNG.hm <- function(my.df = return.ts.ng.studyend.GNG.hm.df()){

  grid.arrange(
ggplot() +
  geom_ribbon(data=my.df, aes(x=xbar.c,ymin=xbar.go, ymax=Inf), color="darkgreen", fill="lightgreen", alpha=.5)+
  geom_ribbon(data=my.df, aes(x=xbar.c,ymax=xbar.ng, ymin=-Inf), color="red", fill="red", alpha=.5)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Observed control mean", y="Observed treatment Mean", color="Decision", title="Observed Control and Treatment means required for Study-end Go and No-Go"),
  ggplot() +
    geom_ribbon(data=my.df, aes(x=xbar.c, ymin=xbar.go - xbar.c, ymax=Inf), color="darkgreen", fill="lightgreen", alpha=.5)+
    geom_ribbon(data=my.df, aes(x=xbar.c, ymax=xbar.ng - xbar.c, ymin=-Inf), color="red", fill="red", alpha=.5)+
    labs(x="Observed control mean", y="Observed treatment effect", color="Decision", title="Observed treatment effect required for Study-end Go and No-Go") +
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0)))

}

