#' @title Make single sample binary Operating characteristic function vs. treatment effect
#'
#' @param results output from get.ss.bin.trt.int.oc.df
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples \donttest{
#' make.ss.bin.trt.int.oc()
#' }
#' @author Greg Cicconetti
make.ss.bin.trt.int.oc <- function(results = get.ss.bin.trt.int.oc.df()){
ggplot(data= results[[1]], aes(x=p.success, y= prob, color=decision)) +
  geom_line(size=.75) + facet_wrap(~interim.n.t) +
  geom_line(data=results[[2]], aes(x=Delta, y=p), linetype=2) +
  labs(x="Underlying proportion of treatment reponsers", y="Probability",
       title="Operating characteristics a function of treatment effect",
       subtitle="Operating characteristcs associated with study-end overlaid with dashed lines",
       color="decision")+
  theme(legend.position = "bottom")
}
