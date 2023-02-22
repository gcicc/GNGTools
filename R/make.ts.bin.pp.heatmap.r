#' @title Make two-sample binary predictive probability heatmap
#'
#' @param my.data output from return.ts.bin.pp.table
#' @param go.thresh go threshold
#' @param ng.thresh no-go threshold
#'
#' @return A ggplot object is returned.
#' @export
#'
#' @examples \donttest{
#' pp.table <- return.ts.bin.int.predprob.table(goparallel=FALSE)
#' make.ts.bin.pp.heatmap(my.data = pp.table, go.thresh=0.8, ng.thresh=0.8)
#' }
make.ts.bin.pp.heatmap <- function(my.data = pp.table, go.thresh=0.8, ng.thresh=0.8) {

        # In general we would bucket group into Go/NO Go and consider and plot those colors

        p1 <- my.data %>% ggplot(aes(x=x.int.con, y=x.int.trt, fill=P.Go)) +
                geom_tile(color="grey80") +
                labs(x="Number of control successes at interim",
                     y="Number of response successes at interim",
                     title="Predictive probabilities of Go and No-Go",
                     subtitle="Posterior predictive probability that study end Go is met given ")

        p2 <- my.data %>% ggplot(aes(x=x.int.con, y=x.int.trt, fill=P.NG)) +
                geom_tile(color="grey80") +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0))+
                labs(x="Number of control successes at interim",
                     y="Number of response successes at interim",
                     title="Predictive probabilities of Go and No-Go",
                     subtitle="Posterior predictive probability that study end No-go is met given control data at interim")

        my.data <- my.data %>% mutate(result = case_when(
                P.Go > go.thresh ~ "Go",
                P.NG > ng.thresh ~ "No-Go",
                TRUE ~ "Consider"))

        p3 <- my.data %>% ggplot(aes(x=x.int.con, y=x.int.trt, fill=result)) +
                geom_tile(color="grey80") +
                labs(x="Number of control successes at interim",
                     y="Number of response successes at interim",
                     title="Predictive probabilities of Go and No-Go",
                     subtitle="Posterior predictive probability that study end Go is met given ") +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0))

        grid.arrange(p1,p2, p3)
}
