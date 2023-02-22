#' @title Pooled standard deviation
#' @param s expression to be plotted
#' @param n number of points to plot
#' @return A vector holding the pooled standard deviation is returned.
#' @description Computes to the pooled standard deviation.
#' @author Greg Cicconetti
#' @examples
#' pooled.sd()
pooled.sd <- function(s = c(4, 5), n = c(14, 20)){
  return(sqrt(sum(s  ^  2 * (n - 1))/(sum(n) - (length(n) - 1))))}


