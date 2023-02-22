#' @title Get normal single sample posterior parameters
#'
#' @param mu.0 prior mean
#' @param n.0 prior effective sample size
#' @param alpha.0 prior alpha parameter
#' @param beta.0 prior beta parameter
#' @param xbar observed sampled mean
#' @param s observed sample standard deviation
#' @param n sample size
#' @param group text string for group label
#'
#' @return Returns normal single sample posterior parameters
#' @export
#' @author Greg Cicconetti
#' @examples
#' my.normal.ss.post <- get.normal.ss.post(mu.0=0, n.0=10, alpha.0=.25,
#' beta.0=1, xbar=4, s=3, n=15, group="Placebo")
#' my.normal.ss.post

get.normal.ss.post <- function(mu.0=0, n.0=10, alpha.0=.25, beta.0=1, xbar=4, s=3, n=15, group="Placebo"){
  expand.grid(mu.0=mu.0, n.0=n.0, alpha.0=alpha.0, beta.0=beta.0, xbar=xbar, s=s, n=n, group=group) %>%
mutate(
    # Posterior parameters
    alpha.n=alpha.0 + n/2,
    beta.n=beta.0 + (n-1)/(2)*s^2 + (n*n.0)/(2*(n+n.0))*(xbar-mu.0)^2,
    mean.tau.n = (alpha.n)/beta.n,
    mean.mu.n=(n*mean.tau.n)/(n*mean.tau.n + n.0*mean.tau.n)*xbar + (n.0*mean.tau.n)/(n*mean.tau.n + n.0*mean.tau.n)*mu.0,
    n.n = n.0 + n,
    # Identifying the mode for use in deriving posterior mean parameter
    # If we were to take a sample from the gamma distribution we'd need to set a seed to ensure reproducibility
    post.mu.tau=n*mean.tau.n + n.0*mean.tau.n)

}
