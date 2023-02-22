#' @title Get single sample normal (with known variance) posterior distribution parameters
#' @description Returns the parameters of the posterior distribution in the case of normal prior with known variance and normal data
#' @param prior.mean prior mean
#' @param prior.var prior variance
#' @param sample.n sample size
#' @param sample.xbar sample mean
#' @param sample.var sample variance
#' @return returns a data.frame hold posterior parameters, sample x.bar, sample variance and posterior mean and variances
#' @author Greg Cicconetti
#' @examples
#'  my.ss.normal.post <- get.ss.normal.post()
#'  head(my.ss.normal.post)
get.ss.normal.post <- function(prior.mean=0, prior.var=1000000, sample.n=10, sample.xbar=seq(-1,1,.01), sample.var=1){
  expand.grid(sample.xbar=sample.xbar, prior.mean=prior.mean, prior.var=prior.var, sample.n=sample.n,sample.var=sample.var) %>%
    mutate(post.var = prior.var*sample.var/(sample.n*prior.var+sample.var),
           post.mean = post.var*(prior.mean/prior.var + sample.n*sample.xbar/sample.var)) %>%
    select(prior.mean, prior.var, sample.n, sample.xbar, sample.var, post.mean, post.var)
}

