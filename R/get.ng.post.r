#' @title Normal-Gamma Posterior Updating
#' @param mu.0 prior mean
#' @param n.0 prior effective sample size
#' @param alpha.0 prior alpha parameter
#' @param beta.0 prior beta parameter
#' @param xbar observed sampled mean
#' @param s observed sample standard deviation
#' @param n sample size
#' @param group text string for group label
#' @return Returns a data.frame with prior, data, and posterior parameters.
#' @examples
#' my.ng.post <- get.ng.post(mu.0 = 0, n.0 = 10, alpha.0 = .25, beta.0 = 1,
#'  xbar = .25, s = 3, n = 15, group = "Control")
#' my.ng.post
#' @author Greg Cicconetti
#' @description Normal-Gamma Posterior Updating

get.ng.post <- function(mu.0 = 0, n.0 = 10, alpha.0 = .25, beta.0 = 1,
                        xbar = .25, s = 3, n = 15, group = "Control"){

  expand.grid(mu.0 = mu.0, n.0 = n.0, alpha.0 = alpha.0, beta.0 = beta.0,
              xbar = xbar, s = s, n = n, group =group) %>%
    mutate(
      # Prior NG parameters
      mu.0 = mu.0,
      n.0 = n.0,
      alpha.0 = alpha.0,
      beta.0 = beta.0,
      # Marginal parameters of t-distribution describing mu.0
      tdf.0 = 2 * alpha.0,
      sigma.0 = beta.0/(alpha.0 * n.0),
      # data
      xbar = xbar,
      s = s,
      n = n,
      # Posterior parameters
      mu.n = (n.0 * mu.0 + n * xbar)/(n.0 + n),
      n.n = n.0 + n,
      alpha.n = alpha.0 + n / 2,
      beta.n = beta.0 + (n - 1)/2 * s ^ 2 + n.0 * n * (xbar - mu.0) ^ 2/(2 * (n.0 + n)),
      # Marginal parameters of t-distribution describing mu.n
      tdf.n = 2 * alpha.n,
      sigma.n = beta.n/(alpha.n * n.n)
    )
}

