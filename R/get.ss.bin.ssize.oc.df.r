#' Get single sample binary sample size operating charactersitics data.frame
#'
#' @param a.trt prior alpha parameter
#' @param b.trt prior beta parameter
#' @param n.trt observed sample size
#' @param Delta.lrv TPP Lower Reference Value aka Min TPP
#' @param Delta.tv TPP Target Value aka Base TPP
#' @param Delta.user User-specified rate
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param SS.OC.N.LB Sample size OC Curve's Lower sample size bound
#' @param SS.OC.N.UB Sample size OC Curve's Upper sample size bound
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.ss.bin.ssize.oc.df <- get.ss.bin.ssize.oc.df()
#' head(my.ss.bin.ssize.oc.df)

get.ss.bin.ssize.oc.df <- function(a.trt = 1, b.trt = 1, n.trt = 40,
         Delta.lrv = .35, Delta.tv=.35, Delta.user = .40,
         tau.tv = 0.10, tau.lrv = 0.8, tau.ng = 0.65,
         SS.OC.N.LB = 20, SS.OC.N.UB = 60){

  # Treatment effect as entered, under null, Base, and Min
  # For a spectrum of sample sizes combination grab the Go/NoGo probabilities
  # For each row we'll run

  specs <- expand.grid(n.trt=SS.OC.N.LB:SS.OC.N.UB) %>%
    mutate(a.trt = a.trt, b.trt = b.trt,
           Delta.lrv = Delta.lrv, Delta.tv=Delta.tv, Delta.user=Delta.user,
           tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)

  # These are custom functions that randomly sample the treatment arm's data
  # Placebo prior and data will be AS ENTERED BY USER
  # Treatment data is sampled from binomial distributons with priors provided by
  # user and prob = DATA ENTERED BY USER

  check <- bind_rows(apply(X = matrix(1:nrow(specs)), MARGIN = 1, FUN = function(x) {
    get.ss.bin.df(a.trt = specs$a.trt[x], b.trt = specs$b.trt[x], n.trt = specs$n.trt[x],
                  x.trt = 0:specs$n.trt[x], Delta.tv = specs$Delta.tv[x],
                  Delta.lrv = specs$Delta.lrv[x],
                  tau.tv = specs$tau.tv[x], tau.lrv = specs$tau.lrv[x],
                  tau.ng = specs$tau.ng[x])})) %>% mutate(
                    p.x.trt.null = dbinom(x = x.trt, size = n.trt, prob = 0),
                    p.x.trt.lrv = dbinom(x = x.trt, size = n.trt, prob = Delta.lrv),
                    p.x.trt.tv = dbinom(x = x.trt, size = n.trt, prob = Delta.tv),
                    p.x.trt.user = dbinom(x = x.trt, size = n.trt, prob = Delta.user))

  check <- check %>% group_by(n.trt, result) %>% summarize(
    s.null = sum(p.x.trt.null),
    s.lrv = sum(p.x.trt.lrv),
    s.tv = sum(p.x.trt.tv),
    s.user = sum(p.x.trt.user)) %>%
    gather(value = value, key=key, -c(n.trt,result), factor_key = T) %>%
    dplyr::filter(result!="Consider") %>%
    mutate(value = ifelse(result=="No-Go", 1-value, value))

  levels(check$key) <- c(
    TeX(paste0("$\\Delta$ = Null = ", round(0,2)*100,'%')),
    TeX(paste0("$\\Delta$ = Min TPP = ", round(Delta.lrv,2)*100,'%')),
    TeX(paste0("$\\Delta$ = Base TPP = ", round(Delta.tv,2)*100,'%')),
    TeX(paste0("$\\Delta$ = User defined = ", round(Delta.user,2)*100,'%')))
  check$key <- factor(check$key, levels(check$key)[order(c(0, Delta.lrv, Delta.tv,
                                                           Delta.user))])

  check <- check %>% mutate(a.trt = a.trt, b.trt = b.trt, n.trt = n.trt,
                            Delta.lrv = Delta.lrv, Delta.tv=Delta.tv, Delta.user = Delta.user,
                            tau.tv = tau.tv, tau.lrv = tau.lrv, tau.ng = tau.ng)
  return(check)
}


