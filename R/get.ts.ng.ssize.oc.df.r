#' @title Get two sample normal-gamma sample size OC data.frame
#'
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size for control group
#' @param s.c sample sd for control group
#' @param group.c group label for control group
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size for treatment group
#' @param s.t sample sd for treatment group
#' @param group.t group label for treatment group
#' @param ARatio randomization ratio
#' @param N Total sample size
#' @param SS.OC.N.LB lower bound for OC curve
#' @param SS.OC.N.UB Upper bound for OC Curve
#' @param Delta.tv Base TPP
#' @param Delta.lrv Min TPP
#' @param Delta.user User's delta
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param npoints number of points
#' @param n.MC n for MC sampling
#' @param goparallel logical to use parallel programming
#' @param cl cluster
#' @param seed random seed
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' my.ts.ng.ssize.oc.df <- get.ts.ng.ssize.oc.df()
#' my.ts.ng.ssize.oc.df
#' }

get.ts.ng.ssize.oc.df <- function(
  mu.0.c = 0, n.0.c = 0.0001, alpha.0.c=.25, beta.0.c = 1, s.c = 4,  group.c="Control",
  mu.0.t = 0, n.0.t = 0.0001, alpha.0.t=.25, beta.0.t = 1, s.t = 4,  group.t="Treatment",
  ARatio=1,
  SS.OC.N.LB = 20, SS.OC.N.UB=200,
  Delta.tv = 1.5, Delta.lrv = 1, Delta.user = 1.40,
  tau.tv = 0, tau.lrv = 0.80, tau.ng = 0.75,
  npoints=1, n.MC = 500,
  seed = 1234, goparallel = FALSE, cl=cl){

  ## We need to do this for a variety of sample sizes
  my.specs <- expand.grid(N=floor(seq(SS.OC.N.LB, SS.OC.N.UB, length.out=npoints))) %>%
    mutate(mu.0.c = mu.0.c,
           n.0.c = n.0.c,
           alpha.0.c=alpha.0.c,
           beta.0.c = beta.0.c,
           s.c = s.c,
           n.c = N - floor(N*(ARatio/(1+ARatio))),
           mu.0.t = mu.0.t,
           n.0.t = n.0.t,
           alpha.0.t=alpha.0.t,
           beta.0.t = beta.0.t,
           s.t = s.t,
           n.t = floor(N*(ARatio/(1+ARatio))),
           Delta.tv = Delta.tv,
           Delta.lrv = Delta.lrv,
           Delta.user=Delta.user,
           tau.tv = tau.tv,
           tau.lrv = tau.lrv,
           tau.ng = tau.ng,
           n.MC = n.MC,
           npoints = npoints,
           seed = seed + 0:(n()-1))

  my.FUN <- function(x){
    set.seed(my.specs$seed[x])
    # Get x-bar samples pulled from the x'th row of my.specs
    # Reminder: Across rows of specs it is n.c and n.t that are changing!
    my.control.null <- rnorm(n = my.specs$n.MC[x], mean = my.specs$mu.0.c[x], sd=my.specs$s.c[x]/sqrt(my.specs$n.c[x]))
    my.treatment.null <- rnorm(n = my.specs$n.MC[x], mean = my.specs$mu.0.c[x], sd=my.specs$s.t[x]/sqrt(my.specs$n.t[x]))
    my.treatment.tv <- rnorm(n = my.specs$n.MC[x], mean = my.specs$mu.0.c[x] + my.specs$Delta.tv[x], sd=my.specs$s.t[x]/sqrt(my.specs$n.t[x]))
    my.treatment.lrv <- rnorm(n = my.specs$n.MC[x], mean = my.specs$mu.0.c[x] + my.specs$Delta.lrv[x], sd=my.specs$s.t[x]/sqrt(my.specs$n.t[x]))
    my.treatment.user <- rnorm(n = my.specs$n.MC[x], mean = my.specs$mu.0.c[x] + my.specs$Delta.user[x], sd=my.specs$s.t[x]/sqrt(my.specs$n.t[x]))

    # Merge with specs - these hold the sample means observed over n.MC trials
    for.apply.null <- cbind(my.control.null, my.treatment.null, my.specs[x,])
    for.apply.tv <- cbind(my.control.null, my.treatment.tv, my.specs[x,])
    for.apply.lrv <- cbind(my.control.null, my.treatment.lrv, my.specs[x,])
    for.apply.user <- cbind(my.control.null, my.treatment.user, my.specs[x,])

    # Now evaluate decision outcomes based on the random sampling
    # These will be tallied later
    message("\n Get results begin\n")
    bind_rows(apply(X = matrix(1:nrow(for.apply.lrv)), MARGIN = 1, FUN = function(y) {

      # Now get.ts.ng.mc will:
      # take these passed parameters and compute posterior parameters.
      # Then it will sample from marginal t-distributions, compute n.MC differences to create estimates for
      # P.R1 = mean(Z.n >= Delta.lrv) and
      # P.R3 = mean(Z.n >= Delta.tv).
      # Then: ifelse(P.R1 > tau.lrv & P.R3 > tau.tv, "Go",  ifelse(P.R3 <= tau.tv & P.R1 <= tau.ng, "No-Go", "Consider"))
      for.return.null <- get.ts.ng.mc(mu.0.c = for.apply.null$mu.0.c[y], n.0.c = for.apply.null$n.0.c[y],
                                      alpha.0.c=for.apply.null$alpha.0.c[y], beta.0.c = for.apply.null$beta.0.c[y],
                                      xbar.c = for.apply.null$my.control.null[y], s.c = for.apply.null$s.c[y], n.c = for.apply.null$n.c[y], group.c="Control",
                                      mu.0.t = for.apply.null$mu.0.t[y], n.0.t = for.apply.null$n.0.t[y], alpha.0.t=for.apply.null$alpha.0.t[y], beta.0.t = for.apply.null$beta.0.t[y],
                                      xbar.t =  for.apply.null$my.treatment.null[y], s.t = for.apply.null$s.t[y], n.t = for.apply.null$n.t[y], group.t="Treatment",
                                      Delta.tv = for.apply.null$Delta.tv[y], Delta.lrv = for.apply.null$Delta.lrv[y],
                                      tau.tv = for.apply.null$tau.tv[y], tau.lrv = for.apply.null$tau.lrv[y],
                                      tau.ng = for.apply.null$tau.ng[y],
                                      n.MC = for.apply.null$n.MC[y]) %>%
        mutate(treatment.effect = "null",
               mu.0.c = for.apply.null$mu.0.c[y], n.0.c = for.apply.null$n.0.c[y],
               alpha.0.c=for.apply.null$alpha.0.c[y], beta.0.c = for.apply.null$beta.0.c[y],
               xbar.c = for.apply.null$my.control.null[y], s.c = for.apply.null$s.c[y], n.c = for.apply.null$n.c[y], group.c="Control",
               mu.0.t = for.apply.null$mu.0.t[y], n.0.t = for.apply.null$n.0.t[y], alpha.0.t=for.apply.null$alpha.0.t[y], beta.0.t = for.apply.null$beta.0.t[y],
               xbar.t =  for.apply.null$my.treatment.null[y], s.t = for.apply.null$s.t[y], n.t = for.apply.null$n.t[y], group.t="Treatment",
               Delta.tv = for.apply.null$Delta.tv[y], Delta.lrv = for.apply.null$Delta.lrv[y],
               tau.tv = for.apply.null$tau.tv[y], tau.lrv = for.apply.null$tau.lrv[y],
               tau.ng = for.apply.null$tau.ng[y],
               n.MC = for.apply.null$n.MC[y])

      for.return.tv <- get.ts.ng.mc(mu.0.c = for.apply.tv$mu.0.c[y], n.0.c = for.apply.tv$n.0.c[y],
                                    alpha.0.c=for.apply.tv$alpha.0.c[y], beta.0.c = for.apply.tv$beta.0.c[y],
                                    xbar.c = for.apply.tv$my.control.null[y], s.c = for.apply.tv$s.c[y], n.c = for.apply.tv$n.c[y], group.c="Control",
                                    mu.0.t = for.apply.tv$mu.0.t[y], n.0.t = for.apply.tv$n.0.t[y], alpha.0.t=for.apply.tv$alpha.0.t[y], beta.0.t = for.apply.tv$beta.0.t[y],
                                    xbar.t =  for.apply.tv$my.treatment.tv[y], s.t = for.apply.tv$s.t[y], n.t = for.apply.tv$n.t[y], group.t="Treatment",
                                    Delta.tv = for.apply.tv$Delta.tv[y], Delta.lrv = for.apply.tv$Delta.lrv[y],
                                    tau.tv = for.apply.tv$tau.tv[y], tau.lrv = for.apply.tv$tau.lrv[y],
                                    tau.ng = for.apply.tv$tau.ng[y],
                                    n.MC = for.apply.tv$n.MC[y]) %>%
        mutate(treatment.effect = "tv",
               mu.0.c = for.apply.tv$mu.0.c[y], n.0.c = for.apply.tv$n.0.c[y],
               alpha.0.c=for.apply.tv$alpha.0.c[y], beta.0.c = for.apply.tv$beta.0.c[y],
               xbar.c = for.apply.tv$my.control.null[y], s.c = for.apply.tv$s.c[y], n.c = for.apply.tv$n.c[y], group.c="Control",
               mu.0.t = for.apply.tv$mu.0.t[y], n.0.t = for.apply.tv$n.0.t[y], alpha.0.t=for.apply.tv$alpha.0.t[y], beta.0.t = for.apply.tv$beta.0.t[y],
               xbar.t =  for.apply.tv$my.treatment.tv[y], s.t = for.apply.tv$s.t[y], n.t = for.apply.tv$n.t[y], group.t="Treatment",
               Delta.tv = for.apply.tv$Delta.tv[y], Delta.lrv = for.apply.tv$Delta.lrv[y],
               tau.tv = for.apply.tv$tau.tv[y], tau.lrv = for.apply.tv$tau.lrv[y],
               tau.ng = for.apply.tv$tau.ng[y],
               n.MC = for.apply.tv$n.MC[y])

      for.return.lrv <- get.ts.ng.mc(mu.0.c = for.apply.lrv$mu.0.c[y], n.0.c = for.apply.lrv$n.0.c[y],
                                     alpha.0.c=for.apply.lrv$alpha.0.c[y], beta.0.c = for.apply.lrv$beta.0.c[y],
                                     xbar.c = for.apply.lrv$my.control.null[y], s.c = for.apply.lrv$s.c[y], n.c = for.apply.lrv$n.c[y], group.c="Control",
                                     mu.0.t = for.apply.lrv$mu.0.t[y], n.0.t = for.apply.lrv$n.0.t[y], alpha.0.t=for.apply.lrv$alpha.0.t[y], beta.0.t = for.apply.lrv$beta.0.t[y],
                                     xbar.t =  for.apply.lrv$my.treatment.lrv[y], s.t = for.apply.lrv$s.t[y], n.t = for.apply.lrv$n.t[y], group.t="Treatment",
                                     Delta.tv = for.apply.lrv$Delta.tv[y], Delta.lrv = for.apply.lrv$Delta.lrv[y],
                                     tau.tv = for.apply.lrv$tau.tv[y], tau.lrv = for.apply.lrv$tau.lrv[y],
                                     tau.ng = for.apply.lrv$tau.ng[y],
                                     n.MC = for.apply.lrv$n.MC[y]) %>%
        mutate(treatment.effect = "lrv",
               mu.0.c = for.apply.lrv$mu.0.c[y], n.0.c = for.apply.lrv$n.0.c[y],
               alpha.0.c=for.apply.lrv$alpha.0.c[y], beta.0.c = for.apply.lrv$beta.0.c[y],
               xbar.c = for.apply.lrv$my.control.null[y], s.c = for.apply.lrv$s.c[y], n.c = for.apply.lrv$n.c[y], group.c="Control",
               mu.0.t = for.apply.lrv$mu.0.t[y], n.0.t = for.apply.lrv$n.0.t[y], alpha.0.t=for.apply.lrv$alpha.0.t[y], beta.0.t = for.apply.lrv$beta.0.t[y],
               xbar.t =  for.apply.lrv$my.treatment.lrv[y], s.t = for.apply.lrv$s.t[y], n.t = for.apply.lrv$n.t[y], group.t="Treatment",
               Delta.tv = for.apply.lrv$Delta.tv[y], Delta.lrv = for.apply.lrv$Delta.lrv[y],
               tau.tv = for.apply.lrv$tau.tv[y], tau.lrv = for.apply.lrv$tau.lrv[y],
               tau.ng = for.apply.lrv$tau.ng[y],
               n.MC = for.apply.lrv$n.MC[y])

      for.return.user <- get.ts.ng.mc(mu.0.c = for.apply.user$mu.0.c[y], n.0.c = for.apply.user$n.0.c[y],
                                      alpha.0.c=for.apply.user$alpha.0.c[y], beta.0.c = for.apply.user$beta.0.c[y],
                                      xbar.c = for.apply.user$my.control.null[y], s.c = for.apply.user$s.c[y], n.c = for.apply.user$n.c[y], group.c="Control",
                                      mu.0.t = for.apply.user$mu.0.t[y], n.0.t = for.apply.user$n.0.t[y], alpha.0.t=for.apply.user$alpha.0.t[y], beta.0.t = for.apply.user$beta.0.t[y],
                                      xbar.t =  for.apply.user$my.treatment.user[y], s.t = for.apply.user$s.t[y], n.t = for.apply.user$n.t[y], group.t="Treatment",
                                      Delta.tv = for.apply.user$Delta.tv[y], Delta.lrv = for.apply.user$Delta.lrv[y],
                                      tau.tv = for.apply.user$tau.tv[y], tau.lrv = for.apply.user$tau.lrv[y],
                                      tau.ng = for.apply.user$tau.ng[y],
                                      n.MC = for.apply.user$n.MC[y]) %>%
        mutate(treatment.effect = "user",
               mu.0.c = for.apply.user$mu.0.c[y], n.0.c = for.apply.user$n.0.c[y],
               alpha.0.c=for.apply.user$alpha.0.c[y], beta.0.c = for.apply.user$beta.0.c[y],
               xbar.c = for.apply.user$my.control.null[y], s.c = for.apply.user$s.c[y], n.c = for.apply.user$n.c[y], group.c="Control",
               mu.0.t = for.apply.user$mu.0.t[y], n.0.t = for.apply.user$n.0.t[y], alpha.0.t=for.apply.user$alpha.0.t[y], beta.0.t = for.apply.user$beta.0.t[y],
               xbar.t =  for.apply.user$my.treatment.user[y], s.t = for.apply.user$s.t[y], n.t = for.apply.user$n.t[y], group.t="Treatment",
               Delta.tv = for.apply.user$Delta.tv[y], Delta.lrv = for.apply.user$Delta.lrv[y],
               tau.tv = for.apply.user$tau.tv[y], tau.lrv = for.apply.user$tau.lrv[y],
               tau.ng = for.apply.user$tau.ng[y],
               n.MC = for.apply.user$n.MC[y])
      for.return <- rbind(for.return.null, for.return.lrv, for.return.tv, for.return.user)

      return(for.return)
    }))
  }

  # Now for each row of my.specs
  # Call my.FUN:
  # gpackage::goparallel()
  if(goparallel==TRUE){
  parallel::clusterExport(cl=cl, varlist=list("get.ts.ng.mc", "my.specs", "my.FUN","get.ts.ng.mc.df", "get.ng.post", "rt_ls"), envir = environment())
  parallel::clusterEvalQ(cl = cl, expr = {

          requireNamespace("dplyr");
          requireNamespace("tidyr")})
  get.results <- bind_rows(parallel::parApply(cl=cl, X = matrix(1:nrow(my.specs)), MARGIN = 1, FUN = my.FUN))} else{
          get.results <- bind_rows(apply(X = matrix(1:nrow(my.specs)), MARGIN = 1, FUN = my.FUN))
  }

  (get.results %>% dplyr::filter(treatment.effect=="tv") %>% .$P.R1)
  (get.results %>% dplyr::filter(treatment.effect=="tv") %>% .$P.R3)
  get.results %>% mutate(result1 = ifelse(P.R1 > tau.lrv & P.R3 > tau.tv, "Go",
                         ifelse(P.R3 <= tau.tv & P.R1 <= tau.ng, "No-Go", "Consider")),
                         result2 = ifelse(P.R1 > tau.lrv & P.R3 >= tau.tv, "Go",
                                          ifelse(P.R3 <= tau.tv & P.R1 <= tau.ng, "No-Go", "Consider"))
                         ) %>% summarize(go1 = mean(result1 == "Go"),
                                         go2 = mean(result2 == "Go"),
                                         no1 = mean(result1 =="No-Go"),
                                         no2 = mean(result2 =="No-Go"))


  get.results %>% dplyr::filter(treatment.effect=="tv") %>% ggplot(aes(x=P.R1, y=P.R3, color=result)) +geom_point()

  get.results %>% dplyr::filter(treatment.effect=="tv") %>% dplyr::select(xbar.t, xbar.c) %>% pivot_longer(cols=1:2) %>% ggplot(aes(x=value,fill=name)) +geom_density()+
    labs(title="Densities of observed sample means sampled from marginal t-distributions")
  get.results %>% dplyr::filter(treatment.effect=="tv") %>% mutate(obs.diff = xbar.t - xbar.c) %>% ggplot(aes(x=obs.diff)) + geom_density()+labs(title="Density of observed differences in sample means")
  get.results %>% dplyr::filter(treatment.effect=="tv") %>% mutate(obs.diff = xbar.t - xbar.c) %>% ggplot(aes(x=obs.diff, fill=result, color=result)) + geom_density(alpha=.25)+geom_rug()+
    labs(title="Density of observed differences grouped by outcome of rule")
  get.results %>% dplyr::filter(treatment.effect=="tv") %>% ggplot(aes(x=xbar.t, y=xbar.c, color=result))+geom_point()+labs(title="Scatter plot of simulated trial results, colored by outcome.")
  get.results %>% dplyr::filter(result=="Consider", treatment.effect=="tv") %>% mutate(obs.diff = xbar.t - xbar.c) %>% ggplot(aes(x=obs.diff)) + geom_density()+labs(title="Density of observed differences for leading to Consider.")



  message("\n get.results finished\n")
  for.plot <- get.results %>% group_by(treatment.effect, mu.0.c, n.0.c, alpha.0.c, beta.0.c, s.c, n.c, group.c, mu.0.t, n.0.t, alpha.0.t, beta.0.t, s.t, n.t, group.t, Delta.tv, Delta.lrv, tau.tv, tau.lrv, tau.ng, n.MC) %>%
    summarize(NoGo=1 - sum(result=="No-Go")/my.specs$n.MC[1],
              Go = sum(result=="Go")/my.specs$n.MC[1]) %>% gather(value=value, key=result, Go, NoGo) %>% mutate(n.total=n.c+n.t)
  for.plot$treatment.effect <- factor(for.plot$treatment.effect, c("null", "lrv", "tv", "user"))
  levels(for.plot$treatment.effect) <-  c(TeX("$\\Delta$ = NULL = 0%"),
                                          TeX(paste("$\\Delta$ = Min TPP = $", Delta.lrv, "$ ")),
                                          TeX(paste("$\\Delta$ = Base TPP = $", Delta.tv, "$ ")),
                                          TeX(paste("$\\Delta$ = User defined = $", Delta.user, "$")))
  for.plot$treatment.effect <- factor(for.plot$treatment.effect, levels(for.plot$treatment.effect)[order(c(0, Delta.lrv, Delta.tv, Delta.user))])
  for.plot$treatment.effect2 <- for.plot$treatment.effect
  levels(for.plot$treatment.effect2)
  message("\nDone with for.plot for make.ts.ng.ssize.oc\n")
  return(for.plot)

}
