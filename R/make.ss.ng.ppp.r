#' @title Make singple sample normal-gamma prior posterior plot
#'
#' @param mu.0.t prior mean
#' @param n.0.t prior effective sample size
#' @param alpha.0.t prior alpha parameter
#' @param beta.0.t prior beta parameter
#' @param xbar.t sample mean for treatment group
#' @param s.t sample sd for treatment group
#' @param n.t sample size for treatment group
#' @param limits upper limits used for visualizations
#' @param gamma.IG.sd.limits limits used for Precision, Variance and standard deviation visualizers
#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' my.ss.ng.ppp <- make.ss.ng.ppp()
#' my.ss.ng.ppp[[1]][[1]]
#' my.ss.ng.ppp[[1]][[2]]
#' gridExtra::grid.arrange(my.ss.ng.ppp[[1]][[1]], my.ss.ng.ppp[[1]][[2]], ncol=2)
#' my.ss.ng.ppp[[2]]
#' my.ss.ng.ppp[[3]]
#' my.ss.ng.ppp[[4]]
#' my.ss.ng.ppp[[5]]
#' my.ss.ng.ppp[[6]]
#' my.ss.ng.ppp[[7]]
make.ss.ng.ppp <- function(mu.0.t = 0, n.0.t = 1, alpha.0.t = .25, beta.0.t = 1,
                           xbar.t = 1.75, s.t = 2, n.t = 50, gamma.IG.sd.limits=c(as.numeric(trimws(unlist(strsplit("5, 25, 25",",")))))) {

  pp <-   get.ng.post(mu.0 = mu.0.t, n.0 = n.0.t, alpha.0 = alpha.0.t, beta.0 = beta.0.t,
                      xbar = xbar.t, s = s.t, n = n.t, group = "Treatment") %>% dplyr::select(group, everything()) %>%
    mutate(prior.ng.mean = paste0("(",round(mu.0,4), ", ", round(alpha.0/beta.0,4), ")"),
           prior.mode = paste0("(",round(mu.0,4), ", ", ifelse(alpha.0 < 0.5, NA, round((alpha.0-1/2)/beta.0,4)), ")"),
           prior.ng.var = paste0("(", round(beta.0/(n.0.t*(alpha.0-1)),4), ", ", round(alpha.0/beta.0^2,4), ")"),
           post.ng.mean = paste0("(", round(mu.n,4), ", ", round(alpha.n/beta.n,4), ")"),
           post.mode = paste0("(",round(mu.n,4), ", ", ifelse(alpha.n < 0.5, NA, round((alpha.n-1/2)/beta.n,4)), ")"),
           post.ng.var = paste0("(", round(beta.n/(n.n*(alpha.n-1)),4), ", ", round(alpha.n/beta.n^2,4), ")"))%>% rename(Group=group)

  tau.limits.prior <- c(pmax(.1, qgamma(p = .0005,shape = pp$alpha.0, rate = pp$beta.0)),
                        qgamma(p = .9995, shape = pp$alpha.0, rate = pp$beta.0))
  tau.limits.post <- c(pmax(.1, qgamma(p = .0005,shape = pp$alpha.n, rate = pp$beta.n)),
                        qgamma(p = .9995,shape = pp$alpha.n, rate = pp$beta.n))

  tau.limits <- c(min(pmax(.1, qgamma(p = .0005,shape = alpha.0.t, rate = beta.0.t)),
                      pmax(.1, qgamma(p = .0005,shape = pp$alpha.n, rate = pp$beta.n))),
                  max(qgamma(p = .9995,shape = alpha.0.t, rate = beta.0.t),
                      qgamma(p = .9995,shape = pp$alpha.n, rate = pp$beta.n)))
  mu.limits.prior <- c(qnorm(p = .0005, mean = pp$mu.0, sd = (pp$alpha.0 / pp$beta.0) ^ (-1)/sqrt(pp$n.0)),
                      qnorm(p = .9995, mean = pp$mu.0, sd = (pp$alpha.0 / pp$beta.0) ^ (-1)/sqrt(pp$n.0)))
  mu.limits.post <- c(qnorm(p = .0005, mean = pp$mu.n, sd = (pp$alpha.n / pp$beta.n) ^ (-1)/sqrt(pp$n.n)),
                      qnorm(p = .9995, mean = pp$mu.n, sd = (pp$alpha.n / pp$beta.n) ^ (-1)/sqrt(pp$n.n)))

  mu.limits <- c(min(qnorm(p = .0005, mean = mu.0.t,
                           sd = (alpha.0.t / beta.0.t) ^ (-1)/sqrt(n.0.t)),
                     qnorm(p = .0005, mean = pp$mu.n,
                           sd = (pp$alpha.n / pp$beta.n) ^ (-1)/sqrt(pp$n.n))),
                 max(qnorm(p = .9995, mean = mu.0.t,
                           sd = (alpha.0.t / beta.0.t) ^ (-1)/sqrt(n.0.t)),
                     qnorm(p = .9995, mean = pp$mu.n,
                           sd = (pp$alpha.n / pp$beta.n) ^ (-1)/sqrt(pp$n.n))))

  prior <- expand.grid(
    tau = seq(tau.limits.prior[1],tau.limits.prior[2], length.out = 100),
    mu = seq(mu.limits.prior[1],mu.limits.prior[2], length.out = 100)) %>%
    mutate(n0 = n.0.t, a0 = alpha.0.t, b0 = beta.0.t, category ="Treatment",
           dens = dnorgam(mu = mu, tau = tau, mu0 = mu.0.t, n0 = n.0.t,
                          a0 = alpha.0.t, b0 = beta.0.t)) %>%
    mutate(
           color = as.numeric(cut((dens),100)),
           group = "Prior",
           density = "Treatment Prior")

  post <- expand.grid(
    tau = seq(tau.limits.post[1],tau.limits.post[2], length.out = 100),
    mu = seq(mu.limits.post[1],mu.limits.post[2], length.out = 100)) %>%
    mutate(n0 = pp$n.n, a0 = pp$alpha.n, b0 = pp$beta.n, category ="Treatment",
           dens = dnorgam(mu = mu, tau = tau, mu0 = pp$mu.n, n0 = pp$n.n,
                          a0 = pp$alpha.n, b0 = pp$beta.n),
           color = as.numeric(cut((dens),100)),
           group = "Posterior",
           density="Treatment Posterior")

  my.df <- rbind(prior,post) %>% mutate(density=factor(density, c("Treatment Prior",
                                                                  "Treatment Posterior")))
  levels(my.df$density) <- c(
    paste0("Treatment Prior: NG(", round(mu.0.t, 2), ", ", round(n.0.t, 2), ", ",
           round(alpha.0.t, 2), ", ", round(beta.0.t, 2), ")"),
    paste0("Treatment Posterior NG(", round(pp$mu.n, 2),", ", round(pp$n.n, 2), ", ",
           round(pp$alpha.n, 2), ", ", round(pp$beta.n, 2), ")" )  )

  my.df$color2 <- as.numeric(cut((my.df$dens),100))

  my.colors <- colorRampPalette(c("black", "red",  "yellow"))(100)
  joint.pdf.prior <- ggplot(data= my.df %>% dplyr::filter(density==levels(density)[1]),
                            aes(x = mu, y = tau, fill = factor(color)))+
    geom_tile() +facet_wrap(~density, scales="free")+guides(fill="none")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0),
                       breaks = seq(0,max(my.df$tau),.25),
                       labels= round((seq(0,max(my.df$tau),.25)) ^ -.5,2),
                       sec.axis = sec_axis(~., name = NULL,
                                           breaks= seq(0,max(my.df$tau),.25)

                                           ))+
    scale_fill_manual(values= my.colors)+
    labs(x="Mean treatment response", y = TeX("standard deviation, $\\sigma$"), caption="")

  joint.pdf.post <- ggplot(data= my.df %>% dplyr::filter(density==levels(density)[2]),
                           aes(x = mu, y = tau, fill = factor(color)))+
    geom_tile() +facet_wrap(~density, scales="free")+guides(fill="none")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0),
                       breaks = seq(0,max(my.df %>% dplyr::filter(density==levels(density)[2]) %>% .$tau),length.out=10),
                       labels= round(seq(0,max(my.df %>% dplyr::filter(density==levels(density)[2]) %>% .$tau), length.out=10) ^ -.5,2),
                       sec.axis = sec_axis(~., name = TeX("precision parameter, $\\tau$"),
                                           breaks=seq(0,max(my.df %>% dplyr::filter(density==levels(density)[2]) %>% .$tau),length.out=10),
                                           labels=round(seq(0,max(my.df %>% dplyr::filter(density==levels(density)[2]) %>% .$tau),length.out=10),2)
                                           ))+
    scale_fill_manual(values= my.colors)+labs(y=NULL, x="Mean treatment response",
                                              caption=TeX(paste0("Treatment data (n, $\\bar{x}$, s): (",
                                                                          round(n.t,2), ", ", round(xbar.t,2), ", ", round(s.t,2),")")))


  marginal.pdf <- rbind(
    gcurve(dt_ls(x,df = 2 * alpha.0.t, mu = mu.0.t,
                 sigma = beta.0.t/(alpha.0.t * n.0.t)) ,
           from = mu.limits[1],
           to = mu.limits[2], n = 1001, category = "Treatment") %>%
      mutate(group="Prior", density="Treatment Prior"),
    gcurve(dt_ls(x,df = 2 * pp$alpha.n, mu = pp$mu.n,
                 sigma = pp$beta.n/(pp$alpha.n * pp$n.n )),
           from = pp$mu.n - 5*pp$beta.n/(pp$alpha.n * pp$n.n ),
           to = pp$mu.n + 5*pp$beta.n/(pp$alpha.n * pp$n.n ), n = 1001, category = "Treatment") %>%
      mutate(group="Posterior", density="Treatment Posterior")) %>%
    mutate(density = factor(density, c("Treatment Prior", "Treatment Posterior")))


  precision.pdf <- rbind(
    gcurve(expr = dgamma(x, shape = alpha.0.t, rate = beta.0.t), from = tau.limits[1], to = tau.limits[2], n = 1001, category = "Treatment" ) %>%
      mutate(group="Prior", density=factor("Treatment precision prior")),
    gcurve(expr = dgamma(x, shape = pp$alpha.n, rate = pp$beta.n), from = tau.limits[1], to = tau.limits[2], n = 1001, category = "Treatment" ) %>%
      mutate(group="Posterior", density=factor("Treatment precision posterior"))

  )


  var.pdf <- rbind(
    gcurve(expr = extraDistr::dinvgamma(x, alpha = alpha.0.t, beta = 1/beta.0.t), from = 1/tau.limits[1], to = 1/tau.limits[2], n = 1001, category = "Treatment" ) %>%
      mutate(group="Prior", density=factor("Treatment variance prior")),
    gcurve(expr = extraDistr::dinvgamma(x, alpha = pp$alpha.n, beta = 1/pp$beta.n), from = 1/tau.limits.post[2], to = 1/tau.limits.post[1], n = 1001, category = "Treatment" ) %>%
      mutate(group="Posterior", density=factor("Treatment variance posterior"))
  )

  var.plot <-    ggplot(data = var.pdf,
         aes(x = x,y = y, color = category))+
          geom_line(size=.75)+ facet_wrap(~density, scales="free")+guides(color="none")+
          theme(legend.position = "bottom")+
          labs(title="Marginal distribution for variance",color = NULL,
               x = TeX("Variance"), y="",
               subtitle = "Variance ~ Inverse-Gamma(shape = alpha, scale = 1/beta)",
               caption="If Z ~ Inverse-Gamma(shape = alpha, scale = 1/beta) then E(Z) = 1/(alpha*beta); Var(Z) ~ 1/((alpha-1)^2(alpha - 2)beta^2)")



  for.labs <- pretty(seq(tau.limits[1],tau.limits[2], length.out = 5), eps.correct=1, n = 8)
  for.labs <- for.labs[for.labs > 0]
  precision.plot <-  ggplot(data = precision.pdf,
                           aes(x = x,y = y, color = category))+
    geom_line(size=.75)+ facet_wrap(~density, scales="free_y")+guides(color="none")+
    theme(legend.position = "bottom")+
    labs(title="Marginal distribution for precision",color = NULL,
         x = TeX("Precision"), y="",
         subtitle="Precision ~ Gamma(alpha, beta)",
         caption="If Z ~ Gamma(alpha,beta) then E(Z) = alpha/beta; Var(Z) ~ alpha/beta^2")

  levels(marginal.pdf$density) <- c(
    paste0("Treatment Prior: t(", 2 * round(alpha.0.t, 2), ", ", round(mu.0.t, 2),
           ", ", round(beta.0.t/(alpha.0.t * n.0.t), 2), ")"),
    paste0("Treatment Posterior: t(", 2 * round(pp$alpha.n), ", ", round(pp$mu.n, 2),
           ", ", round(pp$beta.n/(pp$alpha.n * pp$n.n ), 2), ")")  )

  marginal.plot <-  ggplot(data = marginal.pdf,
                           aes(x = x,y = y, color = category))+
    geom_line(size=.75)+ facet_wrap(~density, scales="free")+guides(color="none")+
    theme(legend.position = "bottom")+
    labs(title="Marginal distribution for the mean",color = NULL,
         x = TeX("Mean treatment response"), y="")



  visualizer.a <-
    bind_rows(
      gcurve(
        dgamma(x, shape = alpha.0.t, rate = beta.0.t), from=0, to=gamma.IG.sd.limits[1], category = paste0("Precision ~ Gamma(shape=",alpha.0.t, ", rate=", beta.0.t, ")")),
      gcurve(
        extraDistr::dinvgamma(x, alpha = alpha.0.t, beta = beta.0.t), from=0, to=gamma.IG.sd.limits[2], category=paste0("Variance ~ IG(shape=", alpha.0.t, ", scale=", beta.0.t, ")")
      )) %>% ggplot(aes(x=x,y=y))+geom_line() + facet_wrap(~category, scales="free") + labs(title="Precision and Variance"
      )


    # sd parameter: sqrt of inverse gamma(alpha, beta)
    visualizer.b <- data.frame(x=sqrt(extraDistr::rinvgamma(n = 10000, alpha = alpha.0.t, beta = beta.0.t))) %>% ggplot(aes(x=x)) +geom_density() + xlim(0,gamma.IG.sd.limits[3]) +
      labs(title="Empirical Density of Standard deviation", subtitle=paste0("Obtained from 10000 draws from IG(shape=", alpha.0.t, ", scale=", beta.0.t, ")")
      )


 list(list(joint.pdf.prior, joint.pdf.post), marginal.plot, precision.plot, var.plot, visualizer.a, visualizer.b, pp[1,c(1:14)], pp[1,c(1,17:22)])

}


