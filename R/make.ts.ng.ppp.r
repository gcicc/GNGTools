#' @title Make two-sample normal-gamma prior/posterior plot
#'
#' @param mu.0.c prior mean for control group
#' @param alpha.0.c prior alpha parameter for control group
#' @param beta.0.c prior beta parameter for control group
#' @param n.0.c prior effective sample size parameter for control group
#' @param mu.0.t prior mean for treatment group
#' @param alpha.0.t prior alpha parameter for treatment group
#' @param beta.0.t prior beta parameter for treatment group
#' @param n.0.t prior effective sample size parameter for treatment group
#' @param xbar.c control mean
#' @param s.c control sd
#' @param n.c control sample size
#' @param xbar.t treatment mean
#' @param s.t treatment sd
#' @param n.t treatment sample size
#' @param limits limits for visualizing precision, variance, standard deviation
#'
#' @return a ggplot object is returned
#' @export
#'
#' @examples
#' my.ts.ng.ppp <- make.ts.ng.ppp()
#' my.ts.ng.ppp[[1]][[1]]
#' my.ts.ng.ppp[[1]][[2]]
#' my.ts.ng.ppp[[1]][[3]]
#' my.ts.ng.ppp[[1]][[4]]
#' gridExtra::grid.arrange(my.ts.ng.ppp[[1]][[1]], my.ts.ng.ppp[[1]][[3]],
#' my.ts.ng.ppp[[1]][[2]], my.ts.ng.ppp[[1]][[4]])
#' my.ts.ng.ppp[[2]]
#' my.ts.ng.ppp[[3]]
#' my.ts.ng.ppp[[4]]
#' my.ts.ng.ppp[[5]]
#' my.ts.ng.ppp[[6]]
#' my.ts.ng.ppp[[7]]

make.ts.ng.ppp <- function(mu.0.c = 0, alpha.0.c = 2.5, beta.0.c = 10, n.0.c = 10,
         mu.0.t = 0, alpha.0.t = .25, beta.0.t = 1, n.0.t = .0001,
         xbar.c = .25, s.c = 1.5, n.c = 55,
         xbar.t = 2.5, s.t = 1.5, n.t = 55,
         limits=c(5, 25, 10)
) {


        pp.t <-   get.ng.post(mu.0 = mu.0.t, n.0 = n.0.t, alpha.0 = alpha.0.t, beta.0 = beta.0.t,
                              xbar = xbar.t, s = s.t, n = n.t, group = "Treatment") %>% dplyr::select(group, everything()) %>%
                mutate(prior.ng.mean = paste0("(",round(mu.0,4), ", ", round(alpha.0/beta.0,4), ")"),
                       prior.mode = paste0("(",round(mu.0,4), ", ", ifelse(alpha.0 < 0.5, NA, round((alpha.0-1/2)/beta.0,4)), ")"),
                       prior.ng.var = paste0("(", round(beta.0/(n.0.t*(alpha.0-1)),4), ", ", round(alpha.0/beta.0^2,4), ")"),
                       post.ng.mean = paste0("(", round(mu.n,4), ", ", round(alpha.n/beta.n,4), ")"),
                       post.mode = paste0("(",round(mu.n,4), ", ", ifelse(alpha.n < 0.5, NA, round((alpha.n-1/2)/beta.n,4)), ")"),
                       post.ng.var = paste0("(", round(beta.n/(n.n*(alpha.n-1)),4), ", ", round(alpha.n/beta.n^2,4), ")"))

        pp.c <-   get.ng.post(mu.0 = mu.0.c, n.0 = n.0.c, alpha.0 = alpha.0.c, beta.0 = beta.0.c,
                              xbar = xbar.c, s = s.c, n = n.c, group = "Control") %>% dplyr::select(group, everything()) %>%
                mutate(prior.ng.mean = paste0("(",round(mu.0,4), ", ", round(alpha.0/beta.0,4), ")"),
                       prior.mode = paste0("(",round(mu.0,4), ", ", ifelse(alpha.0 < 0.5, NA, round((alpha.0-1/2)/beta.0,4)), ")"),
                       prior.ng.var = paste0("(", round(beta.0/(n.0.t*(alpha.0-1)),4), ", ", round(alpha.0/beta.0^2,4), ")"),
                       post.ng.mean = paste0("(", round(mu.n,4), ", ", round(alpha.n/beta.n,4), ")"),
                       post.mode = paste0("(",round(mu.n,4), ", ", ifelse(alpha.n < 0.5, NA, round((alpha.n-1/2)/beta.n,4)), ")"),
                       post.ng.var = paste0("(", round(beta.n/(n.n*(alpha.n-1)),4), ", ", round(alpha.n/beta.n^2,4), ")"))


        pp <-   rbind(pp.c, pp.t) %>% dplyr::select(group, everything()) %>% rename(Group=group)


        tau.limits <- c(min(pmax(.1, qgamma(p = .005, shape = alpha.0.t, rate = beta.0.t)),
                            pmax(.1, qgamma(p = .005, shape = pp.t$alpha.n, rate = pp.t$beta.n)),
                            pmax(.1, qgamma(p = .005, shape = alpha.0.c, rate = beta.0.c)),
                            pmax(.1, qgamma(p = .005, shape = pp.c$alpha.n, rate = pp.c$beta.n))),
                        max(qgamma(p = .9995,shape = alpha.0.t, rate = beta.0.t),
                            qgamma(p = .9995,shape = pp.t$alpha.n, rate = pp.t$beta.n),
                            qgamma(p = .9995,shape = alpha.0.c, rate = beta.0.c),
                            qgamma(p = .9995,shape = pp.c$alpha.n, rate = pp.c$beta.n)
                        ))

        tau.limits.prior <- c(min(pmax(.1, qgamma(p = .005, shape = alpha.0.t, rate = beta.0.t)),
                                  pmax(.1, qgamma(p = .005, shape = alpha.0.c, rate = beta.0.c))),
                              max(qgamma(p = .9995,shape = alpha.0.t, rate = beta.0.t),
                                  qgamma(p = .9995,shape = alpha.0.c, rate = beta.0.c))
        )
        tau.limits.post <- c(min(
                pmax(.05, qgamma(p = .0005, shape = pp.t$alpha.n, rate = pp.t$beta.n)),
                pmax(.05, qgamma(p = .0005, shape = pp.c$alpha.n, rate = pp.c$beta.n))),
                max(
                        qgamma(p = .9995,shape = pp.t$alpha.n, rate = pp.t$beta.n),
                        qgamma(p = .9995,shape = pp.c$alpha.n, rate = pp.c$beta.n)
                ))

        mu.limits <- c(min(qnorm(p = .005, mean = mu.0.t, sd = (alpha.0.t / beta.0.t) ^ (-1)/sqrt(n.0.t)),
                           qnorm(p = .005, mean = pp.t$mu.n, sd = (pp.t$alpha.n / pp.t$beta.n) ^ (-1)/sqrt(pp.t$n.n)),
                           qnorm(p = .005, mean = mu.0.c, sd = (alpha.0.c / beta.0.c) ^ (-1)/sqrt(n.0.c)),
                           qnorm(p = .005, mean = pp.c$mu.n, sd = (pp.c$alpha.n / pp.c$beta.n) ^ (-1)/sqrt(pp.c$n.n))
        ),
        max(qnorm(p = .99995, mean = mu.0.t, sd = (alpha.0.t / beta.0.t) ^ (-1)/sqrt(n.0.t)),
            qnorm(p = .99995, mean = pp.t$mu.n, sd = (pp.t$alpha.n / pp.t$beta.n) ^ (-1)/sqrt(pp.t$n.n)),
            qnorm(p = .99995, mean = mu.0.c, sd = (alpha.0.c / beta.0.c) ^ (-1)/sqrt(n.0.c)),
            qnorm(p = .99995, mean = pp.c$mu.n, sd = (pp.c$alpha.n / pp.c$beta.n) ^ (-1)/sqrt(pp.c$n.n))
        ))


        mu.limits.prior <- c(min(qnorm(p = .005, mean = mu.0.t, sd = (alpha.0.t / beta.0.t) ^ (-1)/sqrt(n.0.t)),

                                 qnorm(p = .005, mean = mu.0.c, sd = (alpha.0.c / beta.0.c) ^ (-1)/sqrt(n.0.c))),
                             max(qnorm(p = .9995, mean = mu.0.t, sd = (alpha.0.t / beta.0.t) ^ (-1)/sqrt(n.0.t)),
                                 qnorm(p = .9995, mean = mu.0.c, sd = (alpha.0.c / beta.0.c) ^ (-1)/sqrt(n.0.c)))
        )

        mu.limits.prior.control <- c(min(
                qnorm(p = .005, mean = mu.0.c, sd = (alpha.0.c / beta.0.c) ^ (-1)/sqrt(n.0.c))),
                max(                           qnorm(p = .9995, mean = mu.0.c, sd = (alpha.0.c / beta.0.c) ^ (-1)/sqrt(n.0.c)))
        )

        mu.limits.prior.treatment <- c(min(qnorm(p = .005, mean = mu.0.t, sd = (alpha.0.t / beta.0.t) ^ (-1)/sqrt(n.0.t))),
                                       max(qnorm(p = .9995, mean = mu.0.t, sd = (alpha.0.t / beta.0.t) ^ (-1)/sqrt(n.0.t)))
        )

        mu.limits.post <- c(min(
                qnorm(p = .005, mean = pp.t$mu.n, sd = (pp.t$alpha.n / pp.t$beta.n) ^ (-1)/sqrt(pp.t$n.n)),
                qnorm(p = .005, mean = pp.c$mu.n, sd = (pp.c$alpha.n / pp.c$beta.n) ^ (-1)/sqrt(pp.c$n.n))),
                max(
                        qnorm(p = .9995, mean = pp.t$mu.n, sd = (pp.t$alpha.n / pp.t$beta.n) ^ (-1)/sqrt(pp.t$n.n)),
                        qnorm(p = .9995, mean = pp.c$mu.n, sd = (pp.c$alpha.n / pp.c$beta.n) ^ (-1)/sqrt(pp.c$n.n))
                ))

        mu.limits.post.control <- c(min(
                qnorm(p = .005, mean = pp.c$mu.n, sd = (pp.c$alpha.n / pp.c$beta.n) ^ (-1)/sqrt(pp.c$n.n))),
                max(
                        qnorm(p = .9995, mean = pp.c$mu.n, sd = (pp.c$alpha.n / pp.c$beta.n) ^ (-1)/sqrt(pp.c$n.n))
                ))

        mu.limits.post.treatment <- c(min(
                qnorm(p = .005, mean = pp.t$mu.n, sd = (pp.t$alpha.n / pp.t$beta.n) ^ (-1)/sqrt(pp.t$n.n))),
                max(
                        qnorm(p = .9995, mean = pp.t$mu.n, sd = (pp.t$alpha.n / pp.t$beta.n) ^ (-1)/sqrt(pp.t$n.n)))
        )

        CON.prior <- expand.grid(
                tau = seq(tau.limits.prior[1],tau.limits.prior[2], length.out = 100),
                mu = seq(mu.limits.prior.control[1],mu.limits.prior.control[2], length.out = 100)) %>%
                mutate(n0 = pp.c$n.0, a0 = pp.c$alpha.0, b0 = pp.c$beta.0, category="Control", group="Prior", density="Control Prior",
                       dens = dnorgam(mu = mu, tau = tau, mu0 = pp.c$mu.0, n0 = pp.c$n.0, a0 = pp.c$alpha.0, b0 = pp.c$beta.0)) %>%
                mutate(color= as.numeric(cut((dens),150)))

        CON.post <- expand.grid(
                tau = seq(tau.limits.post[1],tau.limits.post[2], length.out = 100),
                mu = seq(mu.limits.post[1],mu.limits.post[2], length.out = 100)) %>%
                mutate(n0 = pp.c$n.n, a0 = pp.c$alpha.n, b0 = pp.c$beta.n, category="Control", group="Posterior", density="Control Posterior",
                       dens = dnorgam(mu = mu, tau = tau, mu0 = pp.c$mu.n, n0 = pp.c$n.n, a0 = pp.c$alpha.n, b0 = pp.c$beta.n)) %>%
                mutate(color= as.numeric(cut((dens),150)))

        TRT.prior <- expand.grid(
                tau = seq(tau.limits.prior[1],tau.limits.prior[2], length.out = 100),
                mu = seq(mu.limits.prior.treatment[1],mu.limits.prior.treatment[2], length.out = 100)) %>%
                mutate(n0 = pp.t$n.0, a0 = pp.t$alpha.0, b0 = pp.t$beta.0, category="Treatment", group="Prior", density="Treatment Prior",
                       dens = dnorgam(mu = mu, tau = tau, mu0 = pp.t$mu.0, n0 = pp.t$n.0, a0 = pp.t$alpha.0, b0 = pp.t$beta.0))%>%
                mutate(color= as.numeric(cut((dens),150)))

        # tau.limits.post <- c(.2, .05)
        TRT.post <- expand.grid(
                tau = seq(tau.limits.post[1],tau.limits.post[2], length.out = 100),
                mu = seq(mu.limits.post[1],mu.limits.post[2], length.out = 100)) %>%
                mutate(n0 = pp.t$n.n, a0 = pp.t$alpha.n, b0 = pp.t$beta.n, category="Treatment", group="Posterior", density="Treatment Posterior",
                       dens = dnorgam(mu = mu, tau = tau, mu0 = pp.t$mu.n, n0 = pp.t$n.n, a0 = pp.t$alpha.n, b0 = pp.t$beta.n))%>%
                mutate(color= as.numeric(cut((dens),150)))

        my.df <- rbind(CON.prior,CON.post,TRT.prior, TRT.post)
        # my.df$color <-   as.numeric(cut((my.df$dens),150))
        my.colors <- colorRampPalette(c("black", "red",  "yellow"))(150)

        my.df$density.original <- factor(my.df$density, c("Control Prior", "Control Posterior", "Treatment Prior", "Treatment Posterior"))
        my.df$density <- my.df$density.orginal
        levels(my.df$density) <- c(
                paste0("Control Prior: NG(", round(mu.0.c,2), ", ", round(n.0.c,2), ", ", round(alpha.0.c,2), ", ", round(beta.0.c,2), ")"),
                paste0("Control Posterior: NG(", round(pp.c$mu.n), ", ", round(pp.c$n.n), ", ", round(pp.c$alpha.n), ", ", round(pp.c$beta.n), ")"),
                paste0("Treatment Prior: NG(", round(mu.0.t,2), ", ", round(n.0.t,2), ", ", round(alpha.0.t,2), ", ", round(beta.0.t,2), ")"),
                paste0("Treatment Posterior: NG(", round(pp.t$mu.n), ", ", round(pp.t$n.n), ", ", round(pp.t$alpha.n), ", ", round(pp.t$beta.n), ")")
        )




        joint.pdf.prior.1 <- ggplot(data= my.df %>% dplyr::filter(density.original %in% c("Control Prior")),
                                    aes(x = mu, y = tau, fill = factor(color)))+
                geom_tile() + facet_wrap(~density,nrow=2, scales="free")+
                scale_x_continuous(expand = c(0,0))+
                scale_y_continuous(expand = c(0,0),
                                   breaks = pretty(c(min(my.df$tau),max(my.df$tau)), n=10), labels= round(pretty(c(min(my.df$tau),max(my.df$tau)), n=10)^ -.5,2),
                                   sec.axis = sec_axis(~., name = NULL))+
                scale_fill_manual(values= my.colors)+
                guides(fill="none")+
                labs(x = "Mean response",
                     y = TeX("Response standard deviation, $\\sigma$")
                     # title="Prior and posterior joint distributions for mean and precision parameters",
                     # subtitle=TeX(paste0("Control Data (n, $\\bar{x}$, s): (", round(n.c,2), ", ", round(xbar.c,2), ", ", round(s.c,2),"), ",
                     #                     "Treatment Data: (", round(n.t,2), ", ", round(xbar.t,2), ", ", round(s.t,2),"). Observed treatment effect: $",
                     #                     round(xbar.t,2) - round(xbar.c,2), "$.")
                     #              )

                )
        joint.pdf.prior.2 <- ggplot(data= my.df %>% dplyr::filter(density.original %in% c("Treatment Prior")),
                                    aes(x = mu, y = tau, fill = factor(color)))+
                geom_tile() + facet_wrap(~density,nrow=2, scales="free")+
                scale_x_continuous(expand = c(0,0))+
                scale_y_continuous(expand = c(0,0),
                                   breaks = pretty(c(min(my.df$tau),max(my.df$tau)), n=10), labels= round(pretty(c(min(my.df$tau),max(my.df$tau)), n=10)^ -.5,2),
                                   sec.axis = sec_axis(~., name = NULL))+
                scale_fill_manual(values= my.colors)+
                guides(fill="none")+
                labs(x = "Mean response",
                     y = TeX("Response standard deviation, $\\sigma$"),
                     caption=""
                     # title="Prior and posterior joint distributions for mean and precision parameters",
                     # subtitle=TeX(paste0("Control Data (n, $\\bar{x}$, s): (", round(n.c,2), ", ", round(xbar.c,2), ", ", round(s.c,2),"), ",
                     #                     "Treatment Data: (", round(n.t,2), ", ", round(xbar.t,2), ", ", round(s.t,2),"). Observed treatment effect: $",
                     #                     round(xbar.t,2) - round(xbar.c,2), "$.")
                     #              )

                )


        joint.pdf.post.1 <- ggplot(data= my.df %>% dplyr::filter(density.original %in% c("Control Posterior")),
                                   aes(x = mu, y = tau, fill = factor(color)))+
                geom_tile() + facet_wrap(~density,nrow=2)+
                scale_x_continuous(expand = c(0,0))+
                scale_y_continuous(expand = c(0,0),
                                   breaks = seq(from=min(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau),
                                                to=max(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau), length.out=10),
                                   labels= round(seq(from=min(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau),
                                                     to=max(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau), length.out=10)^ -.5,2),
                                   sec.axis = sec_axis(~., name = TeX("Response precision, $\\tau$"),
                                                       breaks=seq(from=min(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau),
                                                                  to=max(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau), length.out=10),
                                                       labels=round(seq(from=min(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau),
                                                                        to=max(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau), length.out=10),2)))+
                scale_fill_manual(values= my.colors)+
                guides(fill="none")+
                labs(x = "Mean response",
                     y = NULL
                     #title="Prior and posterior joint distributions for mean and precision parameters",
                     # caption=TeX(paste0("Control Data (n, $\\bar{x}$, s): (", round(n.c,2), ", ", round(xbar.c,2), ", ", round(s.c,2),"), ",
                     #                    "Treatment Data: (", round(n.t,2), ", ", round(xbar.t,2), ", ", round(s.t,2),"). Observed treatment effect: $",
                     #                    round(xbar.t,2) - round(xbar.c,2), "$."))
                )

        joint.pdf.post.2 <- ggplot(data= my.df %>% dplyr::filter(density.original %in% c("Treatment Posterior")),
                                   aes(x = mu, y = tau, fill = factor(color)))+
                geom_tile() + facet_wrap(~density,nrow=2)+
                scale_x_continuous(expand = c(0,0))+
                scale_y_continuous(expand = c(0,0),
                                   breaks = seq(from=min(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau),
                                                to=max(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau), length.out=10),
                                   labels= round(seq(from=min(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau),
                                                     to=max(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau), length.out=10)^ -.5,2),
                                   sec.axis = sec_axis(~., name = TeX("Response precision, $\\tau$"),
                                                       breaks=seq(from=min(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau),
                                                                  to=max(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau), length.out=10),
                                                       labels=round(seq(from=min(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau),
                                                                        to=max(my.df %>% dplyr::filter(density.original %in% c("Control Posterior", "Treatment Posterior")) %>% .$tau), length.out=10),2)))+
                scale_fill_manual(values= my.colors)+
                guides(fill="none")+
                labs(x = "Mean response",
                     y = NULL,
                     #title="Prior and posterior joint distributions for mean and precision parameters",
                     caption=TeX(paste0("Control Data (n, $\\bar{x}$, s): (", round(n.c,2), ", ", round(xbar.c,2), ", ", round(s.c,2),"), ",
                                        "Treatment Data: (", round(n.t,2), ", ", round(xbar.t,2), ", ", round(s.t,2),"). Observed treatment effect: $",
                                        round(xbar.t,2) - round(xbar.c,2), "$."))
                )





        marginal.densities <- rbind(
                gcurve(dt_ls(x,df = 2 * alpha.0.c, mu = mu.0.c, sigma = beta.0.c/(alpha.0.c * n.0.c)),
                       from = mu.limits.prior.control[1],
                       to = mu.limits.prior.control[2], n = 1001, category = "Control") %>% mutate(group="Prior", density="Control Prior"),
                gcurve(dt_ls(x,df = 2 * alpha.0.t, mu = mu.0.t, sigma = beta.0.t/(alpha.0.t * n.0.t)) ,
                       from = mu.limits.prior.treatment[1],
                       to = mu.limits.prior.treatment[2], n = 1001, category = "Treatment") %>% mutate(group="Prior", density="Treatment Prior"),
                gcurve(dt_ls(x,df = 2 * pp.c$alpha.n, mu = pp.c$mu.n, sigma = pp.c$beta.n/(pp.c$alpha.n * pp.c$ n.n )),
                       from = mu.limits.post[1],
                       to = mu.limits.post[2], n = 1001, category = "Control") %>% mutate(group="Posterior", density="Control Posterior"),
                gcurve(dt_ls(x,df = 2 * pp.t$alpha.n, mu = pp.t$mu.n, sigma = pp.t$beta.n/(pp.t$alpha.n * pp.t$n.n )),
                       from = mu.limits.post[1],
                       to = mu.limits.post[2], n = 1001, category = "Treatment") %>% mutate(group="Posterior", density="Treatment Posterior")) %>%
                mutate(density=factor(density, c("Control Prior", "Control Posterior", "Treatment Prior", "Treatment Posterior")))



        levels(marginal.densities$density) <- c(
                paste0("Control Prior: t(", round(2 * alpha.0.c, 2), ", ", round(mu.0.c, 2), ", ", round(beta.0.c/(alpha.0.c * n.0.c),2), ")"),
                paste0("Control Posterior: t(",  round(2 * pp.c$alpha.n, 2), ", ", round(pp.c$mu.n, 2), ", ", round(pp.c$beta.n/(pp.c$alpha.n * pp.c$ n.n ),2), ")"),
                paste0("Treatment Prior t(",  round(2 * alpha.0.t, 2), ", ", round(mu.0.t, 2), ", ", round(beta.0.t/(alpha.0.t * n.0.t),2), ")"),
                paste0("Treatment Posterior: t(",  round(2 * pp.t$alpha.n, 2), ", ", round(pp.t$mu.n, 2), ", ", round(pp.t$beta.n/(pp.t$alpha.n * pp.t$ n.n ),2), ")")
        )

        marginal.plot <- ggplot(data=marginal.densities, aes(x = x,y = y, color = category))+
                geom_line(size=.75)+
                theme(legend.position = "bottom")+
                labs(title="Marginal distributions for the means",color = NULL,
                     x = "Mean response", y="")+facet_wrap(~density, scales="free_x")+
                scale_y_continuous(breaks=NULL, labels = NULL) + guides(color="none")


        precision.pdf <- rbind(
                gcurve(dgamma(x, shape = pp.c$alpha.0, rate = pp.c$beta.0),
                       from = tau.limits[1],
                       to = tau.limits[2], n = 1001, category = "Control") %>%
                        mutate(group="Prior", density="Control Prior"),
                gcurve(dgamma(x, shape = pp.t$alpha.0, rate = pp.t$beta.0) ,
                       from = tau.limits[1],
                       to = tau.limits[2], n = 1001, category = "Treatment") %>%
                        mutate(group="Prior", density="Treatment Prior"),
                gcurve(dgamma(x, shape = pp.c$alpha.n, rate = pp.c$beta.n),
                       from = tau.limits.post[1],
                       to = tau.limits.post[2], n = 1001, category = "Control") %>%
                        mutate(group="Posterior", density="Control Posterior"),
                gcurve(dgamma(x, shape = pp.t$alpha.n, rate = pp.t$beta.n),
                       from = tau.limits.post[1],
                       to = tau.limits.post[2], n = 1001, category = "Treatment") %>%
                        mutate(group="Posterior", density="Treatment Posterior")) %>%
                mutate(density=factor(density, c("Control Prior", "Control Posterior", "Treatment Prior", "Treatment Posterior")))



        precision.plot <-  ggplot(data = precision.pdf,
                                  aes(x = x,y = y, color = category))+
                geom_line(size=.75)+ facet_wrap(~density, scales="free_y")+guides(color="none")+
                theme(legend.position = "bottom")+
                labs(title="Marginal distribution for precision",color = NULL,
                     x = TeX("Precision"), y="",
                     subtitle="Precision ~ Gamma(alpha, beta)",
                     caption="If Z ~ Gamma(alpha,beta) then E(Z) = alpha/beta; Var(Z) ~ alpha/beta^2")+
                # scale_y_continuous(breaks=NULL)+
                scale_x_continuous(
                        breaks = pretty(seq(tau.limits[1],tau.limits[2], length.out = 5), eps.correct=1, n = 8),
                        sec.axis = sec_axis( ~ . , name = "Standard deviation scale",
                                             breaks=pretty(seq(tau.limits[1],tau.limits[2], length.out = 5), eps.correct=1, n = 8),
                                             labels=round(1/sqrt(pretty(seq(tau.limits[1],tau.limits[2], length.out = 5), eps.correct=1, n = 8)),3)
                        ))




        # from Single sample
        var.pdf <- rbind(
                gcurve(expr = extraDistr::dinvgamma(x, alpha = alpha.0.c, beta = 1/beta.0.c), from = 1/tau.limits[1], to = 1/tau.limits[2], n = 1001, category = "Treatment" ) %>%
                        mutate(group="Prior", density=factor("Control arm variance prior")),
                gcurve(expr = extraDistr::dinvgamma(x, alpha = pp.c$alpha.n, beta = 1/pp.c$beta.n), from = 1/tau.limits.post[2], to = 1/tau.limits.post[1], n = 1001, category = "Treatment" ) %>%
                        mutate(group="Posterior", density=factor("Control arm variance posterior"))
        )

        var.plot <-  ggplot(data = var.pdf,
                            aes(x = x,y = y, color = category))+
                geom_line(size=.75)+ facet_wrap(~density, scales="free")+guides(color="none")+
                theme(legend.position = "bottom")+
                labs(title="Marginal distribution for variance for Control Arm",color = NULL,
                     x = TeX("Variance"), y="",
                     subtitle = "Variance ~ Inverse-Gamma(shape = alpha, scale = 1/beta)",
                     caption="If Z ~ Inverse-Gamma(shape = alpha, scale = 1/beta) then E(Z) = 1/(alpha*beta); Var(Z) ~ 1/((alpha-1)^2(alpha - 2)beta^2)")+
                # scale_y_continuous(breaks=NULL)+
                scale_x_continuous(
                        breaks = pretty(seq(1/tau.limits[1],1/tau.limits[2], length.out = 5), eps.correct=1, n = 8)
                )


        for.labs <- pretty(seq(tau.limits[1],tau.limits[2], length.out = 5), eps.correct=1, n = 8)
        for.labs <- for.labs[for.labs > 0]
        precision.plot <-  ggplot(data = precision.pdf,
                                  aes(x = x,y = y, color = category))+
                geom_line(size=.75)+ facet_wrap(~density, scales="free_y")+guides(color="none")+
                theme(legend.position = "bottom")+
                labs(title="Marginal distribution for precision",color = NULL,
                     x = TeX("Precision"), y="",
                     subtitle="Precision ~ Gamma(alpha, beta)",
                     caption="If Z ~ Gamma(alpha,beta) then E(Z) = alpha/beta; Var(Z) ~ alpha/beta^2")+
                # scale_y_continuous(breaks=NULL)+
                scale_x_continuous(
                        breaks = pretty(seq(tau.limits[1],tau.limits[2], length.out = 5), eps.correct=1, n = 8),
                        sec.axis = sec_axis( ~ . , name = "Standard deviation scale",
                                             breaks=c(tau.limits[1], for.labs),
                                             labels=round(c(1/tau.limits[1], 1/sqrt(for.labs)),3)
                        ))

        # levels(marginal.pdf$density) <- c(
        #   paste0("Treatment Prior: t(", 2 * round(alpha.0.t, 2), ", ", round(mu.0.t, 2),
        #          ", ", round(beta.0.t/(alpha.0.t * n.0.t), 2), ")"),
        #   paste0("Treatment Posterior: t(", 2 * round(pp$alpha.n), ", ", round(pp$mu.n, 2),
        #          ", ", round(pp$beta.n/(pp$alpha.n * pp$n.n ), 2), ")")  )

        # marginal.plot <-  ggplot(data = marginal.pdf,
        #                          aes(x = x,y = y, color = category))+
        #   geom_line(size=.75)+ facet_wrap(~density, scales="free_y")+guides(color="none")+
        #   theme(legend.position = "bottom")+
        #   labs(title="Marginal distribution for the mean",color = NULL,
        #        x = TeX("Mean treatment response"), y="")

        visualizer.a <-
                bind_rows(
                        gcurve(
                                dgamma(x, shape = alpha.0.c, rate = beta.0.c), from=0, to=limits[1], category = paste0("Precision ~ Gamma(shape=", round(alpha.0.c,4), ", rate=", round(beta.0.c,4), ") for Control Arm")),
                        gcurve(
                                extraDistr::dinvgamma(x, alpha = alpha.0.c, beta = beta.0.c), from=0, to=limits[2], category=paste0("Variance ~ IG(shape=", round(alpha.0.c,4), ", scale=", round(beta.0.c,4), ") for Control Arm")
                        )) %>% ggplot(aes(x=x,y=y))+geom_line() + facet_wrap(~category, scales="free") + labs(title="Precision and Variance for Control Arm"
                        )


        # sd parameter: sqrt of inverse gamma(alpha, beta)
        visualizer.b <- data.frame(x=sqrt(extraDistr::rinvgamma(n = 10000, alpha = alpha.0.c, beta = beta.0.c))) %>% ggplot(aes(x=x)) +geom_density() + xlim(0,limits[3]) +
                labs(title="Empirical Density of Standard deviation for Control Arm", subtitle=paste0("Obtained from 10000 draws from IG(shape=",  round(alpha.0.c,4), ", scale=",  round(beta.0.c,4), ")")
                )

        # Old two-sample




        # grid.arrange(joint.pdf.prior, joint.pdf.post, nrow=1)
        list(list(joint.pdf.prior.1, joint.pdf.prior.2, joint.pdf.post.1, joint.pdf.post.2), marginal.plot, precision.plot, var.plot, visualizer.a, visualizer.b,  pp[,c(1:16)], pp[,c(1, 17:22)])





}
