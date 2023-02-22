# source("get.tte.post.param.r")
# source("gcurve.r")
# TwoSampleTTE ------

# m.con.prior = 10
# m.trt.prior = 10
# HR.prior = .845
# ARatio = 1
# HR.obs = .75
# m.obs = 50

#' @title Make Time to event prior/posterior plot
#'
#' @param m.con.prior prior number of control events
#' @param m.trt.prior prior number of treatment events
#' @param HR.prior HR estimate
#' @param ARatio Randomization ratio
#' @param HR.obs HR observed
#' @param m.obs Number of observed events
#'
#' @return A ggplot object is returned
#' @export
#'
#' @examples
#' make.tte.ppp()
make.tte.ppp <- function(m.con.prior = 10, m.trt.prior = 10, HR.prior = .845,
                         ARatio = 1, HR.obs = .75, m.obs = 50){
  ARatio.send <- 1/(1+ARatio)
  post.params = get.tte.post.param(m.con.prior = m.con.prior,
                                   m.trt.prior = m.trt.prior, HR.prior = HR.prior,
                                   ARatio = ARatio.send, HR.obs = HR.obs, m.obs = m.obs)
  pdfs <- rbind(
    gcurve(expr = dnorm(x, mean = log(HR.prior), sd = sqrt(4/(m.con.prior + m.trt.prior))),
           from = log(0.01), to = log(3), n = 1001, category = "Hazard ratio Prior") %>%
      mutate(HR = HR.prior, events = m.con.prior + m.trt.prior, group="Prior"),
    gcurve(expr = dnorm(x, mean = post.params[1,1], sd = post.params[1,2]),
           from = log(0.01), to = log(3), n = 1001, category="Hazard ratio posterior") %>%
      mutate(HR = post.params[1,1], events = m.con.prior + m.trt.prior + m.obs[2],
             group="Posterior")
  ) %>% mutate(group=factor(group, c("Prior", "Posterior")))
  levels(pdfs$group) <- c(paste0("Prior for HR: ", HR.prior, " worth ",
                                 m.trt.prior+ m.con.prior, " events"),
                          paste0("Posterior for HR: ", round(exp(post.params[1]),4),
                                 " worth ", m.trt.prior+ m.con.prior + m.obs, " events"))

  p1 <- ggplot(data = pdfs, aes(x = exp(x),y = y, color = category))+geom_line(size=.75) +
    facet_wrap(~group)+guides(color = "none")+
    labs(x = "Hazard Ratio", y=NULL,
         title="Prior and posterior distribtuions for the hazard ratio",
         subtitle = paste0("Data: Observed hazard ratio: ", round(HR.obs,2)," based on ",
                           m.obs, " events."),
         color="Density") +
    theme(legend.position = "bottom")+
    scale_x_continuous(breaks = seq(0,3,.2))+
    scale_y_continuous(breaks=NULL, labels = NULL)

  t1 <- data.frame("Prior control events" = m.con.prior, 'Prior treatment events' = m.trt.prior, prior.mean =  round(log(HR.prior),4), 'exp(Prior Mean)' = HR.prior,
             'Randomization Ratio' = ARatio, 'Obs HR' = HR.obs, 'Obs Events' = m.obs) %>%
    mutate(
      post.mean = round(post.params$post.mean,4),
      post.var= round(post.params$post.sd^2,4),
      exp.post.mean = round(exp(post.params$post.mean),4),
           exp.post.var = round(exp(post.params$post.sd^2),4))

  list(p1, t1)
}

# holdit <- make.tte.ppp()
# holdit[[1]]
# holdit[[2]]
