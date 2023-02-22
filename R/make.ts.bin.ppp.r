#' @title Make two-sample binary prior/posterior plot
#'
#' @param a.con prior alpha parameter for control group
#' @param b.con prior beta parameter for control group
#' @param n.con number of trials for control group
#' @param x.con number of responses for control group
#' @param a.trt prior alpha parameter for treatment group
#' @param b.trt prior alpha parameter for treatment group
#' @param n.trt number of trials for treatment group
#' @param x.trt number of responses for treatment group
#'
#' @return A list is returned.
#' @export
#'
#' @examples
#' my.ts.bin.ppp <- make.ts.bin.ppp()
#' my.ts.bin.ppp[[1]]
#' my.ts.bin.ppp[[2]]
make.ts.bin.ppp <- function (a.con = 1, b.con = 1, n.con = 40, x.con = 9,
                             a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 17)
{

        ppp.df.c <- data.frame(alpha.prior = a.con, beta.prior = b.con) %>% mutate(
                prior.mean.c=a.con/(a.con+b.con),
                n.trials = n.con, n.responses = x.con,
                alpha.post = a.con + x.con, beta.post = b.con + (n.con-x.con),
                post.mean.c=(a.con + x.con)/(a.con + x.con+b.con + (n.con-x.con)))
        ppp.df <- bind_rows(
                data.frame(alpha.prior = a.con, beta.prior = b.con) %>% mutate(
                        prior.mean=a.con/(a.con+b.con),
                        prior.mode = case_when(a.con > 1 & b.con > 1 ~ as.character(round((a.con - 1)/(a.con + b.con-2), 4)),
                                               a.con==1 & b.con==1 ~ "any value in (0, 1)", # really any value in (0, 1)
                                               a.con < 1 & b.con < 1 ~ "bimodal {0, 1}", # bimodal (0, 1)
                                               a.con <=1 & b.con > 1 ~ "0",
                                               a.con > 1 & b.con <=1 ~ "1"),
                        prior.var = a.con*b.con/((a.con+b.con)^2*(a.con+b.con+1)),
                        n.trials = n.con, n.responses = x.con,
                        alpha.post = a.con + x.con, beta.post = b.con + (n.con-x.con),
                        post.mean=(a.con + x.con)/(a.con + x.con+b.con + (n.con-x.con)),
                        post.mode =   case_when((a.con + x.con) > 1 & (b.con + (n.con-x.con)) > 1 ~ as.character(round(((a.con + x.con) - 1)/((a.con + x.con) + (b.con + (n.con-x.con)-2)), 4)),
                                                (a.con + x.con)==1 & (b.con + (n.con-x.con))==1 ~ "any value in (0, 1)", # really any value in (0, 1)
                                                (a.con + x.con) < 1 & (b.con + (n.con-x.con)) < 1 ~ "bimodal {0, 1}", # bimodal (0, 1)
                                                (a.con + x.con) <=1 & (b.con + (n.con-x.con)) > 1 ~ "0",
                                                (a.con + x.con) > 1 & (b.con + (n.con-x.con)) <=1 ~ "1"),
                        post.var = (a.con + x.con)*(b.con + (n.con-x.con))/(((a.con + x.con)+(b.con + (n.con-x.con)))^2*((a.con + x.con)+(b.con + (n.con-x.con))+1)),
                        Group="Control") ,
                data.frame(alpha.prior = a.trt, beta.prior = b.trt) %>% mutate(
                        prior.mean=a.trt/(a.trt+b.trt),
                        prior.mode = case_when(a.trt > 1 & b.trt > 1 ~ as.character(round((a.trt - 1)/(a.trt + b.trt-2), 4)),
                                               a.trt==1 & b.trt==1 ~ "any value in (0, 1)", # really any value in (0, 1)
                                               a.trt < 1 & b.trt < 1 ~ "bimodal {0, 1}", # bimodal (0, 1)
                                               a.trt <=1 & b.trt > 1 ~ "0",
                                               a.trt > 1 & b.trt <=1 ~ "1"),
                        prior.var = a.trt*b.trt/((a.trt+b.trt)^2*(a.trt+b.trt+1)),
                        n.trials = n.trt, n.responses = x.trt,
                        alpha.post = a.trt + x.trt, beta.post = b.trt + (n.trt-x.trt),
                        post.mean=(a.trt + x.trt)/(a.trt + x.trt+b.trt + (n.trt-x.trt)),
                        post.mode =   case_when((a.trt + x.trt) > 1 & (b.trt + (n.trt-x.trt)) > 1 ~ as.character(round(((a.trt + x.trt) - 1)/((a.trt + x.trt) + (b.trt + (n.trt-x.trt)-2)), 4)),
                                                (a.trt + x.trt)==1 & (b.trt + (n.trt-x.trt))==1 ~ "any value in (0, 1)", # really any value in (0, 1)
                                                (a.trt + x.trt) < 1 & (b.trt + (n.trt-x.trt)) < 1 ~ "bimodal {0, 1}", # bimodal (0, 1)
                                                (a.trt + x.trt) <=1 & (b.trt + (n.trt-x.trt)) > 1 ~ "0",
                                                (a.trt + x.trt) > 1 & (b.trt + (n.trt-x.trt)) <=1 ~ "1"),
                        post.var = (a.trt + x.trt)*(b.trt + (n.trt-x.trt))/(((a.trt + x.trt)+(b.trt + (n.trt-x.trt)))^2*((a.trt + x.trt)+(b.trt + (n.trt-x.trt))+1)),
                        Group="Treatment") ) %>% dplyr::select(Group, everything())


        prior.df <- rbind(gcurve(expr = dbeta(x, shape1 = a.con, shape2 = b.con),
                                 from = 0, to = 1, n = 1001, category = "Control") %>%
                                  mutate(alpha = a.con, beta = b.con, group = "Prior",
                                         density = "Control Prior"),
                          gcurve(expr = dbeta(x, shape1 = a.trt, shape2 = b.trt),
                                 from = 0, to = 1, n = 1001, category = "Treatment") %>%
                                  mutate(alpha = a.trt, beta = b.trt, group = "Prior",
                                         density = "Treatment Prior"),
                          gcurve(expr = dbeta(x, shape1 = a.con + x.con,
                                              shape2 = b.con + n.con - x.con),
                                 from = 0, to = 1, n = 1001, category = "Control") %>%
                                  mutate(alpha = a.con + x.con, beta = b.con + n.con - x.con,
                                         group = "Posterior",
                                         density = "Control Posterior"),
                          gcurve(expr = dbeta(x, shape1 = a.trt + x.trt,
                                              shape2 = b.trt + n.trt - x.trt),
                                 from = 0, to = 1, n = 1001, category = "Treatment") %>%
                                  mutate(alpha = a.trt + x.trt, beta = b.trt + n.trt - x.trt,
                                         group = "Posterior",
                                         density = "Treatment Posterior")) %>%
                mutate(group = factor(group, c("Prior", "Posterior"))) %>%
                mutate(category = factor(category, c("Treatment", "Control"))) %>%
                mutate(density = factor(density, c("Control Prior", "Control Posterior",
                                                   "Treatment Prior", "Treatment Posterior")))
        levels(prior.df$density) <- c(paste0("Control Prior: Beta(",
                                             a.con, ", ", b.con, ")"),
                                      paste0("Control Posterior: Beta(",
                                             a.con + x.con, ", ", b.con + (n.con - x.con), ")"),
                                      paste0("Treatment Prior: Beta(",
                                             a.trt, ", ", b.trt, ")"),
                                      paste0("Treatment Posterior: Beta(",a.trt + x.trt,
                                             ", ", b.trt + (n.trt - x.trt), ")"))
        levels(prior.df$group) <- c(paste0("Control Prior: Beta(", a.con, ", ",
                                           b.con, ")\nTreatment Prior: Beta(", a.trt,
                                           ", ", b.trt, ")"),
                                    paste0("Control Posterior: Beta(", a.con + x.con, ", ",
                                           b.con + (n.con - x.con),
                                           ")\nTreatment Posterior: Beta(", a.trt + x.trt,
                                           ", ", b.trt + (n.trt - x.trt), ")"
                                    )
        )

        p1 <- ggplot(data = prior.df, aes(x = x, y = y, color = category, linetype=category)) +
                geom_line(size = 0.75) + facet_wrap(~group) +
                labs(x = "Response Rate (%)",
                     y = NULL,
                     color="Group",
                     linetype="Group",
                     title = "Prior and posterior distributions for response rates",
                     subtitle = paste0("Control data: ", round(x.con/n.con,2)*100,
                                       "% (", x.con, "/", n.con, "). ",
                                       "Treatment data: ", round(x.trt/n.trt,2)*100,
                                       "% (", x.trt, "/", n.trt, "). ",
                                       "Observed difference: ", (round(x.trt/n.trt,2)
                                                                 - round(x.con/n.con,2))*100,"%"))+
                scale_linetype_manual(values=c("solid", "dashed"))+
                scale_y_continuous(breaks=NULL)+
                scale_x_continuous(labels = scales::percent)

        list(ppp.df, p1)
}
