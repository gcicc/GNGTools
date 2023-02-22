#' @title Make single sample binary prior/posterior plot
#' @param a.trt prior alpha parameter
#' @param b.trt prior beta parameter
#' @param n.trt number of trials
#' @param x.trt number of responses
#'
#' @return A ggplot object is returned.
#' @export
#'
#' @examples
#' my.ss.bin.ppp <- make.ss.bin.ppp(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 20)
#' my.ss.bin.ppp[[1]]
#' my.ss.bin.ppp[[2]]
#' @author Greg Cicconetti
make.ss.bin.ppp <- function(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 20){

  prior.df <- rbind(
    gcurve(expr = dbeta(x, shape1 = a.trt,shape2 = b.trt), from = 0, to = 1, n = 1001,
           category = "Treatment") %>%
      mutate(alpha = a.trt, beta = b.trt, group="Prior", density="Treatment Prior"),
    gcurve(expr = dbeta(x, shape1 = a.trt + x.trt, shape2  = b.trt + n.trt - x.trt),
           from = 0, to = 1, n = 1001, category = "Treatment") %>%
      mutate(alpha = a.trt + x.trt, beta = b.trt + n.trt - x.trt, group="Posterior",
             density="Treatment Posterior")) %>%
    mutate(density = factor(density, c("Treatment Prior", "Treatment Posterior")))
  levels(prior.df$density) <- c(paste0("Treatment Prior: Beta(", a.trt, ", ", b.trt,
                                       ")"),
                                paste0("Treatment Posterior: Beta(", a.trt + x.trt,
                                       ", ",  b.trt + (n.trt - x.trt), ")"))

  # This ggplot call is same for both plots http://127.0.0.1:7505/#tab-1955-1 and downloads too
 p1<- ggplot(data = prior.df, aes(x = x,y = y, color = category))+geom_line(size=.75) +
    labs(x = "Response Rate (%)",
         y = NULL,
         title="Prior and posterior distributions of treatment response rate",
         subtitle = paste0("Treatment data: ",round(x.trt/n.trt,2)*100,
                           "% (",  x.trt, "/", n.trt, ")"),
         color="Group", linetype="Density") +
    facet_wrap(~density)+
    theme(legend.position = "bottom")+guides(color="none")+
    scale_y_continuous(breaks=NULL, labels=NULL)+
    scale_x_continuous(labels = scales::percent)

 df <- data.frame(a.trt.0 = a.trt, b.trt.0 = b.trt, n.trt = n.trt, x.trt = x.trt, a.trt.n = a.trt + x.trt, b.trt.n = b.trt + (n.trt - x.trt))

 df <- df %>% mutate(prior.mean = a.trt.0/(a.trt.0 + b.trt.0),
               prior.mode = case_when(a.trt > 1 & b.trt > 1 ~ as.character(round((a.trt - 1)/(a.trt + b.trt-2), 4)),
                                      a.trt==1 & b.trt==1 ~ "any value in (0, 1)", # really any value in (0, 1)
                                      a.trt < 1 & b.trt < 1 ~ "bimodal {0, 1}", # bimodal (0, 1)
                                      a.trt <=1 & b.trt > 1 ~ "0",
                                      a.trt > 1 & b.trt <=1 ~ "1"),
               prior.var = a.trt.0*b.trt.0/((a.trt.0 + b.trt.0 + 1)*(a.trt.0 + b.trt.0)^2),
               post.mean = a.trt.n/(a.trt.n + b.trt.n),
               post.mode = case_when((a.trt + x.trt) > 1 & (b.trt + (n.trt-x.trt)) > 1 ~ as.character(round(((a.trt + x.trt) - 1)/((a.trt + x.trt) + (b.trt + (n.trt-x.trt)-2)), 4)),
                                     (a.trt + x.trt)==1 & (b.trt + (n.trt-x.trt))==1 ~ "any value in (0, 1)", # really any value in (0, 1)
                                     (a.trt + x.trt) < 1 & (b.trt + (n.trt-x.trt)) < 1 ~ "bimodal {0, 1}", # bimodal (0, 1)
                                     (a.trt + x.trt) <=1 & (b.trt + (n.trt-x.trt)) > 1 ~ "0",
                                     (a.trt + x.trt) > 1 & (b.trt + (n.trt-x.trt)) <=1 ~ "1"),
               post.var = a.trt.n*b.trt.n/((a.trt.n + b.trt.n + 1)*(a.trt.n + b.trt.n)^2)) %>%
   mutate(prior.mean = round(prior.mean, 4),
          prior.var = round(prior.var, 4),
          post.mean = round(post.mean, 4),
          post.var = round(post.var, 4))
  list(df, p1)
}

