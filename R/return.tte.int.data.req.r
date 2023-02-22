#' @title Return time to event preditive probility
#'
#' @param m.con.prior prior number of control events
#' @param m.trt.prior prior number of treatment events
#' @param HR.prior prior estimate for HR
#' @param ARatio randomization ratio
#' @param interim.HR Interim HR
#' @param interim.m Interim events
#' @param final.m final events
#' @param HR.lrv TPP Lower Reference Value aka Max TPP (large HRs lead to No-Go)
#' @param HR.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param HR.ng HR needed for ng; leave null for standard rule
#' @param HR.go HR needed for go; leave null for standard rule
#' @param go.thresh go threshold for predictive probability
#' @param ng.thresh no-go threshold for predictive probability
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples {
#' return.tte.pp()
#' }
#' @author Greg Cicconetti
return.tte.pp <- function(m.con.prior=50, m.trt.prior=50, HR.prior=.845,
                          ARatio=1,
                          interim.HR=seq(.0025, 5,.0025),
                          interim.m=c(428, 750, 1000),
                          final.m=1500,
                          HR.tv=.139, HR.lrv = .139001,
                          tau.tv=.1, tau.lrv=.8, tau.ng=.65,
                          HR.ng = NULL, HR.go=NULL,
                          go.thresh= 0.8, ng.thresh =0.8){

  # If either are null... Then compute based on Rule in Action
  # This part only needs final.m value passed to m.obs to run
  if(is.null(HR.ng) | is.null(HR.go)){
    results <- get.tte.df(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior,
                          HR.prior = HR.prior, ARatio = ARatio, m.obs=final.m,
                          HR.obs=interim.HR, HR.tv=HR.tv, HR.lrv = HR.lrv,
                          tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng)
    result.go <- results %>% dplyr::filter(result== "Go") %>% slice(dplyr::n())
    # Determine max number of TRT responders for No-Go
    result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(1)
    HR.go <- ifelse(nrow(result.go) == 0, 0, result.go$HR.obs)
    HR.ng <- ifelse(nrow(result.ng) == 0, Inf, result.ng$HR.obs)




    }

  # get.tte.post uses expand.grid, so sending vector for inteirm.m should not pose problem
  # m.obs holds the values of interm.m
  my.params <- get.tte.post(m.con.prior=m.con.prior, m.trt.prior=m.trt.prior, HR.prior=HR.prior,
                            ARatio=ARatio, HR.obs=interim.HR, m.obs=interim.m) %>%
    dplyr::rename(interim.m = m.obs)  %>%
    mutate(
           HR.ng=HR.ng,
           HR.go=HR.go,
           final.m=final.m,
           complement.m=final.m - interim.m) %>%
    mutate(q=log(HR.go)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
           sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m),
           # Go when HR values are small, NoGo when HR values are large
           Go = pnorm(q = log(HR.go)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
                      mean = post.mean,
                      sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m)),
           NoGo = 1 - pnorm(q = log(HR.ng)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
                            mean = post.mean,
                            sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m)),
           Consider =1 - Go - NoGo,
           decision=case_when(Go > go.thresh ~ "Go",
                              NoGo > ng.thresh ~ "No-Go",
                              TRUE ~ "Consider"))
  return(my.params)
}


# This version uses data.frame instead of expand grid via get.tte.df.fg
#' Title
#'
#' @param m.con.prior prior number of control events
#' @param m.trt.prior prior number of treatment events
#' @param HR.prior prior estimate for HR
#' @param ARatio randomization ratio
#' @param interim.HR Interim HR
#' @param interim.m Interim events
#' @param final.m final events
#' @param HR.lrv TPP Lower Reference Value aka Max TPP (large HRs lead to No-Go)
#' @param HR.tv TPP Target Value aka Base TPP
#' @param tau.tv threshold associated with Base TPP
#' @param tau.lrv threshold associated with Min TPP
#' @param tau.ng threshold associated with No-Go
#' @param HR.ng HR needed for ng; leave null for standard rule
#' @param HR.go HR needed for go; leave null for standard rule
#' @param go.thresh go threshold for predictive probability
#' @param ng.thresh no-go threshold for predictive probability
#' @param include_nogo logical
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples \donttest{
#' holdit <- return.tte.int.data.req()
#' head(holdit)
#' }
return.tte.int.data.req <- function(m.con.prior=.001, m.trt.prior=.001, HR.prior=1,
                             ARatio=1,
                             interim.HR=c(.8, .7, .65, .68, .7), # must be vector of same length as interim.m
                             interim.m=c(100, 200, 300, 400, 750),
                             final.m=1000,
                             HR.tv=0.7, HR.lrv = 0.9,
                             tau.tv=.1, tau.lrv=.8, tau.ng=.65,
                             HR.ng = NULL, HR.go=NULL,
                             go.thresh= 0.8, ng.thresh =0.8, include_nogo=TRUE){

  # If either are null... Then compute based on Rule in Action
  # This part only needs final.m value passed to m.obs to run
  if(is.null(HR.ng) | is.null(HR.go)){
    results <- get.tte.df(m.con.prior = m.con.prior, m.trt.prior = m.trt.prior,
                          HR.prior = HR.prior, ARatio = ARatio, m.obs=final.m,
                          HR.obs=interim.HR, HR.tv=HR.tv, HR.lrv = HR.lrv,
                          tau.tv=tau.tv, tau.lrv=tau.lrv, tau.ng=tau.ng)
    result.go <- results %>% dplyr::filter(result== "Go") %>% slice(dplyr::n())
    # Determine max number of TRT responders for No-Go
    result.ng <- results %>% dplyr::filter(result== "No-Go") %>% slice(1)
    HR.go <- ifelse(nrow(result.go) == 0, 0, result.go$HR.obs)
    HR.ng <- ifelse(nrow(result.ng) == 0, Inf, result.ng$HR.obs)
  }

  # get.tte.post uses expand.grid, so sending vector for inteirm.m should not pose problem
  # m.obs holds the values of interm.m

  if(include_nogo==T){
  my.params <- get.tte.post.df(m.con.prior=m.con.prior, m.trt.prior=m.trt.prior, HR.prior=HR.prior,
                               ARatio=ARatio, HR.obs=interim.HR, m.obs=interim.m) %>%
    dplyr::rename(interim.m = m.obs)  %>%
    mutate(
           HR.ng=HR.ng,
           HR.go=HR.go,
           complement.m=final.m - interim.m) %>%
    mutate(q=log(HR.go)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
           sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m),
           # Go when HR values are small, NoGo when HR values are large
           Go = pnorm(q = log(HR.go)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
                      mean = post.mean,
                      sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m)),
           NoGo = 1 - pnorm(q = log(HR.ng)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
                            mean = post.mean,
                            sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m)),
           Consider =1 - Go - NoGo,
           decision=case_when(Go > go.thresh ~ "Go",
                              NoGo > ng.thresh ~ "No-Go",
                              TRUE ~ "Consider"))
  }

  if(include_nogo==F){
    my.params <- get.tte.post.df(m.con.prior=m.con.prior, m.trt.prior=m.trt.prior, HR.prior=HR.prior,
                                 ARatio=ARatio, HR.obs=interim.HR, m.obs=interim.m) %>%
      dplyr::rename(interim.m = m.obs)  %>%
      mutate(
        HR.ng=HR.ng,
        HR.go=HR.go,
        complement.m=final.m - interim.m) %>%
      mutate(q=log(HR.go)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
             sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m),
             # Go when HR values are small, NoGo when HR values are large
             Go = pnorm(q = log(HR.go)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
                        mean = post.mean,
                        sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m)),
             NoGo = 1 - pnorm(q = log(HR.ng)*(interim.m+complement.m)/complement.m - interim.m/complement.m*log(HR.obs),
                              mean = post.mean,
                              sd = sqrt(post.sd^2 + 1/(interim.m * ARatio*(1 - ARatio))/complement.m)),
             Consider =1 - Go - NoGo,
             decision=case_when(Go > go.thresh ~ "Go",
                                TRUE ~ "Consider"))
  }
  return(my.params)
}
