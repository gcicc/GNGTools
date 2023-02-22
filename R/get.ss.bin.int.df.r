#' get.ss.bin.interim.df
#' @title  Get single sample binary interim data.frame
#' @param ss.bin.studyend.GNG The output from get.ss.bin.studyend.GNG
#' @param goThreshold predictive posterior probability of a final Go threshold for an interim Go
#' @param nogoThreshold predictive posterior probability of a final No Go threshold for an interim NoGo
#' @param include_nogo logical
#' @return A dataframe holding interim pred prob of Go and No-Go along with interim decision
#' @export
#'
get.ss.bin.int.df <- function(
  ss.bin.studyend.GNG = get.ss.bin.studyend.GNG(a.trt = 1, b.trt = 1, n.trt = 50, x.trt = 9,
                                                Delta.lrv = .62, Delta.tv = .7,
                                                tau.tv = 0.10, tau.lrv = .80, tau.ng = .65),
  goThreshold = .8,
  nogoThreshold = .8,
  include_nogo =FALSE){

  alpha  = ss.bin.studyend.GNG$result.go$a.prior
  beta = ss.bin.studyend.GNG$result.go$b.prior
  N = ss.bin.studyend.GNG$result.go$n.trt
  goFinal = ss.bin.studyend.GNG$result.go$x.trt
  nogoFinal = ss.bin.studyend.GNG$result.ng$x.trt

  if(include_nogo==T & nogoThreshold < 1){
 for.return <- expand.grid(n = 1:N, x = 0:N) %>%
    dplyr::filter(x<=n) %>%
         dplyr::mutate(
                 Prob.Go = 1-pbetabinom_c(goFinal-x-1,N-n,(alpha+x)/(alpha+beta+n),alpha+beta+n),
                 Prob.NoGo = pbetabinom_c(nogoFinal-x,N-n,(alpha+x)/(alpha+beta+n),alpha+beta+n),
                 Threshold.Go = goThreshold,
                 threshold.NoGo = nogoThreshold,
                 Decision = dplyr::case_when(
                         Prob.Go >= goThreshold ~ 'Go',
                         Prob.NoGo >= nogoThreshold ~ 'No Go',
                         T ~ 'Consider') %>%
                         factor(levels = c('Go','Consider','No Go')),
                 include_nogo = include_nogo
                 )
 } else if(include_nogo==F | nogoThreshold>=1 | length(nogoThreshold)==0){
    for.return <-  expand.grid(n = 1:N, x = 0:N) %>%
            dplyr::filter(x<=n) %>%
            dplyr::mutate(
                    Prob.Go = 1-pbetabinom_c(goFinal-x-1,N-n,(alpha+x)/(alpha+beta+n),alpha+beta+n),
                    Prob.NoGo = pbetabinom_c(nogoFinal-x,N-n,(alpha+x)/(alpha+beta+n),alpha+beta+n),
                    Threshold.Go = goThreshold,
                    threshold.NoGo = 0.8,
                    Decision = dplyr::case_when(
                            Prob.Go >= goThreshold ~ 'Go',
                            T ~ 'Consider') %>%
                            factor(levels = c('Go','Consider')),
                    include_nogo = F # overrideing input!!
                    )
  }
  return(for.return)
}



