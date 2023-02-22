#' @title Get single sample binary interim operating charactersitics simulations
#' @param InterimGNG GNG decision value outputs from get.ss.bin.int.GNG
#' @param TrueRate assumed true response rate. Can be a vector
#' @param runs number of simulation runs
#' @examples \donttest{
#' my.ss.bin.int.df <- get.ss.bin.int.df(ss.bin.studyend.GNG =
#' get.ss.bin.studyend.GNG(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
#'                         Delta.lrv = .2, Delta.tv = .35,
#'                         tau.tv = 0.10, tau.lrv = .80, tau.ng = .65),
#'                    goThreshold = .8,
#'                    nogoThreshold = 1.2,
#'                    include_nogo = TRUE)
#' my.ss.bin.studyend.GNG <- get.ss.bin.studyend.GNG(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
#'                                                   Delta.lrv = .2, Delta.tv = .35,
#'                                                   tau.tv = 0.10, tau.lrv = .80, tau.ng = .65)
#' my.ss.bin.int.GNG <- get.ss.bin.int.GNG( ss.bin.int.df = my.ss.bin.int.df,
#'                                         Interims = 20,
#'                                         ss.bin.studyend.GNG = my.ss.bin.studyend.GNG)
#' my.ss.bin.int.GNG
#' }
#' @return Returns results of a simulation.
#' @export

get.ss.bin.int.oc.sim <- function(InterimGNG,
                                  TrueRate,
                                  runs = 500){

  Sims <- crossing( run = 1:runs, Assessment = 1:nrow(InterimGNG),Rate = TrueRate ) %>%
    mutate(Successes = rbinom(n=n(),size = InterimGNG$Diff[Assessment],prob = Rate)) %>%
    group_by(run,Rate) %>%
    mutate(Total = cumsum(Successes)) %>%
    ungroup() %>%
    mutate(
      Decision = case_when(
        Total >= InterimGNG$MinGo[Assessment] ~ 'Go',
        Total <= InterimGNG$MaxNoGo[Assessment] ~ 'No Go',
        TRUE ~'Consider'
      ) %>% factor(levels = c('Go','Consider','No Go')))

  return(Sims)
}
