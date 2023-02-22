#' @title Get single sample binary interim Go-NoGo
#' @param ss.bin.int.df output from get.ss.bin.interim.df
#' @param Interims vector of sample sizes at which interims will occur
#' @param ss.bin.studyend.GNG output from call to get.ss.bin.studyend.GNG
#' @return A tibble is returned holding the min and max needed for Go and No-Go resp.
#' @examples {
#' my.ss.bin.int.df <- get.ss.bin.int.df(ss.bin.studyend.GNG =
#' get.ss.bin.studyend.GNG(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
#' Delta.lrv = .2, Delta.tv = .35, tau.tv = 0.10, tau.lrv = .80, tau.ng = .65),
#' goThreshold = .8, nogoThreshold = 1.2, include_nogo = TRUE)
#' my.ss.bin.studyend.GNG <-  get.ss.bin.studyend.GNG(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
#'                                                    Delta.lrv = .2, Delta.tv = .35,
#'                                                    tau.tv = 0.10, tau.lrv = .80, tau.ng = .65)
#'    my.ss.bin.int.GNG <- get.ss.bin.int.GNG(
#'    ss.bin.int.df = my.ss.bin.int.df,
#'    Interims = 20,
#'    ss.bin.studyend.GNG = my.ss.bin.studyend.GNG)
#'  }
#'
#' @export
#'
get.ss.bin.int.GNG <- function(
                ss.bin.int.df,
                Interims = 20,
                ss.bin.studyend.GNG
){
            if( !all(Interims %in% ss.bin.int.df$n)) stop("All values of Interims must be in 0:N")

        InterimGNG <- ss.bin.int.df %>%
                dplyr::filter(n %in% Interims) %>%
                group_by(n) %>%
                summarise(
                        MinGo = min(x[Decision == 'Go']),
                        MaxNoGo = max(x[Decision == 'No Go'])
                )

        InterimGNG <- rbind(InterimGNG,
                            tibble(n = ss.bin.studyend.GNG$result.go$n.trt,MinGo = ss.bin.studyend.GNG$result.go$x.trt,MaxNoGo = ss.bin.studyend.GNG$result.ng$x.trt)) %>%
                mutate(Diff = n-lag(n,default = 0))

        return(InterimGNG)
}
