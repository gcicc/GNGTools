#' @title Get single sample binary interim operating characteristics simulation
#' @param ss.bin.int.df call to get.ss.bin.interim.df
#' @param ss.bin.int.GNG call to get.ss.bin.interim.GNG
#' @param lower lower bound
#' @param upper upper bound
#' @param step stepsize
#' @param include_nogo logical
#' @examples \donttest{
#' my.ss.bin.int.df = get.ss.bin.int.df(ss.bin.studyend.GNG =
#' get.ss.bin.studyend.GNG(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
#' Delta.lrv = .2, Delta.tv = .35,
#' tau.tv = 0.10, tau.lrv = .80, tau.ng = .65),
#' goThreshold = .8,
#' nogoThreshold = 1.2,
#' include_nogo = TRUE)
#' my.ss.bin.studyend.GNG <- get.ss.bin.studyend.GNG(a.trt = 1, b.trt = 1, n.trt = 40, x.trt = 9,
#'                                                   Delta.lrv = .2, Delta.tv = .35,
#'      tau.tv = 0.10, tau.lrv = .80, tau.ng = .65)
#'my.ss.bin.int.GNG <- get.ss.bin.int.GNG( ss.bin.int.df = my.ss.bin.int.df,
#'                                         Interims = 20,
#' ss.bin.studyend.GNG = my.ss.bin.studyend.GNG)
#' my.ss.bin.int.GNG
#' }
#'
#' @return Returns a data.frame holding results at interim and final
#' @export
#'

get.ss.bin.int.oc <- function(ss.bin.int.df,
                              ss.bin.int.GNG,
                              lower = 0,upper = 1,step = .025, include_nogo=TRUE){
        if(include_nogo==T){
                purrr::map_df(seq(from= lower,to = upper,by = step), function(x) ocTable_multi(DesignTable = ss.bin.int.GNG,
                                                                                        TargetRate = x)) %>%
                        mutate(
                                Go_upper = `First Go`,
                                Go_lower = 0,
                                Go_line = `First Go`,
                                Consider_upper = `First Go` + `Grey Area`,
                                Consider_lower = `First Go`,
                                Consider_line = `Grey Area`,
                                `No Go_upper`  = 1,
                                `No Go_lower` = 1-`First No Go`,
                                `No Go_line` = `No Go_lower`
                                )%>%
                        select(Analysis,targetRate,ends_with(c('upper','lower','line'))) %>%
                        pivot_longer(ends_with(c('upper','lower','line')),names_to = c("Decision", ".value"),names_pattern = "(.*)_(.*)") %>%
                        mutate(
                                Type = factor(Analysis,levels = c(paste('Analysis',1:nrow(ss.bin.int.GNG)),'Final-solo','Overall'), labels = c(paste('Analysis',1:nrow(ss.bin.int.GNG)),'Final','Any')),
                                Decision = factor(Decision,levels = c('Go','Consider','No Go'))
                                )
                } else {
                        map_df(seq(from= lower,to = upper,by = step),function(x) ocTable_multi(DesignTable = ss.bin.int.GNG,TargetRate = x)) %>%
                                mutate(
                                        Go_upper = `First Go`,
                                        Go_lower = 0,
                                        Go_line = `First Go`,
                                        Consider_upper = `First Go` + `Grey Area`,
                                        Consider_lower = `First Go`,
                                        Consider_line = `Grey Area`,
                                        `No Go_upper`  = 1,
                                        `No Go_lower` = 1-`First No Go`,
                                        `No Go_line` = `No Go_lower`
                                        )%>%
                                select(Analysis,targetRate,ends_with(c('upper','lower','line'))) %>%
                                pivot_longer(ends_with(c('upper','lower','line')),names_to = c("Decision", ".value"),names_pattern = "(.*)_(.*)") %>%
                                mutate(
                                        Type = factor(Analysis,levels = c(paste('Analysis',1:nrow(ss.bin.int.GNG)),'Final-solo','Overall'), labels = c(paste('Analysis',1:nrow(ss.bin.int.GNG)),'Final','Any')),
                                        Decision = factor(Decision,levels = c('Go','Consider','No Go'))
                                )
                }
}



