#' @title Make two-sample binary interim decision table
#'
#' @param DecisionTable the results from get.ts.bin.int.dec
#'
#' @return a simplified table showing the go/nogo thresholds for each number  of control successes.
#' @export
#'
#' @examples
#' my.ts.bin.int.dec <- get.ts.bin.int.dec()
#' my.ts.bin.int.dec.df <- make.ts.bin.int.dec.df(DecisionTable = my.ts.bin.int.dec)
#' head(my.ts.bin.int.dec.df)

make.ts.bin.int.dec.df <- function(DecisionTable){
        DecisionTable %>%
                group_by(Interim,IntermR_C) %>%
                summarise(MinGo = min(c(Inf,IntermR_T[Decision=='Go'])),
                          maxNoGo = max(c(-Inf,IntermR_T[Decision=='No Go'])))
}
