#' @title Make single sample binary interim Group-sequential style table
#' @param InterimDF call to get.ss.bin.interim.df
#' @param goThreshold  predictive probability threshold
#' @param nogoThreshold predictive probability threshold
#'
#' @return A data.frame is returned
#' @export
#'
make.ss.bin.int.GStable <- function(InterimDF, goThreshold = .8, nogoThreshold = .8){
        InterimDF %>%
                group_by(n) %>%
                summarise(
                        minGo = suppressWarnings(min(x[Decision =='Go'])),
                        MaxNoGo = suppressWarnings(max(x[Decision =='No Go'])),
                        Go = minGo/n,
                        `No Go` = MaxNoGo/n
                        ) %>%
                group_by(n) %>% summarize(min.required.for.Go=min(minGo), max.required.for.NoGo=max(MaxNoGo)) %>% mutate(Go.Prop = round(min.required.for.Go/n,4), NoGo.Prop=round(max.required.for.NoGo/n,4)) %>%
                select(n, max.required.for.NoGo, min.required.for.Go, NoGo.Prop, Go.Prop) %>%
                rename('Number of Subjects' = n, 'Min required for Go' = min.required.for.Go, 'Max required for No-Go' = max.required.for.NoGo, "Go Proportion" = Go.Prop, "No-Go Proportion" = NoGo.Prop)
}

