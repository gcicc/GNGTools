#' @title Make two-sample binary interim OC table
#'
#' @param OCsims results from get.ts.bin.int.oc
#' @param Diff the difference of interest for the table, must be a difference value included in OCsims
#' @param interim either a specific interim analysis or 'any'
#'
#' @return returns a 3x3 table showing the probabilities of each outcome for the interim and final analysis.
#' If interm is 'any' the interim analysis is the first non-consider result observed at any interim.
#' @export

make.ts.bin.int.oc.table <- function(OCsims=get.ts.bin.int.oc(),Diff=.1,interim='any'){

        FF <- max(OCsims$assessment)
        runs <- max(OCsims$run)

        if(interim == 'any'){
                Frm <- OCsims %>%
                        filter(near(Effect,Diff )) %>%
                        group_by(run) %>%
                        summarise(
                                Final = decision[assessment == FF],
                                Interim = case_when(
                                        all(decision[assessment != FF] == 'Consider') ~ 'Consider',
                                        T ~ as.character(decision[assessment != FF & decision != 'Consider'][1])
                                )) %>%
                        mutate(
                                Interim = factor(Interim,levels = c('Go','Consider','NoGo'))
                        )

        }else{
                Frm <- OCsims %>%
                        filter(near(Effect,Diff )) %>%
                        group_by(run) %>%
                        summarise(
                                Final = decision[assessment == FF],
                                Interim = decision[assessment == interim]
                        ) %>%
                        mutate(
                                Interim = factor(Interim,levels = c('Go','Consider','NoGo'))
                        )
        }

        holdit <- Frm %>%
                group_by(Final,Interim) %>%
                summarise(Rate = n()/runs) %>%
                ungroup() %>%
                complete(Final,Interim,fill = list(Rate = 0)) %>%
                pivot_wider(names_from ='Final', values_from = 'Rate',values_fill = list(Rate=0)) %>%
                janitor::adorn_totals(where = c('row','col'))
        holdit$Interim <- factor(holdit$Interim)
        holdit$Interim <- recode(holdit$Interim, Consider = "Wait", Go = "Accelerate", 'No Go'="Do not Accelerate", Total="Total" )
        holdit
}
