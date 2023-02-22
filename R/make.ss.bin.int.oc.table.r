#' @title Make single sample binary interim OC table
#' @param Sims call to get.ss.bin.int.oc.sim
#' @param assessment assessment
#' @param rate rate
#'
#' @return A data.frame is returned
#' @export

make.ss.bin.int.oc.table <- function(Sims, assessment = 'all', rate = .3){

  lastAssessment = max(Sims$Assessment)
  runs = max(Sims$run)

  Tbl <- Sims %>%
    filter(Rate > rate -.0001, Rate < rate + .0001)

  if(assessment != 'all'){
    Tbl <- Tbl %>%
      filter(Assessment %in% c(assessment,lastAssessment))
  }

  Tbl %>%
    group_by(run) %>%
    summarise(
      Interim = case_when(
        any(Decision[-lastAssessment] != 'Consider') ~ as.character(Decision[Decision != 'Consider'][1]),
        T ~ 'Consider'
      ),
      Final = as.character(Decision[Assessment == lastAssessment])
    ) %>%
    count(Interim,Final) %>%
    mutate(across(c('Final','Interim'),~ factor(.x,levels = c('Go','Consider','No Go','Total')))) %>%
    complete(Interim,Final) %>%
    filter(Interim !='Total',Final !='Total') %>%
    mutate(n = ifelse(is.na(n),0,n)) %>%
    rbind(.,
          group_by(.,Interim) %>% summarise(.,n = sum(n)) %>% mutate(Final = 'Total')) %>%
    rbind(.,
          group_by(.,Final) %>% summarise(.,n = sum(n)) %>% mutate(Interim = 'Total')) %>%
    mutate(
      Prob = n/runs,
      lab = paste0(formatC(Prob*100,1,format = 'f'),'%')
    ) %>%
    select(-n,-Prob) %>%
    arrange(Interim,Final) %>%
    pivot_wider(names_from = 'Final',values_from = 'lab')
}


