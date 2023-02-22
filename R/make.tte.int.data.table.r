#' @title Make TTE predictive probability table
#'
#' @param my.params output from return.tte.int.data.req
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.tte.int.data.req <- return.tte.int.data.req()
#' my.tte.int.data.table <- make.tte.int.data.table(my.params = my.tte.int.data.req)
#' my.tte.int.data.table

make.tte.int.data.table <- function(my.params = return.tte.int.data.req()){
  my.table <- my.params %>%
    group_by(interim.m) %>%
    dplyr::filter(decision=="Go") %>%
    slice(n()) %>%
    dplyr::select(m.con.prior, m.trt.prior, HR.prior, ARatio, interim.m, complement.m, HR.obs) %>%
    dplyr::rename(HR.Go=HR.obs) %>% left_join(
      my.params %>%
        group_by(interim.m) %>%
        dplyr::filter(decision=="No-Go") %>%
        slice(1) %>%
        dplyr::select(m.con.prior, m.trt.prior, HR.prior, ARatio, interim.m, complement.m, HR.obs) %>%
        dplyr::rename(HR.NG=HR.obs)) %>%
    mutate(study.end.Go = my.params$HR.go[1], study.end.NG = my.params$HR.ng[1])
  return(my.table)
}
