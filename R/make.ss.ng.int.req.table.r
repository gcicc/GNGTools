#' @title Make single sample normal-gamma interim data requirements table
#'
#' @param my.params output from return.ss.ng.int.req
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.ss.ng.int.req.table <- make.ss.ng.int.req.table()
#' head(my.ss.ng.int.req.table)
make.ss.ng.int.req.table <- function(my.params = return.ss.ng.int.req()){
  my.table <- my.params %>%
    group_by(interim.n.t) %>%
    dplyr::filter(decision=="Go") %>%
    slice(1) %>%
    dplyr::select(mu.0, n.0, alpha.0, beta.0, tdf.0, sigma.0, xbar, s,interim.n.t, xbar_ng, xbar_go, complement.n.t) %>%
    rename(xbar.go = xbar) %>%
    left_join(
      my.params %>%
        group_by(interim.n.t) %>%
        dplyr::filter(decision=="No-Go") %>%
        slice(n()) %>%
        dplyr::select(mu.0, n.0, alpha.0, beta.0, tdf.0, sigma.0, xbar, s,interim.n.t, xbar_ng, xbar_go, complement.n.t) %>%
        rename(xbar.ng = xbar)
    ) %>%
    mutate(study.end.Go = my.params$xbar_go[1], study.end.NG = my.params$xbar_ng[1]) %>%
    dplyr::select(mu.0, n.0, alpha.0, beta.0, tdf.0, sigma.0, interim.n.t, complement.n.t,  xbar.go, xbar.ng, s, study.end.Go, study.end.NG)
  return(my.table)
}

