#' @title Make single sample binary pred prob table
#'
#' @param my.df output from return.ss.bin.int.data.req
#'
#' @return A data.frame is returned
#' @export
#'
#' @examples
#' my.ss.bin.int.req.table <- make.ss.bin.int.req.table()
#' head(my.ss.bin.int.req.table)
make.ss.bin.int.req.table <- function(my.df = return.ss.bin.int.req()){
# So for each interim we ask: What's the smallest number of successes we need at this point to ensure Go.
# What's the largest number of success we need to ensure a no-go.
my.df <- my.df %>% group_by(interim.n.t) %>% dplyr::filter(decision=="Go") %>% arrange(interim.n.t, observed) %>% slice(1) %>%
  mutate(prop.needed = observed/interim.n.t) %>% full_join(

    my.df %>% group_by(interim.n.t) %>% dplyr::filter(decision=="No-Go") %>% arrange(interim.n.t, desc(observed)) %>% slice(1) %>%
      mutate(prop.needed = observed/interim.n.t))

  return(my.df)
}


