#' Goparallel - initiate a parallel computing environment
#'
#' @param ncores number of cores
#'
#' @return This function initiates a parallel computing environment based on the parallel package
#' @export

goparallel <- function(ncores=7){
        message(paste("\nCurrent Connections: ", dim(showConnections())[1], "\n"))
        message("\nClosing any open connections...\n")
        closeAllConnections()
        if(exists("cl")) remove(cl)
        message(paste("\nCurrent Connections: ", dim(showConnections())[1], "\n"))
        message(paste("\nStarting new cluster with", ncores, "cores...\n"))
        assign("cl", parallel::makeCluster(spec = ncores, type="PSOCK"), envir=environment())
        #cl <<- parallel::makeCluster(spec = ncores, type="PSOCK")
        message(" Cluster initiation complete\n")
        message(paste("\nCurrent Connections: ", dim(showConnections())[1], "\n"))
        message(paste("\n", exists("cl"), "\n"))

        parallel::clusterEvalQ(cl=cl, expr = {requireNamespace("tidyverse")})
        message("\n\n***\nThe tidyverse pacakge has been sent to each core.\nDo you need other parallel::clusterEvalQ or parallel::clusterExport calls before running your code?\n****\n")
}
