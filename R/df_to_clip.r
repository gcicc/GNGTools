#' Copy a data.frame to the clipboard
#'
#' @param df A data.frame
#' @param digits number of digits for rounding numerical columns
#'
#' @return the data.frame is copied to the clipboard allowing for easy pasting into excel (and subsequently to ppt, etc.)
#' @export
#'
df_to_clip <- function(df, digits=4) {
        nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
        df[,nums] <- round(df[,nums], digits = digits)
        write.table(x=df, file="clipboard", sep="\t", row.names = F)
}
