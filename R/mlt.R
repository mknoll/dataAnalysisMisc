#' @title Melt and merge data w/ pheno info
#'
#' @export
mlt <- function(data, pheno) {
    df <- list()
    for (i in 1:length(data[,1])) {
	df[[length(df)+1]] <- data.frame(VAL=unlist(data[i,]), VAR=rownames(data)[i], pheno)
    }
    return(do.call(rbind, df))
}
