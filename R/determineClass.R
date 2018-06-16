#' @title Tried to determine the class of data.frame columns
#' @export
determineClass <- function(data) {
    ## adjust level
    prg <- apply(data, 2, function(x) { if (any(grepl("[a-zA-Z\\-]", x))) { "factor" } else { "numeric" }  })
    tmp <- sapply(data[,], function(x) as.numeric(as.character(x)))
    toNum <- colnames(tmp)[which(apply(tmp, 2, function(x) !all(is.na(x))))]
    toFact <- colnames(tmp)[which(!colnames(tmp) %in% toNum)]
    data[,which(colnames(data) %in% toFact)] <- sapply(data[,which(colnames(data) %in% toFact)], function(x) as.factor(x))
    data[,which(colnames(data) %in% toNum)] <- sapply(data[,which(colnames(data) %in% toNum)], function(x) as.numeric(as.character(x)))

    return(data)
}
