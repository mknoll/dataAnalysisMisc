#' @title Convert gRNA library screen data per sample
#' 
#' @param data count matrix, guides in rows, samples in cols
#' @param genename vector with gene names, corresponding to order
#' of rows in count matrix
#' @param funGene aggregation per gene within sample
#'
#' @export
convGuideMatrix <- function(data, genename, fun=max) {
    if (length(genename) != length(data[,1])) {
	stop("Dimensions not matching!")
    }

    # aggregate per gene
    print("Aggregate ... ")
    f <- split(data.frame(data), f=genename)
    # TODO: check NA values
    f <- lapply(f, function(x) apply(x, 2, fun))

    # convert to ranks per sample
    print("Convert to ranks .. ")
    f2 <- do.call(rbind, f)
    rownames(f2) <- names(f)
    f3 <- apply(f2, 2, function(x) rank(x))

    # set duplicated min ranks (no guide detected) to NA
    print("Set non-detected guides to NA...")
    f4 <- apply(f3,2, function(x) {
		    x[which((duplicated(x) | duplicated(x, fromLast=T)) & x == min(x) )] <- NA
		    x
})
    f4
}


#' @title Convert preproc gRNA libraray screen data per group
#' @export
convGuidePerGroup <- function(data, group) {
    if (length(data[1,]) != length(group)) {
	stop("Non-matching dimensions!")
    }
    ## split by group
    f <- split(data.frame(t(data)), f=group)
    f2 <- lapply(f, function(x) {
		     y <- apply(x, 2, function(y) max(y, na.rm=T))
		     names(y) <- rownames(data)
		     y
})
    data.frame(t(do.call(rbind, f2)), check.names=F)
}
