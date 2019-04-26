#' @title Calculate localtion sensitive similarity, 450k
#' 
#' @param type pos, neg, both
#' @param n number of neighboring CpGs
#' @param corM type of correlation
#' @param cutPosCor 
#' @param cutNegCor 
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @export
simLoc450k <- function(data, 
		       type="both",
		       n=3,
		       corM="pearson",
		       cutPosCor=0.9,
		       cutNegCor=-0.9) {
    #### TODO: add EPIC annotation
    print("Get annotation ...")
    anno <- minfi::getAnnotation(minfiData::RGsetEx)
    anno <- anno[order(anno$chr, anno$pos),,drop=F]

    ### go parallel
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    ## FIXME
    if (is.na(no_cores)) { no_cores <- 4 }
    doParallel::registerDoParallel(no_cores)
    
    ### match data
    data <- data.frame(data[match(rownames(anno), rownames(data)),])
    ### split by chromosome 
    f <- split(data, f=anno$chr)
    chr <- unique(anno$chr)

    coll <- list()
    for (i in 1:length(f)) {
	print(paste("Processing ", chr[i], "..."))
	sub <- f[[i]]

	#### indices 
	from <- n+1
	to <- length(sub[,1])
	#### Correlation pos
	print("Positive correlation")
	pos <- foreach(j=from:to) %dopar% {
	    if (j %% 100 == 0) { cat(paste("\r  ", round(j/to*100, 2), "%     ")) }
	    m <- cor(sub[(from-n):to,], method=corM)
	    m[which(m < cutPosCor)] <- 0
	    m
	}
	posA <- array(dim=c(length(pos[[1]][1,]),
			    length(pos[[1]][1,]),
			    length(pos)))
	for (j in 1:length(pos)) {
	    posA[,,j] <- pos[[j]]
	}
	pos <- apply(posA, c(1,2), function(x) sum(x, na.rm=T))
	rm(posA)


	#### Correlation neg
	print("Negative correlation")
	neg <- foreach(j=from:to) %dopar% {
	    m <- cor(sub[(from-n):to,], method=corM)
	    m[which(m > cutNegCor)] <- 0
	    diag(m) <- -1
	    m
	}
	negA <- array(dim=c(length(neg[[1]][1,]),
			    length(neg[[1]][1,]),
			    length(neg)))
	for (j in 1:length(neg)) {
	    negA[,,j] <- neg[[j]]
	}
	neg <- apply(negA, c(1,2), function(x) sum(x, na.rm=T))
	rm(negA)

	coll[[length(coll)+1]] <- list(pos=pos, neg=neg, chr=chr[i])
    }

    doParallel::stopImplicitCluster()

    return(coll)
 }
