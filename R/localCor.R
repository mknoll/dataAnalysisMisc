#' @title Calculate genimc local correlation
#' 
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import pheatmap
#' @export
localCor <- function(data, len=3, cutLow=-0.9, cutHigh=0.9, corMethod="pearson") {
    ##use parallel computing
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)

    ## similarity / dissimilary matrices
    stackP <- matrix(ncol=length(data[1,]), nrow=length(data[1,]), 0)
    stackN <- stackP

    ## similariy matrix
    to <- length(data[,1])-len
    tmp <- foreach(i=1:to) %dopar% {
	if (i %% 100 == 0) {
	    cat(paste("\r    ", round(i/to*100), "%", sep=""))
	}
	## select relevatn data
	d <- data[i:(i+len),]
	cM <- cor(d)
	## flatter nmatrix
	m <- as.vector(unlist(cM))
	mN <- mP <- m

	#negative
	mN[which(mN > cutLow | is.na(mN))] <- 0
	#positive
	mP[which(mP < cutHigh | is.na(mP))] <- 0

	list(mN, mP)
    }

    ##aggregate stacks 
    print("\nAggregate pos stack ...")
    pos <- lapply(tmp, function(x) x[[1]])
    posS <- rep(0, length(pos[[1]]))
    posS <- foreach(i=1:length(posS)) %dopar% {
	if (i %% 100 == 0) {
	    cat(paste("\r    ", round(i/to*100), "%", sep=""))
	}
	sum(unlist(lapply(pos, function(x) x[i])), na.rm=T)
    }

    print("\nAggregate neg stack ...")
    neg <- lapply(tmp, function(x) x[[1]])
    negS <- rep(0, length(neg[[1]]))
    negS <- foreach(i=1:length(negS)) %dopar% {
	if (i %% 100 == 0) {
	    cat(paste("\r    ", round(i/to*100), "%", sep=""))
	}
	sum(unlist(lapply(neg, function(x) x[i])), na.rm=T)
    }

    print("Reformat data .. ")
    #restore matrix format
    posM <- matrix(nrow=length(data[1,]), ncol=length(data[1,]), posS)
    negM <- matrix(nrow=length(data[1,]), ncol=length(data[1,]), negS)
    
    doParallel::stopImplicitCluster()

    return(list(pos=posM, neg=negM))
}
