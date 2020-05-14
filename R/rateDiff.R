#' @title Create grid
#'
#' @description Create grid of rates or counts from two dim data+
#' e.g. as output from tSNE or UMAP
#' 
#' @param dat data.frame with two vectors (e.g. umap/tSNE) as first two 
#' colums
#' @param trea treatment vector, same order dat rows
#' @param bin number of bins of the resulting grid
#' @param rate calulate rates (divide by number of rows for the 
#' respective treatment)
#'
#' @export
createGrids <- function(dat, treat, nbin=50, rate=T) {
    ## build base matrix
    m0 <- matrix(nrow=nbin, ncol=nbin,0)
    sqY <-seq(from=min(dat[,2]), to=max(dat[,2]), len=nbin)
    dY <- sqY[2]-sqY[1]
    sqX <-seq(from=min(dat[,1]), to=max(dat[,1]), len=nbin)
    dX <- sqX[2]-sqX[1]

    minX <- min(dat[,1])
    minY <- min(dat[,2])

    ## get per treatment 
    mA <- list()
    for (tr in unique(treat)) {
	m <- m0
	sub <- dat[which(treat == tr),]
	l <- length(sub[,1])
	for (i in 1:length(sub[,1])) {
	    x <- floor((sub[i,1]-minX) / dX)
	    y <- floor((sub[i,2]-minY) / dY)
	    m[y,x] <- m[y,x]+1
	}
	if (rate) {
	    m <- m/l
	}
	mA[[length(mA)+1]] <- m
    }
    names(mA) <- unique(treat)
    mA
}



#' @title Get Quantiles from grid
#' 
#' @param dat output from createGrids()
#' 
#' @export 
getQuantiles<- function(dat, q=c(0.025, 0.975)) {
    allDiff <- list()
    for (i in 1:length(dat)) {
	for (j in 1:length(dat)) {
	    dM <- dat[[i]]-dat[[j]]
	    if (i > j) {
		allDiff[[length(allDiff)+1]] <- dM
	    }
	}
    }
    quantile(unlist(allDiff), q)
}


#' @title Plot pairwise comparisons of grids
#'
#' @param dat output from createGrids()
#' 
#' @export 
plotDiff <- function(dat, q=q) {
    n<- length(dat)
    print(paste("GRID: ",n,"x",n,sep=""))
    par(mfrow=c(n,n), mar=c(0,0,0,0))
    for (i in 1:length(dat)) {
	for (j in 1:length(dat)) {
	    if (i == j) { next }

	    #diff
	    dM <- dat[[i]]-dat[[j]]

	    #### gray background
	    dM2 <- dM
	    dM2[which(dM2 == 0)] <- NA
	    if (i <= j) {
		image(dM2, col=gray.colors(100)[-c(40:60)], axes=F)
	    } else {
		## placeholder
		image(dM2, col="white", axes=F)
	    }

	    ### colored: quantile
	    dM3 <- dM2
	    dM3[which(dM2 > q[1] & dM2 < q[2])] <- NA
	    if (i <= j) { 
		image(dM3, col=redgreen(100)[-c(40:60)], axes=F, add=T)
	    } else {
		### placeholder
		image(dM3, col="white", axes=F, add=T)
	    }
	}
    }
}


#' @title Quantify differences from pairwise comparisons
#' 
#' @param dat output from createGrids()
#' 
#' @export 
quantifyDiff <- function(dat) {
    m <- matrix(nrow=length(dat), ncol=length(dat), 0)
    for (i in 1:length(dat)) {
	for (j in 1:length(dat)) {
	    if (i == j) { next }
	    #diff
	    dM <- dat[[i]]-dat[[j]]
	    if (i <= j) {
		m[i,j] <- sum(dM^2)
	    }
	}
    }
    rownames(m) <- names(dat)
    colnames(m) <-names(dat)
    m
}
