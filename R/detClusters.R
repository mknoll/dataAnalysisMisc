#' @title Determine number of clusters
#' @import pheatmap
#' @import parallel
#' @import foreach
#' @import doParallel
#' @export
detClusters <- function(data, nClust=NULL, fun=mean, lossfun="squared") {
    if (is.null(nClust)) {
	nClust <- c(1, length(data[,1]))
    }

    pm <- pheatmap(data, silent=T)

    ##use parallel computing
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)

    out <- foreach (i=nClust[1]:nClust[2]) %dopar% {
	pm.clust <- cutree(pm$tree_col, i) 
	
	hetTotal <- 0
	ret <- NULL
	tryCatch({
	    for (j in levels(factor(pm.clust))) {	
		if (length(which(pm.clust == j)) == 1) { next } 
		sub <- data[,which(pm.clust == j),drop=F]
		val <- apply(sub, 1, fun)
		if (lossfun == "squared") {
		    het <- sum((sub-val)^2)
		}
		het <- het/length(which(pm.clust == j))
		hetTotal <- hetTotal + het
	    }
	    ret <- data.frame(nClust=i, heterogeneity=het)
	}, error=function(e) { print(e) })
	ret
    }
    
    doParallel::stopImplicitCluster()

    out <- do.call(rbind, out)

    plot(out[,2], t="l", ylab="Heterogeneity", xlab="#Clusters")
    points(out[,2], cex=0.4, pch=19)

    return(out)
}

