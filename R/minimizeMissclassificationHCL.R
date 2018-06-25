#' @title Minimize missclassification by step-wise increase
#' of evaluated candidates
#' 
#' @param data data.frame containing the data to be analyzed
#' @param res rownames of data, ordere by importance
#' @param clst_method clustering method to use for  
#' @param class given classification vector
#' @param test chisq or fisher
#' 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import pheatmap
#' 
#' @export 
#' 
minimizeMissclassificationHCL <- function(data, res, class,
                                          clst_method=c("complete", "ward.D2"),
                                          test="chisq", by=1, complete.cases=T) {
    ## Chekc for NA,Inf
    if (!complete.cases && any(is.na(data) | is.infinite(data))) {
	stop("NA/Inf found. Set complete.cases to T")
    } else {
	data <- data[complete.cases(data*0),]
    }


    ##use parallel computing
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)
    
    total <- list()
    
    cuts <- seq(from=2, to=length(res), by=by)
    
    for (cm in clst_method) {
        print(paste("Clustering method:", cm))
        out <- foreach(cut=cuts) %dopar% {
	    i <- which(cut == cuts)
	    cat(paste("\r    ", round(i/length(cuts)*100), "%  ", sep=""))
	    vec <- NULL
	    tryCatch({
		pm <- pheatmap::pheatmap(data[which(rownames(data) %in% res[1:cut]),], silent=T)
		pm.clust <- cutree(pm$tree_col, 2)
		p.val <- NA
		if (test == "chisq") {
		    p.val <- chisq.test(pm.clust, class)$p.value
		} else {
		    p.val <- fisher.test(table(pm.clust, class))$p.value
		}
		vec <- data.frame(cutoff=cut, p.val=p.val, cl_meth=cm)
	    }, error=function(e) {})
	    vec
        }
        tmp <- do.call(rbind, out)
        tmp <- tmp[order(tmp$p.val),]
        total[[length(total) + 1]] <- tmp 
    }
    
    doParallel::stopImplicitCluster()
    
    return (total)
}
