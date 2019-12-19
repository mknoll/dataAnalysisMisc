#' @title Barnards test (parallel)
#' @param data 0, 1
#'     
#' @import foreach    
#' @import doParallel    
#' @import parallel    
#' @export
barnP <- function(data, grp, side="two", cut=1) {
    ## TODO: sanitize inputs
 
    n <- apply(data, 1, sum)
    sub <- data[which(n >= cut),,drop=F]
    w <- ifelse(side=="two", 1, 2)

    ##use parallel computing
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    ## FIXME
    if (is.na(no_cores)) { no_cores <- 4 }
    doParallel::registerDoParallel(no_cores)

    coll <- foreach(i =1:length(sub[,1]))  %dopar% {
	tbl <- table(sub[i,], grp)
	ret <- NULL
	if (length(tbl) == 4) {
	    p <- barnard.test(tbl[1,1], tbl[1,2], tbl[2,1], tbl[2,2])$p.value[[w]]
	    ret  <- data.frame(p=p, side=side)
	}
	ret
    }
    names(coll) <- rownames(sub)
    coll <- Filter(length, coll)

    doParallel::stopImplicitCluster()

    return(coll)
}
