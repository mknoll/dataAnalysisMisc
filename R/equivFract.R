#' @title Test for euqivalence of proportions
#' 
#' @param data data.frame with 0,1, features in rows and samples in cols
#' @param grp classification, 2 levels
#'
#' @import foreach    
#' @import doParallel    
#' @import parallel    
#' @export
equivFract <- function(data, grp, delta=0.1, z=1.65) {
    ## TODO: parallel
    ## TODO: sanitize inputs
    eMU <- abs(delta)
    eML <- -abs(delta)

    ### parallel
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    if (is.na(no_cores)) { no_cores <- 4 }
    doParallel::registerDoParallel(no_cores)

    coll <- list()
    coll <- foreach(i=1:length(data[,1])) %dopar% {
	vec <- unlist(data[i,])

	tbl <- table(vec, grp)
	p1 <- tbl[2,1]/sum(tbl[,1])
	p2 <- tbl[2,2]/sum(tbl[,2])
	dP <- p1-p2
	s <- sqrt(p1*(1-p1)/sum(tbl[,1]) + p2*(1-p2)/sum(tbl[,2]))
	lwr <- dP-1.65*s
	upr <- dP+1.65*s
	ret <- NULL
	if (lwr > eML && upr < eMU) {     
	    ret <- rownames(mm)[i]
	}
	ret
    }

    doParallel::stopImplicitCluster()

    return(unlist(coll))
}
