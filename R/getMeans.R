#' @title Get Means
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#' @export
getMeans <- function(vec, eps=10^-5, eps2=-0.1, n=1000, frac=0.005, log=T, bw=0.1) {    
    if (log) {    
        v <- log(vec[which(vec != 0)])    
    }    
    v <- v[which(!is.na(v))]    

    ##use parallel computing
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    ## FIXME
    if (is.na(no_cores)) { no_cores <- 4 }
    doParallel::registerDoParallel(no_cores)

    len <- length(v)
    coll <- foreach(i=1:n) %dopar% {    
	ret <- NULL
	tryCatch({
	    mn <- v[sample(1:len, frac*len, replace=T)]    
	    d <- density(mn)    
	    d0 <- diff(d$y)    
	    d00 <- diff(density(d0)$y)    
	    w <- which((d0[-1] < -eps & d0[-length(d0)] > eps) & d00[-c(1:2)] < eps2 )
	    ret <- d$x[w]    
	}, error=function(e) {})
	ret
    }    
    mus <- unlist(coll)    

    d <- density(mus, bw=bw)    
    d0 <- diff(d$y)    
    w <- which((d0[-1] < -eps & d0[-length(d0)] > eps))

    doParallel::stopImplicitCluster()
    
    return(d$x[w])    
}    


#' @title Get Means (MT)
#'
#' @import mixtools
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#' @export 
getMeansMT <- function(vec, maxTry=10, log=T) {    
    v <- vec
    if (log) {    
        v <- log(v[which(v != 0)])    
    }    
    v <- v[which(!is.na(v))]    
    cnt <- 0
    ret <- NULL
    while (cnt < maxTry && is.null(ret)) {
	cnt <- cnt+1
	ret <- NULL
        tryCatch({    
            ret <-normalmixEM(v)
        }, error=function(e) {})    
    }    
    
    return(ret)    
}    

