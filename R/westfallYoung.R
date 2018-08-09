#' @title Westfall Young p-value correction
#' 
#' @description Performs a Westfall-Young p-value adjustment 
#' using a Welch's t-test and a Satterthwaite approximation.
#' Might crash with a C stack overflow, if this happens, 
#' increase the size with ulimit -s 65535
#' 
#' @param grp grouping, exactly two levels required
#' @param data data 
#' @param B Number of permutations
#' 
#' @useDynLib dataAnalysisMisc
#' 
#' @export 
wy <- function(data, grp, B=100) {
    ##check if group has only two levels
    if (length(levels(factor(grp))) != 2) {
        stop("Exactly two group levels required!")
    }
    
    ## Order data frame
    data <- data[,order(grp)]
    grp <- grp[order(grp)]
    
    grp <- as.integer(as.factor(grp))
    
    ###################
    #dyn.load("wy.so")
    ret <- .C("tTest_unpaired_unequalVar", 
              data=as.numeric(unlist(t(data))), 
              lenRows=length(data[,1]), 
              grp=grp, lenGrp=length(grp), 
              B=as.integer(B),
              pval=numeric(length(data[,1])))
    
    return(ret$pval)
}
