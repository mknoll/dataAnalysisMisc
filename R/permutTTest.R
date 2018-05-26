#' @title Permutation t-test
#' 
#' @description Permutation t-test (Welch)
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
permTTest <- function(data, grp, B=10000) {
    ##check if group has only two levels
    if (length(levels(factor(grp))) != 2) {
        stop("Exactly two group levels required!")
    }
    
    ## Order data frame
    data <- data[order(grp)]
    grp <- grp[order(grp)]
    
    grp <- as.integer(as.factor(grp))
    
    ###################
    #dyn.load("~/phd/projects/dataAnalysisMisc/dataAnalysisMisc/src/permTTest.so")
    ret <- .C("perm_tTest_unpaired_unequalVar", 
       data=as.numeric(data), 
       lenRows=length(data), 
       grp=grp, 
       B=as.integer(B),
       pval=numeric(1))
    
    return (ret$pval)
}
