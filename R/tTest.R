#' @title Row-wise t-Test
#' 
#' @description Welch's t-test and a Satterthwaite approximation.
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
tTest <- function(data, grp) {
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
    ret <- .C("tTest_unpaired_unequalVar_simple", 
              data=as.numeric(unlist(t(data))), 
              lenRows=length(data[,1]), 
              grp=grp, lenGrp=length(grp), 
              df=numeric(length(data[,1])),
              statistic=numeric(length(data[,1])))
    
    ### calculate p val
    pval <- 2*pt(-abs(ret$statistic), df=ret$df)

    return(list(pval=pval, statistic=ret$statistic, df=ret$df))
}
