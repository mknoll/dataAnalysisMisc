#' @title Detect the best cutoff to match a categorical classification
#' 
#' @description
#' Separates a continuous variable with different cutoffs 
#' into two groups and performs a chisquared or Fisher's 
#' exact test
#' 
#' @param data vector with the continuous variables
#' @param class classification (two factors) to test
#' against
#' @param delta step size for alterations of 
findOptCutoffClass<- function(data, class, delta=0.1, test="chisq") {
    coll <- NULL
    cutoff <- max(data)
    
    step = abs(max(data, na.rm=T)-min(data, na.rm=T))*delta
    if (is.infinite(step)) {
        return(NULL)
    }
    
    while (cutoff  > min(data)) {
        cat(".")
        grp <- ifelse(data > cutoff, 2, 1)
        cutoff <- cutoff-step
        if (length(levels(factor(grp))) < 2 ) {
            next
        }
        
        if (test == "chisq") {
            p.val <- chisq.test(grp, class)$p.value
        }
        
        vec <- data.frame(cutoff=cutoff+step, p.val=p.val, n1=length(grp[grp == 1]), n2=length(grp[grp == 2]), delta=delta)
        coll <- rbind.fill(coll, vec)
    }
    
    coll <- coll[order(coll$p.val),]
    return(coll)
}