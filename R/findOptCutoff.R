#' @title Find two differening groups in survival data
#' 
#' @description  Tests different cutoffs on a continuous variable 
#' and calculates log-rank tests / survival differences
#' 
#' @param data vector with value to analyze
#' @param srv Surv object with the corresponding data
#' @param delta proportion in which group splitting should 
#' be performed 
#' @param minGrpSize minimum group size for cutoffs
#' 
#' @import foreach
#' @import doParallel
#' @import parallel
#' 
#' @export
#' 
#' @return data.frame with the calculated survival values 
#' (p.value, median survival, group sizes, cutoffsize, 
#' stepsize)
findOptCutoff <- function(data, srv, delta=0.1, minGrpSize=1) {
    ##remove NA SRV data 
    keep <- which(!is.na(srv))
    srv <- srv[keep]
    data <- data[keep]
    
    cutoff <- max(data)
    step = abs(max(data, na.rm=T)-min(data, na.rm=T))*delta
    if (is.infinite(step)) {
        return(NULL)
    }
    
    ##calculate all cutoffs to test
    sq <- seq(from=min(data,na.rm=T), to=max(data,na.rm=T), by=step)
    sq <- sq[which(apply(data.frame(sq), 1, function(x) length(which(data > x)) > minGrpSize  & length(which(data < x)) > minGrpSize ))]
    if(length(sq) == 0) { 
        ## min group size not reached
        return (NULL) 
    }
    
    ##use parallel computing
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)
    
    out <- foreach(cutoff=sq) %dopar% {
        print(cutoff)
        grp <- ifelse(data > cutoff, "high", "low")
        cutoff <- cutoff-step
        if (length(levels(factor(grp))) == 2 ) {
            fit <- survfit(srv~grp)
            sdf <- survdiff(srv~grp)
            p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
            
            ##univariate coxph
            grp <- factor(grp)
            grp <- relevel(grp, "low")
            fit2 <- coxph(srv~factor(grp))
            
            data.frame(cutoff=cutoff+step, 
                       p.val=p.val,
                       nHigh=length(grp[grp == "high"]), 
                       nLow=length(grp[grp == "low"]), 
                       medianHigh=summary(fit)$table["grp=high","median"],
                       medianLow=summary(fit)$table["grp=low","median"],
                       delta=delta,
                       coxph_name=rownames(summary(fit2)$conf.int)[1],
                       coxph_hr=summary(fit2)$conf.int[1,1],
                       coxph_lower=summary(fit2)$conf.int[1,3],
                       coxph_upper=summary(fit2)$conf.int[1,4],
                       coxph_p=summary(fit2)$coeff[1,5],
                       coxph_lrt=summary(fit0)$logtest[[3]],
                       coxph_score=summary(fit0)$sctest[[3]],
                       coxph_wald=summary(fit0)$waldtest[[3]])
        }
    }
    coll <- do.call(rbind, out)
    coll <- coll[order(coll$p.val),]
    
    doParallel::stopImplicitCluster()
    
    return(coll)
}