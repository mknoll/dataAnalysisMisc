#' @title Model based analysis with random effects
#' 
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import nlme
#' @import lme4
#'
#' @export
randEffAnalysis <- function(data, pheno, 
			    frm0=as.formula(VAL~1), 
			    frm=as.formula(VAL~GRP),
			    rand=NULL, nCores=NULL,
			    reCalcREML=T) {
    if (is.null(nCores)) {
	nCores <- parallel::detectCores() - 1
	nCores <- ifelse(nCores == 0, 1, nCores)
	doParallel::registerDoParallel(nCores)
    }

    out <- NULL
    if (is.null(rand)) {
	# Use nlme
	print("Using nlme")
	out <- foreach(i=1:length(data[,1])) %dopar% {
	    ## Obtain Model p-value
	    ret <- NULL
	    tryCatch({
		fit0 <- lme(frm0, data=data, rand=rand, method="ML")
		fit <- lme(frm, data=data, rand=rand, method="ML")
		aP <- anova(fit0, fit)[2,9]
		if (reCalcREML) {
		    fit <- lme(frm, data=data, rand=rand)
		}
		ret <- data.frame(summary(fit)$tTable[-1,,drop=F], anovaP=aP, i=i, ID=rownames(i)[i])
	    }, error=function(e) { })
	    ret
	}
	out <- do.call(rbind, out)
    }
    
    doParallel::stopImplicitCluster()

    return(out)
}
