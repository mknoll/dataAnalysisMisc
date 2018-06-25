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
			    rand=as.formula(~1|ID), nCores=NULL,
			    reCalcREML=T, complete.cases=T) {
    ## Check for missing data
    if (!complete.cases && any(is.na(data) || is.infinite(data))) {
	stop("NAs or Inf values found!. Set complete.cases to T")
    } else {
	data <- data[complete.cases(data*0),,drop=F]
    }

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
