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
			    frm0=as.formula(VAL~1 + (1|ID)), 
			    frm=as.formula(VAL~GRP + (1|ID)),
			    rand=NULL, nCores=NULL,
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
    }
    doParallel::registerDoParallel(nCores)

    out <- NULL
    if (!is.null(rand)) {
	stop("currently DEFUNCT!")
	# Use nlme
	print("Using nlme")
	out <- foreach(i=1:length(data[,1])) %dopar% {
	    ## Obtain Model p-value
	    ret <- NULL
	    tryCatch({
		df <- data.frame(VAL=data[i,], pheno)
		fit0 <- nlme::lme(frm0, data=df, rand=rand, method="ML")
		fit <- nlme::lme(frm, data=df, rand=rand, method="ML")
		### BREAKS THE FUNCTION
		aP <- anova(fit, fit0)[2,9]
		if (reCalcREML) {
		}
		ret <- data.frame(summary(fit)$tTable[-1,,drop=F], anovaP=aP, i=i, ID=rownames(data)[i])
	    }, error=function(e) { print(e) })
	}
	out <- do.call(rbind, out)
    } else {
	# Use lme4
	print("Using lme4")
	out <- foreach(i=1:length(data[,1])) %dopar% {
	    ## Obtain Model p-value
	    ret <- NULL
	    tryCatch({
		df <- data.frame(VAL=data[i,], pheno)
		fit0 <- lmer(frm0, data=df, REML=F)
		fit <- lmer(frm, data=df,  REML=F)
		aP <- anova(fit, fit0)[2,8]
		if (reCalcREML) {
		    fit <- lmer(frm, data=df,  REML=T)
		}
		ret <- data.frame(summary(fit)$coef[-1,,drop=F], anovaP=aP, i=i, ID=rownames(data)[i])
	    }, error=function(e) { print(e) })
	}
	out <- do.call(rbind, out)
	##add two-sided p.value
	out$p.value <- 2*pnorm(-abs(out$t.value))
    }
    
    doParallel::stopImplicitCluster()

    return(out)
}
