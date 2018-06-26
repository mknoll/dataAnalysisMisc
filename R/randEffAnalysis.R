#' @title Model based analysis with random effects
#' 
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import nlme
#' @import lme4
#' 
#' @param data data.frame containg a sample per column and 
#' one gene per row
#' @param pheno Pheno data.frame with columns representing variables
#' and rows patients. Must be matched. 
#' @param frm0 Null-model formula
#' @param frm Full-model formula
#' @param type type of analysis, can be "lm" (lmer), 
#' "nb" (glmer.nb), "qp" (glm, family=quasipoisson(link=log))
#' @param rand random formula required in nlme, if provided,
#' lmer will be used for "lm"
#' @param nCores number of cores to use
#' @param complete.cases automatically remove incomplete data (NA, Inf)
#'
#' @export
randEffAnalysis <- function(data, pheno, 
			    frm0=as.formula(VAL~1 + (1|ID)), 
			    frm=as.formula(VAL~GRP + (1|ID)),
			    type="lm",
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
	    cat(paste("\r  ", round(i/length(data[,1])*100), "%       ", sep=""))
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
	    cat(paste("\r  ", round(i/length(data[,1])*100), "%       ", sep=""))
	    ## Obtain Model p-value
	    ret <- NULL
	    tryCatch({
		df <- data.frame(VAL=data[i,], pheno)

		if (type == "lm") {
		    ##normal distribution
		    fit0 <- lmer(frm0, data=df, REML=F)
		    fit <- lmer(frm, data=df,  REML=F)
		    aP <- anova(fit, fit0)[2,8]
		    if (reCalcREML) {
			fit <- lmer(frm, data=df,  REML=T)
		    }
		    ret <- data.frame(summary(fit)$coef[-1,,drop=F], anovaP=aP, i=i, ID=rownames(data)[i])
		} else if (type == "nb") {
		    ##negative binomial
		    fit0 <- glmer.nb(frm0, data=df)
		    fit <- glmer.nb(frm, data=df)
		    a <- anova(fit0, fit, test="LRT")
		    aP <- a[2,8]
		    ret <- data.frame(summary(fit)$coef[-1,,drop=F], i=i, aP=aP)
		} else if (type == "qp") {
		    ##quasipoisson
		    fit0  <- glm(frm0, family=quasipoisson(link=log), data=df)
		    fit  <- glm(frm, family=quasipoisson(link=log), data=df)
		    a <- anova(fit0, fit, test="LRT")
		    aP <- a[2,5]
		    ret <- data.frame(summary(fit)$coef[-1,,drop=F], i=i,  aP=aP)
		}

	    }, error=function(e) { print(e) })
	}
	out <- do.call(rbind, out)

	if (!"p.value" %in% colnames(out) && "t.value" %in% colnames(out)) {
	    out$p.value <- 2*pnorm(-abs(out$t.value))
	}
    }
    
    doParallel::stopImplicitCluster()

    return(out)
}
