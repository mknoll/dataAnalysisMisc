#' @title Model based analysis with fixed effects
#' 
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import MASS
#' 
#' @param data data.frame containg a sample per column and 
#' one gene per row
#' @param pheno Pheno data.frame with columns representing variables
#' and rows patients. Must be matched. 
#' @param frm0 Null-model formula
#' @param frm Full-model formula
#' @param type type of analysis, can be "lm" (linear), "nb" (negativ binomial),
#' "bin" (logistic regression)
#' @param nCores number of cores to use
#' @param complete.cases automatically remove incomplete data (NA, Inf)
#'
#' @export
fixedEffAnalysis <- function(data, pheno, 
			    frm0=as.formula(VAL~1 ), 
			    frm=as.formula(VAL~GRP ),
			    type="lm",
			    nCores=NULL,
			    complete.cases=F,
			    padj="BH") {
    rownames(data) <- as.character(1:length(data[,1]))

    ## Check for missing data - FIXME
    #if (!complete.cases && (any(is.na(data) || any(is.infinite(data))))) {
#	stop("NAs or Inf values found!. Set complete.cases to T")
#    } else {
#	data <- data[complete.cases(data*0),,drop=F]
#    }

    if (is.null(nCores)) {
	nCores <- parallel::detectCores() - 1
	nCores <- ifelse(nCores == 0, 1, nCores)
    }
    doParallel::registerDoParallel(nCores)

    out <- foreach(i=1:length(data[,1])) %dopar% {
	cat(paste("\r   ", round(i/length(data[,1])*100), "%      ",sep=""))
	## Obtain Model p-value
	ret <- NULL
	df <- data.frame(VAL=unlist(data[i,]), pheno)

	if (type == "lm") {
	    tryCatch({
		##normal distribution
		fit0 <- lm(frm0, data=df)
		fit <- lm(frm, data=df)
		aP <- anova(fit, fit0)[2,6]
		ret <- data.frame(summary(fit)$coef[-1,,drop=F], anovaP=aP, i=i, RN=rownames(data)[i])
	    }, error=function(e) { print(e) })
	} else if (type == "nb") {
	    tryCatch({
		##normal distribution
		fit0 <- glm.nb(frm0, data=df)
		fit <- glm.nb(frm, data=df)
		aP <- anova(fit, fit0)[2,8]
		ret <- data.frame(summary(fit)$coef[-1,,drop=F], anovaP=aP, i=i, RN=rownames(data)[i])
	    }, error=function(e) { print(e) })
	} else if (type == "bin") {
	    tryCatch({
		##normal distribution
		fit0 <- glm(frm0, data=df, family=binomial(link=logit))
		fit <- glm(frm, data=df, family=binomial(link=logit))
		aP <- anova(fit, fit0, test="LRT")[2,5]
		ret <- data.frame(summary(fit)$coef[-1,,drop=F], anovaP=aP, i=i, RN=rownames(data)[i])
	    }, error=function(e) { print(e) })

	} else {
	    stop("Not implemented yet!")
	    return(NULL)
	}
    }
    ## FIXME
    if (F) {
	out <- do.call(rbind, out)

	if (!is.null(padj)) {
	    tmp <- out[which(!duplicated(out$i)),]
	    tmp$padjModell <- p.adjust(tmp$anovaP, padj)
	    out$padjModellLRT <- tmp$padjModell[match(out$i, tmp$i)]
	}
    }
    
    doParallel::stopImplicitCluster()

    return(out)
}
