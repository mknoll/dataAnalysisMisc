#' @title Model based analysis with random effects
#' 
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import lme4
#' @import glmmTMB
#' @import lmtest
#' 
#' @param data data.frame containg a sample per column and 
#' one gene per row
#' @param pheno Pheno data.frame with columns representing variables
#' and rows patients. Must be matched. 
#' @param frm0 Null-model formula
#' @param frm Full-model formula
#' @param type type of analysis, can be "lm" (lmer), 
#' "nb" (glmer.nb), "p" (glmer, family=poisson(link=log)),
#' "l" (glmer, family=binomial(link=logit)), "b" (glmmTMB, family=beta)
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
			    reCalcREML=T, complete.cases=T,
			    transf=F) {
    ## Check for missing data
    #if (!complete.cases && any(is.na(data))) {
#	stop("NAs or Inf values found!. Set complete.cases to T")
#    } else {
#	data <- data[complete.cases(data*0),,drop=F]
 #   }

    if (is.null(nCores)) {
	nCores <- parallel::detectCores() - 1
	nCores <- ifelse(nCores == 0, 1, nCores)
    }
    doParallel::registerDoParallel(nCores)

    out <- NULL
    if (!is.null(rand)) {
	stop("currently DEFUNCT!")
    } else {
	out <- foreach(i=1:length(data[,1])) %dopar% {
	    cat(paste("\r  ", round(i/length(data[,1])*100), "%       ", sep=""))
	    ## Obtain Model p-value
	    ret <- NULL
	    df <- data.frame(VAL=unlist(data[i,]), pheno)

	    if (type == "lm") {
		##normal distribution
		tryCatch({
		    fit0 <- lmer(frm0, data=df, REML=F)
		    fit <- lmer(frm, data=df,  REML=F)
		    aP <- anova(fit, fit0)[2,8]
		    if (reCalcREML) {
			fit <- lmer(frm, data=df,  REML=T)
		    }
		    ret <- data.frame(summary(fit)$coef[-1,,drop=F], anovaP=aP, i=i, ID=rownames(data)[i])
		}, error=function(e) { })
	    } else if (type == "nb") {
		##negative binomial
		tryCatch({
		    fit0 <- glmer.nb(frm0, data=df)
		    fit <- glmer.nb(frm, data=df)
		    a <- anova(fit0, fit, test="LRT")
		    aP <- a[2,8]
		    ret <- data.frame(summary(fit)$coef[-1,,drop=F], i=i, aP=aP)
		}, error=function(e) { })
	    } else if (type == "p") {
		##poisson
		tryCatch({
		    fit0  <- glmer(frm0, family=poisson(link=log), data=df)
		    fit  <- glmer(frm, family=poisson(link=log), data=df)
		    a <- anova(fit0, fit, test="LRT")
		    aP <- a[2,8]
		    ret <- data.frame(summary(fit)$coef[-1,,drop=F], i=i,  aP=aP)
		}, error=function(e) { })
	    } else if (type == "l") {
		##binomial / logisitc regression
		tryCatch({
		    fit0  <- glmer(frm0, family=binomial(link=logit), data=df)
		    fit  <- glmer(frm, family=binomial(link=logit), data=df)
		    a <- anova(fit0, fit, test="LRT")
		    aP <- a[2,8]
		    ret <- data.frame(summary(fit)$coef[-1,,drop=F], i=i,  aP=aP)
		}, error=function(e) { })
	    } else if (type == "b") {
		## beta regression
		tryCatch({
		    ## FIXME
		    df <- df[which(!is.na(df$VAL)),]
		    if (transf) {
			df$VAL <- (df$VAL*(length(df$VAL)-1)+0.5)/length(df$VAL)
		    }
		    fit0 <- glmmTMB(frm0, data=df, family=list(family="beta", link="logit"))
		    fit <- glmmTMB(frm, data=df, family=list(family="beta", link="logit"))
		    pLRT <- lrtest(fit, fit0)[2,5]
		    ret <- data.frame(t(summary(fit)$coef$cond[-1,] ), i=i, ap=pLRT)
		    print(ret)
		}, error=function(e) { })
	    }
	}

	if (!"p.value" %in% colnames(out) && "t.value" %in% colnames(out)) {
	    out$p.value <- 2*pnorm(-abs(out$t.value))
	}
    }

    doParallel::stopImplicitCluster()

    return(out)
}


