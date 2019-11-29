#' @title Fit model
#' 
#' @description Fit a model with the selected features
#'
#' @param obj fitModel instance 
#' @param interval predict or confidence 
#'
#' @export
fitM <- function(obj, interval="predict",  ...) {
    ## max number of features < number of obs
    to <- ifelse(length(obj@cv[,1]) > length(obj@meta[,1]), length(obj@meta[,1]), length(obj@cv[,1])) 
    if (length(obj@cv[,1]) > length(obj@meta[,1])) {
	warning(paste("Using only first", to, "features! (nSample < nSignFeat)"))
    }

    ##FIXME
    dat2 <- data.frame(t(obj@dataSign), obj@meta)
    mss <- obj@cv
    frm <- as.formula(paste(obj@var, "~", paste(obj@cv[1:to,2], collapse="+")))
    if (obj@type == "lm") {
	fit <- lm(frm, data=dat2)
    } else if (obj@type == "lr") {
	fit <- glm(frm, data=dat2, family="binomial")
    }
    off <- 0
    while (off < 10 & any(is.na(summary(fit)$coef[,2]))) {
	off <- off+1
	to <- min(length(dat2[,1]), length(mss[,1]))-off
	frm <- as.formula(paste(obj@var, "~", paste(mss[1:to,2], collapse="+")))
	if (obj@type == "lm") {
	    fit <- lm(frm, data=dat2)
	} else if (obj@type == "lr") {
	    fit <- glm(frm, data=dat2, family="binomial")
	}
    }

    #FIXME
    if (any(is.na(summary(fit)$coef[,2]))) { 
	warning("Model with NA estimates!")
    }

    fitS <- step(fit, ...)
    obj@model <- fitS

    ### prediction
    if (obj@type == "lm") {
	prd <- predict(fitS, newdata=data.frame(t(obj@dataSign)),  interval=interval, se.fit=T)
	df <- data.frame(prd$fit)
	obj@pred <- df
    } else if (obj@type == "lr") {
	prd <- predict(fitS, newdata=data.frame(t(obj@dataSign)), interval=interval, se.fit=T, ...)
	df <- data.frame(fit=prd$fit, lwr=(prd$fit-prd$se.fit), upr=(prd$fit+prd$se.fit))
	obj@pred <- df
    }

    return(obj)
}

#' @title Test significance
#' 
#' @description Select significant features
#' 
#' @param obj fitModel object
#' @param pAdj p-value adjustment methods, as in p.adjust.methods
#' @param pCut p-value cutoff
#' @param frm0 Null model formula
#' @param frm Full model formula 
#' @param singlePVar name of the p-value column for single level/
#' feature tests (when using adjusted tests)
#'
#' @export 
testSign <- function(obj, 
		     pAdj="none",
		     pCut=0.05,
		     frm0=NULL,
		     frm=NULL,
		     singlePVar="Pr...t..",
		     ...) {
    print("Calculate models..")

    if (is.null(frm)) { frm <- as.formula(paste("VAL~", obj@var))  }
    if (is.null(frm0)) { frm0 <- as.formula(paste("VAL~1"))  }

    res<- fixedEffAnalysis(obj@data, obj@meta, frm0=frm0, frm=frm, ...)
    obj@signTest <- res 

    #FIXME: factor /character
    obj@signID <-  as.character(unlist(lapply(res, function(x) x[1,"ID"])))
    ## TODO: perform sanity checks
    obj@signP <- unlist(lapply(res, function(x) x[1,"anovaP"])  )
    obj@signPAdj <- p.adjust(obj@signP, pAdj)
    ### get P value of single var
    if (!is.null(res[[1]])) { 
	w <- which(grepl(obj@var, rownames(res[[1]])))
    } else {
	i <- 1
	while (is.null(res[[i]]) && i < 100) {
	    if (!is.null(res[[i]])) {
		w <- which(grepl(obj@var, rownames(res[[i]])))
	    }
	    i <- i+1
	}
    }
    obj@signP2 <- unlist(lapply(res, function(x) x[w,singlePVar])  )

    len <- length(which(obj@signPAdj  < pCut))
    if (len > 0) {
	### TODO: add adjust, pcut etc.
	print(paste("Identified ", len, " sign. features."))
	obj@dataSign <- obj@data[which(rownames(obj@data) %in% obj@signID[which(obj@signPAdj < pCut & obj@signP2 < pCut )]),,drop=F]
    } else {
	warning("No sign. features found!")
    }
    return(obj)
}

#' @title Crossvalidate linear model
#' 
#' @description Crossvalidation, using features as dependent 
#' and the variable of interest as independent variables. Allows
#' for adjustment of covariates 
#'
#' @param obj fitModel instance
#' @param frm right side of the model formula as character (e.g. "~COVAR1+VAR")
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import DAAG
#'
#' @export
cv <- function(obj, frm=NULL, ...) {
    dat2 <- data.frame(t(obj@dataSign), obj@meta)

    ### start parallel processing
    no_cores <- parallel::detectCores() - 1    
    no_cores <- ifelse(no_cores == 0, 1, no_cores)    
    doParallel::registerDoParallel(no_cores)    

    ## TODO: tryCatch!
    ## TODO:_ covariate frm
    if (is.null(frm)) {
	frm <- paste("~", obj@var)
    }
    mss <- foreach(i=1:length(obj@dataSign[,1])) %dopar% {
	frm <- as.formula(paste(colnames(dat2)[i], frm))
	fit <- lm(frm, data=dat2)
	cv <- cv.lm(dat2, fit, plotit=F)
	ms <- attr(cv, "ms")
	data.frame(ms, colnames(dat2)[i])
    }
    mss <- data.frame(do.call(rbind, mss), check.names=F)
    colnames(mss)[2] <- "ID"
    #TODO: check empty
    mss <- mss[order(mss[,1]),,drop=F]
    obj@cv <- mss
    
    doParallel::stopImplicitCluster()
    return(obj)
}


#' @title Plot pred / trained 
#' 
#' @description Plot predictions / trained data
#' 
#' @param obj fitModel instance
#' 
#' @import ggplot2
#' @export
plotModel <- function(obj, ...) {
    df <- data.frame(obj@pred, GRP=obj@meta[,obj@var])
    df$x <- rank(df$fit, ties.method="first")
    if (obj@type == "lr") {
	g <- ggplot(df, aes(y=fit, x=x, col=factor(GRP))) + geom_point()  + geom_errorbar(aes(ymin=lwr, ymax=upr))
    } else {
	g <- ggplot(df, aes(y=fit, x=x, col=factor(GRP))) + geom_point()  + geom_errorbar(aes(ymin=lwr, ymax=upr))
	#g <- ggplot(df[order(df$GRP),,drop=F], aes(y=fit, x=GRP)) + geom_point()  + geom_errorbar(aes(ymin=lwr, ymax=upr))
    }
    return(g)
}
