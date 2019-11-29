#' @title Find two differening groups in survival data
#' 
#' @description  Tests different cutoffs on a continuous variable 
#' and calculates log-rank tests / survival differences
#' 
#' @param data data.frame, with variables in cols
#' @param srv Surv object with the corresponding data
#' @param delta proportion in which group splitting should 
#' be performed 
#' @param minGrpSize minimum group size for cutoffs
#' @param subject subject IDs if multiple measurements per subjects 
#' are to be tested
#' 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import survival
#'   
#' @export
#' 
#' @return data.frame with the calculated survival values 
#' (p.value, median survival, group sizes, cutoffsize, 
#' stepsize)
findOptCutoff_Model <- function(data, 
			  srv, 
			  frm0="srv~1",
			  frm="srv~GRP",
			  stdVar="VAL",
			  delta=0.1, minGrpSize=1,
			  subjVar=NULL,
			  dist=NULL) {
    ## TODO: test from GRP var!
    #### Environment problems i passed as formula!
    ## FIXME
    frm0 <- as.formula(frm0)
    frm <- as.formula(frm)

    ##remove NA SRV data 
    if (any(is.na(srv))) { 
	warning("Removing NA srv items!")
    }
    ## romve additional rows with na in vars
    tmp <- unlist(strsplit(as.character(frm), "~", fixed=T))
    tmp <- unlist(strsplit(tmp, "+", fixed=T))
    tmp <- trimws(unlist(strsplit(tmp, " ")))
    tmp <- tmp[which(!tmp %in% c("srv", "", "1", "+"))]
    identVar <- c(subjVar, tmp, stdVar)
    if (!is.null(subjVar)) { identVar <- c(identVar, subjVar) }
    identVar <- identVar[which(identVar != "GRP")]
    print(identVar)

    keep <- list()
    keep[[length(keep)+1]] <- !is.na(srv)
    for (var in identVar) {
	rm <- which(is.na(data[,var]))
	if(length(rm) > 0) {
	    warning(paste("Found NA in ", var))
	    keep[[length(keep)+1]] <- !is.na(data[,var])
	}
    }
    keep <- do.call(cbind, keep)
    keep <- which(apply(keep, 1, function(x) all(x == T)))

    #####
    srv <- srv[keep]
    data <- data[keep,,drop=F]

    cutoff <- max(data[,stdVar])
    step = abs(max(data[,stdVar], na.rm=T)-min(data[,stdVar], na.rm=T))*delta
    if (is.infinite(step)) {
        return(NULL)
    }

    ##calculate all cutoffs to test
    sq <- seq(from=min(data[,stdVar],na.rm=T), to=max(data[,stdVar],na.rm=T), by=step)
    sq <- sq[which(apply(data.frame(sq), 1, function(x) length(which(data[,stdVar] > x)) > minGrpSize  & length(which(data[,stdVar] < x)) > minGrpSize ))]
    if(length(sq) == 0) { 
        ## min group size not reached
        return (NULL) 
    }

    ##use parallel computing
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    ## FIXME
    if (is.na(no_cores)) { no_cores <- 4 }
    doParallel::registerDoParallel(no_cores)
    
    #.env <- environment()
    #frm0 <- as.formula(frm0, env=.env)
    #frm <- as.formula(frm, env=.env)

    out <- foreach(cutoff=sq) %dopar% {
        grp <- ifelse(as.numeric(as.character(data[,stdVar])) > cutoff, "high", "low")
        cutoff <- cutoff-step
        if (length(levels(factor(grp))) == 2) {
	    ##univariate coxph
	    grp <- factor(grp)
	    grp <- relevel(grp, "low")
	    ### TODO: test if var already in data
	    data$GRP <- grp

	    ### do we need a cluster term?
	    if (!is.null(subjVar)) {
		frm0 <- as.formula(paste(as.character(frm0), "+cluster(", subjVar, ")"))
		frm <- as.formula(paste(as.character(frm), "+cluster(", subjVar, ")"))
	    }
	    if (is.null(dist)) {
		fit0 <- coxph(frm0, data)
		fit <- coxph(frm, data)
	    } else {
		fit0 <- survreg(frm0, data, dist=dist)
		fit <- survreg(frm, data, dist=dist)
	    }
	    a <- anova(fit0, fit)
	    p.val <- a[2,4]
	    if (!is.null(dist)) {
		p.val <- a[2,7]
	    }

	    df <- NULL
	    tryCatch({
		cf <- summary(fit)$coefficients
		if (!is.null(dist)) {
		    cf <- summary(fit)$table[-1,,drop=F]
		}
		cf <- cf[which(grepl("GRP", rownames(cf))),,drop=F]
		df <- data.frame(cutoff=cutoff+step, 
				 p.val=p.val,
				 nHigh=length(grp[grp == "high"]), 
				 nLow=length(grp[grp == "low"]), 
				 cf
				 )
	    }, error=function(e) { print(e) })
	    df
	}
    }
    coll <- NULL
    tryCatch({
	coll <- do.call(rbind, out)
	coll <- coll[order(coll$p.val),]
    }, error=function(e) { print(e)  })
    
    doParallel::stopImplicitCluster()
    
    return(coll)
}
