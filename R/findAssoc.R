#' @title Find associations 
#' @description find associations between a given groupin
#' and a set of variables
#' @import Barnard
#' @import lme4
#' @import stargazer
#' @import MCMCglmm
#' @export
findAssoc <- function(grp, data, test=NULL, kat="Fisher", filename=NULL,
		      contParam=c("mean", "sd"), latex=F, outdir="../reports/",
		      subject=NULL) {
    ###
    grp <- as.character(grp)
    #if (length(unique(grp)) != 2) {
#	warning("Sry, currently only two group tests are implmented!")
#	return(NULL)
 #   }

    out <- list()
    lvs <- as.character(levels(factor(grp)))
    out[[length(out)+1]] <- data.frame(name0="Group", name1=NA, name2=lvs[1], name3=lvs[2], p=NA)
    out[[length(out)+1]] <- data.frame(name0=NA, name1=NA, name2=as.character(length(which(grp == lvs[1]))), 
						 name3=as.character(length(which(grp == lvs[2]))), p=NA)

    ###
    for (i in 1:length(data[1,])) {
	sub <-  list()
	sub[[length(sub)+1]] <- data.frame(name0=colnames(data)[i], name1=NA, name2=NA, name3=NA, p=NA)
	p.val <- NA
	if (class(data[,i]) %in% c("numeric","integer")) {
	    vals <- data[,i]
	    p.val <- NA
	    tryCatch({
		if (is.null(subject)) {
		    ### non paired data
		    p.val <- t.test(vals~grp)$p.value
		} else {
		    df <- data.frame(VAL=vals, GRP=grp, ID=subject)
		    fit0 <- lme(VAL~1 + (1|ID), data=df, method="ML")
		    fit <- lme(VAL~GRP + (1|ID), data=df, method="ML")
		    p.val <- anova(fit0, fit)[2,9]
		}
	    }, error=function(e) {})
	    if ("mean" %in% contParam) {
		sub[[length(sub)+1]] <- data.frame(name0=NA, name1="Mean", name2=as.character(round(mean(vals[which(grp == lvs[1])], na.rm=T)),2),
						   name3=as.character(round(mean(vals[which(grp == lvs[2])], na.rm=T),2)), p=NA)
	    }
	    if ("sd" %in% contParam) {
		sub[[length(sub)+1]] <- data.frame(name0=NA, name1="SD", name2=as.character(round(sd(vals[which(grp == lvs[1])], na.rm=T)),2),
						   name3=as.character(round(sd(vals[which(grp == lvs[2])], na.rm=T),2)), p=NA)
	    }
	    if ("median" %in% contParam) {
		sub[[length(sub)+1]] <- data.frame(name0=NA, name1="Median", name2=as.character(round(median(vals[which(grp == lvs[1])], na.rm=T)),2),
						   name3=as.character(round(median(vals[which(grp == lvs[2])], na.rm=T),2)), p=NA)
	    }
	    if ("min" %in% contParam) {
		sub[[length(sub)+1]] <- data.frame(name0=NA, name1="Min", name2=as.character(round(min(vals[which(grp == lvs[1])], na.rm=T)),2),
						   name3=as.character(round(min(vals[which(grp == lvs[2])], na.rm=T),2)), p=NA)
	    }
	    if ("max" %in% contParam) {
		sub[[length(sub)+1]] <- data.frame(name0=NA, name1="Max", name2=as.character(round(max(vals[which(grp == lvs[1])], na.rm=T)),2),
						   name3=as.character(round(max(vals[which(grp == lvs[2])], na.rm=T),2)), p=NA)
	    }
	} else {
	    ### factors 
	    vals <- droplevels(factor(data[,i]))
	    ##only one factor: skip
	    if (length(unique(vals)) == 1) { next }

	    ## two factors: use barnard for nonpaired data
	    force <- F
	    if (length(unique(vals)) == 2) {
		tbl <- table(vals, grp)
		if (is.null(subject)) {
		    ## non paired
		    #Barnard two sided
		    p.val <- NA
		    tryCatch({
			p.val <- barnard.test(tbl[1,1], tbl[1,2], tbl[2,1], tbl[2,2])$p.value[2]  
		    }, error=function(e) { } )
		} else {
		    force <- T
		}
		sub[[length(sub)+1]] <- data.frame(name0=NA, name1=rownames(tbl)[1], 
						   name2=as.character(tbl[1,1]), name3=as.character(tbl[1,2]),p=NA)
		sub[[length(sub)+1]] <- data.frame(name0=NA, name1=rownames(tbl)[2], 
						   name2=as.character(tbl[2,1]), name3=as.character(tbl[2,2]),p=NA)
	    } 
	    if (length(unique(vals)) > 2 || force) {
		p.val <- NA
		if (!is.null(subject)) {
		    ## paired data
		    df <- data.frame(VAL=vals, GRP=grp, ID=subject)
		    fit <- MCMCglmm(VAL~GRP,  random=~ID, data=df, family="ordinal")
		    p.val <- summary(fit)$solutions[2,5]
		} else {
		    ## non paired data
		    if (kat == "Fisher") {
			tbl <- table(vals, grp)
			tryCatch({
			    p.val <- fisher.test(tbl)$p.value
			}, error=function(e) { })
		    } else {
			tbl <- table(vals, grp)
			p.val <- NA
			tryCatch({
			    p.val <- chisq.test(tbl)$p.value
			}, error=function(e) { })
		    }
		}
		for (j in 1:length(tbl[,1])) {
		    sub[[length(sub)+1]] <- data.frame(name0=NA, name1=rownames(tbl)[j], 
						       name2=as.character(tbl[j,1]), name3=as.character(tbl[j,2]),p=NA)
		}
	    }
	}
	sub[[1]]$p <- round(p.val, 3)
	out <- c(out, sub)
    }

    pat <- do.call(rbind, out)
    colnames(pat) <- c("Feature","", "Level1", "Level2", "p-value")

    if (!latex) {
        return(pat)
    } else {		
	return (list(tbl=pat, pdffile=preparePdf(pat, outdir, col="llccr",filename=filename))) 
    }
}
