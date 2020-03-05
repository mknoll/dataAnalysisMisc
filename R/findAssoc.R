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
    ##########
    grp <- as.character(grp)

    out <- list()
    lvs <- as.character(levels(factor(grp)))

    ##### Heading
    df <- data.frame(t(c("Group", NA, lvs, NA)))
    out[[length(out)+1]] <- df
    df <- data.frame(t(c(NA, NA, table(grp)[match(lvs, names(table(grp)))], NA)))
    out[[length(out)+1]] <- df

    #### Fill Table 
    for (i in 1:length(data[1,])) {
	sub <-  list()
	#sub[[length(sub)+1]] <- data.frame(name0=colnames(data)[i], name1=NA, name2=NA, name3=NA, p=NA)
	sub[[length(sub)+1]] <- data.frame(t(c(colnames(data)[i], NA, rep(NA, length(lvs)), NA)))
	p.val <- NA

	if (class(data[,i]) %in% c("numeric","integer")) {
	    #### Continuous data
	    #####################
	    vals <- data[,i]
	    p.val <- NA
	    tryCatch({
		if (is.null(subject)) {
		    ### non paired data
		    df <- data.frame(VAL=vals, GRP=factor(grp))
		    w <- which(!is.na(df$VAL) & !is.na(df$GRP))
		    fit0 <- lm(VAL~1, data=df[w,])
		    fit <- lm(VAL~GRP, data=df[w,])
		    p.val <- anova(fit0, fit)[2,6]
		} else {
		    df <- data.frame(VAL=vals, GRP=grp, ID=subject)
		    fit0 <- lmer(VAL~1 + (1|ID), data=df, REML=F)
		    fit <- lmer(VAL~GRP + (1|ID), data=df, REML=F)
		    p.val <- anova(fit0, fit)[2,8]
		}
	    }, error=function(e) { print(e) })

	    for (cp in contParam) {
		eval(parse(text=paste("fun=",cp)))
		valsTmp <- c()
		for (lv in lvs) {
		    v <- as.character(round(apply(data.frame(t(vals[which(grp == lv & !is.na(vals))])),1,fun),2))
		    if (cp == "range") {
			valsTmp <- c(valsTmp, paste(" [",v[1],";",v[2],"]", sep=""))
		    } else {
			valsTmp <- c(valsTmp, v)
		    }
		}
		df <- data.frame(t(c(NA, cp, valsTmp, NA)))
		sub[[length(sub)+1]] <- df
	    }

	} else {
	    ##################################
	    ### factors 
	    vals <- droplevels(factor(data[,i]))
	    ##only one factor: skip
	    ## TODO: -> wilson confidence interval
	    if (length(unique(vals)) == 1) {  next }

	    ## two factors: use barnard for nonpaired data
	    force <- F
	    if (length(unique(vals)) == 2 ) {
		if (length(unique(grp)) == 2) {
		    ## assure that we have a 2x2 table
		    tbl <- table(vals, grp)
		    if (is.null(subject)) {
			## non paired
			#Barnard two sided
			p.val <- NA
			tryCatch({
			    p.val <- barnard.test(tbl[1,1], tbl[1,2], tbl[2,1], tbl[2,2])$p.value[2]  
			}, error=function(e) {  } )
		    } else {
			force <- T
		    }
		} else {
		    force <- T
		}
	    } 
	    if (length(unique(vals)) > 2 || force) {
		tbl <- table(vals, grp)
		p.val <- NA
		if (!is.null(subject)) {
		    ## paired data
		    df <- data.frame(VAL=vals, GRP=grp, ID=subject)
		    fit <- MCMCglmm(VAL~GRP,  random=~ID, data=df, family="ordinal")
		    p.val <- summary(fit)$solutions[2,5]
		} else {
		    ## non paired data
		    if (kat == "Fisher") {
			tryCatch({
			    p.val <- fisher.test(tbl)$p.value
			}, error=function(e) { })
		    } else {
			p.val <- NA
			tryCatch({
			    p.val <- chisq.test(tbl)$p.value
			}, error=function(e) { })
		    }
		}
	    }
	    for (j in 1:length(tbl[,1])) {
		tblNM <- c()
		for (k in 1:length(tbl[1,])) {
		    #tblNM <- c(tblNM, as.character(tbl[j,k]))
		    add <- paste(tbl[j,k], " (", round(tbl[j,k]/sum(tbl[,k])*100), ")",sep="")
		    tblNM <- c(tblNM, add)
		}
		sub[[length(sub)+1]] <- data.frame(t(c(name0=NA, name1=rownames(tbl)[j], 
						   tblNM, 
						   p=NA)))
	    }
	}
	sub[[1]][length(sub[[1]][1,])] <- as.character(round(p.val, 3))
	out <- c(out, sub)
    }

    cn <- c("Feature","", paste("Level", lvs,sep=""), "p-value")
    out <- lapply(out, function(x) {
	       colnames(x) <- cn
	       x <- apply(x, 2, function(x) as.character(x))
	       x
			})
    pat <- do.call(rbind, out)
    #colnames(pat) <- c("Feature","", "Level1", "Level2", "p-value")

    if (!latex) {
        return(data.frame(pat))
    } else {		
	return (list(tbl=pat, pdffile=preparePdf(pat, outdir, col=paste("ll",paste(rep("c", length(lvs)),collapse=""),"r", sep=""),filename=filename))) 
    }
}
