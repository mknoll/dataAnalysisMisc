#' @title Calculated univariate analysis and creates a forest plot
#'
#' @description The function creates a forest plot for a given 
#' number of variables, expect a srv object and a data.frame containing 
#' the selected variables as columns. Univariate Cox PH models 
#' are fitted. A subject vector can be specified to allow for the 
#' analysis of multiple observations per patient (e.g. paired samples),
#' by using marginal model [cluster(subject)]. Errors might occur if the graphic
#' devices dimension is too small (foresplot() fails).
#'
#' @param srv Survival object as created by survival::Surv() function,
#' each observation is linked to one row of the data parameter
#' @param data data.frame containing all variables which will be analyzed.
#' The class of each column determined the type of analysis: numeric cols 
#' will be treated as continous variable, factor and character as factors.
#' @param subject vector identifying independent subjects
#' @param title Plot title
#' @param col Color vector as expected by the forestplot() function
#' @param invalCut Cutoff to set HR, CI and p-values to empty values
#' if HR exceeds the provided cutoff (e.g. if models do not converge)
#' @param removeInval Retain as invalid identified levels (invalCut)
#' 
#' @import forestplot
#' @import survival
#'
#' @export
plotForestParam <- function(srv, data, subject=NULL, title="", col=c("royalblue", "darkblue", "royalblue"), 
		       invalCut=100, removeInval=F, dist="weibull") {
    uv <- list()
    for (i in 1:length(data[1,])) {
	print(colnames(data)[i])
	# Add variable
	# factors
	if (class(data[,i]) %in% c("factor", "character")) {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA,
					      N=NA)
	    if (is.null(subject)) {
		w <- which(!is.na(data[,i]) & !is.na(srv))
		fit <- survreg(srv[w]~factor(data[w,i]), dist=dist)
		tbl <- data.frame(int(fit, dist=dist), 
				  N=summary(fit)$n, 
				  summary(fit)$table[-c(1,length(summary(fit)$table[,1])),,drop=F])
	    } else {
		if (any(as.numeric(srv)[1:length(srv)] ==0)) { 
			warning("Found 0 time in survival object. Removing!")
		}
		w <- which(!is.na(data[,i]) & !is.na(srv) & as.numeric(srv)[1:length(srv)] > 0)
		fit <- survreg(srv[w]~factor(data[w,i])+cluster(subject[w]), dist=dist)
		rmI <- c(1, length(summary(fit)$table[,1]))    
		tbl <- cbind(int(fit), 
			     N=summary(fit)$n, 
			     summary(fit)$table[-rmI,-3,drop=F])    
	    }
	    rownames(tbl) <- substr(rownames(tbl), 19, nchar(rownames(tbl)))
	    for (j in 1:length(tbl[,1])) {
		uv [[length(uv)+1]] <- data.frame(name1=NA, 
						  name2=rownames(tbl)[j],
						  HR=tbl[j,2],
						  LOW=tbl[j, 1],
						  UP=tbl[j, 3],
						  PVAL=tbl[j, 8],
						  N=tbl[j,4])
	    }
	} else if (class(data[,i]) %in% c("numeric", "integer")) {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA,
					      N=NA)
	    if (is.null(subject)) {
		w <- which(!is.na(data[,i]) & !is.na(srv))
		fit <- survreg(srv[w]~data[w,i], dist=dist)
		tbl <- data.frame(int(fit, dist=dist), 
				  N=summary(fit)$n, 
				  summary(fit)$table[-c(1,length(summary(fit)$table[,1])),,drop=F])
	    } else {
		if (any(as.numeric(srv)[1:length(srv)] ==0)) { 
			warning("Found 0 time in survival object. Removing!")
		}
		w <- which(!is.na(data[,i]) & !is.na(srv) & as.numeric(srv)[1:length(srv)] > 0)
		fit <- survreg(srv[w]~data[w,i]+cluster(subject[w]), dist=dist)
		rmI <- c(1, length(summary(fit)$table[,1]))    
		tbl <- cbind(int(fit), 
			     N=summary(fit)$n, 
			     summary(fit)$table[-rmI,-3,drop=F])    
	    }
	    j<-1
	    uv [[length(uv)+1]] <- data.frame(name1=NA, 
					      name2=NA,
					      HR=tbl[j,2],
					      LOW=tbl[j, 1],
					      UP=tbl[j, 3],
					      PVAL=tbl[j, 8],
					      N=tbl[j,4])
	} else {
	    warning(paste("Could not process ", colnames(data)[i]))
	}
    }

    uv <- do.call(rbind, uv)
    ## set invaldi data to NA
    if (F) {
	w <- which(uv[,3] > invalCut)
	if (length(w) >0) {
	    uv[w,c(3:6)] <- NA
	    if (removeInval) {
		uv <- uv[-w,,drop=F]
	    }
	}
    }

    if (dist =="weibull") {
	dN <- "Hazard Ratio"
    } else if (dist == "loglog") {
	dN <- "Odds Ratio"
    }
    tabletext<-cbind(c("", as.character(uv[,1])),
		     c("", as.character(uv[,2])),
		     c(dN, round(uv[,3],2)),
		     c("95% CI", ifelse(uv[,4] == "", "", 
					paste(format(round(uv[,4],2), nsmall=2), "-", format(round(uv[,5],2), nsmall=2), sep=""))),
		     c("p-value", ifelse(round(uv[,6],3) == 0, "<0.001",round(uv[,6],3) )),
		     c("n", paste(uv[,7], sep=""))
		     )
    ## n/nevent
    tabletext[2:(length(tabletext[,1])-1),6] <- tabletext[3:length(tabletext[,1]),6]

    for (i in 1:length(tabletext[,1])) {
	tabletext[i,1] <- paste(tabletext[i,1], tabletext[i,2], collapse="   ")
	tabletext[i,1] <- gsub("NA", "", tabletext[i,1])
	tabletext[i,6] <- gsub("NA/NA", "", tabletext[i,6])
	if (!is.na(tabletext[i,5]) && i > 1) { tabletext[i,6] <- "" }
    }
    tabletext <- tabletext[,-2]
    tabletext[,3] <- gsub("NA-NA", "", tabletext[,3])

    ### boldprint 
    bp <- list()
    for (i in 1:length(tabletext[,1])) {
	bp[[i]] <-list()
	for (j in 1:length(tabletext[1,])) {
	    if (j == 4) {
		if (!is.na(as.numeric(tabletext[i,j])) && (as.numeric(tabletext[i,j]) < 0.05)) {
		    bp[[i]][[j]] <- gpar(fontface="bold")
		} else if (!is.na(tabletext[i,j]) && tabletext[i,j] == "<0.001") {
		    bp[[i]][[j]] <- gpar(fontface="bold")
		} else {
		    bp[[i]][[j]] <- gpar(fontface="plain")
		}
	    } else {
		bp[[i]][[j]] <- gpar(fontface="plain")
		#bp[[i]][[j]] <- gpar(fontface="bold")
	    }
	}
    }

    forestplot(tabletext,
	       txt_gp=fpTxtGp(label=bp),
	       mean  = c(NA, as.numeric(as.character(uv[,3]))),
	       lower = c(NA, as.numeric(as.character(uv[,4]))),
	       upper = c(NA, as.numeric(as.character(uv[,5]))),
	       new_page = TRUE,
	       title=title,
	       is.summary=c(rep(FALSE,length(tabletext[,1]))),
	       clip=c(0.1,3.2),
	       xlog=F,
	       col=fpColors(box=col[1],line=col[2], summary=col[3]),
	       align=1,
	       zero=1)

    return(uv)
}


#' @title Intervals for parametric dist
#' @import SurvRegCensCov
#' @export
int <- function(fit, dist="loglog", z=1.96) {
    ret <- NULL
    if (dist == "loglog") {
	alpha1 <- summary(fit)$table[,1]
	lwr <- alpha1 - summary(fit)$table[,2]*z
	upr <- alpha1 + summary(fit)$table[,2]*z
	alpha1 <- fit$coefficient[]
	p <- 1/fit$scale
	v1 <- exp(-alpha1*p)
	v2 <- exp(-lwr*p)
	v3 <- exp(-upr*p)

	len <- length(v1)
	len <- ifelse(length(v2) < len, length(v2), len)
	len <- ifelse(length(v3) < len, length(v3), len)

	ret <- data.frame(v1[1:len], v2[1:len],v3[1:len])
	o <- order(ret[1,])
	ret <- ret[,o]
	colnames(ret) <- c("LOW","EST","UP")
	ret <- ret[-which(grepl("Intercept", rownames(ret))),,drop=F]
    } else if (dist == "weibull") {
	ret <- ConvertWeibull(fit)$HR
	o <- order(ret[1,])
	ret <- ret[,o,drop=F]
    }
    return(ret)
}
