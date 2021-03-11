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
#' @param MDPI adhere to MDPI requirements
#' 
#' @import forestplot
#' @import survival
#' @import grid
#'
#' @export
plotForest <- function(srv, data, subject=NULL, title="", col=c("royalblue", "darkblue", "royalblue"), 
		       invalCut=100, removeInval=F, MDPI=F) {
    uv <- list()
    for (i in 1:length(data[1,])) {
	# Add variable
	# factors
	if (class(data[,i]) %in% c("factor", "character")) {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA,
					      N=NA,
					      NEVENT=NA)
	    if (is.null(subject)) {
		w <- which(!is.na(data[,i]) & !is.na(srv))
		fit <- coxph(srv[w]~factor(data[w,i]))
		tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int, fit$n, fit$nevent)
	    } else {
		w <- which(!is.na(data[,i]) & !is.na(srv))
		fit <- coxph(srv[w]~factor(data[w,i])+cluster(subject[w]))
		tbl <- cbind(summary(fit)$coef[], summary(fit)$conf.int, fit$n, fit$nevent)
		tbl <- tbl[,-4,drop=F]
	    }
	    rownames(tbl) <- substr(rownames(tbl), 19, nchar(rownames(tbl)))
	    for (j in 1:length(tbl[,1])) {
		uv [[length(uv)+1]] <- data.frame(name1=NA, 
						  name2=rownames(tbl)[j],
						  HR=tbl[j,2],
						  LOW=tbl[j, 8],
						  UP=tbl[j, 9],
						  PVAL=tbl[j, 5],
						  N=tbl[j,10],
						  NEVENT=tbl[j,11])
	    }
	} else if (class(data[,i]) %in% c("numeric", "integer")) {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA,
					      N=NA,
					      NEVENT=NA)
	    if (is.null(subject)) {
		w <- which(!is.na(data[,i]) & !is.na(srv))
		fit <- coxph(srv[w]~data[w,i])
		tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int, fit$n, fit$nevent)
	    } else {
		w <- which(!is.na(data[,i]) & !is.na(srv))
		fit <- coxph(srv[w]~data[w,i]+cluster(subject[w]))
		tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int, fit$n, fit$nevent)
		tbl <- tbl[,-4,drop=F]
	    }
	    j<-1
	    uv [[length(uv)+1]] <- data.frame(name1=NA, 
					      name2=NA,
					      HR=tbl[j,2],
					      LOW=tbl[j, 8],
					      UP=tbl[j, 9],
					      PVAL=tbl[j, 5],
					      N=tbl[j,10],
					      NEVENT=tbl[j,11])
	} else {
	    warning(paste("Could not process ", colnames(data)[i]))
	}
    }
    uv <- do.call(rbind, uv)
    ## set invaldi data to NA
    w <- which(uv[,3] > invalCut)
    if (length(w) >0) {
	uv[w,c(3:6)] <- NA
	if (removeInval) {
	    uv <- uv[-w,,drop=F]
	}
    }

    dash <- ifelse(MDPI, "–", "-")
    tabletext<-cbind(c("", as.character(uv[,1])),
		     c("", as.character(uv[,2])),
		     c("Hazard Ratio", round(uv[,3],2)),
		     c("95% CI", ifelse(uv[,4] == "", "", 
					paste(format(round(uv[,4],2), nsmall=2), dash, 
					      format(round(uv[,5],2), nsmall=2), sep=""))),
		     c("p-value", ifelse(round(uv[,6],3) == 0, "<0.001", round(uv[,6], 3))),
		     c("n/nevent", paste(uv[,7], "/", uv[,8], sep=""))
		     )
    ## n/nevent
    tabletext[2:(length(tabletext[,1])-1),6] <- tabletext[3:length(tabletext[,1]),6]

    for (i in 1:length(tabletext[,1])) {
	tabletext[i,1] <- paste(tabletext[i,1], tabletext[i,2], collapse="   ")
	tabletext[i,1] <- gsub(" NA", "", tabletext[i,1])
	#tabletext[i,1] <- gsub("NA", "", tabletext[i,1])
	if (substr(tabletext[i,1], 1, 3) == "NA ") {
	    tabletext[i,1] <- substr(tabletext[i,1], 3, nchar(tabletext[i,1]))
	}
	tabletext[i,1] <- ifelse(tabletext[i,1] == "NA", "", tabletext[i,1])
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
