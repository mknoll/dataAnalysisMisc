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
plotForestSubgroup <- function(srv, data, ref=NULL, title="", 
			       col=c("royalblue", "darkblue", "royalblue")) {
    uv <- list()

    ### all
    w <- which(!is.na(srv))
    fit <- coxph(srv[w]~data[w,ref])
    tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int, fit$n, fit$nevent)
    j <- 1
    uv [[length(uv)+1]] <- data.frame(name1="All", 
				      name2=NA,
				      HR=tbl[j,2],
				      LOW=tbl[j, 8],
				      UP=tbl[j, 9],
				      PVAL=tbl[j, 5],
				      N=tbl[j,10],
				      NEVENT=tbl[j,11])

    for (i in 1:length(data[1,])) {
	if (colnames(data)[i] == ref) { next }
	uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					  name2=NA,
					  HR=NA, 
					  LOW=NA,
					  UP=NA, 
					  PVAL=NA,
					  N=NA,
					  NEVENT=NA)
	vals <- unique(data[,i])
	if (length(vals) ==0 ) { next }
	vals <-vals[which(!is.na(vals))]
	print(vals)
	for (v in vals) {
	    tryCatch({
		w <- which(data[,i] == v & !is.na(srv))
		fit <- coxph(srv[w]~data[w,ref])
		tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int, fit$n, fit$nevent)

		j <- 1
		uv [[length(uv)+1]] <- data.frame(name1=NA, 
						  name2=v,
						  HR=tbl[j,2],
						  LOW=tbl[j, 8],
						  UP=tbl[j, 9],
						  PVAL=tbl[j, 5],
						  N=tbl[j,10],
						  NEVENT=tbl[j,11])
	    }, error=function(e) {})
	}
    }
    uv <- do.call(rbind, uv)

    tabletext<-cbind(c("", as.character(uv[,1])),
		     c("", as.character(uv[,2])),
		     c("Hazard Ratio", round(uv[,3],2)),
		     c("95% CI", ifelse(uv[,4] == "", "", 
					paste(format(round(uv[,4],2), nsmall=2), "-", 
					      format(round(uv[,5],2), nsmall=2), sep=""))),
		     c("p-value", ifelse(round(uv[,6],3) == 0, "<0.001", round(uv[,6], 3))),
		     c("n/nevent", paste(uv[,7], "/", uv[,8], sep=""))
		     )
    ## n/nevent
    #tabletext[2:(length(tabletext[,1])-1),6] <- tabletext[3:length(tabletext[,1]),6]

    for (i in 1:length(tabletext[,1])) {
	tabletext[i,1] <- paste(tabletext[i,1], tabletext[i,2], collapse="   ")
	tabletext[i,1] <- gsub("NA", "", tabletext[i,1])
	tabletext[i,6] <- gsub("NA/NA", "", tabletext[i,6])
	#if (!is.na(tabletext[i,5]) && i > 1) { tabletext[i,6] <- "" }
    }
    tabletext <- tabletext[,-2]
    tabletext[,3] <- gsub("NA-NA", "", tabletext[,3])

    forestplot(tabletext,
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
