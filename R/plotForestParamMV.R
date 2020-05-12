#' @title Calculates multivariate analysis and creates a forest plot
#' 
#' @description Calculates a multivariate analysis, either utilizing all
#' provided variables supplied in the data data.frame (default), using 
#' the step function for model selection to obtain variables (selection 
#' takes values expected by the direction parameter of step) or retaining 
#' only variables which show a likelhood ratio p-value below a supplied 
#' cutoff in the univariate analysis. Errors might occur if the graphic
#' devices dimension is too small (foresplot() fails).
#'
#' @param srv Survival object as created by survival::Surv() function,
#' each observation is linked to one row of the data parameter
#' @param data data.frame containing all variables which will be analyzed.
#' The class of each column determined the type of analysis: numeric cols 
#' will be treated as continous variable, factor and character as factors.
#' @param subject vector identifying independent subjects. Does not work 
#' with automatic model selection.
#' @param title Plot title
#' @param col Color vector as expected by the forestplot() function
#' @param recalc caluclate LRT per factor
#' 
#' @import forestplot
#' @import survival
#'
#' @export
#'
#' @examples 
#' require(survival)
#' times <- c(100, 87, 96, 20)
#' status <- c(1,1,0,1)
#' srv <- Surv(times, status)
#' data <- data.frame(Surgery=c("yes","yes","no","no"),
#'			Drug=c("no","yes","yes","yes"),
#'			Sex=c("M","F","F","F"),
#'			Age=c(60,65,50,75))
#' subjectIDs <- c(1,2,3,3)
#' 
#' #Use all variables
#' #plotForestMV(srv,data)
#' 
#' #Univariate Cutoff of p-val < 0.2
#' #plotForestMV(srv, data, selection=0.2)
#'
#' #Automatic modell selection
#' #plotForestMV(srv, data, selection="both")
#' 
#' #Observatons from the same individual
#' #plotForestMV(srv, data, subject=subjectIDs)
plotForestParamMV <- function(srv, data, subject=NULL, selection=F, title="",  
			      col=c("royalblue", "darkblue", "royalblue"), 
			      dist="weibull", recalc=F) {
    uv <- list()

    ## any srv time <= 0?
    if (any(as.numeric(srv)[1:length(srv)] <= 0)) {
	warning("srv <0 0! Removing!")
	w <- which(as.numeric(srv)[1:length(srv)] <= 0)
	srv <-srv[-w]
	data <- data[-w,,drop=F]
	if (!is.null(subject)) { 
	    subject <- subject[-w]
	}
    }
    
    #preserve level names
    for (i in 1:length(data[1,])) {
	if (class(data[,i]) == "factor") { 
	    #lvs <- levels(factor(data[,i]))
	    #data[,i] <- as.character(data[,i])
	    #data[,i] <- factor(data[,i], levels=lvs)
	}
    }

    if (is.null(subject)) {
	fit <- survreg(srv~., data=data, dist=dist)
    } else {
	frm <- as.formula(paste("srv~", paste(colnames(data), collapse="+"), "+cluster(subject)"))
	fit <- survreg(frm, data=data, dist=dist)
    }
    if (selection != F) {
	stop("DEFUNC")
    }
    if (is.null(subject)) {
	tbl <- data.frame(int(fit, dist=dist), 
                                   N=summary(fit)$n, 
                                   summary(fit)$table[-c(1,length(summary(fit)$table[,1])),,drop=F])
    } else {
	rmI <- c(1, length(summary(fit)$table[,1]))
	tbl <- cbind(int(fit), summary(fit)$table[-rmI,,drop=F])
    }
    rownames(tbl) <- gsub("`", "", rownames(tbl))

    for (i in 1:length(data[1,])) {
	if (!any(grepl(colnames(data)[i], rownames(tbl)))) { next }

	### do we need to recalculate p-values?
	recalcP <- NULL
	if (recalc) {
	    if (is.null(subject)) {
		frm0 <- as.formula(paste("srv~", colnames(data)[-which(colnames(data) == colnames(data)[i])], collapse="+"))
		fit0 <- survreg(frm0, data=data, dist=dist)
		fit1 <- survreg(srv~., data=data, dist=dist)
		a <- anova(fit0, fit1)
		recalcP <- a[2,7]
	    } else {
		frm <- as.formula(paste("srv~", paste(colnames(data)[-which(colnames(data) == colnames(data)[i])],
						      collapse="+"), "+cluster(subject)"))
		fit0 <- survreg(frm, data=data, dist=dist)
		frm <- as.formula(paste("srv~", paste(colnames(data), collapse="+"), "+cluster(subject)"))
		fit1 <- survreg(frm, data=data, dist=dist)
		a <- anova(fit0, fit1)
		recalcP <- a[2,7]
	    }
	}

	if (class(data[,i]) %in% c("factor", "character")) {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA)
	    w <- which(substr(rownames(tbl), 1, nchar(colnames(data)[i])) == colnames(data)[i])
	    sub <- tbl[w,,drop=F]

	    for (j in 1:length(sub[,1])) {
		var <- substr(rownames(sub)[j], nchar(colnames(data)[i])+1, nchar(rownames(sub)[j]))

		uv [[length(uv)+1]] <- data.frame(name1=NA,
						  name2=var,
						  HR=sub[j,2], 
						  LOW=sub[j,1],
						  UP=sub[j,3], 
						  PVAL=ifelse(is.null(recalcP), sub[j,8],
							      ifelse(j==1,recalcP, NA)))
	    }
	} else if (class(data[,i]) == "numeric") {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA)
	    w <- which(substr(rownames(tbl), 1, nchar(colnames(data)[i])) == colnames(data)[i])
	    sub <- tbl[w,,drop=F]

	    j<-1
	    uv [[length(uv)+1]] <- data.frame(name1=NA,
					      name2=NA,
					      HR=sub[j,2], 
					      LOW=sub[j,1],
					      UP=sub[j,3], 
					      PVAL=ifelse(is.null(recalcP), sub[j,8], recalcP))
	}
    }
    uv <- do.call(rbind, uv)

    if (dist =="weibull") {    
	dN <- "Hazard Ratio"    
    } else if (dist == "loglog") {    
	dN <- "Odds Ratio"    
    }    

    tabletext<-cbind(c(paste(summary(fit)$n), as.character(uv[,1])),
		     c("", as.character(uv[,2])),
		     c(dN, round(uv[,3],2)),
		     c("95% CI", ifelse(uv[,4] == "", "", 
					paste(format(round(uv[,4],2), nsmall=2), "-", 
					      format(round(uv[,5],2), nsmall=2), sep=""))),
		     c("p-value", ifelse(round(uv[,6],3) == 0, "<0.001", round(uv[,6],3)))
		     )
    for (i in 1:length(tabletext[,1])) {
	tabletext[i,1] <- paste(tabletext[i,1], tabletext[i,2], collapse="   ")
	tabletext[i,1] <- gsub("NA", "", tabletext[i,1])
    }
    tabletext <- tabletext[,-2]
    tabletext[,3] <- gsub("NA-NA", "", tabletext[,3])

    forestplot::forestplot(tabletext,
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
