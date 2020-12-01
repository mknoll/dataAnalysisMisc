#' @title Create forest plot
#'
#' @import forestplot
#'
#' @export
plotForestSplit <- function(srv, data, subject=NULL, title="", col=c("royalblue", "darkblue", "royalblue"), split=NULL) {
    lvSplit <- levels(factor(split))
    nSplit <- length(lvSplit)

    uv <- list()
    for (i in 1:length(data[1,])) {
	# Add variable
	# factors
	if (class(data[,i]) %in% c("factor", "character")) {
	    uv [[length(uv)+1]] <- data.frame(name0=colnames(data)[i],
					      name1="",
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA)
	    for (s in 1:nSplit) {
		sel <- which(split ==  lvSplit[s])
		if (is.null(subject)) {
		    fit <- coxph(srv[sel]~factor(data[sel,i]))
		    tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int)
		} else {
		    fit <- coxph(srv[sel]~factor(data[sel,i])+cluster(subject[sel]))
		    tbl <- cbind(summary(fit)$coef[], summary(fit)$conf.int)
		    tbl <- tbl[,-4,drop=F]
		}
		rownames(tbl) <- substr(rownames(tbl), 21, nchar(rownames(tbl)))
		for (j in 1:length(tbl[,1])) {
		    uv [[length(uv)+1]] <- data.frame(name0=NA, 
						      name1=lvSplit[s],
						      name2=rownames(tbl)[j],
						      HR=tbl[j,2],
						      LOW=tbl[j, 8],
						      UP=tbl[j, 9],
						      PVAL=tbl[j, 5])
		}
	    }
	} else if (class(data[,i]) %in% c("numeric", "integer")) {
	    uv [[length(uv)+1]] <- data.frame(name0=colnames(data)[i],
					      name1=NA,
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA)
	    for (s in 1:nSplit) {
		sel <- which(split ==  lvSplit[s])
		if (is.null(subject)) {
		    fit <- coxph(srv[sel]~data[sel,i])
		    tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int)
		} else {
		    fit <- coxph(srv[sel]~data[sel,i]+cluster(subject[sel]))
		    tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int)
		    tbl <- tbl[,-4,drop=F]
		}
		j<-1
		uv [[length(uv)+1]] <- data.frame(name0=NA,
						  name1=lvSplit[s],
						  name2=NA,
						  HR=tbl[j,2],
						  LOW=tbl[j, 8],
						  UP=tbl[j, 9],
						  PVAL=tbl[j, 5])
	    }
	} else {
	    warning(paste("Could not process ", colnames(data)[i]))
	}
    }
    uv <- do.call(rbind, uv)

    tabletext<-cbind(c("", as.character(uv[,1])),
		     c("", as.character(uv[,2])),
		     c("", as.character(uv[,3])),
		     c("Hazard Ratio", round(uv[,4],2)),
		     c("95% CI", paste(round(uv[,5],2), "-", round(uv[,6],2), sep="")),
		     c("p-value", ifelse(round(uv[,7],3) ==0, "<0.001", round(uv[,7],3)))
		     )
    for (i in 1:length(tabletext[,1])) {
	tabletext[i,1] <- paste(tabletext[i,1], tabletext[i,2], tabletext[i,3], collapse="   ")
	#tabletext[i,1] <- gsub("NA", "", tabletext[i,1])
	tabletext[i,1] <- gsub(" NA", "", tabletext[i,1])
	tabletext[i,1] <- ifelse(tabletext[i,1] == "NA", "", tabletext[i,1])
	tabletext[i,1] <- gsub("NA ", "", tabletext[i,1])
    }
    tabletext <- tabletext[,-c(2:3)]
    tabletext[,3] <- gsub("NA-NA", "", tabletext[,3])

    forestplot(tabletext,
	       mean  = c(NA, as.numeric(as.character(uv[,4]))),
	       lower = c(NA, as.numeric(as.character(uv[,5]))),
	       upper = c(NA, as.numeric(as.character(uv[,6]))),
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
