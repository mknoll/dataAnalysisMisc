#' @title Create forest plot for multivariate analysis
#'
#' @import forestplot
#' @import survival
#'
#' @export
plotForestMV <- function(srv, data, selection=F) {
    uv <- list()

    fit <- survival::coxph(srv~., data=data)
    if (selection != F) {
        if (class(selection)=="numeric") {
            selVar <- c()
            for (i in 1:length(data[1,])) {
                fit <- coxph(srv~data[,i])
                if (summary(fit)$logtest['pvalue'][[1]] < selection) {
                    selVar <- c(selVar, i)
                }
            }
            print("Selected: ")
            print(colnames(data)[selVar])
            fit <- coxph(srv~., data=data[,selVar,drop=F])
        } else {
            fit <- step(fit, direction=selection)
        }
    }
    tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int)

    for (i in 1:length(data[1,])) {
	if (!any(grepl(colnames(data)[i], rownames(tbl)))) { next }
	if (class(data[,i]) %in% c("factor", "character")) {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA)
	    sub <- tbl[which(grepl(colnames(data)[i], rownames(tbl))),,drop=F]
	    for (j in 1:length(sub[,1])) {
		var  <- gsub(colnames(data)[i], "", rownames(sub)[j])
		uv [[length(uv)+1]] <- data.frame(name1=NA,
						  name2=var,
						  HR=sub[j,2], 
						  LOW=sub[j,8],
						  UP=sub[j,9], 
						  PVAL=sub[j,5])
	    }
	} else if (class(data[,i]) == "numeric") {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA)
	    sub <- tbl[which(grepl(colnames(data)[i], rownames(tbl))),,drop=F]
	    j<-1
	    uv [[length(uv)+1]] <- data.frame(name1=NA,
					      name2=NA,
					      HR=sub[j,2], 
					      LOW=sub[j,8],
					      UP=sub[j,9], 
					      PVAL=sub[j,5])
	}
    }
    uv <- do.call(rbind, uv)
    tabletext<-cbind(c("", as.character(uv[,1])),
		     c("", as.character(uv[,2])),
		     c("Hazard Ratio", round(uv[,3],2)),
		     c("95% CI", paste(round(uv[,4],2), "-", round(uv[,5],2), sep="")),
		     c("p-value", round(uv[,6],3))
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
      title="",
      is.summary=c(rep(FALSE,length(tabletext[,1]))),
      clip=c(0.1,3.2),
      xlog=F,
      col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
      align=1,
      zero=1)

    return(uv)
}
