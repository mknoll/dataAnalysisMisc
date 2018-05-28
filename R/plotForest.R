#' @title Create forest plot
#'
#' @import forestplot
#'
#' @export
plotForest <- function(srv, data) {
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
					      PVAL=NA)
	    fit <- coxph(srv~factor(data[,i]))
	    tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int)
	    rownames(tbl) <- substr(rownames(tbl), 18, nchar(rownames(tbl)))
	    for (j in 1:length(tbl[,1])) {
		uv [[length(uv)+1]] <- data.frame(name1=NA, 
						  name2=rownames(tbl)[j],
						  HR=tbl[j,2],
						  LOW=tbl[j, 8],
						  UP=tbl[j, 9],
						  PVAL=tbl[j, 5])

	    }
	} else if (class(data[,i]) == "numeric") {
	    uv [[length(uv)+1]] <- data.frame(name1=colnames(data)[i],
					      name2=NA,
					      HR=NA, 
					      LOW=NA,
					      UP=NA, 
					      PVAL=NA)
	    fit <- coxph(srv~data[,i])
	    tbl <- cbind(summary(fit)$coef, summary(fit)$conf.int)
	    j<-1
	    uv [[length(uv)+1]] <- data.frame(name1=NA, 
						  name2=NA,
						  HR=tbl[j,2],
						  LOW=tbl[j, 8],
						  UP=tbl[j, 9],
						  PVAL=tbl[j, 5])
	}
    }

    uv <- do.call(rbind, uv)
    print(uv)
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

    forestplot(tabletext,
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
}
