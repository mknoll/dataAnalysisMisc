#' @title Plots a Kaplan Meier surivival curve
#' 
#' @import survival
#' @import pals
#'
#' @export
plotKM <- function(srv, grp, xlim=NULL, col=NULL, xyleg=NULL, offsetNRisk=-0.2, pval=NULL, ...) {
	# number of groups 
	nGrp <- length(levels(factor(grp)))

	if (is.null(xlim)) { xlim <- c(0, max(as.numeric(srv))) }
	if (is.null(col)) { col <- kelly()[-1] }
	if (is.null(xyleg)) { xyleg <- c(xlim[2]*0.7, 0.8) }
    if (is.null(pval)) { pval <- c(xlim[2]*0.7, 0.2) }
	
	# margin for number at risk table
	par(mar=c(7+nGrp*2,4,2,2))

	# plot KM
	plot(survfit(srv~grp), mark.time=T, xlim=xlim, col=col, ...)

	# plot median os 
	fit <- survfit(srv~grp)
    if (length(data.frame(summary(fit)$table)[1,]) == 1) {
	    tbl <- data.frame(t(summary(fit)$table))
    } else {
	    tbl <- data.frame(summary(fit)$table)
    }
	for (i in 1:length(tbl[,1])) {
		segments(0, 0.5, tbl[i,"median"], 0.5, lty=3, col=col[i])
		segments(tbl[i,"median"], 0.5, tbl[i,"median"], -0.1, lty=3, col=col[i])
	}

	# legende
	legend(xyleg[1], xyleg[2], levels(factor(grp)), fill=col, bty='n', cex=0.8)

    # pvalue
    if (nGrp >= 2) {
        fitC <- coxph(srv~grp)
        p <- summary(fitC)$logtest[3][[1]]
        p <- paste("p=",format(p, scientific=TRUE, digits=3), sep="")
        text(pval[1], pval[2], p, font=2, cex=0.8)
    }

	# create no at risk table 
	par(xpd=T)
	text(-0.25, -0.15+offsetNRisk, "No. at risk", font=2, cex=0.8)
	# times to show numbers at 
	delta <- xlim[2]/4
	pos <- seq(0, xlim[2], delta)
	# get numbers per group
	tbl <- list()
	for (gr in levels(factor(grp))) {
		fit <- survfit(srv[which(grp == gr)]~grp[which(grp == gr)])
		n <- c()
		for (p in pos) {
			nCur <- 0
			tryCatch({ nCur <- fit$n.risk[which(fit$time >= p)][1] }, 
				error=function(e) {} )
			n <- c(n, nCur)
		}
		tbl[[length(tbl)+1]] <- n
	}
	tbl <- do.call(rbind, tbl)
	rownames(tbl) <- levels(factor(grp))

	# plot table 
	for (i in 1:length(tbl[,1])) {
		txt <- c(rownames(tbl)[i], tbl[i,])
		text(y=offsetNRisk-0.15-(i*0.1), x=c(-delta*0.33, pos), txt, cex=0.8)	
	}

	par(xpd=F)
}
