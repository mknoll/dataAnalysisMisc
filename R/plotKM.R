#' @title Plots a Kaplan Meier survival curve with number at risk table
#' 
#' @description Creates Kaplan-Meier survival curves, prints a number at risk
#' table and calculates & plots a Cox-PH LRT p-value. Median event times 
#' are indicted by addiditional lines.
#'
#' @param srv Survival obejct as obtained by survival::Surv()
#' @param grp grouping variable
#' @param xlim xlim range to plot. Default: 0-max time included in the srv object
#' @param xyleg vector contaning x,y coordinates to the legend
#' @param offsetNRisk offset in y dimensions to move the number at risk table
#' @param yDelta y difference between rows of the number at risk table
#' @param nRiskCat Number of breaks of the number a risk table
#' @param pval xy coordinates of where to plot the pvalue
#' @param mar margin parameter, calculated internally (default). Not used if mar 
#' is set to anything different from NULL
#' @param MDPI if set to T, adhere to MDPI requirements
#' 
#' @import survival
#' @import pals
#'
#' @export
#' 
#' @examples
#'	require(survival)
#'	time <- c(10,7,8,4,11,14)
#'	status <- c(1,1,1,0,1,0)
#'	srv <- Surv(time, status)
#' 	grp <- c(rep("A", 3), rep("B", 3))
#' 	plotKM(srv, grp)
plotKM <- function(srv, grp, xlim=NULL, col=NULL, xyleg=NULL, offsetNRisk=-0.2, 
		   yDelta=0.1, nRiskCat=4, pval=NULL, mar=NULL, subject=NULL, 
		   dist=NULL, MDPI=F, plotFit=T, adjAxis=F, ...) {
    # number of groups 
    nGrp <- length(levels(factor(grp)))
    grp <- factor(grp)

    if (is.null(xlim)) { xlim <- c(0, max(as.numeric(srv), na.rm=T)) }
    if (is.null(col)) { col <- kelly()[-1] }
    if (is.null(xyleg)) { xyleg <- c(xlim[2]*0.7, 0.8) }
    if (is.null(pval)) { pval <- c(xlim[2]*0.7, 0.2) }

    # any 0 srv time and dist?
    if (!is.null(dist)){
	if (any(is.na(srv))) {
	    w <- which(is.na(srv))
	    warning("Removing srv NA!")
	    srv <-srv[-w]
	    grp <- grp[-w]
	    if (!is.null(subject)) { subject <-subject[-w] }
	}
	if (any(as.numeric(srv)[1:length(srv)] <= 0)) {
	    w <- which(as.numeric(srv)[1:length(srv)] <= 0)
	    warning("Removing srv <=0!")
	    srv <-srv[-w]
	    grp <- grp[-w]
	    if (!is.null(subject)) { subject <-subject[-w] }
	}
    }

    # margin for number at risk table
    if (is.null(mar)) {
	par(mar=c(7+nGrp*2,4,2,2))
    }

    # plot KM
    if (!is.null(dist) && dist == "loglog") {
	if (adjAxis) {
	    plot(survfit(srv~grp), mark.time=T, xlim=xlim, col=col, xaxt='n', ...)
	} else {
	    plot(survfit(srv~grp), mark.time=T, xlim=xlim, col=col, ...)
	}
    } else {
	plot(survfit(srv~grp), mark.time=T, xlim=xlim, col=col, ...)
    }

    # plot parametric survival curves
    if (!is.null(dist)) {
	warning("repeated measurements not implemented yet! Ignoring subject parameter!")
	fitR <- survreg(srv~grp,dist=dist)
	for (i in 1:nGrp) {
	    y <- predict(fitR, newdata=list(grp=levels(factor(grp))[i]), type="quantile", p=seq(.01,.99,by=.01))
	    x <- seq(.99,.01,by=-.01)
	    if (plotFit) {
		lines(y, x, col=col[i], lty=2, ...)
	    }
	    ww <- which(x == 0.5)
	    print(paste("Median survival time [",levels(factor(grp))[i],"]: ", y[ww], sep=""))
	}
	##adjust xaxis ticks
	if (dist == "loglog" && adjAxis) {
	    fit <- survfit(srv~grp)
	    to <- floor(max(fit$time))
	    if (!missing(xlim)) {
		print(xlim)
		to <-xlim[2]
	    }
	    xat <- c(0, to)                 
	    xat <- seq(from=xat[1], to=xat[2], len=5)
	    print(xat)
	    axis(side=1, at=xat, labels = FALSE)                                          
	    text(x=xat,  par("usr")[3]-0.05,  labels=round(exp(xat)), xpd=T,pos=1, las=2)     
	}
    }

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
	if (is.null(subject)) {
	    fitC <- coxph(srv~grp)
	    p <- summary(fitC)$logtest[3][[1]]
	} else {
	    fitC <- coxph(srv~grp+cluster(subject))
	    print(str(summary(fitC)))
	    p <- summary(fitC)$robscore[3][[1]]
	}
	if (MDPI) {
	    str <- format(p, scientific=TRUE, digits=3)
	    tmp <- strsplit(str, "e")[[1]]
	    expression(tmp[1])
	    res <- paste0("expression(",tmp[1],")")
	    mystr <- paste0("p = ", tmp[1], " x 10")
	    mystr2 <- tmp[2]
	    p <- bquote(.(mystr)^.(mystr2))
	} else {
	    p <- paste("p=",format(p, scientific=TRUE, digits=3), sep="")
	}
	text(pval[1], pval[2], p, font=2, cex=0.8)
    }

    # create no at risk table 
    par(xpd=T)
    text(-0.25, -0.15+offsetNRisk, "No. at risk", font=2, cex=0.8)
    # times to show numbers at 
    delta <- xlim[2]/nRiskCat
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
	text(y=offsetNRisk-0.15-(i*yDelta), x=c(-delta*0.33, pos), txt, cex=0.8)	
    }

    par(xpd=F)
}
