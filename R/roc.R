#' @title Create ROC curve
#' @import pROC
#' @export
rocDetail <- function(val, grp, fit=NULL,xyleg=c(0.5, 0.1)) {
    rc <- roc(grp, val)
    plot(rc)
    text(0.5, 0.2, paste("AUC=",round(rc$auc, 2),sep=""), font=2)
    ## add youden
    ydn <- rc$sensitivities+rc$specificities-1
    points(ydn~seq(0,1,length.out=length(rc$sensitivities)), col="tomato", pch=19, cex=0.2)
    legend(xyleg[1], xyleg[2], "Youden", fill="tomato")
    abline(v=seq(0,1,length.out=length(rc$sensitivities))[which.max(ydn)])
    tr <- rc$thresholds[which.max(ydn)]

    if (!is.null(fit)) {
	prd <- predict(fit, se.fit=T,type="response")
	df <- data.frame(VAL=val, GRP=grp, PRD=prd$fit,
			 LW=prd$fit-prd$se.fit,
			 UP=prd$fit+prd$se.fit)
	df <- df[order(df$PRD),]
	matplot(df[,c("PRD","LW","UP")], type="l", lty=c(1,2,2), 
		col=c("black", "royalblue", "orange"), ylab="Response")
    

	w <- which(!is.na(df$PRD) & !is.infinite(df$PRD))
	mn <- min(df$PRD[w])+0.02*min(df$PRD[w])
	mx <- max(df$PRD[w])-0.02*max(df$PRD[w])
	rg <- range(mn, mx)
	rg <- rg[2]-rg[1]
	
	jit <- rnorm(length(df[,1]), 0, 0.02*rg)
	df$GRP_Y <- ifelse(df$GRP != 0,  mx+jit, mn+jit)
	#print(df)
	points(df$GRP_Y, col=as.numeric(as.factor(df$GRP)), pch=19, cex=0.2)
	abline(h=c(0.5, tr), lty=1:2)
	legend(0.5*length(val), 0.2, c("0.5", "Youden"), lty=1:2)
    }

    return(tr)
}

