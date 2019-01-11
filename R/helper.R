#' @title Draw segments for boxplots
#' @import pals
#' @export
drawSegments <- function(data, y, yDelta=-0.05, xyLeg=NULL, col=NULL,
			 txtBottomOff=-1, segmentDelta=0.2, lwd=10,
			 legText=c(0.55, 0.02)) {
    if (is.null(col)) {
	col <- kelly()[-1]
    }

    par(xpd=T)
    colOff <- 0
    leg <- list()
    for (j in 1:length(data[1,])) {
	val <- data[,j]
	i <- 2
	from <- 1
	cl <- as.numeric(factor(val)) + colOff
	off <- segmentDelta
	while (i <= length(val)) {
	    if (val[i-1] != val[i] || i == length(val)) {
		## draw segmetns
		lst <- ifelse(i == length(val), 1, 0)
		segments(from+off, y, i-off+lst, y, lwd=lwd, col=col[cl[i-1]])
		from <- i
		text(txtBottomOff, y, colnames(data)[j], font=2)
		leg[[length(leg)+1]] <- data.frame(LV=val[i-1], col=cl[i-1], GRP=colnames(data)[j])
	    }
	    if (i < length(val)) {
		i <- i+1
	    } else { break }
	}
	y <- y+yDelta
	colOff <- max(cl)
    }
    leg <- data.frame(do.call(rbind, leg))
    leg <- leg[which(!duplicated(paste(leg$LV, leg$col, leg$GRP))),]
    s <- split(leg, f=leg$GRP)
    for (i in 1:length(s)) {
	text(xyLeg[[i]][1]+legText[1], xyLeg[[i]][2]+legText[2], colnames(data)[i], font=2)
	legend(x=xyLeg[[i]][1], y=xyLeg[[i]][2], s[[i]]$LV,fill=col[s[[i]]$col])
    }
    par(xpd=F)
}

