#' @title Plot some position preserving boxplots
#' 
boxGraph <- function(y, x, xSel=NA, xlim=NULL, ylim=NULL) {
    ##range
    if (is.null(xlim)) {
        xlim <- c(min(x, na.rm=T), max(x, na.rm=T))
    }
    if (is.null(ylim)) {
        ylim <- c(min(y, na.rm=T), max(y, na.rm=T))
    }
    
    par(mar=c(8,5,1,1))
    plot(ylim~xlim, xlim=xlim, ylim=ylim, col="white", las=2, xlab="")
    
    for (xPos in levels(factor(x))) {
        print(xPos)
        yDat <- y[which(x == xPos)]
        width=0.01*(xlim[2]-xlim[1])
        
        cl <- ifelse(xPos==xSel, rgb(0.86, 0.08, 0.24, 0.4), rgb(0.5, 0.5, 0.5, 0.25))
        box(xPos, yDat, width, col=cl)
    }
}


box <- function(xPos, yDat, width, lwd=2, col="white") {
    xPos <- as.numeric(as.character(xPos))
    q <- quantile(as.numeric(as.character(yDat)))
    #1.5x IQR
    upper <- 1.5*(q[4]-q[2])
    lower <- 1.5*(q[2]-q[4])
    
    #colored?
    rect(xPos-width, q[2], xPos+width, q[4], col=col)
    
    ##horizontal
    segments(xPos-width, q[3], xPos+width, lwd=lwd+2)
    
    #upper / lower
    if (any(xPos > upper)) {
        segments(xPos-width, upper, xPos+width, lwd=lwd-1, lty=2)
        segments(xPos, q[4], xPos, upper, lwd=lwd-1, lty=2)
    } else {
        segments(xPos-width, max(yDat), xPos+width, lwd=lwd-1, lty=2)
        segments(xPos, q[4], xPos, max(yDat), lwd=lwd-1, lty=2)
    }
    
    if (any(xPos < lower)) {
        segments(xPos-width, lower, xPos+width, lwd=lwd-1, lty=2)
        segments(xPos, q[2], xPos, lower, lwd=lwd-1, lty=2)
    } else {
        segments(xPos-width, min(yDat), xPos+width, lwd=lwd-1, lty=2)
        segments(xPos, q[2], xPos, min(yDat), lwd=lwd-1, lty=2)
    }
}
