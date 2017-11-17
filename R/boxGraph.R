#' @title Plot some position preserving boxplots
#' @export
boxGraph <- function(y, x, fact=NULL, xSel=NA, xlim=NULL, ylim=NULL, ylab=NA, main="",
                     lines=F, cols=NULL, colsSel=NULL) {
    ##remove nas 
    rm <- which(is.na(y))
    if(length(rm) > 0) {
        warning("NAs removed!")
        y <- y[-rm]
        x <- x[-rm]
    }
    
    ##range
    if (is.null(xlim)) {
        xlim <- c(min(x, na.rm=T), max(x, na.rm=T))
    }
    if (is.null(ylim)) {
        ylim <- c(min(y, na.rm=T), max(y, na.rm=T))
    }
    
    par(mar=c(8,5,3,1))
    plot(ylim~xlim, xlim=xlim, ylim=ylim, col="white", las=2, xlab="", ylab=ylab, main=main)
    
    if (is.null(fact)) {
        fact <- rep(1, length(x))
    }
    
    xBak <- x
    yBak <- y
    u <- 0
    for (lv in levels(factor(fact))) {
        u<-u+1
        x <- xBak[which(fact == lv)]
        y <- yBak[which(fact == lv)]
        
        q4 <- c()
        q3 <- c()
        q2 <- c()
        
        clL <- rgb(0.5, 0.5, 0.5, 0.25)
        clS <- rgb(220/255,20/255,60/255,0.8)
        
        for (xPos in levels(factor(x))) {
            yDat <- y[which(x == xPos)]
            width=0.01*(xlim[2]-xlim[1])
            if (!is.null(cols)) {
                clL <- cols[u]
            }
            if (!is.null(colsSel)) {
                clS <- colsSel[u]
            }
            cl <- ifelse(xPos %in% xSel, clS, clL)
            
            box(xPos, yDat, width, col=cl)
            q3 <- c(q3, quantile(yDat, na.rm=T)[3])
            q2 <- c(q2, quantile(yDat, na.rm=T)[2])
            q4 <- c(q4, quantile(yDat, na.rm=T)[4])
        }
        
        #plot quantile connection lines
        if (lines) {
            lines(q3~levels(factor(x)), col=clS)
            lines(q2~levels(factor(x)), col=clL,lty=2)
            lines(q4~levels(factor(x)), col=clL,lty=2)
        }
    }
}


box <- function(xPos, yDat, width, lwd=1, col="white") {
    xPos <- as.numeric(as.character(xPos))
    q <- quantile(as.numeric(as.character(yDat)))
    #1.5x IQR
    upper <- 1.5*(q[4]-q[2])
    lower <- 1.5*(q[2]-q[4])
    
    #colored?
    rect(xPos-width, q[2], xPos+width, q[4], col=col, lwd=lwd)
    
    ##horizontal
    segments(xPos-width, q[3], xPos+width, lwd=lwd)
    
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
