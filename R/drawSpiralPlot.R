#' @title Draw spiral plot
#' 
#' @description
#' 
#' @import pals
#' @import plotrix
#' 
#' @export
drawSpiralPlot <- function(data, meta, coord, off=1, col=NULL, segm=NULL, leg=F, x.leg=NULL, y.leg=NULL, cex.leg=1, leg.cut=1) {
    ## area
    tmp <- data.frame(do.call(rbind, meta))
    tmp$x <- as.numeric(as.character(tmp$x))
    tmp$y <- as.numeric(as.character(tmp$y))
    #print(tmp)
    
    ## plot + legende
    if (leg) {
        if (is.null(x.leg)) {
            x.leg <- max(coord$x)
        }
        if (is.null(y.leg)) {
            y.leg <- max(coord$y)
        }
        par(mar=c(4,4,4,10))
    }
    plot(coord$y~coord$x, data=tmp,  xlab="", ylab="",  bty='n', col="white", xaxt='n', yaxt='n')
    if (leg) {
        m <- apply(data, 1, max)
        legend(x.leg, y.leg, cex=cex.leg, rownames(data)[which(m > leg.cut)], fill=col[which(m > leg.cut)], bty='n')
    }
    
    
    # points(tmp$y~tmp$x)
    if (!is.null(segm)) {
        for (i in 1:length(segm)) {
            segments(segm[[i]]$xFrom, segm[[i]]$yFrom, segm[[i]]$xTo, segm[[i]]$yTo, lwd=3)
        }
    }
   
    for (i in 1:length(data[1,])) {
        drawCircle(tmp$label[i], centr=c(tmp$x[i], tmp$y[i]),
                   data=data[,i,drop=F], off, scale=tmp$scale[i],
                   phiOff=tmp$phiOff[i], col=col)
    }
}




drawCircle <- function(label, centr=c(0,0), data, off=1, col=NULL, scale=1,phiOff=0) {
    
    # col supplied
    if (is.null(col)) {
        col <- kelly()
    }
    #print(col)
    
    ##central circle with text
    df <- list()
    for (phi in seq(0, 2*pi, 0.01)) { 
        df[[length(df)+1]] <- data.frame(x=centr[1]+off*cos(phi),
                                         y=centr[2]+off*sin(phi),
                                         phi=phi)
    }
    df <- do.call(rbind, df)
    #plot(y~x, data=df, pch=19, ylim=c(-3, 3), xlim=c(-3,3), t="l", lwd=1, 
    #     xaxt='n', yaxt='n', xlab="", ylab="", bty='n')
    lines(y~x, data=df)
    #draw.circle(df$x, df$y, off)
    #print(df)
    text(centr[1], centr[2], label)
    
    
    ##Circles around
    nCirc <- length(data[,1])
    col <- col[order(data[,1])]
    data <- data[order(data[,1]),,drop=F]
    #Area corresponds to fractions
    data$radius <- sqrt(data[,1]/pi)
    
    data$x <- NA
    data$y <- NA
    data$phi <- NA
    
    ####################
    ##calculate phi
    for (i in 1:length(data[,1])) {
        data$phi[i] <- asin(data$radius[i]/(off+data$radius[i]))
    }
    
    ### scale
    scl <- pi/sum(data$phi)*scale
    data$phi <- data$phi*scl
    data$radius <- data$radius*scl
    
    
    for (i in 1:length(data[,1])) {
        sumPhi <- 0
        if (i > 1) {
            sumPhi <- sum(data$phi[1:(i-1)]*2)
        }
        data$x[i] <- centr[1]+(off+data$radius[i])*cos(data$phi[i]+sumPhi+phiOff)
        data$y[i] <- centr[2]+(off+data$radius[i])*sin(data$phi[i]+sumPhi+phiOff)
        
        #circ(data$x[i], data$y[i], data$radius[i])
        draw.circle(data$x[i], data$y[i], data$radius[i], col=col[i])
    }
    
}
