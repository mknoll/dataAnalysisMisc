#' @title Draw spiral plot
#' 
#' @description
#' 
#' @import pals
#' @import plotrix
#' 
#' @export
drawSpiralPlot <- function(data, meta, coord, off=1, col=NULL, segm=NULL, leg=F) {
    ## area
    tmp <- data.frame(do.call(rbind, meta))
    tmp$x <- as.numeric(as.character(tmp$x))
    tmp$y <- as.numeric(as.character(tmp$y))
    #print(tmp)
    
    if (leg) {
        par(mar=c(4,4,4,10))
        par(xpd=T)
        legend(max(coord$x)*1.1, max(coord$y), rownames(data), col=col)
        par(xpd=F)
    }
    
    plot(coord$y~coord$x, data=tmp,  xlab="", ylab="",  bty='n', col="white", xaxt='n', yaxt='n')
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
    
    ##solve to get new r
    #data$phi * scl = c* sin (r / (r+ off))
    a0 <- 0.01
    for (i in 1:length(data[,1])) {
        r <<- data$radius[i]
        k <<- scl
        off <<- off
        phi <<- data$phi[i]
        o <- optim(a0, fnToFind, method="BFGS", hessian=T)
        print(paste("o: ", o$par, " k: ", k, " rad: ", data$radius[i], " radNew: ", abs(o$par)*data$radius[i]))
        data$radius[i] <- abs(o$par)*data$radius[i]
    }
        
    data$phi <- data$phi*scl
    #data$radius <- data$radius*scl
    
    
    
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


fnToFind <- function(a) {
    return(sin(a*r/(a*r+off))-k*phi)
}