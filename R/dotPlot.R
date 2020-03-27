#' @title Dot plot for cor + p-val
#' @import ggplot2
#' @export
dotPlot <- function(m, grp=NULL, na.rm=T, pCut=0.05) {
    ## TODO: sanitize inpu
    ## TODO: p-calc method, cor method

    coll <- list()
    if (is.null(grp)) {
	for (i in 1:length(m[1,])) {
	    for (j in 1:length(m[1,])) {
		#if (j >= i) { next }
		tmp <- m[,c(i,j)]
		if (na.rm) {
		    tmp <- tmp[complete.cases(tmp*0),,drop=F]
		}
		cr <- cor(tmp[,1],tmp[,2])
		p <- cor.test(tmp[,1],tmp[,2])[[3]]
		if (j >= i) { 
		    coll[[length(coll)+1]] <- data.frame(i=i, j=j, cr=NA, p=p, nI=colnames(m)[i], nJ=colnames(m)[j])
		} else {
		    coll[[length(coll)+1]] <- data.frame(i=i, j=j, cr=cr, p=p, nI=colnames(m)[i], nJ=colnames(m)[j])
		}
	    }
	}
    } else {
	if (length(unique(grp)) != 2) { 
	    stop("Exactly two group levels required!")
	}
	### TODO remove redundancy!
	w <- which(grp == unique(grp)[1])
	for (i in 1:length(m[1,])) {
	    for (j in 1:length(m[1,])) {
		#if (j >= i) { next }
		tmp <- m[w,c(i,j)]
		if (na.rm) {
		    tmp <- tmp[complete.cases(tmp*0),,drop=F]
		}
		cr <- cor(tmp[,1],tmp[,2])
		p <- cor.test(tmp[,1],tmp[,2])[[3]]
		if (j > i) { 
		    coll[[length(coll)+1]] <- data.frame(i=i, j=j, cr=cr, p=p, nI=colnames(m)[i], nJ=colnames(m)[j])
		} else if (j == i) { 
		    coll[[length(coll)+1]] <- data.frame(i=i, j=j, cr=NA, p=p, nI=colnames(m)[i], nJ=colnames(m)[j])
		}
	    }
	}
	w <- which(grp == unique(grp)[2])
	for (i in 1:length(m[1,])) {
	    for (j in 1:length(m[1,])) {
		#if (j >= i) { next }
		tmp <- m[w,c(i,j)]
		if (na.rm) {
		    tmp <- tmp[complete.cases(tmp*0),,drop=F]
		}
		cr <- cor(tmp[,1],tmp[,2])
		p <- cor.test(tmp[,1],tmp[,2])[[3]]
		if (i < j) { 
		    coll[[length(coll)+1]] <- data.frame(i=j, j=i, cr=cr, p=p, nI=colnames(m)[j], nJ=colnames(m)[i])
		} 
	    }
	}
    }
    coll <- data.frame(do.call(rbind,coll))

    ### plot
    coll$i <- factor(coll$i, levels=c(1:length(m[,1])))
    coll$j <- factor(coll$j, levels=c(1:length(m[,1])))
    coll$p[which(!is.na(coll$p))] <- ifelse(coll$p[which(!is.na(coll$p))] < pCut, paste0("<",pCut), paste0("≥",pCut))
    coll$p <- as.character(coll$p)
    g <- ggplot(coll, aes(x=i, y=j, size=abs(cr), fill=cr, color=p))+
	geom_point(shape = 21, stroke=2) + #,  fill = "b", size = 5, stroke = 5) + 
	ylim(1, length(m[1,])) + 
	xlim(1, length(m[1,])) +
	labs(color="p", size="correlation", fill="cor", alpha="p-val") +
	scale_alpha(range = c(1, 0.5)) +
	scale_size_continuous(range = c(3,7)) + 
	scale_fill_gradient2(
			      low = "royalblue",
			      high = "tomato", #"#00FF00",
			      space = "Lab",
			      na.value = "grey50",
			      guide = "colourbar",
			      aesthetics = "fill") +
	guides(size=F) + 
	scale_x_discrete(breaks=c(1:length(m[1,])),
			 labels=c(colnames(m))) + 
	scale_y_discrete(breaks=c(1:length(m[1,])),
			 labels=c(colnames(m))) +
	theme(legend.position="top") + 
	scale_color_brewer(palette="Set1")
    if (is.null(grp)) {
	g <- g + labs(x="", y="")
    } else {
	g <- g + labs(x=unique(grp)[2], y=unique(grp)[1])
	g <- g + geom_segment(x=1, y=1, xend=length(m[1,]), yend=length(m[1,]), size=1, col="black")
    }

    print(g)

    return(list(g, coll))
}
