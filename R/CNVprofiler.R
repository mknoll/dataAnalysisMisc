#' @title Find Clusters in CNV data
#'
#' @param data data.frame with samples in columns and segments 
#' in rows. Allowed values: -1, 0, 1. rownames: "1;100;400" 
#' corresponding to chr 1, startPos: 100, endPos: 400
#'
#' @import factoextra
#'
#' @export
findCNVClusters <- function(data, splitChr=";", nSpl=3) {
    rownames(data) <- gsub("chr", "", rownames(data))

    ### annotation
    anCNV <- data.frame(do.call(rbind, strsplit(rownames(data), splitChr)))
    
    colnames(anCNV) <- c("chr", "start", "end")[1:nSpl]
    print(anCNV[1,,drop=F])
    anCNV <- data.frame(apply(anCNV, 2, function(x) as.numeric(as.character(x))), check.names=F)

    ### sort
    print("Sort annotation...")
    rownames(anCNV) <- rownames(data)
    if (nSpl == 3) {
	anCNV <- anCNV[order(anCNV[,1], anCNV[,2], anCNV[,3]),]
    } else if (nSpl == 2) {
	anCNV <- anCNV[order(anCNV[,1], anCNV[,2]),]
    } else {
	stop("Invalid nSpl! Can be 2 (chr:start) or 3 (chr:start:stop)")
    }

    data <- data[match(rownames(anCNV), rownames(data)),]

    ### split data and check order
    f <- split(data.frame(data, check.names=F), f=anCNV$chr)

    f <- lapply(f, function(x) {
		    print(rownames(x)[1])
		    print(paste("DIM: ", paste(dim(x), collapse="x")))
		    res <- NA
		    tryCatch({
			res <- procChr(x)
		    }, error=function(e) {warning(e)})

		    list(cluster=res, data=x)
    })
    
    return(f)
}


#' @title Plot different groups of samples based 
#' on CNV alterations
#' 
#' @import EBImage
#' @export
plotCNVChr <- function(dat) {
    cl <- dat$cluster
    data <- data.matrix(dat$data)
    data[which(data == 0)] <- 0.5
    data[which(data == -1)] <- 0

    f <- split(data.frame(t(data)), f=cl)
    par(mfrow=c(1, length(unique(cl))), mar=c(0.5,0,0,0.5), bty="n")
    for (i in 1:length(f)) {
	image(data.matrix(f[[i]]), xaxt="n", yaxt="n")
    }
}

procChr <- function(dataChr) {
    ### calculate dissimilarity matrix
    mU <- dissimM(dataChr)
    mU[which(is.na(mU))] <- 0
    mL <- t(mU)
    mT <- mU+mL

    ### find clusters
    print("Select # of clusters ... ")
    nc <- fviz_nbclust(t(mT), kmeans, method="silhouette")
    nCl <- which.max(nc$data$y)
    km <- kmeans(t(mT), nCl)

    return(km$cluster)
}


#' @title Compute dissim matrix
#' 
#' @import foreach
#' @import doParallel
#' @import parallel
#' 
#' @export
dissimM <- function(dataChr) {
    ##use parallel computing    
    no_cores <- parallel::detectCores() - 1    
    no_cores <- ifelse(no_cores == 0, 1, no_cores)    
    ## FIXME    
    if (is.na(no_cores)) { no_cores <- 4 }    
    doParallel::registerDoParallel(no_cores)    

    print("Calculate dissimilarity...")
    m <- foreach(i=1:(length(dataChr[1,])-1))%dopar% {
	if (i %% 4 == 0) { cat(paste("\r  ", round(i/length(dataChr[1,])*100, 2), "%"    )) }
	vec <- rep(NA, i)
	for (j in (i+1):length(dataChr[1,])) {
	    n <- length(which(dataChr[,i] != dataChr[,j]))
	    vec <- c(vec, n)
	}
	vec
    }
    doParallel::stopImplicitCluster()


    m[[length(m)+1]] <- rep(NA, length(dataChr[1,]))
    m <- do.call(rbind, m)
    return(m)
}
