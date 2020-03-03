#' @title Prepare data for CNV overlap comparison
#' 
#' @param pos list with positions of CNV segments
#' @param id names of cohorts
#' @param sel list with selected indices 
#' @param diff list with differences
#' 
#' @export
cnvOverlap <- function(pos, nm, sel, diff) {
    # TODO: sanitize input data
    #check if sel in pos!

    # all positions and order
    allP <- unique(unlist(pos))
    bak <- allP
    allP <- data.frame(do.call(rbind, strsplit(allP, ":")))
    rownames(allP) <- bak
    colnames(allP) <- c("chr", "pos")
    allP$chr <- as.numeric(gsub("chr", "", allP[,1]))
    allP$pos <- as.numeric(as.character(allP$pos))
    allP <-allP[order(allP$chr, allP$pos),]

    # split by chr
    f <- list()
    for (cc in unique(allP$chr)) {
	w0 <-which(allP$chr == cc)
	#FIXME
	if (length(w0) < 1) { next }
	print(paste("CHROM", cc))

	coll <- list()
	collD <- list()
	for (j in 1:length(sel)) {
	    vec <- pos[[j]][match(rownames(allP)[w0], pos[[j]])]
	    vecD <- diff[[j]][match(rownames(allP)[w0], names(diff[[j]]))]
	    w <-which(!is.na(vec))
	    vec[w] <- ifelse(vec[w] %in% sel[[j]], 1, 0)
	    w <- which(!is.na(vec))
	    if (length(w) > 0 && w[1] != 1 && length(vec) > 1) {
		vec[1:(w[1]-1)] <- 0
		vecD[1:(w[1]-1)] <- 0
	    }

	    if (length(w0) > 1) {
		for (z in 2:length(vec)) {
		    if (is.na(vec[z])) {
			vec[z] <- vec[z-1]
			vecD[z] <- vecD[z-1]

		    }
		}
	    }
	    print(vec)
	    coll[[length(coll)+1]] <- vec
	    collD[[length(collD)+1]] <- vecD

	}
	if (length(w0) == 1) {
	    s <- data.frame(t(unlist(coll)))
	    sD <- data.frame(t(unlist(collD)))
	} else {
	    s <- apply(do.call(cbind, coll), 2, as.numeric)
	    sD <- do.call(cbind, collD)
	}
	rownames(s) <- rownames(allP)[w0]
	rownames(sD) <- rownames(allP)[w0]
	colnames(s) <-colnames(sD) <- nm

	f[[length(f)+1]] <- list(sel=s, eff=sD)
    }

    f <- lapply(f, function(x) if (any(is.na(x$sel))) { NULL } else { x })

    return(f)

}

#' @title Evaluate overlaps between segments
#' @export
cnvEval <- function(f, cut=NULL) {
    f <- Filter(length, f)
    f <- lapply(f,  function(x) {
		    x$sel <- t(apply(x$sel, 1, function(u) as.numeric(as.character(u))))
		    colnames(x$sel) <-colnames(x$eff)
		    s <-   apply(data.frame(x$sel), 1, sum)
		    print(s)
		    cut0 <- ifelse(is.null(cut), length(x$sel[1,]), cut)
		    if (any(s >= cut0)) {
			#### FIXME: 
			dir <- apply(data.frame(x$eff), 1, function(z) (all(z >0 ) | all(z < 0)))
			x$sel <-data.frame(x$sel, DIR=ifelse(dir == T & s>= cut0,  1,0), SUM=s)
			x
		    } else {
			####FIXME
			NULL
		    }
})
    Filter(length,f)
}


#' @title Augments segment breaks (plottign of steps)
#' @export
cnvAug <- function(f) {
    ### FIXME: chromosomal boundaries!
    anno <- data.frame(minfi::getAnnotation(minfiData::RGsetEx))
    anSub <-anno[,c("chr", "pos")]
    f0 <- split(anSub, f=anSub$chr)
    maxPos <- lapply(f0, function(x) x$pos[which.max(x$pos)])

    f2 <- lapply(f, function(y) {
		     chr0<- do.call(rbind, strsplit(as.character(rownames(y$sel)), ":"))[,1]
		     eff1 <- list()
		     sel1 <- list()
		     chr1 <- list()
		     pos1 <- list()
		     grp1<- list()

		     for (cc in unique(chr0)) {
			 x <- list(sel=y$sel[which(chr0 == cc),],
				   eff=y$eff[which(chr0 == cc),])

			 #### add end of chromosome
			 posM <- max(as.numeric(do.call(rbind, strsplit(as.character(rownames(x$sel)), ":"))[,2]))
			 chrC<- do.call(rbind, strsplit(as.character(rownames(x$sel)), ":"))[,1][1]
			 mP0 <- maxPos[which(names(maxPos) ==  chrC)][[1]]
			 print(paste("chrC:", chrC, " maxPos", mP0, " avail", posM))

			 if (length(mP0)  > 0 && posM != mP0) {
			     mP <- maxPos[which(names(maxPos) ==  chrC)][[1]]

			     print("ADD")
			     #add row
			     tmp <- x$sel[length(x$sel[,1]),,drop=F]
			     tmp2 <- x$eff[length(x$eff[,1]),,drop=F]
			     rownames(tmp) <- paste(chrC, mP, sep=":")
			     rownames(tmp2) <- paste(chrC, mP, sep=":")
			     x$sel <- rbind(data.frame(x$sel), tmp)
			     x$eff <- rbind(data.frame(x$eff), tmp2)
			 }


			 for (i in 1:length(x$eff[1,])) {
			     pos <- as.numeric(do.call(rbind, strsplit(as.character(rownames(x$sel)), ":"))[,2])
			     chr<- do.call(rbind, strsplit(as.character(rownames(x$sel)), ":"))[,1]
			     x$sel <- x$sel[order(pos),]
			     x$eff <- x$eff[order(pos),]


			     vec <- x$sel[,i]
    			     vecD <- x$eff[,i]

			     print("STARTED")
			     print(x$sel)
			     print(paste("VEC", vec))
			     print(paste("POS", pos))
		    
			     f3 <- list()
			     f3D <- list()
			     fP <- list()
			     fC <- list()
			     for (j in 1:(length(vec)-1)) {
				 f3[[length(f3)+1]] <- vec[j]
				 f3D[[length(f3D)+1]] <- vecD[j]
				 fP[[length(fP)+1]] <- pos[j]
				 fC[[length(fC)+1]] <- chr[j]
				 if ((vec[j] == 0 && vec[j+1] == 1)) {
				     f3[[length(f3)+1]] <- vec[j]
				     f3D[[length(f3D)+1]] <- vecD[j]
				     fC[[length(fC)+1]] <- chr[j]
				     fP[[length(fP)+1]] <- pos[j+1]-1
				 }
				 ####
				 if (vec[j] == 1 && vec[j+1] == 1 && vecD[j]*vecD[j+1] < 0) {
				     f3[[length(f3)+1]] <- vec[j]
				     f3D[[length(f3D)+1]] <- vecD[j]
				     fC[[length(fC)+1]] <- chr[j]
				     fP[[length(fP)+1]] <- pos[j+1]-1
				 }

			     }
			     f3[[length(f3)+1]] <- vec[length(vec)]
			     f3D[[length(f3D)+1]] <- vecD[length(vecD)]
			     fC[[length(fC)+1]] <- chr[length(chr)]
			     fP[[length(fP)+1]] <- pos[length(pos)]



			     ####
			     vec <- unlist(f3)
			     vecD <- unlist(f3D)
			     chr <- unlist(fC)
			     pos <- unlist(fP)
			    print(data.frame(vec, pos))

			     f3 <- list()
			     f3D <- list()
			     fP <- list()
			     fC <- list()
			     for (j in 1:(length(vec)-1)) {
				 f3[[length(f3)+1]] <- vec[j]
				 f3D[[length(f3D)+1]] <- vecD[j]
				 fP[[length(fP)+1]] <- pos[j]
				 fC[[length(fC)+1]] <- chr[j]
				 if (vec[j] == 1 && vec[j+1] == 0) {
				     print("ADD")
				     print(paste(vec[j], " - " , pos[j]))
				     print(paste("+1 ->", vec[j+1], " - " , pos[j+1]))
				     f3[[length(f3)+1]] <- vec[j]
				     f3D[[length(f3D)+1]] <- vecD[j]
				     fC[[length(fC)+1]] <- chr[j]
				     fP[[length(fP)+1]] <- pos[j+1]-1
				 }
			     }
			     f3[[length(f3)+1]] <- vec[length(vec)]
			     f3D[[length(f3D)+1]] <- vecD[length(vecD)]
			     fC[[length(fC)+1]] <- chr[length(chr)]
			     fP[[length(fP)+1]] <- pos[length(pos)]

			    print(data.frame(unlist(f3), unlist(fP)))

			     ##
			     eff1[[length(eff1)+1]]<- unlist(f3D)
			     sel1[[length(sel1)+1]]<- unlist(f3)
			     chr1[[length(chr1)+1]] <- unlist(fC)
			     pos1[[length(pos1)+1]] <- unlist(fP)
			     grp1[[length(grp1)+1]] <- rep(colnames(x$sel)[i], length(unlist(f3D)))

			 }
		     }

		     ####
		     ###
		     ###combine to data.frame
		     df <- data.frame(eff=unlist(eff1), sel=unlist(sel1), chr=unlist(chr1), pos=unlist(pos1), grp=unlist(grp1))
})
    f2
}

#' @title Merge segments
#' @export
segMerge <- function(tmp) {
    go <- T
    while(go) {
	id <- apply(tmp, 1, function(x) paste(x, collapse=";"))
	if (length(tmp[,1]) == 1) { 
	    go <- F
	} else {
	    for (i in 2:(length(tmp[,1])-1)) {
		if (id[i-1] == id[i] && id[i] == id[i+1]) {
		    tmp <- tmp[-i,]
		    break
		}
		if (i == length(tmp[,1])-1) {
		    go <- F
		}
	    }
	    if (length(tmp[,1]) < 3) {
		go <- F
	    }
	}
    }
    return(tmp)
}

#' @title Get genes in overlapping segments 
#' 
#' @description Currently only for human
#' 
#' @export
cnvSegGene <- function(f, cut=NULL, type="hs",mergeSeg=T) {
    if (type == "hs") {
	sel <- as.data.frame(org.Hs.egSYMBOL)
    } else {
	stop("Only hs implemented yet!")
    }

    ###FIXME: borders!
    f2 <- lapply(f, function(x) {
		     ### merge intervals
		     if (mergeSeg) {
			 x$sel <- segMerge(x$sel)
			 x$eff <- x$eff[match(rownames(x$sel), rownames(x$eff)),]
		     }

		     cut0 <- ifelse(is.null(cut), length(x$eff[1,]), cut)
		     print(cut0)
		     w <- which(x$sel$DIR == 1 & x$sel$SUM >= cut0)
		     ret <- list()
		     if (length(w) > 0) {
			 myc <-data.frame(do.call(rbind, strsplit(rownames(x$sel), ":")))
			 colnames(myc) <- c("chr", "pos")
			 myc$pos <- as.numeric(as.character(myc$pos))

			 w2 <-w+1
			 wPos <- myc$pos[w]
			 if (max(w2) > length(x$sel[,1])) {
			     w2[which.max(w2)] <- length(x$sel[,1])
			     wPos2 <- myc$pos[w2]-1
			     wPos2[length(wPos2)] <- wPos2[length(wPos2)]+1
			 } else {
			     wPos2 <- myc$pos[w2]-1
			 }

			 dfT <- data.frame(chrom=myc$chr[w], start=wPos, end=wPos2)

			 if (length(dfT[,1]) >= 1) {
			     for (z in 1:length(dfT[,1])) {
				 print(dfT[z,])
				 myc.gr <- makeGRangesFromDataFrame(dfT[z,])
				 cands <- data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), myc.gr))
				 ret[[length(ret)+1]] <- sel[which(sel$gene_id %in% cands$gene_id),2]
				 names(ret)[length(ret)] <- paste(dfT[z,1], dfT[z,2], dfT[z,3], sep=":")
			     }
			 }

		     }
		     ret

})
}

