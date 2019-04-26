#' @title Assign class for test data based
#' based on correlation
#'
#' @param ref Reference data (cols: samples)
#' @param refCl Class assignment for reference data
#' @param test Test dataset
#' @param testCl Class testset 
#' 
#' @export
assignCl <- function(ref, refCl, test, testCl=NULL, alternative="two.sided") {
    ## match test data 
    test <- subValid[match(rownames(ref), rownames(subValid)),]

    coll <- list()
    for (corM in c("spearman", "pearson", "kendall")) {
	## calc correlations
	for (i in 1:length(ref[1,])) {
	    for (j in 1:length(test[1,])) {
		cr <- cor(ref[,i], test[,j], method=corM)
		coll[[length(coll)+1]] <- data.frame(REF=i, TEST=j, COR=cr, CORM=corM)
	    }   
	}
    }
    coll <- do.call(rbind, coll)
    coll$CLUSTER <- refCl[coll$REF]

    ## get class assignments
    out <- list()
    for (corM in unique(coll$CORM)) {
	for (i in unique(coll$TEST)) {
	    sub <- coll[which(coll$TEST == i & coll$CORM == corM),]
	    #pT <- t.test(coll$COR~coll$CLUSTER)$p.value
	    #pW <- wilcox.test(coll$COR~coll$CLUSTER)$p.value
	    sub$RANK <- rank(sub$COR)
	    aggRM <- aggregate(sub$RANK, by=list(CLUSTER=sub$CLUSTER), median)
	    out[[length(out)+1]] <- data.frame(SAMPLE=i, 
					       CLASS=aggRM$CLUSTER[which.max(aggRM$x)],
					       TYPE="RANK_MED",
					       CORM=corM)

	    aggCM <- aggregate(sub$COR, by=list(CLUSTER=sub$CLUSTER), median)
	    out[[length(out)+1]] <- data.frame(SAMPLE=i, 
					       CLASS=aggCM$CLUSTER[which.max(aggCM$x)],
					       TYPE="COR_MED",
					       CORM=corM)
	}
    }
    out <- data.frame(do.call(rbind,out))

    mtch <- list()
    if (!is.null(testCl)) {
	## calculate matches
	for (corM in unique(out$CORM)) {
	    for (type in unique(out$TYPE)) {
		sub <- out[which(out$CORM == corM & out$TYPE == type),]
		tbl <- table(sub$CLASS, testCl)
		pF <- NA
		tryCatch({
		    pF <- fisher.test(tbl, alternative=alternative)$p.value
		}, error=function(e) print(e) )
		mtch[[length(mtch)+1]] <- data.frame(pF=pF, CORM=corM, TYPE=type)
	    }
	}
	mtch <- do.call(rbind, mtch)
	mtch <- mtch[order(mtch[,1]),]
    }

    return(list(agg=mtch, class=out))
}


