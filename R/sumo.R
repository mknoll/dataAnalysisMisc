#' @title Import from SUMO
#' @export 
impSUMO <- function(data) {

}

#' @title Export to SUMO
#' @export
expSUMO <- function(data) {
}


# An S4 Helper class to allow NULL values
setClassUnion("dfOrNULL", c("data.frame", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("matrixOrNULL", c("matrix", "NULL"))

#' An S4 class representing a SUMO matrix experiment
#'
#' @slot sumoData SUMO exported data
#' @name matrixSUMO-class
#' @exportClass matrixSUMO
matrixSUMO <- setClass("matrixSUMO", 
			 slots=c(sumoMatrixFile="character",
				 sumoData="dfOrNULL",
				 expr="matrixOrNULL", #data.matrix would be better
				 pheno="dfOrNULL",
				 anno="dfOrNULL",
				 selRow="numericOrNULL",
				 selCol="numericOrNULL",
				 log="list")
			 )

#' @title Construcotr matrixSUMO
#' @import data.table
#' @export 
setMethod("initialize", "matrixSUMO",
	  function(.Object,
		   sumoMatrixFile=character) {
	      .Object@sumoMatrixFile <- sumoMatrixFile

	      print("Loading data ...")
	      .Object@sumoData <- data.frame(fread(sumoMatrixFile))

	      print("Splitting data ... ")
	      ## First probeID / row split 
	      wR <- which(!is.na(as.numeric(as.character(exp@sumoData$ProbeID))))[1]
	      toR <- length(exp@sumoData[,1])
	      ## First sampleID / col split 
	      wC <- which(colnames(exp@sumoData) == "Ari.mean")+1
	      toC <- length(exp@sumoData[1,])
	      
	      .Object@expr <- as.matrix(data.matrix(exp@sumoData[wR:toR, wC:toC]))
	      .Object@anno <- data.frame(exp@sumoData[wR:toR, 1:(wC-1)])
	      .Object@pheno <- data.frame(t(data.frame(exp@sumoData[1:(wR-1), wC:toC])))
	      colnames(.Object@pheno) <- exp@sumoData$ProbeID[1:(wR-1)]
	      ## TODO_: Check dimensions
	      
	      .Object@selRow <- 1:length(.Object@expr[,1])
	      .Object@selCol <- 1:length(.Object@expr[1,])
	      print("Finished!")

	      l <- list()
	      l[[length(l)+1]] <- data.frame(Sys.time(), "Created instance")
	      .Object@log <- l

	      .Object
	  })

#' @title Aggregate epxression
#' @export 
agg <- function(obj, fun=mean, var="ILMN_Gene") {
    cl <- which(colnames(obj@anno) == var)
    agg <- aggregate(obj@expr, by=list(ORD=obj@anno[,cl]), fun)
    anno <- obj@anno[match(agg$ORD, obj@anno[,cl]),]

    agg <- agg[,-1]
    obj@expr <- as.matrix(data.matrix(agg))
    obj@anno <- anno

    obj@selRow <- 1:length(agg[,1])

    return(obj)
}

#' @title Keep only most variant
#' @export
filter <- function(obj, fun=mad, cutoff=0.5) {
    mads <- apply(obj@expr, 1, fun)
    sel <- which(mads > quantile(mads, na.rm=T, cutoff)[[1]])

    obj@expr <- obj@expr[sel,]
    obj@anno <- obj@anno[sel,]
    return(obj)
}


#' @title Transform expression data
#' @export 
transf <- function(obj, fun=log) {
    obj@expr <- as.matrix(data.frame(apply(obj@expr, 2, fun)))
    return(obj)
}

#' @title Analyze fixed Effects
#' @export
feAnalyze <- function(obj, sel, frm, frm0=as.formula(VAL~1), padj="BH") {
    res <- fixedEffAnalysis(obj@expr[,sel], 
			    obj@pheno[sel,],
			    frm=frm,
			    frm0=frm0,
			    padj=padj)
    res$SYMBOL <- obj@anno$ILMN_Gene[as.numeric(res$RN)]
    return(res)
}

#' @title Plot single gene
#' @import lattice
#' @export	
plotSingle <- function(obj, gene, sel, frm=as.formula(value~GRP|TREAT*TIME)) {
    mD <- melt(obj@expr[which(obj@anno$ILMN_Gene == gene),sel,drop=F])
    mD <- cbind(mD, obj@pheno[match(mD[,2], rownames(obj@pheno)),])
    xyplot(value~GRP|TREAT*TIME , data=mD, type=c("p"), pch=19 , main=gene,
	   panel=function(x,y, ...) {
	       panel.xyplot(x,y,...)
	       panel.linejoin(x,y,horizontal=F, lwd=3, col="royalblue")
	   })
}

