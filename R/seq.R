#' @title Most frequent number per col of integer matrix
#' 
#' @description get large matrix metrics
#'
#' @param dm matrix
#' 
#' @useDynLib dataAnalysisMisc
#'
#' @export
getMF <- function(dm) {
    mode(dm) <- "integer"    
    n <- integer(length(dm[1,]))       
    cons <- integer(length(dm[1,]))       
    len1 <- as.integer(length(dm[,1]))    
    len2 <- as.integer(length(dm[1,]))    
    katMax <- as.integer(max(dm)+1)    

    ret <- .C("cons", matrix=dm, len1=len1, len2=len2, sel=cons, n=n, katMax=katMax)    
    return(list(ret$sel, ret$n))
}


#' @title Match and subst seq
#' 
#' @description ... 
#'
#' @param seq seq
#' @param search pattern
#' 
#' @useDynLib dataAnalysisMisc
#'
#' @export
matchSub <- function(seq, search) {
    seqN <- integer(nchar(seq))    
    len <- as.integer(nchar(seq))    
    posM <- as.integer(1)

    ret <- .C("mtchSub", len=len, pattern=search, string=seq, seqN=seqN, posM=posM)    
    #ret    
    found <- any(ret$seqN != 0)

    if (found) {
	return(list(found=found, num=ret$seqN[-c(1:ret$posM)]))
    } else {
	return(list(found=found, num=ret$seqN))
    }
}


#' @title Match and subst seq rev
#' 
#' @description ... 
#'
#' @param seq seq
#' @param search pattern
#' 
#' @useDynLib dataAnalysisMisc
#'
#' @export
matchSubRev <- function(seq, search) {
    seqN <- integer(nchar(seq))    
    len <- as.integer(nchar(search))    
    posM <- as.integer(1)

    ret <- .C("mtchSubRev", len=len, pattern=search, string=seq, seqN=seqN, posM=posM)    
    #ret    
    found <- any(ret$seqN != 0)

    if (found) {
	return(list(found=found, num=rev(ret$seqN[-c((ret$posM+nchar(search)):length(ret$seqN))])))
    } else {
	warning("Not implmeneted yet!")
	#return(list(found=found, num=ret$seqN))
    }
}


#' @title Test occurence of substring
#' 
#' @description ... 
#'
#' @param seq seq
#' @param search pattern
#' 
#' @useDynLib dataAnalysisMisc
#'
#' @export
matchPattern <- function(search, seq) {
    res <- integer(1)

    ret <- .C("mtch", pattern=search, string=seq, res=res)    
    return (ret$res)
}
