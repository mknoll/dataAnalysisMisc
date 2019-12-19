#' @title Create patient characteristics table 
#'
#' @description Creates a patient characteristics table for a provided
#' data.frame, as used by the plotForest and plotForestMV functions. 
#' Multiple observations per patient are indicated by the subject parameter,
#' defaults to NULL. In this case, each row is considered an independent 
#' subject. If a subject vector is provided and if subjRef is set to "subject", 
#' the reference Value (100%) is set to the number of subject, otherwise the 
#' number of rows is selected. If for any variable, no duplicate analysis
#' makes sense (e.g. sex), the selected indices for this variable can be 
#' provided as names list via the subjVar parameter. 
#'
#' @param latex create .tex files and/or compile it.
#' @param data data.frame containing the variables 
#' @param subject vecotr containing subject ids 
#' @param subjVar list conatining indices of which rows of data to used for 
#' any given variable (provided as name of the list element)
#' @param na.rm remove na?
#' @param outdir directory supplied to preparePdf(). Here, the tex and pdf files 
#' will be stored.
#' @param filename filename of the pdf/tex file
#' 
#' @return if latex is set to False, a data.frame with the patient characterisitcs.
#' Otherwise, a list containing the patient characteristics table and information
#' returned by the preparePdf() function (.tex filename, .pdf filename) 
#'
#' @export 
createPatChar <- function(data, subject=NULL, subjVar=NULL, subjRef=NULL, na.rm=T, latex=T, outdir="../reports/", filename=NULL, fract=F) {
    pat <- list()
    ## All
    if (is.null(subject)) {
        pat[[length(pat)+1]] <- data.frame(name1="All",
                                           name2=NA,
                                           n=length(data[,1]),
                                           per=100)
        nTotal <- length(data[,1])
    } else {
        pat[[length(pat)+1]] <- data.frame(name1="All",
                                           name2=NA,
                                           n=NA,
                                           per=NA)
        nTotal <- ifelse(is.null(subjRef) || subjRef == "subject", length(unique(subject)), length(data[,1]))
        print(nTotal)
        pat[[length(pat)+1]] <- data.frame(name1=NA,
                                           name2="Samples",
                                           n=length(data[,1]),
                                           per=round(length(data[,1])/nTotal*100, 2))
        pat[[length(pat)+1]] <- data.frame(name1=NA,
                                           name2="Subjects",
                                           n=length(unique(subject)),
                                           per=round(length(unique(subject))/nTotal*100, 2))
    }

    for (i in 1:length(data[1,])) {
	print(colnames(data)[i])
        if (class(data[,i]) %in% c("character", "factor")) {
            pat[[length(pat)+1]] <- data.frame(name1=colnames(data)[i],
                                               name2=NA,
                                               n=NA,
                                               per=NA)
            for (lv in levels(factor(data[,i]))) {
                if (colnames(data)[i] %in% names(subjVar)) {
                    sel <- subjVar[which(names(subjVar) == colnames(data)[i])][[1]]
                    pat[[length(pat)+1]] <- data.frame(name1=NA,
                                                       name2=lv,
                                                       n=length(which(data[sel,i] == lv)),
                                                       per=round(length(which(data[sel,i] == lv))/nTotal*100, 2))
                } else {
                    pat[[length(pat)+1]] <- data.frame(name1=NA,
                                                       name2=lv,
                                                       n=length(which(data[,i] == lv)),
                                                       per=round(length(which(data[,i] == lv))/nTotal*100, 2))
                }
            }
        } else if (class(data[,i]) == "numeric") {
            pat[[length(pat)+1]] <- data.frame(name1=colnames(data)[i],
                                               name2=NA,
                                               n=NA,
                                               per=NA)
            if (colnames(data)[i] %in% names(subjVar)) {
                sel <- subjVar[which(names(subjVar) == colnames(data)[i])][[1]]
                pat[[length(pat)+1]] <- data.frame(name1=NA,
                                                   name2=NA,
                                                   n=paste(round(median(data[sel,i], na.rm=na.rm), 2)),
                                                   per=round(length(which(!is.na(data[sel,i])))/nTotal*100, 2))
                pat[[length(pat)+1]] <- data.frame(name1=NA,
                                                   name2=NA,
                                                   n=paste("[",round(min(data[sel,i], na.rm=na.rm),2),
                                                           ";",round(max(data[sel,i], na.rm=na.rm),2), "]",
                                                           sep=""),
                                                   per=NA)
            } else {
                pat[[length(pat)+1]] <- data.frame(name1=NA,
                                                   name2=NA,
                                                   n=paste(round(median(data[,i], na.rm=na.rm), 2)),
                                                   per=round(length(which(!is.na(data[,i])))/nTotal*100, 2))
                pat[[length(pat)+1]] <- data.frame(name1=NA,
                                                   name2=NA,
                                                   n=paste("[",round(min(data[,i], na.rm=na.rm),2),
                                                           ";",round(max(data[,i], na.rm=na.rm),2), "]",
                                                           sep=""),
                                                   per=NA)
            }
        }
    }
    pat <- do.call(rbind, pat)
    colnames(pat) <- c("Feature", "", "n", "%")

    if (!latex) {
        return(pat)
    } else {		
	return (list(tbl=pat, pdffile=preparePdf(pat, outdir, col="llcc", filename=filename))) 
    }
}

