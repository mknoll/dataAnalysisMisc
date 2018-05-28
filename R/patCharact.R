#' @title Create patient characteristics table 
#' @import stargazer
#' @export 
createPatChar <- function(data, subject=NULL, subjVar=NULL, subjRef=NULL, na.rm=T, latex=T, outdir="../reports/") {
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
        sgOut <- stargazer(pat, summary=F, rownames=F)
        for (i in 1:length(sgOut)) {
            if (!(i >= 10 && i != 11 && i < length(sgOut)-3)) { next } 
            tmp <- sgOut[[i]]
            #print(paste("I: ",i, ":", tmp))
            if (substr(tmp,1,1) != " ") {
                #find first & and remove
                feat <- substr(tmp, 1, regexpr("&", tmp, fixed=T)[[1]][1]-1)
                rest <- substr(tmp, regexpr("&", tmp, fixed=T)[[1]][1]+1, nchar(tmp))
                tmp <- paste("\\multicolumn{2}{l}{\\textbf{", feat, "}}", rest, collapse="")
                sgOut[[i]] <- tmp
            }
        }
        sgOut[[7]] <- "\\begin{tabular}{@{\\extracolsep{5pt}} llcc}"

	tex <- "\\documentclass{article} \n \\usepackage[margin=0pt]{geometry} \n \\pagestyle{empty} \\begin{document} "
	tex <- paste(tex, paste(sgOut, collapse="\n"), sep="\n")
	tex <- paste(tex, "\\end{document}", sep="\n")

	tf <- paste(trimws(tempfile()), ".tex", sep="")
	write(tex, tf)
	folder <- paste(outdir, substr(Sys.time(), 1, 10), sep="")
	system(paste("mkdir ", folder, sep=""))
	pdfF <- strsplit(tf, "\\.")
	pdfF <- pdfF[[1]][length(pdfF[[1]])-1]
	pdfF <- strsplit(pdfF, "/")
	pdfF <- pdfF[[1]][length(pdfF[[1]])]
	pdfF <- paste(pdfF, ".pdf", sep="")
	cmd <- paste("cd ",folder," && pdflatex ", tf, sep="")
	system(cmd)
	system(paste("cp ", tf, " ", folder,"/", sep=""))
	#system(paste("mupdf ", folder, "/",pdfF, sep=""))

	return (list(tbl=pat, pdffile=paste(folder,"/",pdfF,sep="")))
    }
}

