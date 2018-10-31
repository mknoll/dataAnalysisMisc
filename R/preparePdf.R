#' @title Create tex/pdf files 
#' @export
preparePdf <- function(pat, outdir, col="llcc", filename=NULL) {
    sgOut <- stargazer(pat, summary=F, rownames=F)
    for (i in 1:length(sgOut)) {
	#print(paste(i, ": ", sgOut[[i]]))
	if (i <= 6)  {
	    sgOut[[i]] <- ""
	}
	if (!(i >= 10 && i != 11 && i < length(sgOut)-3)) { next } 
	tmp <- sgOut[[i]]
	if (substr(tmp,1,1) != " ") {
	    #find first & and remove
	    feat <- substr(tmp, 1, regexpr("&", tmp, fixed=T)[[1]][1]-1)
	    rest <- substr(tmp, regexpr("&", tmp, fixed=T)[[1]][1]+1, nchar(tmp))
	    tmp <- paste("\\multicolumn{2}{l}{\\textbf{", feat, "}}", rest, collapse="")
	    sgOut[[i]] <- tmp
	}
    }
    sgOut[[7]] <- paste("\\begin{longtable}{@{\\extracolsep{5pt}} ", col, "}")
    sgOut[[length(sgOut)-1]] <- ""
    sgOut[[length(sgOut)]] <- "\\end{longtable}"


    tex <- "\\documentclass{article} \n \\usepackage[margin=0.1in]{geometry}  \n \\usepackage{longtable} \n \\pagestyle{empty} \\begin{document} "
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

    retName <- paste(folder, "/", pdfF, sep="")
    if (!is.null(filename)) {
	system(paste(" mv ", tf, " ", folder, "/", filename, ".tex", sep=""))
	system(paste(" mv ", folder, "/", pdfF , " ", folder, "/", filename, ".pdf", sep=""))
	retName <- paste(folder, "/", filename, ".pdf", sep="")
    } 
    #system(paste("mupdf ", folder, "/",pdfF, sep=""))
    return(retName)
}
