#' @title Cantor Encode
#' @export
encCant <- function(text) {
    words <- unlist(strsplit(text, " "))
    out <- c()

    ## TODO: forecah
    for (word in words) {
        wordR <- rev(strsplit(word, "|")[[1]])
        vals <- sapply(wordR, utf8ToInt)-32
        out <- c(out, cant(vals))
    }

    names(out) <- words
    return(out)
}



#' @title Cantor
#' @export 
cant <- function(x) {
    if (length(x) > 2) {
        tmp <- cant(x[-length(x)])
        print(paste("   ", tmp, " len: ", length(x), " x=", paste(x, collapse="|")))
        return (cant(c(tmp, x[length(x)])))
    } else {
        x0 <- x[1]
        y0 <- x[2]
        return (0.5*(x0+y0)*(x0+y0+1)+y0)
    }
}

