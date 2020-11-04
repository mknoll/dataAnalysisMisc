#' @title Create mut matrix
#'
#' @param files maf files 
#' @param excl Variant_Classification to exclude
#' 
#' @import data.table
#' @export
mutMatrix <- function(files, excl=c("Silent")) {
    mafs <- list()    
    for (i in 1:length(fl)) {    
	print(i)    
	tmp <- data.frame(fread(fl[i]))    
	type <- strsplit(fl[i], "\\.")[[1]][3]    
	mafs[[length(mafs)+1]] <- tmp    
	names(mafs)[length(mafs)] <- type    
    }    
    mafs2 <- lapply(mafs, function(x) x[which(!x$Variant_Classification %in% excl),])    

    #TODO: check Sample Barcode -> type (tumor, ...)
    mafs2 <- lapply(mafs, function(x) {    
			x$ID <- substr(x$Tumor_Sample_Barcode, 1,12)    
			x })    
    genes <- unique(unlist(lapply(mafs2, function(x) x$SYMBOL)))    
    pat <- unique(unlist(lapply(mafs2, function(x) x$ID)))    
    m<-list()    
    for (i in 1:length(mafs2)) {    
	print(i)    
	tmp <- mafs2[[i]]    
	subM <- list()    
	for (p in pat) {    
	    tmp2 <- tmp[which(tmp$ID == p),]    
	    subM[[length(subM)+1]] <- ifelse(genes %in% tmp2$SYMBOL, 1,0)    
	}    
	subM <- do.call(cbind, subM)    
	colnames(subM) <- pat    
	m[[length(m)+1]] <- subM    
	names(m)[length(m)] <- names(m)[i]    
    }    
    m2 <- m[[1]]
    for (i in 2:length(m)) {
	m2 <- m2+m[[i]]
    }
    rownames(m2) <- genes    
    colnames(m2) <- pat

    m2
}
