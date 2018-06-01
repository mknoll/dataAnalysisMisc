#' @title Rename vars
#' @export
renameVars <- function(data, rename) {
    cpy <- data
    for (i in 1:length(data[1,])) {
	#varnames
	w <- which(map$OLD == colnames(data)[i] & map$TYPE == "VARNAME")
	if (length(w) == 0) { next }
	colnames(cpy)[i] <- as.character(map$NEW[w])
	
	#levels
	mp <- which(map$VARNAME == colnames(data)[i])
	if (length(mp) == 0) { next }
	sub <- map[mp,,drop=F] 
	cpy[,i] <- as.character(cpy[,i])
	for (j in 1:length(sub[,1])) {
	    cpy[which(cpy[,i] == as.character(sub$OLD[j])),i] <- as.character(sub$NEW[j])
	}
    }

    return(cpy)
}
