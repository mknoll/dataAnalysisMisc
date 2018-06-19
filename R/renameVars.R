#' @title Rename variables in a data frame using a mapping file
#' 
#' @description Allows to easily utilize different version of variable
#' naming, e.g. for diffnerent languages.
#' 
#' @param data data.frame containing each variable in a spearate column
#' @param rename mapping file (csv) to rename vars. Contains the following columns:
#' OLD, NEW, TYPE, VARNAME, REFERENCE. OLD contains the old names, NEW the new ones,
#' TYPE is VARNAME for variables (colnames of data), and LEVEL for the corresponding
#' variable levels. VARNAME contains the old variable name to which each levels 
#' corresponds to. REFERENCE is the ordre for factors, 1 being the reference
#'
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
	sub$OLD <- trimws(sub$OLD)
	cpy[,i] <- as.character(cpy[,i])
	for (j in 1:length(sub[,1])) {
	    cpy[which(cpy[,i] == as.character(sub$OLD[j])),i] <- as.character(sub$NEW[j])
	}
	
	#reorder levels
	lvsNew <- levels(factor(cpy[,i]))
	lvsOrder <- sub$REFERENCE[match(lvsNew, sub$NEW)]
	lvsNew <- lvsNew[order(lvsOrder)]
	tryCatch({
	    cpy[,i] <- factor(cpy[,i], lvsNew) 
	}, error=function(e) {})
    }

    return(cpy)
}
