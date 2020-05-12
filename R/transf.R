#' @title Transform data with zero
#'
#' @description Performs a transformation as proposed 
#' by Smithson & Verkuilen (2006) 
#' 
#' @param data vector or data.frame 
#' @param dir direction, 1: per row (default), 2: per column
#'
#' @export
transf <- function(dat, dir=1) {
    if (class(dat) %in%  c("data.frame", "matrix", "data.matrix")) {
	dat <- apply(dat, dir, function(x) (x*(length(x) -1) + 0.5)/length(x))
    } else {
	dat <- (dat*(length(dat)-1)+0.5)/length(dat)
    }
    dat
}
