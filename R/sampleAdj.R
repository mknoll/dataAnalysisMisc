#' @title Sample clusters from adjacence matrix
#'
#' @description Desc
#'
#' @param adj Adjazenzmatrix
#' @param size Size of clusters (number of links)
#' @param n Number of clusters
#' @param maxTry numbers of random tries to establish clusters
#' @param starTopo create star topology
#' @param type If set to 0, n refers to number of edges, choose
#' 1 for number of nodes
#'
#' @useDynLib dataAnalysisMisc
#'
#' @export 
sampleAdj <- function(adj, size=5, 
		      n=5, maxTry=1000,
		      starTopo=1,
              type=1) {

    m <- c(as.matrix(unlist(adj), nrow=1))
    erg <- c(as.integer(rep(0, length(m))))
    ret <- .C("startAdj", matrix=as.integer(m), len=length(adj[1,]), 
	      n=as.integer(n), size=as.integer(size), erg=erg, 
	      maxTry=as.integer(maxTry), star=as.integer(starTopo),
          type=as.integer(type))
    tbl <- table(ret$erg)
    if (length(tbl)-1 != n) { 
	warning("To few clusters identified! Try increasing the maxTry parameter.")
    }
    print(tbl)

    matrix(ret$erg, ncol=length(adj[1,]), nrow=length(adj[1,]))
}
