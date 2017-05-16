#' @title
#' Find hierarchical cluster-based groups for survival differences.
#' 
#' @description 
#' Expectes data which can be clustered as well as corresponding survival
#' data. Different clusterign methods are then tried on this data.
#' 
#' @export 
#' 
#' @import plyr
findHCLSurv <- function(data, srv, depth=7, minGrSize=3, p.valCutoff=0.1, dist=NULL) {
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)
    
    if (!is.null(dist)) {
        distances <- dist
    } else {
        distances <- c("ward.D", "ward.D2", "single", "complete", "average")#, "mcquitty", "median", "centroid")
    }
    
    if (length(levels(factor(as.matrix(data)))) < 2) {
        return (NULL)
    }
    
    ## collect data
    collectHC <- NULL
    for (disM in distances) {
        for (de in 2:depth) {
            pm <- pheatmap(data, silent = T, clustering_method = disM)
            pm.clust <- cutree(pm$tree_col, de)
            clust <- pm.clust
            
            ##test all single groups against each other
            out <- foreach (i=1:de) %dopar% {
                ret <- NULL
                for (j in (i+1):de) {
                    if (j > de) { break }
                    if (length(clust[clust==i]) <= minGrSize ||
                            length(clust[clust==j]) <= minGrSize) {
                        next
                    }
                    if (all(is.na(srv[clust==i])) ||
                            all(is.na(srv[clust==j]))) {
                        next
                    }
                    sdf <- survdiff(srv[clust == i | clust == j]~clust[clust == i | clust == j])
                    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
                    if (p.val <= p.valCutoff) {
                        vec <- data.frame(dist=disM, n=de, i=i, j=j,
                                          p.val=p.val, 
                                          gr_i=length(clust[clust==i]), gr_j=length(clust[clust==j]))
                        #collectHC <- rbind.fill(collectHC, vec)
                        ret <- rbind.fill(ret, vec)
                    }
                }
                ret
            }
            collectHC <- rbind.fill(collectHC, do.call(rbind, out))
        }
    }
    doParallel::stopImplicitCluster()
    
    if (!is.null(collectHC)) {
        collectHC <- collectHC[order(collectHC$p.val),,drop=F]
    }
    
    return(collectHC)
}


#' @title 
#' Print Kaplan Meier plots 
printHCLCand <- function(data, srv, n, i, j, dist, main="") {
    pm <- pheatmap(data, silent=T, clustering_method = dist)
    pm.clust <- cutree(pm$tree_col, n)
    clust <- pm.clust
    pheatmap(data, annotation_col=data.frame(clust), clustering_method = dist, cluster_rows = F)
    
    par(mfrow=c(2,2))
    plot(survfit(srv~clust), col=1:n, lwd=2, xlim=c(0,30), main=paste(main, ":", n,i,j,dist))
    plot(survfit(srv[clust == i | clust == j]~clust[clust == i | clust == j]), col=c(i,j), 
         lwd=2, xlim=c(0,30), main=paste(main, ":", n,i,j,dist))
    
    sdf <- survdiff(srv[clust == i | clust == j]~clust[clust == i | clust == j])
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    
    text(5, 0.05, round(p.val, 4), font=2)
}


