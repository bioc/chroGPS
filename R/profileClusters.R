## New function to profile clusters
## x: combined dataset as output from diffGPS
## clus: clustering results output from clusGPS
## clusName: the clustering solution we want to profile. If null, first one is used
## normalize: if TRUE, log2 ratio for mark proportion against whole genome is returned. if FALSE, only proportion of genes with that mark in that cluster.
## log2: log2 of ratio/proportion results.
## mc.cores: ...
profileClusters <- function(x,clus,clusName=NULL,normalize=FALSE,mc.cores=1,...)
    {
        x <- as.data.frame(x)
        ## uid <- apply(x,1,paste,collapse=',') # this is much faster via eval parse
        uid <- NULL # to make it visible for check
        txt <- paste("uid <- paste(",paste("x[,",1:ncol(x),"]",collapse=","),", sep=',')",sep="")
        eval(parse(text=txt))
        if (is.null(clusName)) clusName <- clusNames(clus)[1]
        clusID <- clus@clus[[clusName]]$id[uid]
        if (normalize)
            ans <- do.call(rbind,parallel::mclapply(sort(unique(clusID)),function(i) { x.sel <- x[clusID==i,]; log2(1+colMeans(x.sel)) - log2(1+colMeans(x)) },mc.cores=mc.cores))
        else ans <- do.call(rbind,parallel::mclapply(sort(unique(clusID)),function(i) { x.sel <- x[clusID==i,]; log2(1+(colMeans(x.sel))) },mc.cores=mc.cores))
        rownames(ans) <- sort(unique(clusID))
        return(ans)
    }
