## standard Generic
setGeneric("mergeReplicates", function(x,id,mergeBy='any',mc.cores=1) standardGeneric("mergeReplicates"))

## Genes x Factors matrix
setMethod("mergeReplicates",signature(x='matrix'),
          function(x,id,mergeBy='any',mc.cores=1) {
              mergeReplicateMatrix(x=x,id=id,mergeBy=mergeBy,mc.cores=mc.cores)
}
)

## Genes x Factors data frame
setMethod("mergeReplicates",signature(x='data.frame'),
          function(x,id,mergeBy='any',mc.cores=1) {
              mergeReplicateMatrix(x=x,id=id,mergeBy=mergeBy,mc.cores=mc.cores)
}
)

## Genomic Intervals GRangesList
setMethod("mergeReplicates",signature(x='GRangesList'),
          function(x,id,mergeBy='any',mc.cores=1) {
              mergeReplicateList(x=x,id=id,mergeBy=mergeBy,mc.cores=mc.cores)
}
)

## Genomic Intervals List
setMethod("mergeReplicates",signature(x='list'),
          function(x,id,mergeBy='any',mc.cores=1) {
              x <- GRangesList(x)
              mergeReplicates(x=x,id=id,mergeBy=mergeBy,mc.cores=mc.cores)
}
)

## Factors MDS (Genes too, but no much sense doing that...)
setMethod("mergeReplicates",signature(x='mds'),
          function(x,id) {
              mergeReplicateMDS(x=x,id=id)
}
)
## Merge replicates (Genes x Factors)
## Basically collapses columns to proportion of replicates with factor and then updates them before returning depending on selected mergeBy option
mergeReplicateMatrix <- function(x,id,mergeBy='any',mc.cores=mc.cores)
    {
        if (is.numeric(mergeBy))
            if (!(mergeBy<=1 & mergeBy>0)) stop('mergeBy has to be either any, all, or a number between (0,1]')
        if (is.character(mergeBy) & !(mergeBy %in% c('any','all'))) stop('mergeBy has to be either any, all, or a number between (0,1]')
        
        if (ncol(x) != length(id)) stop('Non-conformable data / id')
        
        ## Collapse with selected policy
        x <- do.call(cbind,lapply(sort(unique(id)),function(i) rowMeans(as.data.frame(x[,id==i]))))
        if (mergeBy=='any') x[x>0] <- 1
        else if (mergeBy=='all') x[x<1] <- 0
        else if (mergeBy>0 & mergeBy<=1) x[x<mergeBy] <- 0
        else stop('mergeBy has to be either any, all, or a number between (0,1]')
        colnames(x) <- sort(unique(id))
        return(x)
    }

## Merge Replicates (Genomic Intervals list)
## For each replicate group checks proportion of samples for each replicate set group. Only peaks matching the desired mergeBy option are returned for that replicate
mergeReplicateList <- function(x,id,mergeBy='any',mc.cores=mc.cores)
{
    if (length(x)!=length(id)) stop('Non-conformable arguments')
    ans <- lapply(sort(unique(id)),function(i)
                    {
                        sel <- which(id==i)
                        if (length(sel)==1) {
                            return(x[[sel]])
                        }
                        else if (length(sel)>1)
                            {
                                ## Concatenate sites
                                #txt <- paste("allsites <- c(",paste("x[[",sel,"]]",collapse=','),")",sep='')
                                #eval(parse(text=txt))
                                allsites <- unlist(x[sel])
                                suppressWarnings(n <- do.call(cbind,lapply(as.list(x[sel]),function(y) as.numeric(unlist(allsites %over% y)))))
                                n <- rowSums(n) / ncol(n)
                                if (mergeBy=='any') return(allsites[n>0])
                                else if (mergeBy=='all') return(allsites[n==1])
                                else return(allsites[n>=mergeBy])
                            }
                    })
    names(ans) <- sort(unique(id))
    ans <- GRangesList(ans)
    return(ans)
}

## mergeReplicates (MDS map, meant for factors, will work with anything as long as things are conformable)
mergeReplicateMDS <- function(x,id)
    {
        p <- getPoints(x)
        if (nrow(p) != length(id)) stop('Non-conformable arguments')
        p <- matrix(unlist(by(p,id,colMeans)),ncol=ncol(p),byrow=TRUE)
        rownames(p) <- sort(unique(id))
        x@points <- p
        return(x)
    }
