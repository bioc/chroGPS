## Functions to perform differential gene maps analysis

## ## Merge replicates (columns)
## mergeReplicates <- function(x,id,mergeBy='any')
##     {
##         if (is.numeric(mergeBy))
##             if (!(mergeBy<=1 & mergeBy>0)) stop('mergeBy has to be either any, all, or a number between (0,1]')
##         if (is.character(mergeBy) & !(mergeBy %in% c('any','all'))) stop('mergeBy has to be either any, all, or a number between (0,1]')

##         if (ncol(x) != length(id)) stop('Non-conformable data / id')

##         ## Collapse with selected policy
##         x <- do.call(cbind,lapply(sort(unique(id)),function(i) rowMeans(as.data.frame(x[,id==i]))))
##         if (mergeBy=='any') x[x>0] <- 1
##         if (mergeBy=='all') x[x<1] <- 0
##         if (mergeBy>0 & mergeBy<=1) x[x<mergeBy] <- 0
##         colnames(x) <- sort(unique(id))
##         return(x)
##     }

## distGPS extensions to generate differential maps
## x1, x2 are tables of genesxfactors, with no repeated row or colnames (ie one column per factor, unified by whatever criteria you like before calling this)
combineGenesMatrix <- function(x,y,label.x,label.y,minFactors=10,minGenes=1000)
    ## stacks together genes from x and y, using common factors and removing genes without marks in those
    {
        ## Common factors
        selfactors <- sort(intersect(colnames(x),colnames(y)))
        if (length(selfactors)<minFactors) stop('Not enough common factors between both datasets')
        ## Genes with marks for common factors in both them
        x.genes <- rownames(x)[rowSums(x[,selfactors])>0]
        y.genes <- rownames(y)[rowSums(y[,selfactors])>0]
        allgenes <- sort(intersect(x.genes,y.genes))
        if (length(allgenes)<minGenes) stop('Not enough common genes with information in both datasets')
        ans <- rbind(x[allgenes,selfactors],y[allgenes,selfactors])
        rownames(ans) <- c(paste(allgenes,label.x,sep='.'),paste(allgenes,label.y,sep='.'))
        ans
    }

diffGenes <- function(xy,m,clus,label.x,label.y,clusName=NULL,fdr=TRUE,mc.cores=1)
    ## retrieves cluster ID information for genes in x and y, with cluster id and pp per gene
    ## returns nice table with gene geneid, map position, cluster id and pp for X and Y
    {
        label.xx <- paste('.',label.x,sep='')
        label.yy <- paste('.',label.y,sep='')
        geneid <- as.character(rownames(xy)) # nice geneids without background identifiers
        geneid <- gsub(label.yy,'',gsub(label.xx,'',geneid))
        if (is.null(clusName)) { clusName <- clusNames(clus)[1] }
        uid <- NULL # to make it visible to the cmd check
        txt <- paste("uid <- paste(",paste("xy[,",1:ncol(xy),"]",collapse=","),",sep=',')",sep="")
        eval(parse(text=txt))
        coords <- m@points[uid,]
        clusid <- clus@clus[[clusName]]$id
        propclass <- unlist(parallel::mclapply(names(clusid),function(p) clus@clus[[clusName]]$propclass[[clusid[p]]][p],mc.cores=mc.cores))
        ans <- data.frame(type=c(rep(label.x,nrow(xy)/2),rep(label.y,nrow(xy)/2)),geneid=geneid,uid=uid,X=coords[uid,1],Y=coords[uid,2],ClusID=clusid[uid],ProbClus=propclass[uid])
        sel <- ans$type==label.x
        ans1 <- ans[sel,-1]; 
        ans2 <- ans[!sel,-1:-2]; 
        ## Add FDR if wanted
        if (fdr)
            {
                ans1 <- ans1[order(ans1$ProbClus),]
                ans1$FDR <- find.fdr(ans1$ProbClus,ans1$ProbClus)$fdr
                ans2 <- ans2[order(ans2$ProbClus),]
                ans2$FDR <- find.fdr(ans2$ProbClus,ans2$ProbClus)$fdr
            }
        colnames(ans1)[-1] <- paste(colnames(ans1)[-1],label.x,sep='.')
        colnames(ans2) <- paste(colnames(ans2),label.y,sep='.')
        #ClusIDs <- paste(ans1$ClusID.X,ans2$ClusID.Y,sep='.')
        #cbind(ans1,ans2,ClusIDs)
        cbind(ans1,ans2)
    }

# Add functions to cut PDE using a given FDR threshold (MDA find.threshold)
# I also use this function to compute FDRs for each individual pp
find.fdr <- function(v,threshold) {
    # Input
    #   v: vector with posterior probabilities of differential expression (missing values are removed)
    #   k: vector with thresholds to declare differential expression (defaults to a sequence between min(v),max(v))
    # Output: a list containing the following elements
    #   threshold: the k vector
    #   fdr: estimated FDR
    #   fnr: estimated FNR

    if (missing(v)) stop('The vector of posterior probabilities v must be specified')
    if (sum(is.na(v))>0) v <- v[!is.na(v)]
    if (missing(threshold)) threshold <- seq(min(v),max(v)-diff(range(v))/100,diff(range(v))/100)

    fdr <- fnr <- rep(NA,length(threshold))
    for (i in 1:length(threshold)) {
        fdr[i] <- sum((v>=threshold[i])*(1-v))/sum(v>=threshold[i])
        fnr[i] <- sum((v<threshold[i])*v)/sum(v<threshold[i])
    }
    
    return(list(threshold=threshold,fdr=fdr,fnr=fnr))
}

find.threshold <- function(pde,fdr,eps=.0001) {
    # Input
    #   pde: vector of probabilities of differential expression
    #   fdr: desired FDR
    # Output:
    #   threshold: threshold for probability of differential expression with estimated FDR closest to 'fdr'
    #   fdrest: estimated FDR

    pde <- pde[order(pde)]
    imin <- 1; imax <- length(pde); i <- round(imax/2)
    fdrest <- find.fdr(pde,pde[i])$fdr
    while ((imax-imin>1) & (abs(fdr-fdrest)>eps)) {
        if (fdrest>fdr) { imin <- i } else { imax <- i }
        i <- round(.5*(imin+imax))
        fdrest <- find.fdr(pde,pde[i])$fdr
        print(fdrest)
      }
    return(list(threshold=pde[i],fdrest=fdrest))

}

# Compute cluster centroids
#Compute and plot cluster centroids and plot cluster number
# m is mds object
# clus is clusGPS object
# clusName for cluster configuration in case more than one is within object
plotCenters <- function(m,clus,clusName=NULL,col='#00000090',bg='#FFFFFF90',lwd=4,pch=21,cex=5,text.cex=2,text.col='black',plot=TRUE,...)
{
    if (is.null(clusName)) { clusName <- clusNames(clus)[1]; warning('More than one cluster configuration available, used first element') }
    m@points <- m@points[names(clus@clus[[clusName]]$id),] # ensure same order of elements is used
    centers <- do.call(rbind,by(m@points,clus@clus[[clusName]]$id,colMeans)) # Compute mean of coordinates for every cluster
    if (plot)
        {
            points(centers,pch=pch,cex=cex,col=col,bg=bg,lwd=lwd)
            text(x=centers[,1],y=centers[,2],1:nrow(centers),cex=2,col='black')
        }
    centers
}

plotTransitions <- function(x,m,clus,res,plotCentroids=TRUE,transID,xlim,ylim,...)
    {
        # split data from each background
        fdrpos <- length(grep('FDR',colnames(res)))
        if (fdrpos==0)
            {
                x1 <- res[,1:6]
                x2 <- res[,c(1,7:11)]
            }
        else if (fdrpos==2)
            {
                x1 <- res[,1:7]
                x2 <- res[,c(1,8:13)]
            }
        else stop('Invalid Differential results table')

        ## check
        rownames(x1) <- x1$geneid
        rownames(x2) <- x2$geneid
        
        ## Select genes and plot
        transitions <- res$cc
        genes2plot <- as.character(res$geneid[transitions==transID])
        points(x1[genes2plot,c(3,4)],pch=19,col='blue',xlim=xlim,ylim=ylim)
        points(x2[genes2plot,c(3,4)],pch=19,col='red',xlim=xlim,ylim=ylim)
        segments(x1[genes2plot,3],x1[genes2plot,4],x2[genes2plot,3],x2[genes2plot,4],col=rgb(0,1,0,.35),lwd=2)
        if (plotCentroids & length(genes2plot)>=5)
            {
                points(ellipse::ellipse(x=cov(x=x1[genes2plot,c(3,4)]),centre=colMeans(x1[genes2plot,c(3,4)]),level=.5),type='l',xlim=xlim,ylim=ylim,col='blue',lwd=4)
                points(ellipse(x=cov(x=x2[genes2plot,c(3,4)]),centre=colMeans(x2[genes2plot,c(3,4)]),level=.5),type='l',xlim=xlim,ylim=ylim,col='red',lwd=4)
            }
    }

## Add functions to plot these results...
## x: combined GenesxFactors matrix as output from diffGPS
## m: an MDS obtained from x via distGPS and mds
## clus: a clustering from x and m, obtained via clusGPS
## res: a differential analysis result obtained via diffGPS.clus
plotDiffGenes <- function(x,m,clus,res,label.x='1',label.y='2',clusName=NULL,transitions=NULL,plotIDs=TRUE,fdr1=0.05,fdr2=0.05,minGenes=1,plotCentroids=TRUE,centroid.lwd=2,centroid.lty=1,point.col='lightgrey',probContour=0.75,contour.lwd=2,contour.lty=2,dolegend=TRUE,poslegend='bottomright',xlim,ylim,...)
    {
        ## checks
        if (is.null(clusName)) clusName <- clusNames(clus)[1]
        if (missing(xlim)) xlim <- range(m@points)
        if (missing(ylim)) ylim <- range(m@points)

        ## Generate transition table and prefilter by FDR if desired
        fdr <- res[,grep('FDR',colnames(res))]
        if (!ncol(fdr) %in% c(0,2)) stop('Invalid Differential results table')
        if (ncol(fdr)==2) res <- res[fdr[,1]<fdr1 & fdr[,2]<fdr2,]
        res$cc <- apply(res[,grep('ClusID',colnames(res))],1,paste,collapse='.')

        # Plot selected transitions
        if (is.null(transitions)) transitions <- sort(unique(res$cc))
        for (mycc in transitions)
            {
                ## Main map with cluster contours and IDs
                plot(m,point.col=point.col,xlim=xlim,ylim=ylim,...)
                plot(clus,type='contours',k=clusName,lwd=contour.lwd,probContour=probContour,lty=contour.lty,xlim=xlim,ylim=ylim,...)
                ## Plot Genes changing clusters (plotTransitions)
                if (sum(res$cc==mycc)>=minGenes) plotTransitions(x,m,clus,res,plotCentroids=TRUE,transID=mycc,xlim=xlim,ylim=ylim)
                ## Legend with cluster size
                if (dolegend) {
                    #tclus <- table(clusterID(clus,clusName))
                    #names(tclus) <- clusID
                    tclus <- tabClusters(clus,clusName)
                    legend(poslegend,legend=paste(names(tclus),as.numeric(tclus),sep=':'),col=rainbow(length(tclus)),pch=19,bty='n',cex=1,title='Cluster: # of elements')
                    }
                ## Compute cluster centers and plot cluster number
                #clusID <- 1:max(clusterID(clus,clusName))
                plotCenters(m,clus,clusName=clusName)
                ## Final legend
                legend('topleft',legend=sprintf('r2=%.3f / stress=%.3f',m@R.square,m@stress))
                legend('topright',legend=sprintf('Cluster Transition: %s\nFDR=%.2f/%.2f (n=%d)',mycc,fdr1,fdr2,sum(res$cc==mycc)))
                legend('bottomleft',col=c('blue','red'),lwd=2,legend=c(label.x,label.y))
                ##legend('bottomleft',col=c('blue','red'),lwd=2,legend=c('S2','BG3'),title=sprintf('elliptic hull (50%%)\n(S2 genes lost %d (%.2f %%))',table(genes.clus$u)[as.character(cc)],loseprop[as.character(cc)]*100),bty='n')
            }
    }

enrichedGO.plot <- function()
    {
        ## Plot GO enrichments
        #sel <- goenrich$S2.BG3.Clus==cc
        #if (sum(sel)>0)
        #    {
        #        par(new=TRUE)
        #        bp <- barplot(goenrich$pp[sel],horiz=TRUE,main=sprintf('Enriched GO terms in cluster transition %s',cc),col='white',xlim=c(-8,6),axes=FALSE,border='darkgrey')
        #        abline(v=-log10(c(0.1,0.01,0.001,0.0001,0.000001)),col='lightgrey')
        #        abline(v=-log10(0.05),col='red')
        #        text(x=.1,y=bp,goenrich$name[sel],cex=.75,adj=0)
        #        axis(3,at=1:6,cex.axis=.5)
        #    }
        #else text(x=.7,y=mean(range(m@points)),'No significant enrichments found',cex=1.5,adj=0)
    }


