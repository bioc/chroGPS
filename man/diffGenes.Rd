\name{diffGenes}
\alias{diffGenes}
\alias{find.fdr}
\alias{find.threshold}
\alias{diffGenes-methods}
\alias{plotDiffGenes}
\title{
Performs differential analysis of chroGPS-genes maps based on Bayesian
Correct Classification Rates.
}
\description{
The function uses Bayesian posterior probability of correct
classification (CCR) of chroGPS-genes maps as computed by the
\code{clusGPS} function to identify statistically significant epigenetic
changes between two different backgrounds (technical and/or biological).
The functions \code{find.fdr} and \code{find.threshold} are used
internally to compute approximated FDRs from a vector of input posterior
probabilities of correct classification (CCRs).
}
\usage{
diffGenes(xy,m,clus,label.x,label.y,clusName=NULL,fdr=TRUE,mc.cores=1)
}
\arguments{
  \item{xy}{
    A \code{matrix} or \code{data.frame} object, containing a combined
    matrix with whole-genome epigenetic information for backgrounds X
    and Y, as obtained by the \code{combineGenesMatrix} function.
  }
  \item{m}{
    A \code{mds} object, containing a chroGPS-genes MDS map, obtained
    for the combined matrix xy.
  }
  \item{clus}{
    A \code{clusGPS} object, containing a chroGPS-genes cluster solution, obtained
    for the provided background via the \code{clusGPS} and
    \code{mergeClusters} functions.
  }
  \item{label.x}{
    Name for background X, by default 'mds1'
  }
  \item{label.y}{
    Name for background Y, by default 'mds2'
  }
  \item{clusName}{
    If \code{clus} contains more than one clustering solution
    (i.e. \code{cutree} dendrogram height cut), the name of the selected one. If no name is provided, the first element is used.
  }
  \item{fdr}{
    Set to TRUE (default) to compute False Discovery Rate of cluster
    correct classification from posterior probabilities.
  }
  \item{mc.cores}{
    Number of cores to use for parallel computations.
  }
}
\value{
The function returns a data frame with posterior probabilities and FDRs
of cluster classifications for all genes in backgrounds X and Y, so that
significant changes can be traced and reported easily.
}
\author{
Oscar Reina.
}
\seealso{
See functions \code{distGPS, mds, clusGPS, profileClusters} for
generating chroGPS-genes maps. See function \code{plotTransitions} for
examples on how to visualize differential gene maps results.
}
\examples{

# Summarize factor replicates with method 'any' so that 1 replicate
#  having the mark is enough
# Assuming s2.tab and bg3.tab contain the full datasets for dm3 genome
#  and all factors used in Font-Burgada et al.
# s2.tab <- mergeReplicates(s2.tab,s2names$Factor,'any')
# bg3.tab <- mergeReplicates(bg3.tab,bg3names$Factor,'any')

# Join, use common factors. Then use common genes only from those that
#  have at least one mark in both s2 and bg3
# x <- combineGenesMatrix(s2.tab,bg3.tab,'S2','BG3')

# Build map and cluster as always
# d <- distGPS(x,metric='tanimoto',uniqueRows=TRUE)

# Not run
# m <- mds(d,type='classic',splitMDS=TRUE,split=0.16,mc.cores=4)
# mm <- mds(d,m,type='boostMDS',samplesize=0.005,mc.cores=6)
# unique(rownames(m@points)==rownames(mm@points)) # sanity check
# This should be incorporated in the function code...
#m@points <- m@points[rownames(d@d),]
#mm@points <- mm@points[rownames(d@d),]

# Cluster
# h <- hclust(as.dist(d@d),method='average')
# clus <-
# clusGPS(d,mm,h,k=max(cutree(h,h=0.5)),ngrid=10000,mc.cores=8,recalcDist=FALSE,verbose=FALSE)
# clus.merged <- mergeClusters(clus,brake=0,mc.cores=8)
# clus
# clus.merged

# Use new function to profile clusters
# pc <- profileClusters2(x,clus.merged,normalize=TRUE)
# pheatmap(pc,trace='none',scale='none',col=bluered(100))

# Perform differential analysis
# x.diff <- res <-
# diffGPS.clus(x,mm,clus.merged,label.x='S2',label.y='BG3',clusName=clusNames(clus.merged)[1],fdr=TRUE,mc.cores=8)
# write.csv(x.diff,'diffgenes_fdrest.csv')

# Select genes changing clusters with FDR 0.05
# xx.diff <- x.diff[x.diff$ClusID.S2!=x.diff$ClusID.BG3 &
# x.diff$FDR.S2<0.25 & x.diff$FDR.BG3<0.25,]
# xx.diff$CC <- paste(xx.diff$ClusID.S2,xx.diff$ClusID.BG3,sep='.')
# head(sort(table(xx.diff$CC),decreasing=TRUE))
# write.csv(xx.diff,'kk_fdrest2.csv')

# Perform enrichment test using getEnrichedGO from chippeakanno package
# library(ChIPpeakAnno)
# library(org.Dm.eg.db)

# enriched.GO <- lapply(c('2.9','5.2'),function(cc) {
#        fbid <- as.character(xx.diff$geneid[xx.diff$CC==cc])
#            if (length(fbid)>=25)
#                        ans <-
#        getEnrichedGO(annotatedPeak=fbid,orgAnn='org.Dm.eg.db',maxP=0.05,multiAdjMethod='BH')
#            else ans <- NULL
#            return(ans)
#    })


# names(enriched.GO) <- c('2.9','5.2')
# enriched.GO <- enriched.GO[unlist(lapply(enriched.GO,length))>0]
# enriched.GO <- lapply(enriched.GO,function(x) lapply(x,function(y)
# unique(y[,-ncol(y)])))
# lapply(enriched.GO,head)

# Plot results with diffGPS.plot function
# res.sel <- res[res$ClusID.S2!=res$ClusID.BG3,]

# Plot
# diffGenes.plot(x,mm,clus.merged,res.sel,transitions='10.2',label.x='S2',label.y='BG3',fdr1=0.25,fdr2=0.25)

}
\keyword{ ~mds }
\keyword{ ~clustering }
