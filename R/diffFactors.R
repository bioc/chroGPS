
## Subset MDS to selected points
subsetMDS <- function(m,sel)
    {
        mm <- new('mds',points=getPoints(m)[sel,],R.square=NULL,stress=NULL)
        mm
    }

pause = function()
    {
        if (interactive())
            {
                invisible(readline(prompt = "Press <Enter> to continue..."))
            }
    }

## Function to compare two MDS chroGPS Factor maps and report Procrustes sum of squares differences between maps by factor
## m1 and m2 contain valid MDS class objects with common points, with range normalized to -1,1
## ... gives additional graphical parameters for plotting
diffFactors <- function(m1,m2,name1='mds1',name2='mds2',minPoints=10,poslegend='topleft',plot=TRUE,pointcol.m1='red',pointcol.m2='green',textcol.m1='black',textcol.m2='grey',pch.m1=19,pch.m2=19,cex.m1=1,cex.m2=1,segcol='blue',...)
    {
        ## Perform Procrustes analysis, first adjust M2 to match M1, then report final differences
        p1 <- getPoints(m1)
        p2 <- getPoints(m2)
                
        ## Subset to common elements only, they have to be in rownames
        sel <- intersect(rownames(p1),rownames(p2))
        p1 <- p1[sel,]
        p2 <- p2[sel,]
        if (!unique(rownames(p1)==rownames(p2)) | (nrow(p1)<minPoints)) stop('Not enough common points to perform differential analysis')
        
        # Procrust and generate new MDS with solution
        pp <- vegan::procrustes(p1,p2)
        #pp <- mergeMDS(p1,p2)
        m3 <- new('mds',points=pp$Yrot,Type=m1@Type,Adj=TRUE,R.square=-1,stress=-1)
        p3 <- getPoints(m3)
                
        ## Plot before adjustment
        plot(0,xlim=c(-1,1),ylim=c(-1,1),xlab='',ylab='',xaxt='n',yaxt='n',col='NA',main='Before Adjustment')
        segments(p1[,1],p1[,2],p2[,1],p2[,2],col=segcol)
        par(new=TRUE)
        plot(m1,labels=rownames(p1),point.pch=pch.m1,point.cex=cex.m1,text.cex=0.75,point.col=pointcol.m1,text.col=textcol.m1,xlim=c(-1,1),ylim=c(-1,1),...)
        par(new=TRUE)
        plot(m2,labels=rownames(p2),point.pch=pch.m2,point.cex=cex.m2,text.cex=0.75,point.col=pointcol.m2,text.col=textcol.m2,xlim=c(-1,1),ylim=c(-1,1),...)
        legend(poslegend,col=c(pointcol.m1,pointcol.m2),pch=c(pch.m1,pch.m2),legend=c(name1,name2))
        pause()
        
        ## Plot after adjustment
        plot(0,xlim=c(-1,1),ylim=c(-1,1),xlab='',ylab='',xaxt='n',yaxt='n',col='NA',main='After Adjustment')
        segments(p1[,1],p1[,2],p3[,1],p3[,2],col=segcol)
        par(new=TRUE)
        plot(m1,labels=rownames(p1),point.pch=pch.m1,point.cex=cex.m1,text.cex=0.75,point.col=pointcol.m1,text.col=textcol.m1,xlim=c(-1,1),ylim=c(-1,1),...)
        par(new=TRUE)
        plot(m3,labels=rownames(p3),point.pch=pch.m2,point.cex=cex.m2,text.cex=0.75,point.col=pointcol.m2,text.col=textcol.m2,xlim=c(-1,1),ylim=c(-1,1),...)
        legend(poslegend,col=c(pointcol.m1,pointcol.m2),pch=c(pch.m1,pch.m2),legend=c(name1,name2))        
        pause()
                
        ## Plot Procrustes results
        plot(pp)
        pause()
        
        ## Plot summary of errors by point
        par(las=1,mar=c(4,12,4,4)); barplot(sort(residuals(pp),decreasing=TRUE),horiz=TRUE,xlim=c(0,max(residuals(pp))+.1),col=heat.colors(length(residuals(pp))),main='Procrustes errors')

        ## Return full results list
        return(list(mds1=m1,mds2=m2,mds3=m3,procrustes=pp))
    }
