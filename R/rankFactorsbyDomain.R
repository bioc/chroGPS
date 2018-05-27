5## standard Generic
setGeneric("rankFactorsbyDomain", function(d,sampleinfo,ranktype='domainDist',selName='Color',selValue,k=NULL,mc.cores=1) standardGeneric("rankFactorsbyDomain"))

## Genes x Factors matrix
#setMethod("rankFactorsbyDomain",signature(d="matrix"),
#          function(d,sampleinfo,ranktype='domainDist',selName='Color',selValue,k=NULL,mc.cores=1) {                        
#              rankFactorsbyDomain(d,sampleinfo,ranktype,selName,selValue,k,mc.cores)
#          }
#          )

## distGPS
setMethod("rankFactorsbyDomain",signature(d="distGPS"),
          function(d,sampleinfo,ranktype='domainDist',selName='Color',selValue,k=NULL,mc.cores=1) {
              d <- chroGPS::as.matrix(d)
              rankFactorsbyDomain(d,sampleinfo,ranktype,selName,selValue,k,mc.cores)
          }
          )

## Genes x Factors data frame
setMethod("rankFactorsbyDomain",signature(d="data.frame"),
          function(d,sampleinfo,ranktype='domainDist',selName='Color',selValue,k=NULL,mc.cores=1) {
              d <- as.matrix(d)
              rankFactorsbyDomain(d,sampleinfo,ranktype,selName,selValue,k,mc.cores)
          }
          )

## Rank Factors by Domain, using intra/inter domain distance
setMethod("rankFactorsbyDomain",signature(d="matrix"),
          function(d,sampleinfo,ranktype='domainDist',selName='Color',selValue,k=NULL,mc.cores=1)
          {
              domains <- sampleinfo
              if (!selName %in% colnames(domains)) stop('Invalid sel name in sampleinfo')
              ff <- rownames(domains)[domains[,selName]==selValue]
              ff <- combs(ff,k=k)
              cat(sprintf('\n### Evaluating %s domain ###\n',selValue))
              dd <- parallel::mclapply(1:nrow(ff),function(i) {
                  ##cat(sprintf('\n# %s #\n',paste(ff[i,],collapse='.')))
                  sel <- rownames(domains) %in% ff[i,] | domains[,selName]!=selValue
                  inter <- domainDist(d[sel,sel],gps='factors',domain=domains[,selName][sel],type='inter',plot=FALSE,avg=FALSE)
                  inter <- inter[grep(selValue,names(inter))]
                  inter <- mean(unlist(inter))
                  intra <- mean(domainDist(d[sel,sel],gps='factors',domain=domains[,selName][sel],type='intra',plot=FALSE,avg=TRUE)$Avg)
                  return(list(inter=inter,intra=intra))
              },mc.cores=mc.cores)
              names(dd) <- apply(ff,1,paste,collapse=' + ')
              dd
          }
 )

