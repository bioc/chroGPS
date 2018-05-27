## Rank Factors Globally, using whole map information and glm/lm
rankFactorsbyProfile <- function(x,minFactors=5,ranktype='glm',glm.threshold=0.5,verbose=TRUE,maxIter=ncol(x),mc.cores=1)
    ## removes first those factors that can be more accurately predicted using the rest
    {
        factorTable <- factorTable2 <- as.matrix(uniqueCount(x)[,2:(ncol(x)+1)])
        res <- vector('list',min(ncol(factorTable)-minFactors,maxIter))
        res.idx <- 1
        while(ncol(factorTable)>minFactors & res.idx<=maxIter)
            {
                maxpred <- 0; maxfactor <- 0
                if (ranktype=='glm')
                    {
                        ppred <- unlist(parallel::mclapply(1:ncol(factorTable),function(i)
                                                 {
                                                     glm1 <- glm(factorTable[,i] ~ factorTable[,-c(i)],family=binomial)
                                                     lpred <- predict(glm1) # returns real predictions, have to rescale to binary 1s and 0s
                                                     ppred <- 1/(1+exp(-lpred)) # Probabilities between 0 and 1. If ppred <.5, value is 0, else value is 1
                                                     ppred <- ifelse(ppred>glm.threshold,1,0) # 1 if p>.5, 0 otherwise
                                                     prop <- length(which(ppred==factorTable[,i]))/length(ppred) # Proportion of good predictions
                                                     ##if (prop>maxpred) { maxpred <- prop; maxfactor <- i; factorName <- colnames(factorTable)[maxfactor] }
                                                     ##print(sprintf('# Prediction capability for %s: %.3f',colnames(factorTable)[i],prop))
                                                     prop
                                                 },mc.cores=mc.cores))
                    }
                else if (ranktype=='lm')
                    {
                        ppred <- unlist(parallel::mclapply(1:ncol(factorTable),function(i)
                                                 {
                                                     lm1 <- lm(factorTable[,i] ~ factorTable[,-c(i)])
                                                     prop <- summary(lm1)$r.squared
                                                     ##print(sprintf('Linear regression R2 for %s: %.3f',colnames(factorTable)[i],prop))
                                                     ##if (prop>maxpred) { maxpred <- prop; maxfactor <- i; factorName <- colnames(factorTable)[maxfactor] }
                                                     ##print(sprintf('# Prediction capability for %s: %.3f',colnames(factorTable)[i],prop))
                                                     prop
                                                 },mc.cores=mc.cores))
                    }
                else stop('Wrong predictor specification')
                names(ppred) <- colnames(factorTable)
                factorName <- names(which.max(ppred))
                maxpred <- as.numeric(which.max(ppred))
                if (verbose) print(sprintf('# Removing %s factor with %s prediction score %.3f',factorName,ranktype,ppred[maxpred]))
                ##d2 <- distGPS(factorTable2,metric='tanimoto')
                ##mycor <- cor(d@d,d2@d) # if you want to compute cors
                ##print(sprintf('# Correlation of new gene distances with those from full dataset: %.3f',mycor))
                #print(ppred)
                res[[res.idx]] <- ppred
                names(res)[res.idx] <- factorName
                res.idx <- res.idx + 1
                factorTable <- factorTable[,-c(maxpred)]
            }
        return(res)
    }
