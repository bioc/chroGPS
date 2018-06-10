\name{rankFactorsbyProfile}
\alias{rankFactorsbyProfile}
\title{
Function to help selecting candidate epigenetic factors based on their predictive capabilities.
}
\description{
Function to help selecting candidate epigenetic factors based on their predictive capabilities.
This function evaluates how information from a certain set of epigenetic factors can be used to accurately predict information from the rest using linear and logistic regression.
}
\usage{
rankFactorsbyProfile(x,minFactors=5,ranktype='glm',glm.threshold=0.5,verbose=TRUE,maxIter=ncol(x),mc.cores=1)
}
\arguments{
  \item{x}{
A genes x factors table as used in the distGPS function for generation of epigenetic gene maps.
}
  \item{minFactors}{
Minimum set of 'core' factors to retain. Factors are removed based on how accurately they can be predicted by others.
}
  \item{ranktype}{
Either 'lm' for using linear regression, or 'glm' for logistic.
}
  \item{glm.threshold}{
For 'glm', threshold to round values obtained from the predict function so that they are compared with the real values. Defaults to 0.5.
}
  \item{verbose}{
Informs about progression of the linear or logistic regression steps.
}
  \item{maxIter}{
Deprecated.
}
  \item{mc.cores}{
Cores to use in calls to parallel::mclapply.
}
}
\value{
A data frame with the epigenetic factors in the order they are removed and their prediction accuracy.
}
\author{
Oscar Reina.
}
\seealso{
\code{\link{distGPS}}, \code{\link{domainDist}}, \code{\link{rankFactorsbyDomain}}
}
\examples{
## Not run

data(s2)

# Unknown domains
# Perform computation
glm.rank <- rankFactorsbyProfile(s2.tab,ranktype='glm',glm.threshold=0.75,mc.cores=1)

# Returned objects are lists named by the factor with highest prediction accuracy in each interation, as well as the rest, let's generate a matrix
library(gtools)
glm.rank <- do.call(smartbind,glm.rank)

# Now let's order based on which factor is removed in each iteration
glm.rank <- glm.rank[,rownames(glm.rank)]
glm.rank[1:5,1:5]

# And finally some plots...
boxplot(glm.rank,horizontal=TRUE,lwd=2,ylab='',xlab='Prediction accuracy',names=rownames(glm.rank),col=rainbow(nrow(glm.rank)))

}
