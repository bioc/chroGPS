\name{rankFactorsbyDomain}
\alias{rankFactorsbyDomain}
\alias{rankFactorsbyDomain-methods}
\alias{rankFactorsbyDomain,data.frame-method}
\alias{rankFactorsbyDomain,distGPS-method}
\alias{rankFactorsbyDomain,matrix-method}

\title{
Function to help selecting candidate epigenetic factors based on
Epigenetic Domain cohesion/separation.
}
\description{
Function to help selecting candidate epigenetic factors based on
Epigenetic Domain cohesion/separation. This function ranks all combinations of a certain size of epigenetic factors in the selected domain based on how much they contribute to provide a good cohesion and separation of the points of their epigenetic domain against the rest. This allows selection of the 'best' combination of epigenetic factors in a certain domain to generate informative epigenetic maps.
}
\usage{
rankFactorsbyDomain(d, sampleinfo, ranktype = "domainDist", selName = "Color", selValue, k = NULL, mc.cores = 1)
}
\arguments{
  \item{d}{
A distGPS object with the epigenetic distances, see function distGPS.
}
  \item{sampleinfo}{
A data.frame object with at least the domain information for the epigenetic factors in the provided distance matrix. Sampleinfo rownames need to match row and colnames in the distance matrix.
}
  \item{ranktype}{
Character vector indicating the method to use. Currently only supporting 'domainDist' method.
}
  \item{selName}{
Name of the column containing the epigenetic domain information in the sampleinfo dataframe.
}
  \item{selValue}{
Character indicating the name of the domain to be evaluated.
}
  \item{k}{
Size of the factor combination to be evaluated.
}
  \item{mc.cores}{
Cores to use in calls to parallel::mclapply.
}
}
\value{
Named list with each evaluated combination of factors and:
\item{inter }{Their resulting average inter-domain distance.}
\item{intra }{Their resulting average intra-domain distance.}
}
\author{
Oscar Reina.
}
\seealso{
\code{\link{distGPS}}, \code{\link{domainDist}}
}
\examples{

## Rank Factors by Domain, using intra/inter domain distance

data(s2)
data(toydists)
#d <- distGPS(s2,metric='avgdist',mc.cores=8) # Compute distances
rownames(s2names) <- s2names$ExperimentName

# Known domains
# Call rankFactorsbyDomain for HP1a repression domain, select a combination of 4
# factors
library(caTools)
rank.factors.4 <- rankFactorsbyDomain(d,s2names,ranktype='domainDist',selName='Color',selValue='lightblue',k=3,mc.cores=1) # Test HP1a repression
ddd <- as.data.frame(do.call(rbind,lapply(rank.factors.4,unlist)))
ddd <- ddd[order(ddd$intra,decreasing=FALSE),]
head(ddd)

}
