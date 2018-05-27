\name{combineGenesMatrix}
\alias{combineGenesMatrix}
\title{
Combine two datasets with epigenetic factor profiles at gene level in
order to produce a differential chroGPS-genes map.
}
\description{
Combine two datasets with epigenetic factor profiles at gene level in
order to produce a differential chroGPS-genes map.
}
\usage{
combineGenesMatrix(x, y, label.x, label.y, minFactors = 10, minGenes = 1000)
}
\arguments{
  \item{x}{
Data frame or matrix with genes (rows) and epigenetic factors (columns),
and values of 0/1 to indicate binding of a given epigenetic factor over genes.
}
  \item{y}{
Same as x, containing the second dataset to compare. Column names have
to match between both datasets, with at least a minimum number of common
ones.
}
  \item{label.x}{
Name for the x dataset.
}
  \item{label.y}{
Name for the y dataset.
}
  \item{minFactors}{
Minimum number of common factors between both datasets.
}
  \item{minGenes}{
Minimum number of genes in each dataset with epigenetic information
(that is, minimum number of rows where all values are not zeros).
}
}
\value{
A matrix with the combined epigenetic profile dataset to perform
differential analysis.
}
\author{
Oscar Reina.
}
\seealso{
See also \code{distGPS, mds, diffGenes}.
}
\examples{
## Not run
## See example in diffGenes function.
}

