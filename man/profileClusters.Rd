\name{profileClusters}
\alias{profileClusters}
\title{
Assess epigenetic profiles for genes present in each
cluster as obtained by the \code{clusGPS} function.
}
\description{
Assess epigenetic profiles for genes present in each cluster as obtained by the \code{clusGPS} function.
}
\usage{
profileClusters(x, clus, clusName = NULL, normalize = FALSE, mc.cores = 1, ...)
}
\arguments{
  \item{x}{
  Matrix or data.frame with epigenetic profiles used to generate chroGPS-genes map. See \code{distGPS}.
}
  \item{clus}{
  A valid \code{clusGPS} object from the chroGPS-genes map generated from x.
}
  \item{clusName}{
  If desired, the name of a clustering solution within the \code{clusGPS} object. See \code{clusNames}.
}
  \item{normalize}{
  Normalize epigenetic profile information with respect to whole genome (to assess enrichments/depletions).
}
  \item{mc.cores}{
  Cores to use in call to parallel::mclapply.
}
  \item{\dots}{
  Additional arguments.
}
}
\value{
 Data frame with epigenetic enrichment profiles for all clusters (rows) and epigenetic factors (columns), to be used for downstream assessments or visualization (heatmaps).
}
\author{
Oscar Reina.
}
\seealso{
\code{distGPS}, \code{clusGPS}.
}
\examples{
## Not run
}
