\name{mergeReplicates}
\alias{mergeReplicates}
\alias{mergeReplicates-methods}
\alias{mergeReplicates,GRangesList-method}
\alias{mergeReplicates,data.frame-method}
\alias{mergeReplicates,matrix-method}
\alias{mergeReplicates,list-method}
\alias{mergeReplicates,mds-method}
\alias{mergeReplicateList}
\alias{mergeReplicateMatrix}
\alias{mergeReplicateMDS}
\title{Merges information from epigenetic replicates at factor, gene and chroGPS MDS map level.}
\description{
The function performs effectively merging of epigenetic replicate information, either as genomic intervals or as gene epigenetic profiles. It also allows merging of 
final chroGPS MDS factor maps so that multiple points for the same epigenetic factor (i.e. different antibodies or experimental sources) are returned as a single location in the map (centroid).
}
\section{Methods}{
\describe{
  
  \item{signature(x='GRangesList')}{ Each element in \code{x} is
  assumed to indicate the binding sites for a different sample,
  e.g. epigenetic factor. Typically \code{space(x)} indicates the
  chromosome, \code{start(x)} the start position and \code{end(x)} the
  end position (in bp). Strand information is ignored.  }

\item{signature(x='matrix')}{ Rows in \code{x} contain individuals for
  which we want to compute distances. Columns in \code{x} contain the
  variables, and should only contain either 0's and 1's or \code{FALSE}
  and \code{TRUE}.}

\item{signature(x='mds')}{A \code{MDS} object containing a chroGPS-factors map.}

}}
\usage{
mergeReplicates(x,id,mergeBy='any',mc.cores=1)
}
\arguments{
  \item{x}{Object containing the epigenetic information used to generate either factors or gene maps. Can also be a final chroGPS-factors MDS map object.}
  \item{id}{Identifier for epigenetic factors (list elements for x='GRangesList', column names for x='matrix' or 'data.frame', row names for x='mds'). This information will be used to identify valid replicates to be merged.}
  \item{mergeBy}{Either a character vector with possible values 'any' or 'all' to merge replicate genomic intervals or gene epigenetic profiles based on a logical 'OR' or 'AND' fashion, or a decimal number between (0,1] indicating the minimum proportion of replicates needed to consider.}
  \item{mc.cores}{Number of cores to use in calls to parallel::mclapply}
}
\value{
	An object from the same class as x, containing a data set with merged replicates as specified by the 'mergeBy' argument.
}
\seealso{ \code{\link{mds}} to create MDS-oriented objects }
\examples{

## Not run

#data(s2)
#data(bg3)

#names(s2)
#names(bg3)

# Unify replicates
#mnames <- sort(unique(intersect(s2names$Factor,bg3names$Factor)))
#sel <- s2names$Factor %in% mnames 
#s2.repset <- mergeReplicateList(s2[sel],id=s2names$Factor[sel],mergeBy='any')
#sel <- bg3names$Factor %in% mnames
#bg3.repset <- mergeReplicateList(bg3[sel],id=bg3names$Factor[sel],mergeBy='any')

#names(s2.repset)
#names(bg3.repset)

}
\keyword{multivariate}
\keyword{clustering}
