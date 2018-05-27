\name{bg3}
\alias{bg3}
\alias{bg3names}
\alias{bg3.tab}
\docType{data}
\title{
  Sample binding site and related data from BG3 cell line in Drosophila melanogaster.
}
\description{
  chroGPS example dataset including ChIP-CHIP (modEncode) and 
  ChIP-Seq (NCBI GEO GSE19325) data for Drosophila melanogaster BG3 cell
  line. The object \code{toydists} stores precomputed
  \code{distGPS} objects (called \code{d, d2, d3}) for the epigenetic factors used in the dynamic
  vignette that comes with the package.
}
\usage{data(bg3)}
\source{
  \code{http://www.modencode.org}
  \code{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19325}
}
\references{
  \code{http://www.modencode.org}
  \code{http://www.ncbi.nlm.nih.gov/geo/}
}
\seealso{
  \code{s2} for the analog S2 cell line object. \code{toydists} for
  precomputed distance matrixes.
}
\examples{
data(bg3)
class(bg3)
bg3
bg3names$Factor
# See vignette examples for several uses of these datasets.
}
\keyword{datasets}
\keyword{modEncode}
