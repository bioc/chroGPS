\name{distGPS-class}
\Rdversion{1.1}
\docType{class}
\alias{distGPS-class}
\alias{show,distGPS-method}
\alias{as.matrix,distGPS-method}
\alias{splitDistGPS-class}

\title{Class \code{"distGPS"}}
\description{
  Pairwise distances between elements. Function \code{distGPS} creates
  objects of this class.
  \code{splitDistGPS} in an private class used internally for parallel
  Multidimensional Scaling.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("distGPS", ...)} to
generate chroGPS-compliant distance matrices with user-defined
metrics. The internal distance matrix can be extracted with the function
\code{"as.matrix"} }
\section{Slots}{
  \describe{
    \item{\code{d}:}{Object of class \code{"matrix"} with pairwise
      dissimilarities (distances) between elements. }
    \item{\code{metric}:}{Object of class \code{"character"} indicating
      the metric type used for calculating distances. See function \code{"distGPS"}. }
    \item{\code{type}:}{Object of class \code{"character"}, deprecated. }
  }
}
\author{
  Oscar Reina
}
\examples{
showClass("distGPS")
data(s2)
data(toydists)
d
class(as.matrix(d))
as.matrix(d)[1:5,1:5]
}
\keyword{classes}
