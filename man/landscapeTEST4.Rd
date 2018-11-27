\name{landscapeTEST4}
\alias{landscapeTEST4}
\docType{data}
\title{
  A landscape for testing
}
\description{
  A landscape SpatialPolygons for testing, composed of 153 fields.
}
\usage{landscapeTEST4}
\format{
  The format is formal class 'SpatialPolygons' [package "sp"]
}
\details{
  The landscape structure is simulated using a T-tessellation algorithm (Kiêu et al. 2013) in order to control specific features such as number, area and shape of the fields.
}
\references{
   Kiêu K, Adamczyk-Chauvat K, Monod H, Stoica RS. A completely random T-tessellation model and Gibbsian extensions. Spat. Stat. 2013;6:118-38. doi: http://dx.doi.org/10.1016/j.spasta.2013.09.003.
}
\examples{
  library(sp)
  library(landsepi)
  landscapeTEST4
  plot(landscapeTEST4)
}
\keyword{datasets}
