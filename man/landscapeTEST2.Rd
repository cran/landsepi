\name{landscapeTEST2}
\alias{landscapeTEST2}
\docType{data}
\title{
  A landscape for testing
}
\description{
  A landscape SpatialPolygons for testing, composed of 154 fields.
}
\usage{landscapeTEST2}
\format{
  The format is a formal class 'SpatialPolygons' [package "sp"] 
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
  landscapeTEST2
  plot(landscapeTEST2)
}
\keyword{datasets}
