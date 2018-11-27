\name{dispP_5}
\alias{dispP_5}
\docType{data}
\title{
  Pathogen dispersal data for testing.
}
\description{
  Five vectorised dispersal matrices of the pathogen, associated with landscape 5, composed of 156 fields.
}
\usage{dispP_5}
\format{
  The format is:
 num [1:24336] 8.82e-01 9.04e-10 2.32e-10 6.44e-07 3.48e-09 ...
}
\details{
  The pathogen dispersal matrix gives the probability for a pathogen in a field i (row) to migrate to field i' (column) through dispersal. It is computed based on a dispersal kernel and the euclidian distance between each point in fields i and i', using the CaliFloPP algorithm (Bouvier et al. 2009).
  In this example, the dispersal kernel is an isotropic power-law function of equation:
  f(x)=((b-2)*(b-1)/(2*pi*a^2)) * (1 + x/a)^(-b)
  with a=40 a scale parameter and b=7 related to the weight of the dispersal tail. The expected mean dispersal distance is given by 2*a/(b-3)=20 m.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\references{
  Bouvier A, Kiêu K, Adamczyk K, Monod H. Computation of the integrated flow of particles between polygons. Environ. Model Softw. 2009;24(7):843-9. doi: http://dx.doi.org/10.1016/j.envsoft.2008.11.006.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
dispP_5
summary(dispP_5)
## maybe str(dispP_5) ; plot(dispP_5) ...
}
\keyword{datasets}
