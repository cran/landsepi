\name{dispP_2}
\alias{dispP_2}
\docType{data}
\title{
  Pathogen dispersal data for testing.
}
\description{
  Five vectorised dispersal matrices of the pathogen, associated with landscape 2, composed of 154 fields.
}
\usage{dispP_2}
\format{
  The format is:
 num [1:23716] 8.70e-01 5.87e-08 3.35e-07 8.82e-10 2.55e-08 ...
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
dispP_2
summary(dispP_2)
## maybe str(dispP_2) ; plot(dispP_2) ...
}
\keyword{datasets}
