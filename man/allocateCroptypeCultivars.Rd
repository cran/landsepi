% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods-LandsepiParams.R
\name{allocateCroptypeCultivars}
\alias{allocateCroptypeCultivars}
\title{Allocate cultivars to one croptype}
\usage{
allocateCroptypeCultivars(
  croptypes,
  croptypeName,
  cultivarsInCroptype,
  prop = NULL
)
}
\arguments{
\item{croptypes}{a dataframe containing all croptypes, initialised via
\code{\link{loadCroptypes}}}

\item{croptypeName}{the name of the croptype to be allocated}

\item{cultivarsInCroptype}{name of cultivars composing the croptype}

\item{prop}{vector of proportions of each cultivar in the croptype. Default to
balanced proportions.}
}
\value{
a croptype data.frame updated for the concerned croptype.
}
\description{
Updates a given croptype by allocating cultivars composing it.
}
\examples{
\dontrun{
simul_params <- createSimulParams()
cultivar1 <- loadCultivar(name = "Susceptible", type = "wheat")
cultivar2 <- loadCultivar(name = "Resistant1", type = "wheat")
cultivar3 <- loadCultivar(name = "Resistant2", type = "wheat")
cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
simul_params <- setCultivars(simul_params, cultivars)
croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Mixture"))
croptypes
croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Resistant1", "Resistant2"))
croptypes
}
}
\seealso{
\link{setCroptypes}, \link{setCultivars}
}
