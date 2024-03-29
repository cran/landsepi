% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods-LandsepiParams.R
\name{createSimulParams}
\alias{createSimulParams}
\title{Create a LandsepiParams object.}
\usage{
createSimulParams(outputDir = "./")
}
\arguments{
\item{outputDir}{ouput directory for simulation (default: current directory)}
}
\value{
a LandsepiParams object initialised with the following context:
\itemize{
\item random seed
\item all pathogen parameters fixed at 0
\item no between-polygon dispersal (neither pathogen nor host)
\item no pathogen introduction
\item no resistance gene
\item no chemical treatment
\item no output to generate.
}
}
\description{
Creates a default object of class LandsepiParams.
}
\details{
Create a default object of class LandsepiParams used to store all
simulation parameters. It also creates a subdirectory in \code{outputDir}
using the date; this directory will contain all simulation outputs.
}
\examples{
\dontrun{
createSimulParams()
}
}
