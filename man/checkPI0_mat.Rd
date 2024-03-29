% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods-LandsepiParams.R
\name{checkPI0_mat}
\alias{checkPI0_mat}
\title{Check the array PI0_mat when entered manually in \code{loadInoculum()}.}
\usage{
checkPI0_mat(mat, params)
}
\arguments{
\item{mat}{a 3D array of dimensions (1:Nhost,1:Npatho,1:Npoly)}

\item{params}{a LandsepiParams object.}
}
\value{
the same array at mat, possibly corrected if incompatibility has been detected
}
\description{
Checks validity of the array.
}
