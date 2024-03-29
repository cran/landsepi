% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Class-LandsepiParams.R
\docType{class}
\name{LandsepiParams}
\alias{LandsepiParams}
\alias{LandsepiParams-class}
\title{Class LandsepiParams}
\description{
Landsepi simulation parameters
}
\details{
An object of class LandsepiParams that can be created by calling \code{\link{createSimulParams}}
}
\section{Slots}{

\describe{
\item{\code{Landscape}}{a landscape as sf object. See \code{\link{loadLandscape}}}

\item{\code{Croptypes}}{a dataframe with three columns named 'croptypeID' for croptype index,
'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype.
See \code{\link{loadCroptypes}}, \code{\link{setCroptypes}} and \code{\link{allocateCroptypeCultivars}}}

\item{\code{Cultivars}}{a dataframe of parameters associated with each host genotype (i.e. cultivars, lines)
when cultivated in pure crops.See \code{\link{loadCultivar}} and \code{\link{setCultivars}}}

\item{\code{CultivarsGenes}}{a list containing, for each host genotype, the indices of carried resistance genes.
See \code{\link{allocateCultivarGenes}}}

\item{\code{Genes}}{a data.frame of parameters associated with each resistance gene and with the evolution of
each corresponding pathogenicity gene. See \code{\link{loadGene}} and \code{\link{setGenes}}}

\item{\code{Pathogen}}{a list of i. pathogen aggressiveness parameters on a susceptible host
for a pathogen genotype not adapted to resistance and ii. sexual reproduction parameters. See \code{\link{loadPathogen}} and \code{\link{setPathogen}}}

\item{\code{ReproSexProb}}{a vector of size TimeParam$nTSpY + 1 (end of season) of the probabilities for an infectious host to reproduce via sex rather
than via cloning at each step.}

\item{\code{PI0}}{initial probability for the first host (whose index is 0) to be infectious (i.e. state I)
at the beginning of the simulation. Must be between 0 and 1. See \code{\link{setInoculum}}}

\item{\code{DispHost}}{a vectorized matrix giving the probability of host dispersal
from any field of the landscape to any other field. See \code{\link{loadDispersalHost}} and \code{\link{setDispersalHost}}}

\item{\code{DispPathoClonal}}{a vectorized matrix giving the probability of pathogen dispersal (clonal propagules)
from any field of the landscape to any other field. See \code{\link{loadDispersalPathogen}} and \code{\link{setDispersalPathogen}}}

\item{\code{DispPathoSex}}{a vectorized matrix giving the probability of pathogen dispersal (sexual propagules)
from any field of the landscape to any other field. See \code{\link{loadDispersalPathogen}} and \code{\link{setDispersalPathogen}}}

\item{\code{Treatment}}{a list of parameters to simulate the effect of chemical treatments on the pathogen,
see \code{\link{loadTreatment}} and \code{\link{setTreatment}}}

\item{\code{OutputDir}}{the directory for simulation outputs}

\item{\code{OutputGPKG}}{the name of the output GPKG file containing parameters of the deployment strategy}

\item{\code{Outputs}}{a list of outputs parameters. See \code{\link{loadOutputs}} and \code{\link{setOutputs}}}

\item{\code{TimeParam}}{a list of time parameters. See \code{\link{setTime}}}

\item{\code{Seed}}{an integer used as seed value (for random number generator). See \code{\link{setTime}}}
}}

