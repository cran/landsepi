% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods-LandsepiParams.R
\name{saveDeploymentStrategy}
\alias{saveDeploymentStrategy}
\title{Save landscape and deployment strategy}
\usage{
saveDeploymentStrategy(
  params,
  outputGPKG = "landsepi_landscape.gpkg",
  overwrite = FALSE
)
}
\arguments{
\item{params}{a LandsepiParams Object.}

\item{outputGPKG}{name of the GPKG output (default: "landsepi_landscape.gpkg") to be generated.}

\item{overwrite}{a boolean specifying if existing files can be overwritten (TRUE) or not
(FALSE, default).}
}
\value{
an updated LandsepiParams object.
}
\description{
Generates a GPKG file containing the landscape and all parameters of
the deployment strategy
}
\details{
The function generates a GPKG file in the simulation path.
The GPKG file contains all input parameters needed to restore the landscape (sf object)
and deployment strategy (croptypes, cultivars and genes).
}
\examples{
\dontrun{
## Initialisation
simul_params <- createSimulParams(outputDir = getwd())
## Time parameters
simul_params <- setTime(simul_params, Nyears = 10, nTSpY = 120)
## Landscape
simul_params <- setLandscape(simul_params, loadLandscape(1))
## Genes
gene1 <- loadGene(name = "MG 1", type = "majorGene")
gene2 <- loadGene(name = "MG 2", type = "majorGene")
genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
simul_params <- setGenes(simul_params, genes)
## Cultivars
cultivar1 <- loadCultivar(name = "Susceptible", type = "wheat")
cultivar2 <- loadCultivar(name = "Resistant1", type = "wheat")
cultivar3 <- loadCultivar(name = "Resistant2", type = "wheat")
cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
simul_params <- setCultivars(simul_params, cultivars)
## Allocate genes to cultivars
simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
## Allocate cultivars to croptypes
croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
, "Resistant crop 1"
, "Resistant crop 2"))
croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
simul_params <- setCroptypes(simul_params, croptypes)
## Allocate croptypes to landscape        
rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
rotation_period <- 0 ## same croptypes every years
prop <- c(1 / 3, 1 / 3, 1 / 3) ## croptypes proportions
aggreg <- 10 ## aggregated landscape
simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = rotation_period,
rotation_sequence = rotation_sequence,
rotation_realloc = FALSE, prop = prop, aggreg = aggreg)
## Save into a GPKG file
simul_params <- saveDeploymentStrategy(simul_params)
}
}
