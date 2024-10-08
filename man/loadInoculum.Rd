% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods-LandsepiParams.R
\name{loadInoculum}
\alias{loadInoculum}
\title{Load Inoculum}
\usage{
loadInoculum(
  params,
  pI0_all = NULL,
  pI0_host = NULL,
  pI0_patho = NULL,
  pI0_poly = NULL,
  pI0_mat = NULL
)
}
\arguments{
\item{params}{a LandsepiParams object.}

\item{pI0_all}{a numeric indicating the (same) probability to infect a host for all
pathogen genotypes, all cultivars and in all polygons}

\item{pI0_host}{a vector of length Nhost indicating the probabilities to infect an host,
for each cultivar (for all pathogen genotypes and all polygons).}

\item{pI0_patho}{a vector of length Npatho indicating the probabilities to infect an host,
for each pathogen genotype (for all cultivars and all polygons).}

\item{pI0_poly}{a vector of length Npoly indicating the probabilities to infect an host,
for each polygon (for all pathogen genotypes and all cultivars).}

\item{pI0_mat}{a 3D array of dimensions (1:Nhost,1:Npatho,1:Npoly) indicating the
probability to infect an host, for each cultivar, pathogen genotype and polygon.
Note that \code{pI0_all}, \code{pI0_host}, \code{pI0_patho} and \code{pI0_poly}
are not accounted if \code{pI0_mat} is filled.}
}
\value{
a 3D array of dimensions (1:Nhost,1:Npatho,1:Npoly)
}
\description{
Loads an inoculum for the beginning of the simulation (t=0), with
controlled localisation (polygons), infected cultivars and pathogen genotypes.
Note that landscape, gene, cultivar and croptype parameters must be set before
loading the inoculum.
}
\details{
The different options enable different types of inoculum (localisation,
infected cultivars and pathogen genetic diversity,
see different options in Examples).\cr
Unless the array \code{pI0_mat} is filled, the probability for a host to be infected
at the beginning of the simulation is computed in every polygon (poly), cultivar (host)
and pathogen genotype (patho) with
\code{pI0[host, patho, poly] = pI0_all * pI0_patho[patho] * pI0_host[host] * pI0_poly[poly]}. \cr
Before loading the inoculum, one can use
\code{getMatrixGenePatho()}, \code{getMatrixCultivarPatho()} and \code{getMatrixCroptypePatho()}
to acknowledge which pathogen genotypes are adapted to which genes, cultivars and croptypes.\cr
Once \code{setInoculum()} is used, one can call \code{inoculumToMatrix()} to get the inoculum
as a 3D array (1:Nhost,1:Npatho,1:Npoly)\cr
}
\examples{
\dontrun{
simul_params <- createSimulParams()
simul_params <- setTime(simul_params, Nyears = 1, nTSpY = 80)
basic_patho_param <- loadPathogen(disease = "rust")
simul_params <- setPathogen(simul_params, patho_params = basic_patho_param)
simul_params <- setLandscape(simul_params, loadLandscape(id = 1))
simul_params <- setDispersalPathogen(simul_params, loadDispersalPathogen(id = 1)[[1]])
gene1 <- loadGene(name = "MG 1", type = "majorGene")
gene2 <- loadGene(name = "MG 2", type = "majorGene")
genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
simul_params <- setGenes(simul_params, genes)
cultivar1 <- loadCultivar(name = "Susceptible", type = "wheat")
cultivar2 <- loadCultivar(name = "Resistant", type = "wheat")
cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
simul_params <- setCultivars(simul_params, cultivars)
simul_params <- allocateCultivarGenes(simul_params, "Resistant", c("MG 1", "MG 2"))
croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Resistant crop"))
croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop", c("Resistant"))
simul_params <- setCroptypes(simul_params, croptypes)
simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = 0
, rotation_sequence = croptypes$croptypeID
, prop = c(1/2,1/2), aggreg = 1, graphic = FALSE)

#### Definition of the inoculum ####

### Scenario 1. Only the avirulent pathogen on the susceptible cultivar ###
# In this situation, the susceptible cultivar must be entered
# at the first line of the table cultivars

## Global inoculum (i.e. in the whole landscape)
# Option 1: simply use the default parameterisation
simul_params <- setInoculum(simul_params, 5E-4)

# Option 2: use loadInoculum()
Npatho <- prod(simul_params@Genes$Nlevels_aggressiveness)
Nhost <- nrow(simul_params@Cultivars)
pI0 <- loadInoculum(simul_params,
                    pI0_all=5E-4,
                    pI0_host=c(1,rep(0, Nhost-1)),
                    pI0_patho=c(1,rep(0, Npatho-1)))
simul_params <- setInoculum(simul_params, pI0)
inoculumToMatrix(simul_params)

## Local inoculum (i.e. in some random polygons only)
Npatho <- prod(simul_params@Genes$Nlevels_aggressiveness)
Nhost <- nrow(simul_params@Cultivars)
Npoly <- nrow(simul_params@Landscape)
Npoly_inoc <- 5  ## number of inoculated polygons
## whether the avr pathogen can infect the polygons
compatible_poly <- getMatrixPolyPatho(simul_params)[,1]
## random polygon picked among compatible ones
id_poly <- sample(grep(1, compatible_poly), Npoly_inoc)
pI0_poly <- as.numeric(1:Npoly \%in\% id_poly)  
pI0 <- loadInoculum(simul_params,
                    pI0_all=5E-4,
                    pI0_host=c(1,rep(0, Nhost-1)),
                    pI0_patho=c(1,rep(0, Npatho-1)), 
pI0_poly=pI0_poly)
simul_params <- setInoculum(simul_params, pI0)
inoculumToMatrix(simul_params)

### Scenario 2. Diversity of pathogen genotypes in the inoculum ###
# in this example,  Nhost=2 cultivars, Npatho=4

## Global inoculum (i.e. in all polygons of the landscape)
pI0 <- loadInoculum(simul_params, pI0_patho=c(1E-3,1E-4,1E-4,1E-5), pI0_host=c(1,1))
simul_params <- setInoculum(simul_params, pI0)
inoculumToMatrix(simul_params)[,,1:5]

## Local inoculum (i.e. in some polygons only) ##
Npoly <- nrow(simul_params@Landscape)
Npoly_inoc <- 5  ## number of inoculated polygons 
id_poly <- sample(1:Npoly, Npoly_inoc)  ## random polygon 
pI0_poly <- as.numeric(1:Npoly \%in\% id_poly) 
pI0 <- loadInoculum(simul_params, pI0_patho=c(1E-3,1E-4,1E-4,1E-5),
pI0_host=c(1,1), pI0_poly=pI0_poly)
simul_params <- setInoculum(simul_params, pI0)
inoculumToMatrix(simul_params)
}
}
\seealso{
\link{inoculumToMatrix}, \link{getMatrixGenePatho}, \link{getMatrixCultivarPatho},
\link{getMatrixCroptypePatho}, \link{setInoculum}
}
