# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
#                    Julien Papaix <julien.papaix@inrae.fr>
#                    Jean-François Rey <jean-francois.rey@inrae.fr>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,i
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


# Default data.frame column names
.croptypesColNames <- c("croptypeID", "croptypeName")
.cultivarsColNames <- c("cultivarName", "initial_density", "max_density", "growth_rate"
                        , "reproduction_rate", "yield_H", "yield_L", "yield_I"
                        , "yield_R", "planting_cost", "market_value")
.cultivarsGenesColNames <- c()
.geneColNames <- c("geneName", "efficiency", "age_of_activ_mean", "age_of_activ_var"
                   , "mutation_prob", "Nlevels_aggressiveness", "adaptation_cost"
                   , "relative_advantage"
                   , "tradeoff_strength", "target_trait", "recombination_sd")


#' @title LandsepiParams
#' @description Creates and initialises a LandsepiParams object with default parameters.
#' @param .Object a LandsepiParam object. 
#' @param Landscape a landscape as sf object.
#' @param Croptypes a dataframe with three columns named 'croptypeID' for croptype index,
#' 'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar 
#' within the croptype.
#' @param Cultivars a dataframe of parameters associated with each host genotype 
#' (i.e. cultivars, lines) when cultivated in pure crops.
#' @param CultivarsGenes a list containing, for each host genotype, the indices of 
#' carried resistance genes.
#' @param Genes a data.frame of parameters associated with each resistance gene and with 
#' the evolution of each corresponding pathogenicity gene.
#' @param Pathogen a list of pathogen aggressiveness parameters on a susceptible host
#' for a pathogen genotype not adapted to resistance.
#' @param PI0 vector of length Npoly.Nhost.Npatho filled with the initial probabilities for hosts 
#' to be infectious (i.e. state I), for each pathogen genotype,
#' at the beginning of the simulation.
#' @param DispHost a vectorized matrix giving the probability of host dispersal
#' from any polygon of the landscape to any other polygon
#' @param DispPathoClonal a vectorized matrix giving the probability of pathogen dispersal
#' from any polygon of the landscape to any other polygon.
#' @param DispPathoSex a vectorized matrix giving the probability of pathogen dispersal
#' from any polygon of the landscape to any other polygon (sexual propagule).
#' @param Treatment a list of chemical treatment parameters (indices of treated cultivars, 
#' times of application, efficiency and degradation rate)
#' @param OutputDir the directory for simulation outputs 
#' @param OutputGPKG the name of the output GPKG file containing parameters of the 
#' deployment strategy
#' @param Outputs a list of outputs parameters.
#' @param TimeParam a list of time parameters.
#' @param Seed an integer used as seed value (for random number generator).
#' @param ... more options
#' @rdname initialize-methods
# @aliases LandsepiParams-method
#' @include Class-LandsepiParams.R
#' @include GPKGTools.R tools.R
setMethod(
  "initialize", "LandsepiParams",
  function(.Object,
           Landscape = st_sf(st_sfc()),
           Croptypes = data.frame(),
           Cultivars = data.frame(matrix(ncol = length(.cultivarsColNames)
                                         , nrow = 0, dimnames = list(NULL, .cultivarsColNames))),
           CultivarsGenes = data.frame(),
           Genes = data.frame(matrix(ncol = length(.geneColNames)
                                     , nrow = 0, dimnames = list(NULL, .geneColNames))),
           Pathogen = list(
             name = "no pathogen",
             survival_prob = 0,
             repro_sex_prob = 0,
             infection_rate = 0,
             propagule_prod_rate = 0,
             latent_period_mean = 0,
             latent_period_var = 0,
             infectious_period_mean = 0,
             infectious_period_var = 0,
             sigmoid_kappa = 0,
             sigmoid_sigma = 0,
             sigmoid_plateau = 1,
             sex_propagule_viability_limit  = 0,
             sex_propagule_release_mean = 0,
             clonal_propagule_gradual_release = 0
           ),
           # ReproSexProb = vector(),
           PI0 = 0, #vector(),
           DispHost = vector(),
           DispPathoClonal = vector(),
           DispPathoSex = vector(),
           Treatment = list(
             treatment_degradation_rate = 0.1,
             treatment_efficiency = 0,
             treatment_timesteps = vector(),
             treatment_cultivars = vector(),
             treatment_cost = 0,
             treatment_application_threshold = vector()
           ),
           OutputDir = normalizePath(character(getwd())),
           OutputGPKG = "landsepi_landscape.gpkg",
           Outputs = list(epid_outputs = "", evol_outputs = ""
                          , thres_breakdown = NA, audpc100S = NA), 
           TimeParam = list(Nyears = 0, nTSpY = 0),
           Seed = NULL,
           ...) {
    # .Object <- callNextMethod(...)
    .Object@Landscape <- Landscape
    .Object@Croptypes <- Croptypes
    .Object@Cultivars <- Cultivars
    .Object@CultivarsGenes <- CultivarsGenes
    .Object@Genes <- Genes
    .Object@Pathogen <- Pathogen
    # .Object@ReproSexProb <- ReproSexProb
    .Object@PI0 <- PI0
    .Object@DispHost <- DispHost
    .Object@DispPathoClonal <- DispPathoClonal
    .Object@DispPathoSex <- DispPathoSex
    .Object@Treatment <- Treatment
    .Object@OutputDir <- OutputDir
    .Object@OutputGPKG <- OutputGPKG
    .Object@Outputs <- Outputs
    .Object@TimeParam <- TimeParam
    .Object@Seed <- Seed

    validObject(.Object)
    .Object
  }
)

#' @name print
#' @title print
#' @description Prints a LandsepiParams object.
#' @param x a LandsepiParams object
#' @param ... print options
#' @rdname print-methods
#' @aliases print,LandsepiParams-method
#' @export
setMethod("print", "LandsepiParams", function(x, ...) {
  print("## LandsepiParams values :")
  print("### Landscape")
  print(x@Landscape)
  if (nrow(x@Landscape) != 0) {
    plot(st_geometry(x@Landscape))
  } else {
    print("Nothings to plot")
  }
  print("### Croptypes")
  print(x@Croptypes)
  print("### Cultivars")
  print(x@Cultivars)
  print("### CultivarsGenes")
  print(x@CultivarsGenes)
  print("### Genes")
  print(x@Genes)
  print("### Pathogen")
  print(x@Pathogen)
  # print("### Pathogen ReproSexProb")
  # print(x@ReproSexProb)

  print("### Treatment")
  print(x@Treatment)
  
  print("### Inoculum : ")
  print(inoculumToMatrix(x))
  print("### Nyears : ")
  print(x@TimeParam$Nyears)
  print("### nTSpY Number of step by year : ")
  print(x@TimeParam$nTSpY)
  print("### Output Directory : ")
  print(x@OutputDir)
  print("### Output GPKG : ")
  print(x@OutputGPKG)
  print("### Seed : ")
  print(x@Seed)
  print("### Outputs")
  print(x@Outputs)
})


#' @name summary
#' @title summary
#' @description Prints the summary of a LandsepiParams object.
#' @param object a LandsepiParams object.
#' @rdname summary-methods
#' @aliases summary,LandsepiParams-method
#' @export
setMethod("summary", "LandsepiParams", function(object) {
  message("## LandsepiParam Object slots:\n")

  message("### Landscape : ")
  if (nrow(object@Landscape) == 0) {
    message("\tnot set (see setLandscape method)")
  } else {
    summary(object@Landscape)
  }

  message("### Croptypes (proportions of Cultivars in each croptype) : ")
  if (nrow(object@Croptypes) == 0) {
    message("\tnot set (see setCroptypes method)")
  } else {
    summary(object@Croptypes)
  }

  message("### Cultivars (cultivars parameters) : ")
  if (nrow(object@Croptypes) == 0) {
    message("\tnot set (see setCultivars method)")
  } else {
    summary(object@Cultivars)
  }

  message("### CultivarsGenes (List of Genes by Cultivars) : ")
  if (nrow(object@Croptypes) == 0) {
    message("\tnot set (see setCultivarsGenes method)")
  } else {
    summary(object@CultivarsGenes)
  }

  message("### Genes (Genes parameters) : ")
  if (nrow(object@Genes) == 0) {
    message("\tnot set (see setGenes method)")
  } else {
    summary(object@Genes)
  }

  message("### Pathogen parameters : ")
  if (length(object@Pathogen) == 0) {
    message("\tnot set (see loadPathogen and setPathogen methods)")
  } else {
    summary(object@Pathogen)
  }
  # message("### Reproduction Sex Probabilities :")
  # summary(object@ReproSexProb)

  message("### Pathogen Dispersal Matrix (as vector) : ")
  if (length(object@DispPathoClonal) == 0) {
    message("\tnot set (see loadDispersalPathogen and setDispersalPathogen methods)")
  } else {
    summary(object@DispPathoClonal)
  }
  
  message("### Pathogen Dispersal Repro Sex Matrix (as vector) : ")
  if (length(object@DispPathoSex) == 0) {
    message("\tnot set (see loadDispersalPathogen and setDispersalPathogen methods)")
  } else {
    summary(object@DispPathoSex)
  }

  message("### Host Dispersal Matrix (as vector) : ")
  if (length(object@DispHost) == 0) {
    message("\tnot set (see loadDispersalHost and setDispersalHost methods)")
  } else {
    summary(object@DispHost)
  }

  message("### Inoculum : ", inoculumToMatrix(object))
  
  message("### Treatment")
  summary(object@Treatment)
  
  message("### Nyears : ", object@TimeParam$Nyears)
  message("### nTSpY Number of step by year : ", object@TimeParam$nTSpY)
  message("### Output Directory : ", object@OutputDir)
  message("### Output GPKG : ", object@OutputGPKG)
  message("### Seed : ", object@Seed)
  message("### Outputs")
  message(object@Outputs)
})


#' @name show
#' @title show
#' @description Shows a LandsepiParams object.
#' @param object a LandsepiParams object
#' @rdname show-methods
#' @aliases show,LandsepiParams-method
#' @export
setMethod("show", "LandsepiParams", function(object) {
  print(object)
})



#' @name checkSimulParams
#' @title Check simulation parameters
#' @description Checks validity of a LandsepiParams object.
#' @param params a LandsepiParams Object.
#' @return TRUE if OK for simulation, FALSE otherwise
#' @export
checkSimulParams <- function(params) {
  validity <- TRUE
  validity <- validity && checkLandscape(params)
  validity <- validity && checkCultivars(params)
  validity <- validity && checkCroptypes(params)
  validity <- validity && checkGenes(params)
  validity <- validity && checkCultivarsGenes(params)
  validity <- validity && checkPathogen(params)
  validity <- validity && checkTreatment(params)
  validity <- validity && checkDispersalHost(params)
  validity <- validity && checkDispersalPathogen(params)
  validity <- validity && checkInoculum(params)
  #validity <- validity && checkInoculumLandscape(params)
  validity <- validity && checkTime(params)
  validity <- validity && checkOutputs(params)
  
  return(validity)
}



#' @name createSimulParams
#' @title Create a LandsepiParams object.
#' @description Creates a default object of class LandsepiParams.
#' @param outputDir ouput directory for simulation (default: current directory)
#' @details Create a default object of class LandsepiParams used to store all 
#' simulation parameters. It also creates a subdirectory in \code{outputDir} 
#' using the date; this directory will contain all simulation outputs.
#' @return a LandsepiParams object initialised with the following context:
#' \itemize{
#' \item random seed
#' \item all pathogen parameters fixed at 0
#' \item no between-polygon dispersal (neither pathogen nor host)
#' \item no pathogen introduction
#' \item no resistance gene
#' \item no chemical treatment
#' \item no output to generate.
#' }
#' @examples \dontrun{
#' createSimulParams()
#' }
#' @export
createSimulParams <- function(outputDir = "./") {

  ## Avoid subdirectory creation
  if (length(grep("simul_landsepi_", normalizePath(outputDir))) != 0) {
    outputDir <- dirname(normalizePath(outputDir))
  }

  ## create a subdirectory with time
  timeSimul <- paste(strsplit(as.character(Sys.time()), " ")[[1]], collapse = "_")
  nameDir <- paste(outputDir, "/simul_landsepi_", gsub(":", "-", timeSimul), sep = "")
  dir.create(nameDir)

  message("Created output directory : ", normalizePath(nameDir))

  lp <- new("LandsepiParams",
    OutputDir = normalizePath(nameDir),
    Seed = setSeedValue()
  )

  return(lp)
}


#' @name loadSimulParams
#' @title Load simulation parameters
#' @description Loads a GPKG file from the output of a landsepi simulation.
#' @details See \code{\link{saveDeploymentStrategy}}.
#' @param inputGPKG name of the GPKG file.
#' @return a LandsepiParams object.
#' @export
loadSimulParams <- function(inputGPKG = "") {
  lp <- new("LandsepiParams",
    OutputDir = normalizePath(dirname(inputGPKG)),
    OutputGPKG = basename(inputGPKG)
    # Seed = setSeedValue(seed),
    # TimeParam = list(Nyears = Nyears, nTSpY = nTSpY)
  )

  lp <- setLandscape(lp, st_read(dsn = inputGPKG, layer = "croptypeID"))
  lp <- setCroptypes(lp, CroptypeBDD2Params(inputGPKG))
  lp <- setGenes(lp, GeneBDD2Params(inputGPKG))
  lp <- setCultivars(lp, CultivarBDD2Params(inputGPKG))
  lp@CultivarsGenes <- CultivarGeneBDD2Params(inputGPKG)

  ## TODO get all parameters from GPKG and parameters.txt if exist
  ## TODO: doesn't seem to work with croptypes, cultivars and cultivarGenes

  message("Not yet implemented for DispHost, DispPathoClonal, PI0, pathogen, time, outputs...")

  return(lp)
}


#' @name saveDeploymentStrategy
#' @title Save landscape and deployment strategy 
#' @description Generates a GPKG file containing the landscape and all parameters of 
#' the deployment strategy
#' @details The function generates a GPKG file in the simulation path. 
#'  The GPKG file contains all input parameters needed to restore the landscape (sf object) 
#'  and deployment strategy (croptypes, cultivars and genes).
#' @param params a LandsepiParams Object.
#' @param outputGPKG name of the GPKG output (default: "landsepi_landscape.gpkg") to be generated.
#' @param overwrite a boolean specifying if existing files can be overwritten (TRUE) or not 
#' (FALSE, default).
#' @return an updated LandsepiParams object.
#' @examples
#' \dontrun{
#' ## Initialisation
#' simul_params <- createSimulParams(outputDir = getwd())
#' ## Time parameters
#' simul_params <- setTime(simul_params, Nyears = 10, nTSpY = 120)
#' ## Landscape
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' ## Genes
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' ## Cultivars
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' ## Allocate genes to cultivars
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
#' ## Allocate cultivars to croptypes
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
#' , "Resistant crop 1"
#' , "Resistant crop 2"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
#' simul_params <- setCroptypes(simul_params, croptypes)
#' ## Allocate croptypes to landscape        
#' rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
#' rotation_period <- 0 ## same croptypes every years
#' prop <- c(1 / 3, 1 / 3, 1 / 3) ## croptypes proportions
#' aggreg <- 10 ## aggregated landscape
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = rotation_period,
#' rotation_sequence = rotation_sequence,
#' rotation_realloc = FALSE, prop = prop, aggreg = aggreg)
#' ## Save into a GPKG file
#' simul_params <- saveDeploymentStrategy(simul_params)
#' }
#' @export
saveDeploymentStrategy <- function(params, outputGPKG = "landsepi_landscape.gpkg"
                                   , overwrite = FALSE) {
  
  params@OutputGPKG <- outputGPKG
  
  if (!dir.exists(params@OutputDir)) {
    warning("Directory ", params@OutputDir, " does not exist")
    return(params)
  }

  # setwd(params@OutputDir)
  # message("Move to ", params@OutputDir, " directory for simulation")

  if (file.exists(paste0(params@OutputDir, "/", params@OutputGPKG))) {
    if (overwrite == FALSE) {
      warning(params@OutputGPKG, " already exists, can't overwrite it.")
      warning("use overwrite = TRUE to allow replacement of existing files")
      return(params)
    }
    else {
      message("Will overwrite existing files in ", params@OutputDir)
      file.remove(paste0(params@OutputDir, "/", params@OutputGPKG))
    }
  }

  # try to add one more polygons (year_Nyears+1), if missing
  if (length(grep("^year_", colnames(params@Landscape))) == params@TimeParam$Nyears) {
    params@Landscape[, paste0("year_", params@TimeParam$Nyears + 1)] <- as.data.frame(
      params@Landscape[, paste0("year_", params@TimeParam$Nyears)])[, 1]
    message("Add one more year of simulation (only for simulation model constraints)")
  }
  ## create gpkg file
  message("Create ", paste0(params@OutputDir, "/", params@OutputGPKG), " file")
  ## save only years_ and area ?
  land <- params@Landscape[, grep("^year_", colnames(params@Landscape))]
  gpkgFile <- createLandscapeGPKG(land, paste0(params@OutputDir, "/", params@OutputGPKG))
  ## add data tables
  GPKGAddTables(gpkgFile)
  ## Fill data tables
  if (nrow(params@Cultivars) > 0) GPKGAddInputData(gpkgFile, table = "Cultivar"
                                                   , data = params2CultivarBDD(params)
                                                   , deleteExistingData = TRUE)
  if (nrow(params@Croptypes) > 0) GPKGAddInputData(gpkgFile, table = "CultivarList"
                                                   , data = params2CroptypeBDD(params)
                                                   , deleteExistingData = TRUE)
  if (nrow(params@Genes) > 0) GPKGAddInputData(gpkgFile, table = "Gene"
                                               , data = params2GeneBDD(params)
                                               , deleteExistingData = TRUE)
  if (nrow(params@CultivarsGenes) > 0) GPKGAddInputData(gpkgFile, table = "GeneList"
                                                        , data = params2GeneListBDD(params)
                                                        , deleteExistingData = TRUE)

  return(params)
}


#' @name runSimul
#' @title Run a simulation
#' @description Runs a simulation with landsepi, 
#' a stochastic, spatially-explicit, demo-genetic model simulating the spread and evolution
#' of a pathogen in a heterogeneous landscape and generating a wide range of epidemiological, 
#' evolutionary and economic outputs.
#' @param params a LandsepiParams Object containing all simulation parameters. Must be initialised 
#' with \code{\link{createSimulParams}} and updated using \code{set*()} methods 
#' (see vignettes for details).
#' @param graphic a logical indicating if graphics must be generated (TRUE, default) 
#' or not (FALSE).
#' @param writeTXT a logical indicating if outputs must be written in text files (TRUE, default) 
#' or not (FALSE).
#' @param videoMP4 a logical indicating if a video must be generated (TRUE) or not (FALSE, default).
#' Works only if graphic=TRUE and audpc_rel is computed.
#' @param keepRawResults a logical indicating if binary files must be kept after the end of 
#' the simulation (default=FALSE). Careful, many files may be generated if keepRawResults=TRUE.
#' @details See \code{?landsepi} for details on the model, assumptions and outputs, and our 
#' vignettes for tutorials (\code{browseVignettes("landsepi")}). The function runs the model 
#' simulation using a LandsepiParams object.
#' Briefly, the model is stochastic, spatially explicit (the basic spatial unit is an 
#' individual field or polygon), based on a SEIR (‘susceptible-exposed-infectious-removed’, 
#' renamed HLIR for 'healthy-latent-infectious-removed' to avoid confusions with 'susceptible host') 
#' structure with a discrete time step. It simulates the spread and
#'  evolution (via mutation, recombination through sexual reproduction, selection and drift) 
#'  of a pathogen in a heterogeneous cropping landscape, across cropping seasons split 
#'  by host harvests which impose potential bottlenecks to the pathogen. A wide array of 
#'  resistance deployment strategies (possibly including chemical treatments) 
#'  can be simulated and evaluated using several possible 
#'  outputs to assess the epidemiological, evolutionary and economic performance
#'  of deployment strategies.
#' 
#' @return A list containing all required outputs.
#' A set of text files, graphics and a video showing epidemic dynamics can be generated.
#' If keepRawResults=TRUE, a set of binary files is generated for every year of simulation and
#' every compartment: \itemize{
#'  \item H: healthy hosts,
#'  \item Hjuv: juvenile healthy hosts (for host reproduction),
#'  \item L: latently infected hosts,
#'  \item I: infectious hosts,
#'  \item R: removed hosts,
#'  \item P: propagules.}
#' Each file indicates for every time step the number of individuals in each polygon, and when 
#' appropriate for each host and pathogen genotype. Additionally, a binary file called TFI is 
#' generated and gives the Treatment Frequency Indicator (expressed as the number of treatment 
#' applications per polygon).
#' @seealso \link{demo_landsepi}
#' @examples \dontrun{
#' ### Here is an example of simulation of a mosaic of three cultivars (S + R1 + R2). See our 
#' ## tutorials for more examples.
#' ## Initialisation
#' simul_params <- createSimulParams(outputDir = getwd())
#' ## Seed & Time parameters
#' simul_params <- setSeed(simul_params, seed = 1)
#' simul_params <- setTime(simul_params, Nyears = 10, nTSpY = 120)
#' ## Pathogen parameters
#' simul_params <- setPathogen(simul_params, loadPathogen("rust"))
#' ## Landscape & dispersal
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' simul_params <- setDispersalPathogen(simul_params, loadDispersalPathogen[[1]])
#' ## Genes
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' ## Cultivars
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' ## Allocate genes to cultivars
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
#' ## Allocate cultivars to croptypes
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
#' , "Resistant crop 1"
#' , "Resistant crop 2"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
#' simul_params <- setCroptypes(simul_params, croptypes)
#' ## Allocate croptypes to landscape        
#' rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
#' rotation_period <- 0 ## same croptypes every years
#' prop <- c(1 / 3, 1 / 3, 1 / 3) ## croptypes proportions
#' aggreg <- 10 ## aggregated landscape
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = rotation_period,
#' rotation_sequence = rotation_sequence,
#' rotation_realloc = FALSE, prop = prop, aggreg = aggreg)
#' Set the inoculum
#' simul_params <- setInoculum(simul_params, 5e-4)
#' ## list of outputs to be generated
#' simul_params <- setOutputs(simul_params, loadOutputs())
#' ## Check simulation parameters
#' checkSimulParams(simul_params)
#' ## Save deployment strategy into GPKG file
#' simul_params <- saveDeploymentStrategy(simul_params)
#' ## Run simulation
#' runSimul(simul_params)
#' 
#' ### Simulation of rust epidemics in a 1-km^2 patch cultivated with a susceptible wheat cultivar
#'seed=10
#'Nyears=5
#'disease="rust"
#'hostType="growingHost"
#'simul_params <- createSimulParams(outputDir = getwd())
#'
#'## Seed and time parameters
#'simul_params <- setSeed(simul_params, seed)
#'simul_params <- setTime(simul_params, Nyears, nTSpY=120)
#'
## Pathogen parameters
#'simul_params <- setPathogen(simul_params, loadPathogen(disease))
#'myLand <- Polygons(list(Polygon(matrix(c(0,0,1,1,0,1,1,0)*1000, nrow=4))), "ID1")
#'myLand <- SpatialPolygons(list(myLand))
#'simul_params <- setLandscape(simul_params, myLand)
#'
#'## Simulation, pathogen, landscape and dispersal parameters
#'simul_params <- setDispersalPathogen(simul_params, c(1))
#'
#'## Cultivars
#'simul_params <- setCultivars(simul_params, loadCultivar(name = "Susceptible", type = hostType))
#'
#'## Croptypes
#'croptype <- data.frame(croptypeID = 0, croptypeName = c("Fully susceptible crop")
#', Susceptible = 1)
#'simul_params <- setCroptypes(simul_params, croptype)
#'simul_params <- allocateLandscapeCroptypes(simul_params,
#'rotation_period = 0, rotation_sequence = list(c(0)),
#'rotation_realloc = FALSE, prop = 1, aggreg = 1)
#'
#'## Inoculum
#' simul_params <- setInoculum(simul_params, 5e-4)
#' 
#'## list of outputs to be generated
#'outputlist <- loadOutputs(epid_outputs = "all", evol_outputs = "")
#'simul_params <- setOutputs(simul_params, outputlist)
#'
#'## Check, save and run simulation
#'checkSimulParams(simul_params)
#'runSimul(simul_params, graphic = TRUE)
#' }
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018).
#' Assessing the durability andefficiency of landscape-based strategies to deploy plant 
#' resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#' @export
runSimul <- function(params, graphic=TRUE, writeTXT=TRUE, videoMP4=FALSE, keepRawResults=FALSE) {

  ### !!!!! BE CAREFUL !!!!! ###
  ### croptypes, cultivars, genes and cultivarsGenes have to be ordered all in the same way
  ### ID have to match row index and col index

  initPath <- getwd()
  setwd(params@OutputDir)

  ## remove genes not used from CultivarsGenes and Genes
  if(ncol(params@CultivarsGenes) >= 1) {
    drop_genes <- lapply(1:ncol(params@CultivarsGenes), FUN = function(c) {
      if(sum(params@CultivarsGenes[,c]) == 0) return(c);
    })
    drop_genes <- which(!sapply(drop_genes,is.null))
  } else { drop_genes <- NULL }

  if( !is.null(drop_genes) && length(drop_genes) > 0){
    print(paste0("Genes not affected ",params@Genes[drop_genes,1]))
    cultivarsGenes_tmp <- as.data.frame(params@CultivarsGenes[,-drop_genes])  
    Genes_tmp <- as.data.frame(params@Genes[-drop_genes,])
  } else {
    cultivarsGenes_tmp <- params@CultivarsGenes
    Genes_tmp <- params@Genes 
  }
  cultivars_genes_list <- lapply(1:nrow(params@Cultivars), FUN = function(i) {
    return(which(cultivarsGenes_tmp[i, ] == 1) - 1)
  })

  cdf <- as.data.frame(params@Landscape)
  ncol <- length(grep("^year_", colnames(cdf)) %in% colnames(cdf))
  ## TODO: use value of Nyears in previous line?
  rotation <- as.matrix(cdf[, grep("^year_", colnames(cdf))], ncol = ncol)
  croptypes_cultivars_prop <- params2CroptypeBDD(params)[, c(2, 3, 4)]

  
  
  ## Run the simulation
  outputs <- simul_landsepi(
    seed = params@Seed,
    time_param = params@TimeParam,
    croptype_names = params@Croptypes$croptypeName,
    cultivars = params@Cultivars,
    cultivars_genes_list = cultivars_genes_list,
    genes = Genes_tmp,
    landscape = as_Spatial(st_geometry(params@Landscape)),
    area = as.vector(params@Landscape$area[, 1]),
    rotation = rotation,
    croptypes_cultivars_prop = croptypes_cultivars_prop,
    basic_patho_param = params@Pathogen,
    # repro_sex_prob = params@ReproSexProb,
    disp_patho_clonal = params@DispPathoClonal,
    disp_patho_sex = params@DispPathoSex,
    disp_host = params@DispHost,
    treatment = params@Treatment,
    pI0 = params@PI0,
    epid_outputs = params@Outputs[["epid_outputs"]],
    evol_outputs = params@Outputs[["evol_outputs"]],
    thres_breakdown = params@Outputs[["thres_breakdown"]],
    audpc100S = params@Outputs[["audpc100S"]],
    graphic = graphic,
    writeTXT = writeTXT,
    videoMP4 = videoMP4,
    keepRawResults = keepRawResults
  )

  setwd(initPath)
  return(outputs)
}


#' @name setSeed
#' @title Set the seed
#' @description Updates a LandsepiParams object with a seed value for random number generator
#' @param params a LandsepiParams Object.
#' @param seed an integer used as seed value (for random number generator).
#' @return a LandsepiParams object.
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setSeed(simul_params, 100)
#' simul_params@Seed
#' }
#' @export
setSeed <- function(params, seed) {
  params@Seed <- setSeedValue(seed)
  return(params)
}


#' @name setTime
#' @title Set time parameters
#' @description Updates a LandsepiParams object with time parameters : Nyears and nTSpY
#' @param params a LandsepiParams Object.
#' @param Nyears an integer giving the number of cropping seasons (e.g. years) to simulate.
#' @param nTSpY an integer giving the number of time steps per cropping season (e.g. days).
#' @return a LandsepiParams object.
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setTime(simul_params, Nyears=10, nTSpY=120)
#' simul_params@TimeParam
#' }
#' @export
setTime <- function(params, Nyears, nTSpY) {
  land <- params@Landscape
  
  if (nrow(land) > 0){
    st_geometry(land) <- NULL
    ldf <- as.data.frame(land)
    
    if (length(grep("^year_", colnames(ldf))) < Nyears) {
      message("Landscape croptypes affectation by year have to be regenerated")
    }
  }
  
  params@TimeParam <- list(Nyears = as.numeric(Nyears), nTSpY = as.numeric(nTSpY))
  checkTime(params)
  return(params)
}

#' @name checkTime
#' @title Check time
#' @description Checks time parameters validity
#' @param params a LandsepiParams Object.
#' @return a boolean TRUE if times are setted.
checkTime <- function(params) {
  validity <- TRUE
  if(is.null(params@TimeParam$Nyears) || !is.numeric(params@TimeParam$Nyears) 
     || length(params@TimeParam$Nyears) != 1 ) {
    warning("Invalid nb of years of simulation: use setTime()")
    validity <- FALSE
  } else if (!is.wholenumber(params@TimeParam$Nyears) 
             || !is.strict.positive(params@TimeParam$Nyears) ){
    warning("Nb of years of simulation must be a whole number > 0")
    validity <- FALSE
  }
  
  
  if(is.null(params@TimeParam$nTSpY) || !is.numeric(params@TimeParam$nTSpY) 
     || length(params@TimeParam$nTSpY) > 1 ) {
    warning("Invalid nb of steps per year: use setTime()")
    validity <- FALSE
  } else if (!is.wholenumber(params@TimeParam$nTSpY) 
             || !is.strict.positive(params@TimeParam$nTSpY)){
    warning("Nb of steps per year must be a whole number > 0")
    validity <- FALSE
  }
  
  return(validity)
}


#' @name setLansdcape
#' @title Set the landscape
#' @description Updates a LandsepiParams object with a sp or sf object as landscape.
#' @details The landscape should be a sp or sf object. Built-in landscape are available using 
#' \code{\link{loadLandscape}}. 
#' See our tutorial (vignettes) for details on how to use your own landscape.
#' If the landscape contains only polygons, croptypes can be allocated later using 
#' \code{\link{allocateLandscapeCroptypes}}.
#' Otherwise the landscape has to contain a data.frame specifying for every year, the index 
#' of the croptype cultivated in each polygon.
#' Each features has a field identified by "year_XX" (XX <- seq(1:Nyears+1)) and containing 
#' the croptype ID.
#'
#' | Features/fields | year_1 | year_2 | ... year_Nyears+1 |
#' |---------------- | ------ | ------ | ----------------- |
#' | polygons1       | 13     | 10     | 13                |
#' | polygonsX       | 2      | 1      | 2                 |
#'
#' @param params a LandsepiParams Object.
#' @param land a landscape as sp or sf object
#' @return a LandsepiParams object.
#' @seealso \link{loadLandscape}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' simul_params@Landscape
#' }
#' @importFrom sf st_as_sf
#' @export
setLandscape <- function(params, land) {
  if (class(land)[1] != "sf") {
    params@Landscape <- st_as_sf(land)
  } else {
    params@Landscape <- land
  }

  params@Landscape$area <- data.frame(area = st_area(params@Landscape))
  
  ## Initialise host and pathogen dispersal with diagonal matrices
  if (length(params@DispHost) == 0 |
      length(params@DispHost) != (nrow(params@Landscape))^2){
    disp_host <- loadDispersalHost(params, type = "no")
    params <- setDispersalHost(params, disp_host)
  }
  if (length(params@DispPathoClonal) == 0){
    disp_patho_clonal <- loadDispersalHost(params, type = "no") # (same function)
    disp_patho_sex <- loadDispersalHost(params, type = "no")
    params <- setDispersalPathogen(params, disp_patho_clonal, disp_patho_sex)
  }
  return(params)
}


#' @name loadLandscape
#' @title Load a landscape
#' @description Loads one of the five built-in landscapes simulated using a T-tesselation algorithm 
#' and composed of 155, 154, 152, 153 and 156 polygons, respectively.
#' Each landscape is identified by a numeric from 1 to 5.
#' @param id a landscape ID between 1 to 5 (default = 1)
#' @return a landscape in sp format
#' @seealso \link{landscapeTEST}, \link{setLandscape}
#' @examples
#' land <- loadLandscape(1)
#' length(land)
#' @export
loadLandscape <- function(id = 1) {
  if (id >= 1 && id <= 5) {
    land <- get(paste0("landscapeTEST", id))
  }
  else {
    stop("Indices of available landscapes are 1 to 5")
  }

  return(land)
}


#' @name checkLandscape
#' @title Check the landscape
#' @description Checks landscape validity
#' @param params a LandsepiParams Object.
#' @return TRUE if Ok, FALSE otherwise
checkLandscape <- function(params) {
  ret <- TRUE
  ## TODO : check bbox, proj4string, epsg and geometry as POLYGON

  land <- params@Landscape
  st_geometry(land) <- NULL
  ldf <- as.data.frame(land)

  # check CroptypeID present in Landscape
  if (sum(!unique(as.integer(unlist(ldf[, grep("^year_", colnames(ldf))]))) 
          %in% params@Croptypes$croptypeID) != 0) {
    warning("Croptypes undef in Landscape")
    warning("Croptypes ID : ", params@Croptypes$croptypeID)
    warning("Croptypes ID in Landscape : ", unique(unlist(ldf)))
    ret <- FALSE
  }

  # layer croptypeID need Nyears + 1 fields
  if (length(grep("^year_", colnames(ldf))) != params@TimeParam$Nyears + 1) {
    warning("Landscape Fields 'year_X' numbers [", length(grep("^year_", colnames(ldf)))
            , "] differ from Nyears ", params@TimeParam$Nyears, " +1")
    ret <- FALSE
  }

  return(ret)
}


#' @name setDispersalPathogen
#' @title Set pathogen dispersal
#' @description Updates a LandsepiParams object with a pathogen dispersal matrix.
#' Note that landscape parameters must be set before updating setting dispersal.
#' @details See tutorial (vignettes) on how to 
#' use your own landscape and compute your own pathogen dispersal kernel. 
#' The dispersal matrix a square matrix whose size is the number of polygons in the landscape 
#' and whose elements are, for each line i and each column i' the probability that propagules 
#' migrate from polygon i to polygon i'. 
#' Lines of the matrix can be normalised to sum to 1 (reflective boundaries); 
#' otherwise propagules dispersing outside the landscape are lost (absorbing boundaries).  
#' @param params a LandsepiParams Object.
#' @param mat_clonal a square matrix giving the probability of pathogen dispersal 
#' (clonal propagules) from any polygon of the landscape to any other polygon. 
#' It can be generated manually, or, alternatively, via \code{\link{loadDispersalPathogen}}. 
#' The size of the matrix must match the number of polygons in the landscape, and lines of 
#' the matrix may sum to 1 (reflecting boundaries) or be <1 (absorbing boundaries).
#' @param mat_sex a square matrix giving the probability of pathogen dispersal (sexual propagules) 
#' from any polygon of the landscape to any other polygon (default identity matrix) . 
#' It can be generated manually, or, alternatively, via \code{\link{loadDispersalPathogen}}. 
#' The size of the matrix must match the number of polygons in the landscape, and lines of 
#' the matrix may sum to 1 (reflecting boundaries) or be <1 (absorbing boundaries).
#' @return a LandsepiParam object.
#' @seealso \link{loadDispersalPathogen}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' d <- loadDispersalPathogen(1)
#' simul_params <- setDispersalPathogen(simul_params, d[[1]], d[[2]])
#' simul_params@DispPathoClonal
#' }
#' @export
setDispersalPathogen <- function(params, mat_clonal, mat_sex=NULL) {
  if (class(mat_clonal)[1] == "matrix") {
    params@DispPathoClonal <- as.vector(mat_clonal) # By columns
  } else {
    params@DispPathoClonal <- mat_clonal
  }
  
  if(is.null(mat_sex)) {
    params@DispPathoSex = c(diag(1,sqrt(length(params@DispPathoClonal))
                                 , sqrt(length(params@DispPathoClonal))))
  }
  else {
    if (class(mat_sex)[1] == "matrix") {
      params@DispPathoSex <- as.vector(mat_sex) # By columns
    } else {
      params@DispPathoSex <- mat_sex
    }
  }
  
  checkDispersalPathogen(params)

  return(params)
}


#' @name loadDispersalPathogen
#' @title Load pathogen dispersal matrices
#' @description It loads one of the five built-in vectorised dispersal matrices of rust fungi 
#' associated with the five built-in landscapes. Landscape and DispersalPathogen ID must be 
#' the same. And set a vectorized identity matrix for sexual reproduction dispersal.
#' @param id a matrix ID between 1 to 5 (must match the ID of the landscape loaded with 
#' \code{\link{loadLandscape}}).
#' @details *landsepi* includes built-in dispersal matrices to represent rust dispersal in the 
#' five built-in landscapes. These have been computed from a power-law dispersal kernel: 
#' \eqn{ g(d) = ((b-2)*(b-1) / (2*pi*a^2)) * (1 +  d/a)^{-b} }
#'  with a=40 the scale parameter and b=7 a parameter related to the width of the dispersal kernel. 
#'  The expected mean dispersal distance is given by \eqn{ 2*a /(b-3) = 20 m }.
#' @return a vectorised dispersal matrix representing the dispersal of clonal propagules, 
#' and a vectorised dispersal identity matrix for sexual propagules. All by columns.
#' @seealso \link{dispP}, \link{setDispersalPathogen}
#' @examples
#' d <- loadDispersalPathogen(1)
#' d
#' @export
loadDispersalPathogen <- function(id = 1) {
  if (id >= 1 && id <= 5) {
    disp_clonal <- get(paste0("dispP_", id))
  }
  else {
    warning("Indices of available pathogen dispersal matrices are 1 to 5")
    disp_clonal <- numeric()
  }
  
  disp_sex <- c(diag(1,sqrt(length(disp_clonal)),sqrt(length(disp_clonal))))

  return( list(disp_clonal=disp_clonal, disp_sex=disp_sex) )
}



#' @name checkDispersalPathogen
#' @title Check pathogen dispersal
#' @description Checks pathogen dispersal validity
#' @param params a LandsepiParams Object.
#' @return a boolean TRUE if OK, FALSE otherwise
checkDispersalPathogen <- function(params) {
  ret <- TRUE
  if (length(params@DispPathoClonal) != nrow(params@Landscape) * nrow(params@Landscape)) {
    warning("Size of pathogen dispersal is not landscape features^2")
    ret <- FALSE
  }

  if (sum(params@DispPathoClonal > 1) != 0 || sum(params@DispPathoClonal < 0) != 0) {
    warning("Probabilities of pathogen dispersal must be in [0,1]")
    warning(params@DispPathoClonal[which(params@DispPathoClonal > 1)])
    warning(params@DispPathoClonal[which(params@DispPathoClonal < 0)])
    ret <- FALSE
  }

  if (length(params@DispPathoSex) != nrow(params@Landscape) * nrow(params@Landscape)) {
    warning("Size of pathogen dispersal is not landscape features^2")
    ret <- FALSE
  }
  
  if (sum(params@DispPathoSex > 1) != 0 || sum(params@DispPathoSex < 0) != 0) {
    warning("Probabilities of Sexual pathogen dispersal must be in [0,1]")
    warning(params@DispPathoSex[which(params@DispPathoSex > 1)])
    warning(params@DispPathoSex[which(params@DispPathoSex < 0)])
    ret <- FALSE
  }
  
  return(ret)
}

#' @name setDispersalHost
#' @title Set host dispersal
#' @description Updates a LandsepiParams object with a host dispersal matrix.
#' Note that landscape parameters must be set before updating setting dispersal.
#' @details the dispersal matrix gives the probability for a host individual in a polygon i (row)
#' to migrate to polygon j (column) through dispersal. 
#' If the host is a cultivated plant: seeds are harvested and do not disperse. 
#' Thus the dispersal matrix is the identity matrix.
#' @param params a LandsepiParams Object.
#' @param mat a square matrix giving the probability of host dispersal
#' from any polygon of the landscape to any other polygon. 
#' It can be generated manually, or, alternatively, via \code{\link{loadDispersalHost}}.
#' The size of the matrix must match the number of polygons in the landscape.
#' @return a LandsepiParam object.
#' @seealso \link{loadDispersalHost}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' d <- loadDispersalHost(simul_params)
#' simul_params <- setDispersalHost(simul_params, d)
#' simul_params@DispHost
#' }
#' @export
setDispersalHost <- function(params, mat) {
  if (class(mat)[1] == "matrix") {
    params@DispHost <- as.vector(mat) # By columns
  } else {
    params@DispHost <- mat
  }

  checkDispersalHost(params)

  return(params)
}



#' @name loadDispersalHost
#' @title Load a host dispersal matrix
#' @description It loads a vectorised diagonal matrix to simulate no host dispersal.
#' @details as the size of the matrix depends on the number of polygons in the landscape, 
#' the landscape must be defined before calling \code{loadDispersalHost}. 
#' @param params a LandsepiParams Object.
#' @param type a character string specifying the type of dispersal ("no" for no dispersal)
#' @return a vectorised dispersal matrix.
#' @seealso \link{setDispersalHost}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' d <- loadDispersalHost(simul_params)
#' d
#' }
#' @export
loadDispersalHost <- function(params, type = "no") {
  if (nrow(params@Landscape) == 0) {
    warning("Lanscape has to be set before loading host dispersal matrix")
  }
  else {
    disp <- switch(type,
                   "no" = diag(1, nrow(params@Landscape), nrow(params@Landscape)),
                   diag(1, nrow(params@Landscape), nrow(params@Landscape))
    )# By columns
  }
  
  return(disp)
}


#' @name checkDispersalHost
#' @title Check host dispersal
#' @description Checks host dispersal matrix validity.
#' @param params a LandsepiParams Object.
#' @return a boolean TRUE if OK, FALSE otherwise
checkDispersalHost <- function(params) {
  ret <- TRUE
  if (length(params@DispHost) != nrow(params@Landscape) * nrow(params@Landscape)) {
    warning("Size of vector of host dispersal is not landscape features^2")
    ret <- FALSE
  }

  if (sum(params@DispHost > 1) != 0 || sum(params@DispHost < 0) != 0) {
    warning("Host dispersal probabilities must be in [0,1]")
    warning(params@DispHost[which(params@DispHost > 1)])
    warning(params@DispHost[which(params@DispHost < 0)])
    ret <- FALSE
  }

  return(ret)
}



#' @name allocateLandscapeCroptypes
#' @title Allocate croptypes to the landscape
#' @description Updates the landscape of a LandsepiParams object with croptype allocation in 
#' every polygon of the landscape and every year of simulation. Allocation is based on an algorithm 
#' which controls croptype proportions (in surface) and spatio-temporal aggregation.
#' Note that time, landscape and croptype parameters must be set before allocating 
#' landscape croptypes.
#' @param params a LandsepiParams Object.
#' @param rotation_period number of years before rotation of the landscape. There is no rotation 
#' if rotation_period=0 or rotation_period=Nyears.
#' @param rotation_sequence a list, each element of the list contains indices of croptypes that 
#' are cultivated during a period given by "rotation_period". There is no change in cultivated 
#' croptypes if the list contains only one element (e.g. only one vector c(0,1,2), indicating 
#' cultivation of croptypes 0, 1 and 2).
#' @param rotation_realloc a logical indicating if a new random allocation of croptypes is 
#' performed when the landscape is rotated (FALSE=static allocation, TRUE=dynamic allocation). 
#' Note that if rotation_realloc=FALSE, all elements of the list "rotation_sequence" must have 
#' the same length, and only the first element of the lists "prop" and "aggreg" will be used.
#' @param prop a list of the same size as "rotation_sequence", each element of the list contains 
#' a vector of the proportions (in surface) associated with the croptypes in "rotation_sequence". 
#' A single vector can be given instead of a list if all elements of "rotation_sequence" are 
#' associated with the same proportions.
#' @param aggreg a list of the same size as "rotation_sequence", each element of the list is a 
#' single double indicating the degree of
#' aggregation of the landscape. This double must greater or equal 0; the greater its value, 
#' the higher the degree of spatial aggregation (roughly, aggreg between 0 and 0.1 for fragmented 
#' landscapes, between 0.1 and 0.5 for balanced landscapes, between 0.5 and 3 for aggregated 
#' landscapes, and above 3 for highly aggregated landscapes). A single double can be given 
#' instead of a list if all elements of "rotation_sequence" are associated with the same level 
#' of aggregation.
#' @param algo the algorithm used for the computation of the variance-covariance matrix 
#' of the multivariate normal distribution: "exp" for exponential function, "periodic" 
#' for periodic function, "random" for random draw (see details of function multiN). 
#' If algo="random", the parameter aggreg is not used. 
#' Algorithm "exp" is preferable for big landscapes.
#' @param graphic a logical indicating if graphics must be generated (TRUE) or not (FALSE).
#' @details An algorithm based on latent Gaussian fields is used to allocate two different 
#' croptypes across the simulated landscapes (e.g. a susceptible and a resistant cultivar, 
#' denoted as SC and RC, respectively). This algorithm allows the control of the proportions of 
#' each croptype in terms of surface coverage, and their level of spatial aggregation. 
#' A random vector of values is drawn from a multivariate normal distribution with expectation 0 
#' and a variance-covariance matrix which depends on the pairwise distances between
#' the centroids of the polygons. Next, the croptypes are allocated to different polygons 
#' depending on whether each value drawn from the multivariate normal distribution is above 
#' or below a threshold. The proportion of each cultivar in the landscape is controlled by 
#' the value of this threshold. To allocate more than two croptypes, \code{AgriLand} uses 
#' sequentially this algorithm. For instance, the allocation of three croptypes 
#' (e.g. SC, RC1 and RC2) is performed as follows:
#' \enumerate{
#' \item the allocation algorithm is run once to segregate the polygons where the susceptible 
#' cultivar is grown, and
#' \item the two resistant cultivars (RC1 and RC2) are assigned to the remaining candidate 
#' polygons by re-running the allocation algorithm.
#' }
#' @return a LandsepiParams object with Landscape updated with the layer "croptypeID". 
#' It contains croptype allocation in every polygon of the landscape for all years of simulation.
#' @examples
#' \dontrun{
#' ## Initialisation
#' simul_params <- createSimulParams(outputDir = getwd())
#' ## Time parameters
#' simul_params <- setTime(simul_params, Nyears = 10, nTSpY = 120)
#' ## Landscape
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' ## Cultivars
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' ## Allocate cultivars to croptypes
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
#' , "Resistant crop 1"
#' , "Resistant crop 2"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
#' simul_params <- setCroptypes(simul_params, croptypes)
#' ## Allocate croptypes to landscape        
#' rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
#' rotation_period <- 0 ## same croptypes every years
#' prop <- c(1 / 3, 1 / 3, 1 / 3) ## croptypes proportions
#' aggreg <- 10 ## aggregated landscape
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = rotation_period,
#' rotation_sequence = rotation_sequence,
#' rotation_realloc = FALSE, prop = prop, aggreg = aggreg)
#' simul_params@Landscape
#' }
#' @export
allocateLandscapeCroptypes <- function(params, rotation_period, rotation_sequence
                                       , rotation_realloc = FALSE
                                       , prop, aggreg, algo = "periodic", graphic = TRUE) {
  ### TODO Check validity of params slot before Agriland

  orig_landscape <- params@Landscape

  croptypeSP <- AgriLand(as_Spatial(params@Landscape),
    Nyears = params@TimeParam$Nyears,
    rotation_period = rotation_period,
    rotation_sequence = rotation_sequence,
    rotation_realloc = rotation_realloc,
    prop = prop,
    aggreg = aggreg,
    algo = algo,
    croptype_names = params@Croptypes$croptypeName,
    graphic = graphic,
    outputDir = params@OutputDir
  )
  params@Landscape <- st_as_sf(croptypeSP)
  params@Landscape$area <- data.frame(area = st_area(params@Landscape))
  if (length(orig_landscape$ID) != 0 && length(orig_landscape$Name) != 0) {
    params@Landscape$ID <- orig_landscape$ID
    params@Landscape$NAME <- orig_landscape$NAME
  }

  return(params)
}

#' @name loadTreatment
#' @title Load treatment parameters
#' @description Loads a list of treatment parameters for a specific disease (initialised at 0
#' , i.e. absence of treatments)
#' @param disease a disease name, among "mildew", "sigatoka" and "no pathogen"
#' @details Chemical treatment is applied in a polygon only if disease severity (i.e. I/N) in 
#' this polygon exceeds the threshold given by `treatment_application_threshold`. 
#' Treatment efficiency is maximum (i.e. equal to the parameter treatment_efficiency) 
#' at the time of treatment application (noted \eqn{t*}); then it decreases with time 
#' (i.e. natural pesticide degradation) and host growth (i.e. new biomass is not protected 
#' by treatments):                                                                                                                                                                                     protected by treatments):Efficiency of the treatment at time t after the application date is given by:
#' \eqn{ efficiency(t) = treatment\_efficiency / (1 + exp(a-b*C(t))) }
#' with \eqn{ C(t)= C_1 * C_2}: \itemize{
#' \item{\eqn{C_1 = exp(- treatment\_degradation\_rate * \Delta t) } is the reduction of 
#' fungicide concentration due to time (e.g. natural degradation, volatilization, weathering), 
#' with \eqn{\Delta t = t - t*} the timelag passed since the time of 
#' treatment application.}
#' \item{ \eqn{ C_2 = min(1, N(t*) / N(t)) } is the reduction of fungicide concentration due 
#' to plant growth, since new plant tissue is not covered by fungicide. 
#' \eqn{N(t*)} and \eqn{N(t)} being the number of 
#' host individuals  a the time of treatment \eqn{t*} and at time \eqn{t}, respectively.}
#' \item{ \eqn{a \in [3.5 ; 4.5]} and \eqn{b \in [8 ; 9]} are shape parameters.}
#' } 
#' @return a list of treatment parameters:
#' \itemize{ 
#' \item treatment_degradation_rate = degradation rate (per time step) of chemical concentration,
#' \item treatment_efficiency = maximal efficiency of chemical treatments (i.e. fractional reduction 
#' of pathogen infection rate at the time of application),
#' \item treatment_timesteps = vector of time steps corresponding to treatment application dates,
#' \item treatment_cultivars = vector of indices of the cultivars that receive treatments,
#' \item treatment_cost = cost of a single treatment application (monetary units/ha)
#' \item treatment_application_threshold = vector of thresholds (i.e. disease severity, one 
#' for each treated cultivar) above which the treatment is applied in a polygon.
#' }
#' @seealso \link{setTreatment}
#' @examples
#' treat <- loadTreatment("sigatoka")
#' treat
#' @export
loadTreatment <- function(disease="no pathogen") {
  treatment <- switch(disease,
                      "mildew" = list(treatment_degradation_rate = 0.1,
                                      treatment_efficiency = 0.0,
                                      treatment_timesteps = numeric(),
                                      treatment_cultivars = numeric(),
                                      treatment_cost = 0,
                                      treatment_application_threshold = numeric())
                      
                      , "sigatoka" = list(treatment_degradation_rate = 0.1,
                                          treatment_efficiency = 1.0,
                                          treatment_timesteps = seq(5,182,5),
                                          treatment_cultivars = 0,
                                          treatment_cost = 1,
                                          treatment_application_threshold = 0)
                      
                      , list(treatment_degradation_rate = 0.1,
                             treatment_efficiency = 0.0,
                             treatment_timesteps = numeric(),
                             treatment_cultivars = numeric(),
                             treatment_cost = 0,
                             treatment_application_threshold = numeric())
  )
  
  if (length(treatment) == 0) {
    warning('Unknown type of disease: "', disease
            , '". Currently the only possible types are: "midlew", "sigatoka"')
  }
  return(treatment)
}

#' @name setTreatment
#' @title Set chemical treatments
#' @description Updates a LandsepiParams object with treatment parameters
#' @details Chemical treatment is applied in a polygon only if disease severity (i.e. I/N) in 
#' this polygon exceeds the threshold given by `treatment_application_threshold`. 
#' Treatment efficiency is maximum (i.e. equal to the parameter treatment_efficiency) 
#' at the time of treatment application (noted \eqn{t*}); then it decreases with time 
#' (i.e. natural pesticide degradation) and host growth (i.e. new biomass is not protected 
#' by treatments):                                                                                                                                                                                     protected by treatments):Efficiency of the treatment at time t after the application date is given by:
#' \eqn{ efficiency(t) = treatment\_efficiency / (1 + exp(a-b*C(t))) }
#' with \eqn{ C(t)= C_1 * C_2}: \itemize{
#' \item{\eqn{C_1 = exp(- treatment\_degradation\_rate * \Delta t) } is the reduction of 
#' fungicide concentration due to time (e.g. natural degradation, volatilization, weathering), 
#' with \eqn{\Delta t = t - t*} the timelag passed since the time of 
#' treatment application.}
#' \item{ \eqn{ C_2 = min(1, N(t*) / N(t)) } is the reduction of fungicide concentration due 
#' to plant growth, since new plant tissue is not covered by fungicide. 
#' \eqn{N(t*)} and \eqn{N(t)} being the number of 
#' host individuals  a the time of treatment \eqn{t*} and at time \eqn{t}, respectively.}
#' \item{ \eqn{a \in [3.5 ; 4.5]} and \eqn{b \in [8 ; 9]} are shape parameters.}
#' } 
#' An empty list of treatments (i.e. absence of application) can be loaded using 
#' \code{\link{loadPathogen}}.
#' @param params a LandsepiParams Object.
#' @param treatment_params list of parameters related to pesticide treatments: \itemize{ 
#' \item treatment_degradation_rate = degradation rate (per time step) of chemical concentration,
#' \item treatment_efficiency = maximal efficiency of chemical treatments 
#' (i.e. fractional reduction of pathogen infection rate at the time of application),
#' \item treatment_timesteps = vector of time steps corresponding to treatment application dates,
#' \item treatment_cultivars = vector of indices of the cultivars that receive treatments,
#' \item treatment_cost = cost of a single treatment application (monetary units/ha)
#' \item treatment_application_threshold = vector of thresholds 
#' (i.e. disease severity, one for each treated cultivar) above which the treatment 
#' is applied in a polygon.
#' }
#' @return a LandsepiParams object
#' @seealso \link{loadTreatment}
#' @examples
#' \dontrun{
#' t <- loadTreatment()
#' simul_params <- setTreatment(simul_params, t)
#' simul_params@Treatment
#' }
#' @export
setTreatment <- function(params, treatment_params) {
  params@Treatment <- treatment_params
  checkTreatment(params)
  
  return(params)
}

#' @name checkTreatment
#' @title Check treatment
#' @description Checks treatment validity
#' @param params a LandsepiParams Object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkTreatment <- function(params) {
  
  ret <- TRUE
  
  if (length(params@Treatment$treatment_cultivars) > 0) {
    warning("Simulation with chemical treatment applications")
    
    if (!is.numeric(params@Treatment$treatment_efficiency) ||
        !is.in.01(params@Treatment$treatment_efficiency)){
      warning("Treatment efficiency must be between 0 and 1")
      ret <- FALSE
    }
    if ( !is.numeric(params@Treatment$treatment_degradation_rate) ||
         !is.positive(params@Treatment$treatment_degradation_rate) ){
      warning("Treatment degradation rate must be >=0")
      ret <- FALSE
    }
    if ( !is.numeric(params@Treatment$treatment_cost) ||
         !is.positive(params@Treatment$treatment_cost) ){
      warning("Treatment cost must be >=0")
    }
    if ( !is.numeric(params@Treatment$treatment_timesteps) ||
         sum( !is.strict.positive(params@Treatment$treatment_timesteps) )>0 ){
      warning("Timesteps of treatments must be >0")
      ret <- FALSE
    }
    if (!is.numeric(params@Treatment$treatment_application_threshold) ||
        sum(!is.in.01(params@Treatment$treatment_application_threshold))>0){
      warning("Treatment application threshold must be between 0 and 1")
      ret <- FALSE
    }
  }
  return(ret)
}

#' @name loadPathogen
#' @title Load pathogen parameters
#' @description Loads default pathogen parameters for a specific disease
#' @details Available diseases:
#' * "no pathogen"
#' * "rust" (genus \emph{Puccinia}, e.g. stripe rust, stem rust and leaf rust of wheat and barley)
#' * "mildew" (\emph{Plasmopara viticola}, downy mildew of grapevine)
#' * "sigatoka" (\emph{Pseudocercospora fijiensis}, black sigatoka of banana)
#' Note that when disease = "mildew" a price reduction between 0% and 5% is applied to the 
#' market value according to disease severity. 
#' @param disease a disease name, among "rust" (default), "mildew", "sigatoka" and "no pathogen"
#' @return a list of pathogen parameters on a susceptible host
#' for a pathogen genotype not adapted to resistance
#' @seealso \link{setPathogen}
#' @examples
#' basic_patho_params <- loadPathogen()
#' basic_patho_params
#' @export
loadPathogen <- function(disease = "rust") {
  patho <- switch(disease,
    "rust" = list(
      name = "rust",
      survival_prob = 1e-4,
      repro_sex_prob = 0,
      infection_rate = 0.4,
      propagule_prod_rate = 3.125,
      latent_period_mean = 10,
      latent_period_var = 9,
      infectious_period_mean = 24,
      infectious_period_var = 105,
      sigmoid_kappa = 5.333,
      sigmoid_sigma = 3,
      sigmoid_plateau = 1,
      sex_propagule_viability_limit  = 1,
      sex_propagule_release_mean = 1,
      clonal_propagule_gradual_release = 0
    ),  
    "mildew" = list(name = "mildew",
      survival_prob = 1e-4,
      repro_sex_prob = 0,
      infection_rate = 0.9,
      propagule_prod_rate = 2.0,
      latent_period_mean = 7,
      latent_period_var = 8,
      infectious_period_mean = 14,
      infectious_period_var = 22,
      sigmoid_kappa = 5.333,
      sigmoid_sigma = 3,
      sigmoid_plateau = 1,
      sex_propagule_viability_limit  = 5,
      sex_propagule_release_mean = 1,
      clonal_propagule_gradual_release = 1
    ),
    "sigatoka" = list(name = "sigatoka",
                    survival_prob = 1/2,
                    repro_sex_prob = 0,
                    infection_rate = 0.02,
                    propagule_prod_rate = 90.9,
                    latent_period_mean = 25.5,
                    latent_period_var = 1.5,
                    infectious_period_mean = 22,
                    infectious_period_var = 14,
                    sigmoid_kappa = 5.333,
                    sigmoid_sigma = 3,
                    sigmoid_plateau = 1,
                    sex_propagule_viability_limit  = 1,
                    sex_propagule_release_mean = 1,
                    clonal_propagule_gradual_release = 0
    ),
    list(
      name = "no pathogen",
      survival_prob = 0,
      repro_sex_prob = 0,
      infection_rate = 0,
      propagule_prod_rate = 0,
      latent_period_mean = 0,
      latent_period_var = 0,
      infectious_period_mean = 0,
      infectious_period_var = 0,
      sigmoid_kappa = 0,
      sigmoid_sigma = 0,
      sigmoid_plateau = 1,
      sex_propagule_viability_limit  = 1,
      sex_propagule_release_mean = 1, 
      clonal_propagule_gradual_release = 0)
  )

  if (length(patho) == 0) {
    warning('Unknown type of disease: "', disease
            , '". Currently the only possible types are: "rust", "midlew", "sigatoka"')
  }

  return(patho)
}


#' @name setPathogen
#' @title Set the pathogen
#' @description Updates a LandsepiParams object with pathogen parameters
#' @details a set of parameters representative of rust fungi, downy mildew or black sigatoka 
#' can be loaded via \code{\link{loadPathogen}}.
#' @param params a LandsepiParams Object.
#' @param patho_params a list of pathogen aggressiveness parameters on a susceptible host
#' for a pathogen genotype not adapted to resistance: \itemize{
#' \item infection_rate = maximal expected infection rate of a propagule on a healthy host,
#' \item propagule_prod_rate = maximal expected effective propagule production rate of an 
#' infectious host per time step,
#' \item latent_period_mean = minimal expected duration of the latent period,
#' \item latent_period_var = variance of the latent period duration,
#' \item infectious_period_mean = maximal expected duration of the infectious period,
#' \item infectious_period_var = variance of the infectious period duration,
#' \item survival_prob = probability for a propagule to survive the off-season,
#' \item repro_sex_prob = probability for an infectious host to reproduce via sex rather 
#' than via cloning,
#' \item sigmoid_kappa = kappa parameter of the sigmoid contamination function,
#' \item sigmoid_sigma = sigma parameter of the sigmoid contamination function,
#' \item sigmoid_plateau = plateau parameter of the sigmoid contamination function,
#' \item sex_propagule_viability_limit = maximum number of cropping seasons up to which 
#' a sexual propagule is viable
#' \item sex_propagule_release_mean = average number of seasons after which a sexual 
#' propagule is released.
#' \item clonal_propagule_gradual_release = whether or not clonal propagules surviving 
#' the bottleneck are gradually released along the following cropping season.
#' }
#' It can be generated manually, or, alternatively, via \code{\link{loadPathogen}}.
#' @return a LandsepiParams object
#' @seealso \link{loadPathogen}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setPathogen(simul_params, loadPathogen())
#' simul_params@Pathogen
#' }
#' @export
setPathogen <- function(params, patho_params) {
  params@Pathogen <- patho_params
  
  # if( length(params@Pathogen$repro_sex_prob) != params@TimeParam$nTSpY+1 ) {
  #   params@Pathogen$repro_sex_prob <- rep(patho_params$repro_sex_prob, params@TimeParam$nTSpY +1)
  # }
  checkPathogen(params)

  return(params)
}


#' @name updateReproSexProb
#' @title Update the probability of sexual reproduction
#' @description set the probabilities for an infectious host to reproduce via sex rather 
#' than via cloning at every time step.
#' Note that time parameters must be set before updating sexual reproduction probabilities.
#' @param params a LandsepiParams object
#' @param vec a vector of size TimeParam$nTSpY +1 (season end) with the probabilities 
#' for an infectious host to reproduce via sex rather than via cloning at each time step. 
#' @return a LandsepiParams object updated
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setTime(simul_params, Nyears=10, nTSpY=120)
#' simul_params <- setPathogen(simul_params, loadPathogen("rust"))
#' repro_sex_probs <- c(rep(0.0, 120), 1.0)  
#' simul_params <- updateReproSexProb(simul_params, repro_sex_probs)
#' simul_params@Pathogen
#' }
#' @export
updateReproSexProb <- function(params, vec) {
  if( params@TimeParam$nTSpY+1 != length(vec) ){
    warning("Vector of Probability of sexual reproduction not compatible with the nSTpY value")
  }
  else {
    params@Pathogen$repro_sex_prob <- vec
  }

  return(params)
}


#' @name updateSurvivalProb
#' @title Update pathogen survival probability during the off-season
#' @description update survival probability of the pathogen with a probability value for every 
#' simulated year (number of years = Nyears) and every croptype (number of croptypes = Ncroptypes). 
#' Note that time parameters, pathogen and croptypes must be set before updating 
#' survival probabilities.
#' @param params a LandsepiParams object
#' @param mat_year a vector of size Nyear, giving survival probabilities for every year 
#' (replicated for every croptype).
#' @param mat_croptype a vector of size Ncroptypes, giving survival probabilities for 
#' every croptype (replicated for every year).
#' @param mat a matrix of dimension (Nyears, Ncroptypes) giving survival probabilities 
#' for every year (rows) 
#' and every croptype (columns).
#' @details Unless the matrix \code{mat} is filled, the matrix containing the survival 
#' probability during the offseason 
#' is computed for every year and croptype with 
#' \code{mat[year, croptype] = mat_year[year] * mat_croptype[croptype]}. \cr
#' @return a LandsepiParams object updated.
#' @seealso \link{survivalProbToMatrix}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setTime(simul_params, Nyears=10, nTSpY=120)
#' simul_params <- setPathogen(simul_params, loadPathogen("rust"))
#'
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' 
#' croptypes <- loadCroptypes(simul_params
#' , names = c("Susceptible crop", "Resistant crop", "Mixture"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop", "Resistant")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Susceptible", "Resistant"))
#' simul_params <- setCroptypes(simul_params, croptypes)
#' Ncroptypes <- nrow(simul_params@Croptypes)
#' Nyears <- simul_params@TimeParam$Nyears
#' 
#' ## Same probability in every croptype:
#' simul_params <- updateSurvivalProb(simul_params, mat_year=1:Nyears/100)
#' simul_params@Pathogen
#' ## Same probability every year:
#' simul_params <- updateSurvivalProb(simul_params, mat_croptype=1:Ncroptypes/10)
#' simul_params@Pathogen
#' ## specific probability for different croptypes and years:
#' simul_params <- updateSurvivalProb(simul_params
#' , mat_year=1:Nyears/100, mat_croptype=1:Ncroptypes/10)
#' simul_params@Pathogen
#' ## One probability per year and per croptype:
#' simul_params <- updateSurvivalProb(simul_params
#' , mat=matrix(runif(Nyears*Ncroptypes), ncol=Ncroptypes))
#' simul_params@Pathogen
#' survivalProbToMatrix(simul_params)
#' }
#' @export
updateSurvivalProb <- function(params, mat_year=NULL, mat_croptype=NULL, mat=NULL) {
  Nyears <- params@TimeParam$Nyears
  Ncroptypes <- nrow(params@Croptypes)
  
  if (Nyears==0 | Ncroptypes==0){
    stop("Please set the number of years and a set of croptypes before updating 
         survival probabilities")
  }
  
  if (all(missing(mat), missing(mat_year), missing(mat_croptype)))
    stop("Missing argument: The probability of survival must be defined")
  
  if (!is.null(mat)){
    if(any(dim(mat) != c(Nyears, Ncroptypes))){
      stop("The matrix of survival probabilities ('mat') must be of dimensions 
           (Nyears, Ncroptypes)")
    }else{
      if (any(!is.null(mat_year), !is.null(mat_croptype)))
        warning("'mat_year', 'mat_croptype' are not accounted if 'mat' is filled")
    }
  } else {  ## i.e. mat is null
    
    if (is.null(mat_year)){
      mat_year <- rep(1, Nyears)
    }else{
      if (length(mat_year) != Nyears)
        stop("'mat_year' must have the same length as the number of years")
    }
    if (is.null(mat_croptype)){
      mat_croptype <- rep(1, Ncroptypes)
    }else{
      if (length(mat_croptype) != Ncroptypes)
        stop("'mat_croptype' must have the same length as the number of croptypes")
    } 
    
    ## Computation of the matrix of survival probabilities
    mat <- matrix(0, nrow=Nyears, ncol=Ncroptypes)
    for (year in 1:Nyears){
      for (croptype in 1:Ncroptypes){
        mat[year, croptype] <- mat_year[year] * mat_croptype[croptype]
      }
    }
  } ## else mat is null
  
  params@Pathogen$survival_prob <- as.vector(mat)
  return(params)
}


#' @name survivalProbToMatrix
#' @title Survival probability To Matrix
#' @description Transform the off-season survival probability of the pathogen 
#' (1D vector of length Nyears*Ncroptypes) into a matrix (for visualization purpose)
#' @details After updating the off-season survival probability with \code{updateSurvivalProb()}, 
#' this function returns the probability as a matrix for every year (rows) and croptypes (columns) 
#' as well as, if croptypes have been previously allocated to a landscape, a matrix for every 
#' polygon (rows) and year (columns). 
#' @param params a LandsepiParams object.
#' @return a list containing a matrix of dimensions (Nyears, Ncroptypes) as well as a matrix of 
#' dimensions (Npoly, Nyears)
#' @seealso \link{updateSurvivalProb}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setTime(simul_params, Nyears=10, nTSpY=120)
#' simul_params <- setPathogen(simul_params, loadPathogen("rust"))
#'
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' 
#' croptypes <- loadCroptypes(simul_params
#' , names = c("Susceptible crop", "Resistant crop", "Mixture"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop", "Resistant")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Susceptible", "Resistant"))
#' simul_params <- setCroptypes(simul_params, croptypes)
#' 
#' Ncroptypes <- nrow(simul_params@Croptypes)
#' Nyears <- simul_params@TimeParam$Nyears
#' 
#' landscape <- loadLandscape(1)
#' simul_params <- setLandscape(simul_params, landscape)
#' simul_params <- allocateLandscapeCroptypes(simul_params,
#' rotation_period = 0, rotation_sequence = croptypes$croptypeID,
#' rotation_realloc = FALSE,
#' prop = rep(1/Ncroptypes, Ncroptypes),
#' aggreg = 0.05, graphic = FALSE)
#' 
#' ## One probability per year and per croptype:
#' simul_params <- updateSurvivalProb(simul_params
#' , mat=matrix(runif(Nyears*Ncroptypes), ncol=Ncroptypes))
#' simul_params@Pathogen
#' survivalProbToMatrix(simul_params)
#' }
#' @export
survivalProbToMatrix <- function(params){ 
  Nyears <- params@TimeParam$Nyears
  croptypes <- params@Croptypes
  Ncroptypes <- nrow(croptypes)
  cdf <- as.data.frame(params@Landscape)
  ncol <- length(grep("^year_", colnames(cdf)) %in% colnames(cdf))
  rotation <- as.matrix(cdf[, grep("^year_", colnames(cdf))], ncol = ncol)
  Npoly <- nrow(rotation)
  
  survival_prob <- matrix(params@Pathogen$survival_prob, nrow=Nyears, ncol=Ncroptypes)
  rownames(survival_prob) <- paste0("year_", 1:Nyears)
  colnames(survival_prob) <- croptypes$croptypeName
  
  if(Npoly>0){
    survival_prob_poly <- matrix(NA, nrow=Npoly, ncol=Nyears)
    colnames(survival_prob_poly) <- paste0("year_", 1:Nyears)
    rownames(survival_prob_poly) <- paste0("poly_", 1:Npoly)
    for (poly in 1:Npoly){
      for (year in 1:Nyears){
        survival_prob_poly[poly, year] <- survival_prob[year, rotation[poly, year]+1] 
        ##+1 because of C indices starting at 0
      }
    }
    return(list(survival_prob=survival_prob, survival_prob_poly=survival_prob_poly))
    
  }else{
    return(survival_prob)
  }
}


#' @name checkPathogen
#' @title Check pathogen
#' @description Checks pathogen validity
#' @param params a LandsepiParams Object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkPathogen <- function(params) {
  
  ret <- TRUE
  if (length(params@Pathogen) == 0 ||
      sum( sapply(params@Pathogen[names(params@Pathogen)!="repro_sex_prob" 
                                  & names(params@Pathogen)!="survival_prob"]
                  , length) != rep(1,length(params@Pathogen)-2) ) > 0 ){
    warning("Invalid parameters for Pathogen, use setPathogen()")
    ret <- FALSE
    return(ret)
  }

  if ( !is.numeric(params@Pathogen$infection_rate) ||
       !is.in.01(params@Pathogen$infection_rate) ){
    warning("Infection rate must be between 0 and 1")
    ret <- FALSE
  }
  if ( !is.numeric(params@Pathogen$propagule_prod_rate) ||
       !is.positive(params@Pathogen$propagule_prod_rate) ){
    warning("Propagule production rate must be >= 0")
    ret <- FALSE
  }
  if (!is.numeric(params@Pathogen$survival_prob) ||
      sum(params@Pathogen$survival_prob < 0) > 0 || 
      sum(params@Pathogen$survival_prob > 1) > 0){      
    warning("Survival probability must be between 0 and 1")
    ret <- FALSE
  }
  if (!is.numeric(params@Pathogen$repro_sex_prob) || 
      sum(params@Pathogen$repro_sex_prob < 0) > 0 || 
      sum(params@Pathogen$repro_sex_prob > 1) > 0){
    warning("Probability of sexual reproduction must be between 0 and 1")
    ret <- FALSE
  }

  # if( params@TimeParam$nTSpY+1 != length(params@Pathogen$repro_sex_prob)){
  #   warning("Vector of Probability of sexual reproduction not compatible with the nSTpY value")
  #   ret <- FALSE
  # }
  
  if (!is.numeric(params@Pathogen$sigmoid_plateau) || 
      !is.in.01(params@Pathogen$sigmoid_plateau) ){
    warning("sigmoid_plateau must be between 0 and 1")
    ret <- FALSE
  }
  if (!is.numeric(params@Pathogen$latent_period_mean) ||
      !is.numeric(params@Pathogen$latent_period_var) ||
      !is.numeric(params@Pathogen$infectious_period_mean) ||
      !is.numeric(params@Pathogen$infectious_period_var) ||
      !is.numeric(params@Pathogen$sigmoid_kappa) ||
      !is.numeric(params@Pathogen$sigmoid_sigma) ||
      
      !is.positive(params@Pathogen$latent_period_mean) || 
      !is.positive(params@Pathogen$infectious_period_mean) ||
      !is.positive(params@Pathogen$latent_period_var) ||
      !is.positive(params@Pathogen$infectious_period_var) ||
      !is.positive(params@Pathogen$sigmoid_sigma) || 
      !is.positive(params@Pathogen$sigmoid_kappa) ){
    warning("Latent period, infectious period and sigmoid parameters must be >= 0")
    ret <- FALSE
  }
  
  if (!is.numeric(params@Pathogen$sex_propagule_viability_limit ) ||
      sum(!is.wholenumber(params@Pathogen$sex_propagule_viability_limit) > 0) ||
      !is.strict.positive( params@Pathogen$sex_propagule_viability_limit)){
    warning("sex_propagule_viability_limit  must be a whole number > 0")
    ret <- FALSE
  }
  if (!is.numeric(params@Pathogen$sex_propagule_release_mean) ||
      !is.strict.positive( params@Pathogen$sex_propagule_release_mean)){
    warning("sex_propagule_release_mean  must be > 0 ")
    ret <- FALSE
  }

  return(ret)
}



#' @name loadCroptypes
#' @title Load Croptypes
#' @description Creates a data.frame containing croptype parameters and filled with 0
#' @param params a LandsepiParams Object.
#' @param croptypeIDs a vector of indices of croptypes (must start at 0 and match with 
#' croptype IDs in the landscape)
#' @param names a vector containing the names of all croptypes
#' @details Croptypes need to be later updated with \code{\link{allocateCroptypeCultivars}}.
#' If neither croptypeIDs nor names are given, it will automatically generate
#' 1 croptype per cultivar.
#' @return a data.frame with croptype parameters
#' @seealso \link{setCroptypes}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Mixture"))
#' croptypes
#' }
#' @export
loadCroptypes <- function(params, croptypeIDs = NULL, names = NULL) {
  cultivar_names <- params@Cultivars$cultivarName
  Ncultivars <- length(cultivar_names)

  if (is.null(croptypeIDs) & is.null(names)) {
    croptypeIDs <- 0:(Ncultivars - 1)
  }

  if (is.null(croptypeIDs)) {
    Ncroptypes <- length(names)
    croptypeIDs <- 0:(Ncroptypes - 1)
  } else if (is.null(names)) {
    Ncroptypes <- length(croptypeIDs)
    names <- paste("Crop", 1:Ncroptypes)
  }


  Ncroptypes <- length(croptypeIDs)
  cropt <- list(
    croptypeID = croptypeIDs,
    croptypeName = names
  )
  prop_tmp <- rep(0, Ncultivars)
  names(prop_tmp) <- c(cultivar_names)
  cropt <- c(cropt, prop_tmp)

  cropt <- as.data.frame(cropt, stringsAsFactors = FALSE)
  return(cropt)
}


#' @name allocateCroptypeCultivars
#' @title Allocate cultivars to one croptype
#' @description Updates a given croptype by allocating cultivars composing it.
#' @param croptypes a dataframe containing all croptypes, initialised via 
#' \code{\link{loadCroptypes}}
#' @param croptypeName the name of the croptype to be allocated
#' @param cultivarsInCroptype name of cultivars composing the croptype
#' @param prop vector of proportions of each cultivar in the croptype. Default to 
#' balanced proportions.
#' @return a croptype data.frame updated for the concerned croptype.
#' @seealso \link{setCroptypes}, \link{setCultivars} 
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Mixture"))
#' croptypes
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Resistant1", "Resistant2"))
#' croptypes
#' }
#' @export
allocateCroptypeCultivars <- function(croptypes, croptypeName, cultivarsInCroptype, prop = NULL) {
  n <- length(cultivarsInCroptype) ## number of cultivars composing the croptype
  if (is.null(prop)) {
    prop <- rep(1 / n, n)
  } else if (length(prop) == 1) {
    prop <- rep(prop, n)
  }

  for (k in 1:n) {
    croptypes[croptypes$croptypeName == croptypeName, cultivarsInCroptype[k]] <- prop[k]
  }

  return(croptypes)
}



#' @name setCroptypes
#' @title Set croptypes
#' @description Updates a LandsepiParams object with croptypes and their composition with regard 
#' to cultivar proportions.
#' Note that landscape and cultivar parameters may be required if not all information is 
#' present to set croptypes.
#' @details
#' The data.frame for cultivar allocations into croptypes must take this format (example):
#'
#' | croptypeID | croptypeName  | cultivarName1 | cultivarName2 | ... |
#' | ---------- | ------------- | ------------- | ------------- | --- |
#' | 0          |  "cropt1"     |  1            | 0             | ... |
#' | 1          |  "cropt2"     |  0.5          | 0.5           | ... |
#'
#' croptypeIDs must start at 0 and match with values from landscape "croptypeID" layer with 
#' feature year_X. 
#' Cultivars names have to match cultivar names in the cultivars data.frame.
#'
#' @param params a LandsepiParams Object.
#' @param dfCroptypes a data.frame containing cultivar proportions in each croptype (see details). 
#' It can be generated manually, or initialised with \code{\link{loadCroptypes}} and later 
#' updated with \code{\link{allocateCroptypeCultivars}}.
#' @return a LandsepiParams object
#' @seealso \link{loadCroptypes}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Mixture"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Resistant1", "Resistant2"))
#' simul_params <- setCroptypes(simul_params, croptypes)
#' simul_params@Croptypes
#' }
#' @export
setCroptypes <- function(params, dfCroptypes) {

  # no croptypeID and croptypeName
  if (is.null(dfCroptypes$croptypeID) && is.null(dfCroptypes$croptypeName)) {
    warning("Can't find croptype ID or Name in the data.frame")
    return(params)
  }
  else {
    # Try to add Name
    if (!is.null(dfCroptypes$croptypeID)) {
      rownames(dfCroptypes) <- dfCroptypes$croptypeID
      if (length(params@Landscape$Name) != 0 && length(params@Landscape$ID) != 0) {
        id_name <- as.data.frame(params@Landscape[, c("ID", "Name")], stringsAsFactors = FALSE)
        id_name$geometry <- NULL
        id_name <- unique(id_name)
        id_name <- id_name[which(dfCroptypes$croptypeID == id_name$ID), "Name"]
        dfCroptypes <- data.frame(croptypeName = id_name, dfCroptypes, stringsAsFactors = FALSE)
      }
    }
    else {
      # Try to add ID
      if (length(params@Landscape$ID) != 0 && length(params@Landscape$Name) != 0) {
        id_name <- as.data.frame(params@Landscape[, c("ID", "Name")], stringsAsFactors = FALSE)
        id_name$geometry <- NULL
        id_name <- unique(id_name)
        id_name <- id_name[which(dfCroptypes$croptypeName == id_name$Name), "ID"]
        dfCroptypes <- data.frame(croptypeID = id_name, dfCroptypes, stringsAsFactors = FALSE)
        rownames(dfCroptypes) <- dfCroptypes$croptypeID
      }
      else {
        warning("Can't retrieve croptypeID from croptypeName")
      }
    }
  }

  params@Croptypes <- data.frame(dfCroptypes, stringsAsFactors = FALSE)

  checkCroptypes(params)

  return(params)
}


#' @name checkCroptypes
#' @title Check croptypes
#' @description checks croptypes validity
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkCroptypes <- function(params) {
  ret <- TRUE

  if (nrow(params@Croptypes) == 0) {
    message("Croptypes data.frame undef")
    ret <- FALSE
  }

  # check cultivars proportion by croptypes
  lcultivars <- as.matrix(params@Croptypes[, -which(.croptypesColNames 
                                                    %in% colnames(params@Croptypes))])
  ret_tmp <- apply(params@Croptypes[, -which("croptypeName" == colnames(params@Croptypes))],
    MARGIN = 1,
    FUN = function(l) {
      if (sum(as.numeric(l[-1])) != 0 && sum(as.numeric(l[-1])) != 1.0) {
        message("Croptypes ", l[1], " have a proportion of cultivars not egal to 1")
        return(FALSE)
      }
      else {
        return(TRUE)
      }
    }
  )
  if (sum(!ret_tmp) > 0) ret <- FALSE

  if (nrow(params@Landscape) > 0) {
    lc <- as.data.frame(params@Landscape)[, grep("^year_", colnames(params@Landscape))]
    if (length(lc) > 0 &
        sum(!params@Croptypes$croptypeID %in% unique(unlist(lc))) != 0) {
      ret <- FALSE
      message("croptypeID from Croptypes not found in the landscape")
    }
  }

  if (nrow(params@Cultivars) > 0) {
    if ((ncol(params@Croptypes) - length(which(.croptypesColNames 
                                               %in% colnames(params@Croptypes)))) 
        > nrow(params@Cultivars)) {
      message("Croptypes have more Cultivars than those defined in Cultivars data.frame")
      ret <- FALSE
    }
    if (length(which(colnames(params@Croptypes)[-which(.croptypesColNames 
                                                       %in% colnames(params@Croptypes))] 
                     %in% rownames(params@Cultivars)))
    != (ncol(params@Croptypes) - length(which(.croptypesColNames 
                                              %in% colnames(params@Croptypes))))) {
      message("Cultivars in Croptypes data.frame are undef in Cultivars data.frame")
      ret <- FALSE
    }
  }

  return(ret)
}


#' @name loadCultivar
#' @title Load a cultivar
#' @description create a data.frame containing cultivar parameters depending of his type
#' @param name a character string (without space) specifying the cultivar name.
#' @param type the cultivar type, among: "growingHost" (default), "nongrowingHost", 
#' "grapevine", "banana" or "nonCrop".
#' @details 
#' * "growingHost" is adapted to situations where the infection unit is a piece of leaf 
#' (e.g. where a fungal lesion can develop); the number of available infection units 
#' increasing during the season due to plant growth (as typified by cereal crops). 
#' * "nongrowingHost" corresponds to situations where the infection unit is the whole plant 
#' (e.g. for viral systemic infection); thus the number of infection units is constant. 
#' * "grapevine" corresponds to parameters for grapevine (including host growth).
#' * "banana" corresponds to parameters for banana (including host growth).
#' * "nonCrop" is not planted, does not cost anything and does not yield anything 
#' (e.g. forest, fallow).
#' @return a dataframe of parameters associated with each host genotype 
#' (i.e. cultivars, lines) when cultivated in pure crops.
#' @seealso \link{setCultivars}
#' @include Cultivars_List.R
#' @examples
#' c1 <- loadCultivar("winterWheat", type = "growingHost")
#' c1
#' c2 <- loadCultivar("forest", type = "nonCrop")
#' c2
#' @export
loadCultivar <- function(name, type = "growingHost") {
  culti <- Cultivars_list[[type]]
  culti["cultivarName"] <- name

  culti <- as.data.frame(culti, stringsAsFactors = FALSE)
  if (length(culti) <= 1) {
    warning('Unknown type of host: "', type
            , '". Possible types are: 
            "growingHost", "nongrowingHost", "grapevine", "banana", "nonCrop"')
  } else {
    # To be sure of the columns names
    colnames(culti) <- .cultivarsColNames
  }

  return(culti)
}


#' @name setCultivars
#' @title Set cultivars
#' @description Updates a LandsepiParams object with cultivars parameters
#' @param params a landsepiParams object.
#' @param dfCultivars a data.frame defining the cultivars (see details). It can be generated 
#' manually or, alternatively, via \code{\link{loadCultivar}}.
#'
#' @details dfCultivars is a dataframe of parameters associated with each host genotype 
#' (i.e. cultivars, lines) when cultivated in pure crops. Columns of the dataframe are:\itemize{
#' \item cultivarName: cultivar names (cannot accept space),
#' \item initial_density: host densities (per square meter) at the beginning of the cropping season 
#' as if cultivated in pure crop,
#' \item max_density: maximum host densities (per square meter) at the end of the cropping season 
#' as if cultivated in pure crop,
#' \item growth rate: host growth rates,
#' \item reproduction rate: host reproduction rates,
#' \item yield_H: theoretical yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status H as if cultivated in pure crop,
#' \item yield_L: theoretical yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status L as if cultivated in pure crop,
#' \item yield_I: theoretical yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status I as if cultivated in pure crop,
#' \item yield_R: theoretical yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status R as if cultivated in pure crop,
#' \item planting_cost = planting costs (in monetary units / ha / cropping season) as if 
#' cultivated in pure crop,
#' \item market_value = market values of the production (in monetary units / weight or volume unit).
#' }
#' 
#' The data.frame must be defined as follow (example):
#' 
#' | cultivarName | initial_density | max_density | growth_rate | reproduction_rate | yield_H | yield_L | yield_I |yield_R | planting_cost | market_value |
#' | ------------ | --------------- | ----------- | ----------- | ----------------- | ------- | ------- | ------- | ------ | ------------- | ------------ |
#' | Susceptible  | 0.1             |  2.0        | 0.1         | 0.0               | 2.5     | 0.0     | 0.0     | 0.0    | 225           | 200          |
#' | Resistant1   | 0.1             |  2.0        | 0.1         | 0.0               | 2.5     | 0.0     | 0.0     | 0.0    | 225           | 200          |
#' | Resistant2   | 0.1             |  2.0        | 0.1         | 0.0               | 2.5     | 0.0     | 0.0     | 0.0    | 225           | 200          |
#'
#' @return a LandsepiParams object
#' @seealso \link{loadCultivar}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params@Cultivars
#' }
#' @export
setCultivars <- function(params, dfCultivars) {
  if (!is.null(dfCultivars$cultivarName)) {
    rownames(dfCultivars) <- dfCultivars$cultivarName
  }

  params@Cultivars <- data.frame(dfCultivars[, .cultivarsColNames], stringsAsFactors = FALSE)

  checkCultivars(params)

  return(params)
}


#' @name checkCultivars
#' @title Check cultivars
#' @description check cultivars validity
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkCultivars <- function(params) {

  ret <- TRUE

  if (is.null(params@Cultivars) || nrow(params@Cultivars) == 0) {
    warning("Cultivars is NULL, use setCultivars()")
    ret <- FALSE
    return(ret)
  }

  if (sum(.cultivarsColNames %in% colnames(params@Cultivars)) != length(.cultivarsColNames)) {
    warning("Missing columns in Cultivars data.frame : ", .cultivarsColNames)
    ret <- FALSE
    return(ret)
  }
  
  if (!is.character(params@Cultivars$cultivarName) ||
      sum(grepl(" ", params@Cultivars$cultivarName)) > 0){
    warning("Cultivar names must be character strings without spaces")
    ret <- FALSE
  }
  
  if (!is.numeric(params@Cultivars$growth_rate) ||
      !is.numeric(params@Cultivars$reproduction_rate) ||
      # !is.numeric(params@Cultivars$death_rate) ||
      
      sum(!is.in.01(params@Cultivars$growth_rate) > 0) || 
      sum(!is.in.01(params@Cultivars$reproduction_rate) > 0) 
      # || sum(!is.in.01(params@Cultivars$death_rate) > 0) 
      ){
    warning("growth and reproduction rates must be between 0 and 1")
    ret <- FALSE
  } 
  
  if(!is.numeric(params@Cultivars$initial_density) ||
     !is.numeric(params@Cultivars$yield_H) ||
     !is.numeric(params@Cultivars$planting_cost) ||
     !is.numeric(params@Cultivars$market_value) ||
     
     sum( !is.positive(params@Cultivars$initial_density) ) > 0 ||
     sum( !is.positive(params@Cultivars$yield_H) ) > 0 ||
     sum( !is.positive( params@Cultivars$planting_cost) ) > 0 ||
     sum( !is.positive(params@Cultivars$market_value) ) > 0 ){
    warning("initial_density, yield_H, planting_cost and market_value must be >= 0")
    ret <- FALSE
  }
  
  if(!is.numeric(params@Cultivars$yield_L) ||
     !is.numeric(params@Cultivars$yield_I) ||
     !is.numeric(params@Cultivars$yield_R) ||

     sum( !is.positive(params@Cultivars$yield_L) ) > 0 ||
     sum( !is.positive(params@Cultivars$yield_I) ) > 0 ||
     sum( !is.positive(params@Cultivars$yield_R) ) > 0 ||

     sum( params@Cultivars$yield_L > params@Cultivars$yield_H ) > 0 ||
     sum( params@Cultivars$yield_I > params@Cultivars$yield_H ) > 0 ||
     sum( params@Cultivars$yield_R > params@Cultivars$yield_H ) > 0 ){

    warning("yield_L, yield_I and yield_R must be >= 0 and smaller or equal to yield_H")
    ret <- FALSE
  }
  
  if(!is.numeric(params@Cultivars$max_density) ||
     sum( !is.strict.positive(params@Cultivars$max_density) ) > 0 ||
     sum( params@Cultivars$max_density < params@Cultivars$initial_density ) > 0 ){
    warning("Maximal density must be strictly positive and greater or equal to initial density")
    ret <- FALSE
  }
  

  if (ncol(params@Croptypes) > 0) {
    if (nrow(params@Cultivars) 
        < (length(params@Croptypes[, -which(.croptypesColNames 
                                            %in% colnames(params@Croptypes))]) - 1)) {
      warning("Cultivars number is less than thoses defined in Croptypes data.frame")
      ret <- FALSE
    }
    if (length(which(rownames(params@Cultivars) 
                     %in% colnames(params@Croptypes)[
                       -which(.croptypesColNames %in% colnames(params@Croptypes))])) 
        != (ncol(params@Croptypes) - length(which(.croptypesColNames
                                                  %in% colnames(params@Croptypes))))) {
      warning("Cultivars are undef in Cultivars data.frame compared to croptypes")
      ret <- FALSE
    }
  }

  return(ret)
}


#' @name loadGene
#' @title Load a gene
#' @description Creates a data.frame containing parameters of a gene depending of his type
#' @param name name of the gene
#' @param type type of the gene: "majorGene", "APR", "QTL" or "immunity" (default = "majorGene")
#' @details 
#' * "majorGene" means a completely efficient gene that can be broken down via a single 
#' pathogen mutation
#' * "APR" means a major gene that is active only after a delay of 30 days after planting
#' * "QTL" means a partial resistance (50% efficiency) that requires several pathogen mutations 
#' to be completely eroded
#' * "immunity" means a completely efficient resistance that the pathogen has no way to adapt 
#' (i.e. the cultivar is non-host).  
#' 
#' For different scenarios, the data.frame can be manually updated later.
#' @return a data.frame with gene parameters
#' @seealso \link{setGenes}
#' @examples
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene1
#' gene2 <- loadGene(name = "Lr34", type = "APR")
#' gene2
#' @export
loadGene <- function(name, type = "majorGene") {
  gene <- switch(type,
    "majorGene" = list(
      "geneName" = name,
      "efficiency" = 1.0,
      "age_of_activ_mean" = 0.0,
      "age_of_activ_var" = 0.0,
      "mutation_prob" = 0.0000001,
      "Nlevels_aggressiveness" = 2,
      "adaptation_cost" = 0.5,
      "relative_advantage" = 0.5,
      "tradeoff_strength" = 1.0,
      "target_trait" = "IR",
      "recombination_sd" = 1.0
    ),
    "APR" = list(
      "geneName" = name,
      "efficiency" = 1.0,
      "age_of_activ_mean" = 30.0,
      "age_of_activ_var" = 30.0,
      "mutation_prob" = 0.0000001,
      "Nlevels_aggressiveness" = 2,
      "adaptation_cost" = 0.5,
      "relative_advantage" = 0.5,
      "tradeoff_strength" = 1.0,
      "target_trait" = "IR",
      "recombination_sd" = 1.0
    ),
    "QTL" = list(
      "geneName" = name,
      "efficiency" = 0.5,
      "age_of_activ_mean" = 0.0,
      "age_of_activ_var" = 0.0,
      "mutation_prob" = 0.0001,
      "Nlevels_aggressiveness" = 6,
      "adaptation_cost" = 0.5,
      "relative_advantage" = 0.5,
      "tradeoff_strength" = 1.0,
      "target_trait" = "IR",
      "recombination_sd" = 0.27
    ),
    "immunity" = list(
      "geneName" = name,
      "efficiency" = 1.0,
      "age_of_activ_mean" = 0.0,
      "age_of_activ_var" = 0.0,
      "mutation_prob" = 0,
      "Nlevels_aggressiveness" = 1,
      "adaptation_cost" = 0,
      "relative_advantage" = 0,
      "tradeoff_strength" = 1,
      "target_trait" = "IR",
      "recombination_sd" = 1.0
    ),
    list()
  )

  gene <- as.data.frame(gene, stringsAsFactors = FALSE)
  if (length(gene) == 0) {
    warning('Unknown type of gene: "', type
            , '". Possible types are: "majorGene", "APR", "QTL", "immunity')
  } else {
    # To be sure of the columns names
    colnames(gene) <- .geneColNames
  }
  return(gene)
}


#' @name setGenes
#' @title Set genes
#' @description Updates a LandsepiParams object with parameters associated with resistance genes
#' and pathogen adaptation.
#' @details dfGenes is a data.frame of parameters associated with each resistance gene and 
#' with the evolution of each corresponding pathogenicity gene. Columns of the dataframe are:
#' \itemize{
#' \item geneName: names of resistance genes,
#' \item target_trait: aggressiveness components ("IR", "LAT", "IP", or "PR") targeted by 
#' resistance genes,
#' \item efficiency: resistance gene efficiencies, i.e. the percentage of reduction of the targeted 
#' aggressiveness component (IR, 1/LAT, IP and PR),
#' \item age_of_activ_mean: expected delays to resistance activation (for APRs),
#' \item age_of_activ_var: variances of the delay to resistance activation (for APRs),
#' \item mutation_prob: mutation probabilities for pathogenicity genes (each of them 
#' corresponding to a resistance gene),
#' \item Nlevels_aggressiveness: number of adaptation levels related to each resistance gene 
#' (i.e. 1 + number of required mutations for a pathogenicity gene to fully adapt to the 
#' corresponding resistance gene),
#' \item adaptation_cost: fitness penalties paid by pathogen genotypes 
#' fully adapted to the considered resistance genes on all hosts, 
#' \item relative_advantage: fitness advantages of pathogen genotypes fully adapted to the 
#' resistance genes on hosts carrying these genes, relative to those that do not carry these genes,
#' \item tradeoff_strength: strengths of the trade-off relationships between the 
#' level of aggressiveness on hosts that do and do not carry the resistance genes.
#' \item recombination_sd: standard deviation of the normal distribution used for recombination 
#' of quantitative traits during sexual reproduction (infinitesimal model)
#' }
#' 
#' The data.frame must be defined as follow (example):
#'
#' | geneName | efficiency | age_of_activ_mean | age_of_activ_var | mutation_prob | Nlevels_agressiveness | adaptation_cost | relative advantage | tradeoff_strength | target_trait | recombination_sd |
#' | -------- | ---------- | ----------------- | ----------------- | ------------- | --------------------- | ------------ | --------------- | ----------------- | ------------ | ------------------------ |
#' | MG1      |  1         |  0                | 0                 | 1e-07         | 2                     | 0.5          | 0.5             | 1                 | IR           | 0.27                     |
#' | QTL1     | 0.5        |  0                | 0                 | 0.0001        | 10                    | 0.74         | 0.74            | 1                 | LAT          | 0.27                     |
#'
#' @param params a LandsepiParams object
#' @param dfGenes a data.frame containing gene parameters. It can be defined manually, or, 
#' alternatively, with \code{\link{loadGene}}.
#' @return a LandsepiParams object.
#' @seealso \link{loadGene}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' simul_params@Genes
#' }
#' @export
setGenes <- function(params, dfGenes) {
  if (!is.null(dfGenes$geneName)) {
    rownames(dfGenes) <- dfGenes$geneName
  }
  params@Genes <- dfGenes[, .geneColNames]

  checkGenes(params)

  return(params)
}


#' @name checkGenes
#' @title Check genes
#' @description checks Genes data.frame validity
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkGenes <- function(params) {
  ret <- TRUE

  if (nrow(params@Genes) == 0) {
    warning("Simulation with no resistance gene")
    return(ret)
  }

  if (is.null(params@Genes$geneName)) warning("missing 'geneName' column into genes data.frame")

  if (sum(!.geneColNames %in% colnames(params@Genes)) > 0) {
    warning("Genes data.frame column(s) missing")
    warning("Genes colnames are ", .geneColNames)
    ret <- FALSE
    return(ret)
  }
  
  validTraits <- c("IR","LAT","PR","IP")
  if(!is.character(params@Genes$target_trait) ||
     is.na( sum(match(params@Genes$target_trait, validTraits)) ) ){
    warning( "Error: valid target traits are:", paste(validTraits, collapse = ", ") )
    ret <- FALSE
  }

  if (!is.numeric(params@Genes$efficiency) ||
      !is.numeric(params@Genes$mutation_prob) ||
      !is.numeric(params@Genes$adaptation_cost) ||
      !is.numeric(params@Genes$relative_advantage) ||
      
      sum(!is.in.01(params@Genes$efficiency) > 0) || 
      sum(!is.in.01(params@Genes$mutation_prob) > 0) || 
      sum(!is.in.01(params@Genes$adaptation_cost) > 0) ||
      sum(!is.in.01(params@Genes$relative_advantage) > 0) ){
    warning("efficiencies, mutation probabilities, adaptation costs and relative advantage 
            must be between 0 and 1")
    ret <- FALSE
  } 
  
  if(!is.numeric(params@Genes$age_of_activ_mean) ||
     !is.numeric(params@Genes$age_of_activ_var) ||
     
     sum(!is.positive(params@Genes$age_of_activ_mean) > 0) ||
     sum(!is.positive(params@Genes$age_of_activ_var) > 0) ){
    warning("Expectation and variance of the times to resistance activation must be >= 0")
    ret <- FALSE
  }
  
  if(!is.numeric(params@Genes$tradeoff_strength) ||
     sum(!is.strict.positive(params@Genes$tradeoff_strength) > 0) ){
    warning("tradeoff strengths must be > 0")
    ret <- FALSE
  }
  
  if(!is.numeric(params@Genes$Nlevels_aggressiveness) ||
     sum( !is.wholenumber(params@Genes$Nlevels_aggressiveness) > 0) ||
     sum( !is.strict.positive(params@Genes$Nlevels_aggressiveness) > 0) ){
    warning("Number of levels of aggressiveness must be a whole number >= 1")
    ret <- FALSE
  }
  
  if (!is.numeric(params@Genes$recombination_sd) ||
      sum(!is.strict.positive(params@Genes$recombination_sd) > 0) ){
    warning("recombination_sd must be > 0")
    ret <- FALSE
  } 
  
  return(ret)
}


#' @name allocateCultivarGenes
#' @title Allocate genes to a cultivar
#' @description Updates a LandsepiParams object with, for a given cultivar, the list of genes 
#' it carries
#' @param params a LandsepiParams object.
#' @param cultivarName the name of the cultivar to be allocated.
#' @param listGenesNames the names of the genes the cultivar carries
#' @param force.clean force to clean previous allocated genes to all cultivars
#' @return a LandsepiParams object
#' @seealso \link{setGenes}, \link{setCultivars}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant", c("MG 1", "MG 2"))
#' simul_params@CultivarsGenes
#' }
#' @export
allocateCultivarGenes <- function(params, cultivarName, listGenesNames = c("")
                                  , force.clean=FALSE) {
  if (isTRUE(force.clean) || length(params@CultivarsGenes) == 0 
      || nrow(params@CultivarsGenes) != nrow(params@Cultivars) 
      || nrow(params@Genes) != ncol(params@CultivarsGenes)) {
    params@CultivarsGenes <- data.frame(matrix(rep(0, nrow(params@Genes) * nrow(params@Cultivars))
                                               , nrow = nrow(params@Cultivars))
                                        , row.names = rownames(params@Cultivars))
    colnames(params@CultivarsGenes) <- params@Genes$geneName
  }

  if (cultivarName %in% params@Cultivars$cultivarName &&
    sum(listGenesNames %in% params@Genes$geneName) == length(listGenesNames)) {
    params@CultivarsGenes[which(rownames(params@CultivarsGenes) 
                                == cultivarName), listGenesNames] <- 1
    #params@CultivarsGenes[which(rownames(params@CultivarsGenes) 
    #                            != cultivarName), listGenesNames] <- 0
  } else {
    stop("Can't find cultivarName or geneName from data.frame")
  }
  return(params)
}


#' @name resetCultivarsGenes
#' @title Reset cultivars genes
#' @description Resets the lists of genes carried by all cultivars
#' @param params a LandsepiParams object.
#' @return a LandsepiParams object
resetCultivarsGenes <- function(params) {
  ## (used in shinyApp)
  params@CultivarsGenes <- data.frame()
  return(params)
}


#' @name checkCultivarsGenes
#' @title Check cultivars genes
#' @description Checks CultivarsGene data.frame validity
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkCultivarsGenes <- function(params) {
  ret <- TRUE
  if (length(params@CultivarsGenes) > 0 & 
      ( nrow(params@CultivarsGenes) != nrow(params@Cultivars) || 
       nrow(params@Genes) != ncol(params@CultivarsGenes)) ) {
    warning("Cultivars Genes undef (some genes are not allocated to any cultivar)")
    ret <- FALSE
  }
  return(ret)
}


#' @name loadInoculum
#' @title Load Inoculum
#' @description Loads an inoculum for the beginning of the simulation (t=0), with 
#' controlled localisation (polygons), infected cultivars and pathogen genotypes.
#' Note that landscape, gene, cultivar and croptype parameters must be set before 
#' loading the inoculum.
#' @param params a LandsepiParams object.
#' @param pI0_all a numeric indicating the (same) probability to infect a host for all 
#' pathogen genotypes, all cultivars and in all polygons
#' @param pI0_host a vector of length Nhost indicating the probabilities to infect an host, 
#' for each cultivar (for all pathogen genotypes and all polygons).
#' @param pI0_patho a vector of length Npatho indicating the probabilities to infect an host, 
#' for each pathogen genotype (for all cultivars and all polygons).
#' @param pI0_poly a vector of length Npoly indicating the probabilities to infect an host, 
#' for each polygon (for all pathogen genotypes and all cultivars).
#' @param pI0_mat a 3D array of dimensions (1:Nhost,1:Npatho,1:Npoly) indicating the 
#' probability to infect an host, for each cultivar, pathogen genotype and polygon. 
#' Note that \code{pI0_all}, \code{pI0_host}, \code{pI0_patho} and \code{pI0_poly} 
#' are not accounted if \code{pI0_mat} is filled.
#' @details The different options enable different types of inoculum (localisation, 
#' infected cultivars and pathogen genetic diversity, 
#' see different options in Examples).\cr
#' Unless the array \code{pI0_mat} is filled, the probability for a host to be infected 
#' at the beginning of the simulation is computed in every polygon (poly), cultivar (host) 
#' and pathogen genotype (patho) with 
#' \code{pI0[host, patho, poly] = pI0_all * pI0_patho[patho] * pI0_host[host] * pI0_poly[poly]}. \cr
#' Before loading the inoculum, one can use 
#' \code{getMatrixGenePatho()}, \code{getMatrixCultivarPatho()} and \code{getMatrixCroptypePatho()} 
#' to acknowledge which pathogen genotypes are adapted to which genes, cultivars and croptypes.\cr
#' Once \code{setInoculum()} is used, one can call \code{inoculumToMatrix()} to get the inoculum 
#' as a 3D array (1:Nhost,1:Npatho,1:Npoly)\cr
#' @return a 3D array of dimensions (1:Nhost,1:Npatho,1:Npoly)
#' @seealso \link{inoculumToMatrix}, \link{getMatrixGenePatho}, \link{getMatrixCultivarPatho}, 
#' \link{getMatrixCroptypePatho}, \link{setInoculum}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setTime(simul_params, Nyears = 1, nTSpY = 80)
#' basic_patho_param <- loadPathogen(disease = "rust")
#' simul_params <- setPathogen(simul_params, patho_params = basic_patho_param)
#' simul_params <- setLandscape(simul_params, loadLandscape(id = 1))
#' simul_params <- setDispersalPathogen(simul_params, loadDispersalPathogen(id = 1)[[1]])
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant", c("MG 1", "MG 2"))
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Resistant crop"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop", c("Resistant"))
#' simul_params <- setCroptypes(simul_params, croptypes)
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = 0
#' , rotation_sequence = croptypes$croptypeID
#' , prop = c(1/2,1/2), aggreg = 1, graphic = FALSE)
#' 
#' #### Definition of the inoculum ####
#' 
#' ### Scenario 1. Only the avirulent pathogen on the susceptible cultivar ###
#' # In this situation, the susceptible cultivar must be entered
#' # at the first line of the table cultivars
#' 
#' ## Global inoculum (i.e. in the whole landscape)
#' # Option 1: simply use the default parameterisation
#' simul_params <- setInoculum(simul_params, 5E-4)
#' 
#' # Option 2: use loadInoculum()
#' Npatho <- prod(simul_params@Genes$Nlevels_aggressiveness)
#' Nhost <- nrow(simul_params@Cultivars)
#' pI0 <- loadInoculum(simul_params,
#'                     pI0_all=5E-4,
#'                     pI0_host=c(1,rep(0, Nhost-1)),
#'                     pI0_patho=c(1,rep(0, Npatho-1)))
#' simul_params <- setInoculum(simul_params, pI0)
#' inoculumToMatrix(simul_params)
#' 
#' ## Local inoculum (i.e. in some random polygons only)
#' Npatho <- prod(simul_params@Genes$Nlevels_aggressiveness)
#' Nhost <- nrow(simul_params@Cultivars)
#' Npoly <- nrow(simul_params@Landscape)
#' Npoly_inoc <- 5  ## number of inoculated polygons
#' ## whether the avr pathogen can infect the polygons
#' compatible_poly <- getMatrixPolyPatho(simul_params)[,1]
#' ## random polygon picked among compatible ones
#' id_poly <- sample(grep(1, compatible_poly), Npoly_inoc)
#' pI0_poly <- as.numeric(1:Npoly %in% id_poly)  
#' pI0 <- loadInoculum(simul_params,
#'                     pI0_all=5E-4,
#'                     pI0_host=c(1,rep(0, Nhost-1)),
#'                     pI0_patho=c(1,rep(0, Npatho-1)), 
#' pI0_poly=pI0_poly)
#' simul_params <- setInoculum(simul_params, pI0)
#' inoculumToMatrix(simul_params)
#' 
#' ### Scenario 2. Diversity of pathogen genotypes in the inoculum ###
#' # in this example,  Nhost=2 cultivars, Npatho=4
#' 
#' ## Global inoculum (i.e. in all polygons of the landscape)
#' pI0 <- loadInoculum(simul_params, pI0_patho=c(1E-3,1E-4,1E-4,1E-5), pI0_host=c(1,1))
#' simul_params <- setInoculum(simul_params, pI0)
#' inoculumToMatrix(simul_params)[,,1:5]
#' 
#' ## Local inoculum (i.e. in some polygons only) ##
#' Npoly <- nrow(simul_params@Landscape)
#' Npoly_inoc <- 5  ## number of inoculated polygons 
#' id_poly <- sample(1:Npoly, Npoly_inoc)  ## random polygon 
#' pI0_poly <- as.numeric(1:Npoly %in% id_poly) 
#' pI0 <- loadInoculum(simul_params, pI0_patho=c(1E-3,1E-4,1E-4,1E-5),
#' pI0_host=c(1,1), pI0_poly=pI0_poly)
#' simul_params <- setInoculum(simul_params, pI0)
#' inoculumToMatrix(simul_params)
#' }
#' @export
loadInoculum <- function(params, pI0_all=NULL, pI0_host=NULL, pI0_patho=NULL, pI0_poly=NULL
                         , pI0_mat=NULL){ 
  Nhost <- nrow(params@Cultivars)
  Npatho <-prod(params@Genes$Nlevels_aggressiveness)
  Npoly <- nrow(params@Landscape)
  mat_croptypes <- params@Croptypes
  cultivarNames <- params@Cultivars$cultivarName
  rotation_0 <- params@Landscape$year_1
  mat_cultivar_patho <- getMatrixCultivarPatho(params)

  ## Security check
  if (Nhost==0 | Npoly==0) {
    stop("Please define landscape, cultivars and pathogen first")
  }
  if (all(missing(pI0_mat), missing(pI0_all), missing(pI0_patho), missing(pI0_host)
          , missing(pI0_poly)))
    stop("Missing argument: The probability of infection must be defined")

  if (!is.null(pI0_mat)){
    if(any(dim(pI0_mat) != c(Nhost, Npatho, Npoly))){
      stop("'pI0_mat' must be of dimensions (Nhost, Npatho, Npoly)")
    }else{
      if (any(!is.null(pI0_all), !is.null(pI0_patho), !is.null(pI0_host), !is.null(pI0_poly)))
        warning("'pI0_all', 'pI0_patho', 'pI0_host' and 'pI0_patho' are not accounted if 
                'pI0_mat' is filled")
      pI0 <- checkPI0_mat(pI0_mat, params)
      
    }
  } else {
    
    if (is.null(pI0_all)){
      pI0_all <- 1
    }else{
      if (length(pI0_all) != 1)
        stop("'pI0_all' must be of length 1")
    }
    if (is.null(pI0_patho)){
      pI0_patho <- rep(1, Npatho)
    }else{
      if (length(pI0_patho) != Npatho)
        stop("'pI0_patho' must have the same length as the number of pathogen genotypes")
    }
    if (is.null(pI0_host)){
      pI0_host <- rep(1, Nhost)
    }else{
      if (length(pI0_host) != Nhost)
        stop("'pI0_host' must have the same length as the number of cultivars")
    } 
    if (is.null(pI0_poly)){
      pI0_poly <- rep(1, Npoly)
    }else{
      if (length(pI0_poly) != Npoly)
        stop("'pI0_poly' must have the same length as the number of polygons")
    } 
    
    ## Computation of pI0
    pI0 <- array(data = 0, dim = c(Nhost, Npatho, Npoly))
    for (host in 1:Nhost){
      for (patho in 1:Npatho){
        for (poly in 1:Npoly){
          id_croptype <- rotation_0[poly]+1  ## +1 because of C code
          host_present <- as.numeric(mat_croptypes[id_croptype, cultivarNames[host]] > 0)  
          ## i.e. is the host present in the croptype of poly ?
          host_compatible <- as.numeric(mat_cultivar_patho[host, patho])  
          ## i.e. is the patho able to infect the host ?
          
          pI0[host, patho, poly] <- host_present * host_compatible * 
            pI0_all * pI0_host[host] * pI0_patho[patho] * pI0_poly[poly]
        }
      }
    }
  } ## else pI0_mat is null
  
    return(pI0)
}


#' @name setInoculum
#' @title Set inoculum
#' @description Updates a LandsepiParams object with the initial probability for an individual host 
#' to be infectious (i.e. state I) at the beginning of the simulation (i.e. t=0). 
#' @param params a LandsepiParams object.
#' @param val a numeric value (default = 5e-4) indicating the probability for the first cultivar 
#' to be infected by the first pathogen genotype in all polygons of the landscape 
#' (must be between 0 and 1). 
#' The parameter can also be entered as a 3D array of dimensions (1:Nhost,1:Npatho,1:Npoly) 
#' indicating the initial probability to be infectious, for each cultivar, pathogen genotype and 
#' polygon (independently from the possible presence of cultivars carrying resistance genes). 
#' It can be generated manually or, alternatively, via \code{\link{loadInoculum}}.
#' @details Before setting the inoculum, one can use \code{getMatrixGenePatho()}, 
#' \code{getMatrixCultivarPatho()}, 
#' \code{getMatrixCroptypePatho()} and \code{getMatrixPolyPatho()} to acknowledge which 
#' pathogen genotypes are compatible to which genes, cultivars, croptypes and polygons.\cr
#' Once \code{setInoculum()} is used, one can call \code{inoculumToMatrix()} to get 
#' the inoculum as a 3D array (1:Nhost,1:Npatho,1:Npoly)\cr
#' @return a LandsepiParams object
#' @seealso \link{inoculumToMatrix}, \link{loadInoculum}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setInoculum(simul_params, 1E-3)
#' simul_params@PI0
#' }
#' @export
setInoculum <- function(params, val = 5e-4) {
  params@PI0 <- as.vector(val)
  checkInoculum(params)
  
  return(params)
}


#' @name inoculumToMatrix
#' @title Inoculum To Matrix
#' @description Transform the inoculum pI0 (1D vector of length Nhost*Npatho*Npoly) into 
#' a 3D array (for visualization purpose)
#' @details After defining the inoculum with \code{setInoculum()}, this function returns 
#' the inoculum as a 3D array. 
#' @param params a LandsepiParams object.
#' @return a 3D array of structure (1:Nhost,1:Npatho,1:Npoly)
#' @seealso \link{setInoculum}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setTime(simul_params, Nyears = 1, nTSpY = 80)
#' simul_params <- setPathogen(simul_params, loadPathogen(disease = "rust"))
#' simul_params <- setLandscape(simul_params, loadLandscape(id = 1))
#' simul_params <- setDispersalPathogen(simul_params, loadDispersalPathogen(id = 1)[[1]])
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant", c("MG 1", "MG 2"))
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Resistant crop"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop", c("Resistant"))
#' simul_params <- setCroptypes(simul_params, croptypes)
#' simul_params@Croptypes
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = 0
#' , rotation_sequence = croptypes$croptypeID
#' , prop = c(1/2,1/2), aggreg = 1, graphic = FALSE)
#' pI0 <- loadInoculum(simul_params, pI0_patho=c(1E-3,1E-4,1E-4,1E-5), pI0_host=c(1,1))
#' simul_params <- setInoculum(simul_params, pI0)
#' inoculumToMatrix(simul_params)[,,1:5]
#' }
#' @export
inoculumToMatrix <- function(params){ 
  pI0_vect <- params@PI0
  Nhost <- nrow(params@Cultivars)
  Npatho <- prod(params@Genes$Nlevels_aggressiveness)
  Npoly <- nrow(params@Landscape)
  mat_cultivar_patho <- getMatrixCultivarPatho(params)
  pI0 <- array(data = 0, dim = c(Nhost, Npatho, Npoly))
  i=1
  for(poly in 1:Npoly)  {
    for(patho in 1:Npatho)    {
      for(host in 1:Nhost)      {
        pI0[host, patho, poly] <- pI0_vect[i]
        i=i+1
      }
    }
    rownames(pI0) <- rownames(mat_cultivar_patho)
  }
  return(pI0)
}


#' @name getMatrixGenePatho
#' @title Get the "resistance gene/pathogen genotype" compatibility matrix.
#' @description Build the matrix indicating if infection is possible at the beginning of the season 
#' for every combination of plant resistance gene (rows) and pathogen genotype (columns).
#' @details For hosts carrying each resistance gene, there is either possibility of infection 
#' by the pathogen genotype (value of 1), either complete protection (value of 0). 
#' Complete protection only occurs if the resistance gene targets the infection rate, 
#' has a complete efficiency, and is expressed from the beginning of the cropping season 
#' (i.e. this is not an APR).  
#' @param params a LandsepiParams object.
#' @return an interaction matrix composed of 0 and 1 values.
#' @seealso \link{getMatrixCultivarPatho}, \link{getMatrixCroptypePatho}, \link{getMatrixPolyPatho}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' getMatrixGenePatho(simul_params)
#' }
#' @export
getMatrixGenePatho <- function(params){
  Ngenes <- length(params@Genes$geneName)
  Nlevels_aggressiveness <- params@Genes$Nlevels_aggressiveness
  Npatho <- prod(Nlevels_aggressiveness)
  if(Ngenes==0){
    message("No genes are set")
  }
  mat_gene_patho <- array(data = 0, dim = c(Ngenes, Npatho))
  for (patho in 1:Npatho) {
    mat_gene_patho[,patho] <- switch_patho_to_aggr(patho-1, Ngenes, Nlevels_aggressiveness)  
    ## (function in output.R)
  }
  ## Infection is possible for:
  indices_infection <- which(params@Genes$age_of_activ_mean != 0 | 
                               ## Adult Plant Resistance (delayed)
                             params@Genes$efficiency < 1 |   ## Partial resistance
                             params@Genes$target_trait != "IR")  ## IR is not targeted 
  mat_gene_patho[indices_infection,] <- 1
  rownames(mat_gene_patho) <- params@Genes$geneName
  return(mat_gene_patho)
}


#' @name getMatrixCultivarPatho
#' @title Get the "cultivar/pathogen genotype" compatibility matrix.
#' @description Build the matrix indicating if infection is possible at the beginning of the season 
#' for every combination of cultivar (rows) and pathogen genotype (columns).
#' @details For each cultivar, there is either possibility of infection by the 
#' pathogen genotype (value of 1), or complete protection (value of 0).
#' @param params a LandsepiParams object.
#' @return an interaction matrix composed of 0 and 1 values.
#' @seealso \link{getMatrixGenePatho}, \link{getMatrixCroptypePatho}, \link{getMatrixPolyPatho}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "monoResistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "monoResistant2", type = "growingHost")
#' cultivar4 <- loadCultivar(name = "Pyramid", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3, cultivar4)
#' , stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params <- allocateCultivarGenes(simul_params, "monoResistant1", c("MG 1"))
#' simul_params <- allocateCultivarGenes(simul_params, "monoResistant2", c("MG 2"))
#' simul_params <- allocateCultivarGenes(simul_params, "Pyramid", c("MG 1", "MG 2"))
#' getMatrixCultivarPatho(simul_params)
#' }
#' @export
getMatrixCultivarPatho <- function(params){
  Nhost <- nrow(params@Cultivars)
  Npatho <- prod(params@Genes$Nlevels_aggressiveness)
  Ngenes <- length(params@Genes$geneName)
  mat_cultivar_patho <- array(data = 0, dim = c(Nhost, Npatho))
  rownames(mat_cultivar_patho) <- params@Cultivars$cultivarName
  mat_cultivar_gene <- params@CultivarsGenes
  
  if(Nhost==0){
    stop("No cultivars are set")
  }else {
    
    if(Ngenes==0 | length(mat_cultivar_gene)==0){
      message("No genes are set, or they are not allocated")
      mat_cultivar_patho[,] <- 1
    }else{
      mat_gene_patho <- getMatrixGenePatho(params)
      for (host in 1:Nhost){
        indices_genes <- which(mat_cultivar_gene[host,] == 1)  
        ## indices of genes carried by the host
        if (ncol(mat_gene_patho)==1){
          mat_cultivar_patho[host,] <- prod(mat_gene_patho[indices_genes,])
        }else if (length(indices_genes)==1){
          
          mat_cultivar_patho[host,] <- mat_gene_patho[indices_genes,]
        }else{  ## (in case length(indices_genes)==0, the result is 1)
          mat_cultivar_patho[host,] <- apply(mat_gene_patho[indices_genes,], 2, prod)  
          ## product of the compatibilities for each carried gene
        }
      } ## for host
      
      indices_nonCrop <- which(params@Cultivars$initial_density == 0)  
      ## indices of cultivars that are not planted
      mat_cultivar_patho[indices_nonCrop,] <- 0
      
    } ## else Nhost > 0
    
    return (mat_cultivar_patho)
  }
}

#' @name getMatrixCroptypePatho
#' @title Get the "croptype/pathogen genotype" compatibility matrix.
#' @description Build the matrix indicating if infection is possible at the beginning of the season 
#' for every combination of croptype (rows) and pathogen genotype (columns).
#' @details For each croptype, there is either possibility of infection by the pathogen genotype 
#' (value of 1), either complete protection (value of 0) 
#' @param params a LandsepiParams object.
#' @return an interaction matrix composed of 0 and 1 values.
#' @seealso \link{getMatrixGenePatho}, \link{getMatrixCultivarPatho}, \link{getMatrixPolyPatho}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivar4 <- loadCultivar(name = "Pyramid", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3, cultivar4)
#' , stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
#' simul_params <- allocateCultivarGenes(simul_params, "Pyramid", c("MG 1", "MG 2"))
#' croptypes <- loadCroptypes(simul_params,
#'                            names = c("Susceptible crop",
#'                                      "Resistant crop 1",
#'                                      "Mixture S+R",
#'                                      "Mixture R1+R2",
#'                                      "Pyramid crop"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture S+R", c("Susceptible", "Resistant1"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture R1+R2", c("Resistant1", "Resistant2"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Pyramid crop", c("Pyramid"))
#' simul_params <- setCroptypes(simul_params, croptypes)
#' getMatrixCroptypePatho(simul_params)
#' }
#' @export
getMatrixCroptypePatho <- function(params){
  Npatho <- prod(params@Genes$Nlevels_aggressiveness)
  mat_croptypes <- params@Croptypes
  Ncroptype <- nrow(mat_croptypes)
  mat_croptype_patho <- array(data = 0, dim = c(Ncroptype, Npatho))
  rownames(mat_croptype_patho)<- mat_croptypes$croptypeName  
  
  if(Ncroptype==0){
    stop("No croptypes are set")
  }else{
    
    mat_cultivar_patho <- getMatrixCultivarPatho(params)
    for (crop in 1:Ncroptype){
      indices_hosts <- which(mat_croptypes[crop, params@Cultivars$cultivarName] > 0)  
      ## indices of cultivars composing the croptype
      if (ncol(mat_croptype_patho)==1){
        mat_croptype_patho[crop,] <- sum(mat_cultivar_patho[indices_hosts,]) > 0  
        ## infection is possible if any of the cultivars can be infected
      }else if (length(indices_hosts)==1){
        mat_croptype_patho[crop,] <- mat_cultivar_patho[indices_hosts,]
      }else{
        mat_croptype_patho[crop,] <- apply(mat_cultivar_patho[indices_hosts,], 2
                                           , function(x) sum(x==1)>0)  
        ## sum of the compatibilities for each composing cultivar
        # mat_croptype_patho[crop,] > 0 <- 1
      }
    } ## for crop
    
    return(mat_croptype_patho)
  }
}


#' @name getMatrixPolyPatho
#' @title Get the "polygon/pathogen genotype" compatibility matrix.
#' @description Build the matrix indicating if infection is possible at the beginning of the season 
#' for every combination of polygon (rows) and pathogen genotype (columns).
#' @details For each polygon, there is either possibility of infection by the pathogen genotype 
#' (value of 1), either complete protection (value of 0) 
#' @param params a LandsepiParams object.
#' @return an interaction matrix composed of 0 and 1 values.
#' @seealso \link{getMatrixGenePatho}, \link{getMatrixCultivarPatho}, \link{getMatrixCroptypePatho}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setTime(simul_params, Nyears = 1, nTSpY = 80)
#' simul_params <- setLandscape(simul_params, loadLandscape(id = 1))
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivar4 <- loadCultivar(name = "Pyramid", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3, cultivar4)
#' , stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
#' simul_params <- allocateCultivarGenes(simul_params, "Pyramid", c("MG 1", "MG 2"))
#' croptypes <- loadCroptypes(simul_params,
#'                            names = c("Susceptible crop",
#'                                      "Resistant crop 1",
#'                                      "Mixture S+R",
#'                                      "Mixture R1+R2",
#'                                      "Pyramid crop"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture S+R", c("Susceptible", "Resistant1"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture R1+R2", c("Resistant1", "Resistant2"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Pyramid crop", c("Pyramid"))
#' simul_params <- setCroptypes(simul_params, croptypes)
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = 0,
#' prop=rep(1/5,5), aggreg=3 , rotation_sequence = croptypes$croptypeID)
#' getMatrixPolyPatho(simul_params)
#' }
#' @export
getMatrixPolyPatho <- function(params){
  Npatho <- prod(params@Genes$Nlevels_aggressiveness)
  Npoly <- nrow(params@Landscape)
  rotation_0 <- params@Landscape$year_1
  mat_poly_patho <- array(data = 0, dim = c(Npoly, Npatho))
  
  if(Npoly==0){
    stop("Landscape is not defined")
  }else{
    
    mat_croptype_patho <- getMatrixCroptypePatho(params)
    for (poly in 1:Npoly){
      id_croptype <- rotation_0[poly]+1  ## +1 because of C code
      mat_poly_patho[poly,] <- mat_croptype_patho[id_croptype,]
    } ## for poly
    
    return(mat_poly_patho)
  }
}


#' @name checkInoculum
#' @title Check inoculum
#' @description Checks inoculum validity.
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkInoculum <- function(params) {
  ret <- TRUE
  
  if(length(params@PI0)==0){
    warning("The inoculum has not been defined")
  }
  if (all(params@PI0 == 0)) {
    warning("The vector PI0 is only filled with 0")
  }
  if( any(!is.numeric(params@PI0)) || any(!is.in.01(params@PI0)) ){
    warning("Invalid inoculum value: must be in [0,1]")
    ret <- FALSE
  }
  if (length(params@PI0) > 1){
    Nhost <- nrow(params@Cultivars)
    Npatho <-prod(params@Genes$Nlevels_aggressiveness)
    Npoly <- nrow(params@Landscape)
    if (length(params@PI0) != Nhost*Npatho*Npoly){
      warning("The vector PI0 doesn't have the correct dimensions")
      ret <- FALSE
    }
  }
  
  return(ret)
}


#' @name checkPI0_mat
#' @title Check the array PI0_mat when entered manually in \code{loadInoculum()}.
#' @description Checks validity of the array.
#' @param mat a 3D array of dimensions (1:Nhost,1:Npatho,1:Npoly)
#' @param params a LandsepiParams object.
#' @return the same array at mat, possibly corrected if incompatibility has been detected
checkPI0_mat <- function(mat, params){
  warn = 0
  rotation_0 <- params@Landscape$year_1
  mat_croptypes <- params@Croptypes
  cultivarNames <- params@Cultivars$cultivarName
  mat_cultivar_patho <- getMatrixCultivarPatho(params)
  
  for (host in 1:dim(mat)[1]){
    for (patho in 1:dim(mat)[2]){
      for (poly in 1:dim(mat)[3]){
        id_croptype <- rotation_0[poly]+1  ## +1 because of C code
        host_present <- as.numeric(mat_croptypes[id_croptype, cultivarNames[host]] > 0)  
        ## i.e. is the host present in the croptype of poly ?
        host_compatible <- as.numeric(mat_cultivar_patho[host, patho])  
        ## i.e. is the patho able to infect the host ?
        
        if (host_present*host_compatible == 0 & mat[host, patho, poly] > 0){
          warn <- 1
          mat[host, patho, poly] <- host_present * host_compatible * mat[host, patho, poly]
        }
      }
    }
  }
  
  if (warn){
    warning("Incompatibility between the inoculum and the landscape (cultivar not present or 
            incapacity of pathogen genotype unable to infect the cultivar)")
    warning("The inoculum has been corrected to account for the landscape")  
  }
  
  return(mat)
}


#' @name compute_audpc100S
#' @title Compute AUDPC in a single 100% susceptible field
#' @description Compute AUDPC in a single field cultivated with a susceptible cultivar.
#' @param disease a disease name, among "rust" (default), "mildew" and "sigatoka"
#' @param hostType cultivar type, among: "growingHost" (default), "nongrowingHost", "grapevine".
#' @param nTSpY number to time steps per cropping season
#' @param area area of the field (must be in square meters).
#' @param seed an integer used as seed value (for random number generator).
#' @details audpc100S is the average AUDPC computed in a non-spatial simulation. 
#' @return The AUDPC value (numeric)
#' @examples 
#' \dontrun{
#' compute_audpc100S("rust", "growingHost", area=1E6)
#' compute_audpc100S("mildew", "grapevine", area=1E6)
#' compute_audpc100S("sigatoka", "banana", area=1E6, nTSpY=182)
#' }
#' @seealso \link{loadOutputs}
#' @export
compute_audpc100S <- function(disease="rust", hostType="growingHost"
                              , nTSpY=120, area=1E6, seed=12345){
  message(paste("Computing audpc100S for", disease, "in a single susceptible field of"
                , area, "m^2 during", nTSpY, "time steps"))
  
  res=simul_landsepi(seed=seed
                     , time_param = list(Nyears = 5, nTSpY = nTSpY)
                     , area=area
                     , basic_patho_param=loadPathogen(disease)
                     , cultivars=loadCultivar(name="Susceptible", type=hostType)
                     , epid_outputs = c("audpc") #, "audpc_rel", "gla", "gla_rel")
                     , evol_outputs = "", writeTXT = FALSE, graphic = FALSE)
  audpc100S <- mean(res$epid_outputs$audpc$total[2:5])
  return(audpc100S)    
}


#' @name loadOutputs
#' @title Load outputs
#' @description Creates an output list
#' @param epid_outputs a character string (or a vector of character strings if several outputs 
#' are to be computed) specifying the type of epidemiological and economic outputs to generate 
#' (see details):\itemize{
#' \item "audpc" : Area Under Disease Progress Curve (average number of diseased host individuals
#' per time step and square meter) 
#' \item "audpc_rel" : Relative Area Under Disease Progress Curve (average proportion of 
#' diseased host individuals relative to the total number of existing hosts)
#' \item "gla" : Green Leaf Area (average number of healthy host individuals per time step and 
#' square meter)
#' \item "gla_rel" : Relative Green Leaf Area (average proportion of healthy host individuals 
#' relative to the total number of existing hosts)
#' \item "eco_yield" : total crop yield (in weight or volume units per ha) 
#' \item "eco_cost" : operational crop costs (in monetary units per ha) 
#' \item "eco_product" : total crop products (in monetary units per ha) 
#' \item "eco_margin" : Margin (products - operational costs, in monetary units per ha) 
#' \item "contrib": contribution of pathogen genotypes to LIR dynamics
#' \item "HLIR_dynamics", "H_dynamics", "L_dynamics", "IR_dynamics", "HLI_dynamics", etc.: 
#' Epidemic dynamics related to the specified sanitary status (H, L, I or R and all their 
#' combinations). Graphics only, works only if graphic=TRUE.
#' \item "all" : compute all these outputs (default)
#' \item "" : none of these outputs will be generated.
#' }
#' @param evol_outputs a character string (or a vector of character strings if several 
#' outputs are to be computed) specifying the type of evolutionary outputs to generate :\itemize{
#' \item "evol_patho": Dynamics of pathogen genotype frequencies
#' \item "evol_aggr": Evolution of pathogen aggressiveness
#' \item "durability": Durability of resistance genes
#' \item "all": compute all these outputs (default)
#' \item "": none of these outputs will be generated.
#' }
#' @return a list of outputs and parameters for output generation
#' @seealso \link{setOutputs}, \link{compute_audpc100S}
#' @examples 
#' outputList <- loadOutputs(epid_outputs = "audpc", evol_outputs = "durability")
#' outputList
#' @export
loadOutputs <- function(epid_outputs = "all", evol_outputs = "all"){
  outputList <- list(epid_outputs = epid_outputs
                     , evol_outputs = evol_outputs
                     , thres_breakdown = 50000
                     , audpc100S = 0.76) ## audpc100S = 8.48 for grapevine mildew
  return(outputList)
}


#' @name setOutputs
#' @title Set outputs
#' @description Updates a LandsepiParams object with a list of output parameters.
#' @param params a LandsepiParams object.
#' @param output_list a list of outputs to be generated and parameters for output generation. 
#' It can be generated manually or, alternatively, via \code{\link{loadOutputs}}. This list 
#' is composed of:\itemize{
#' \item epid_outputs = epidemiological outputs to compute (see details)
#' \item evol_outputs = evolutionary outputs to compute (see details)
#' \item thres_breakdown = an integer (or vector of integers) giving the threshold 
#' (i.e. number of infections) above which a pathogen genotype is unlikely to go extinct, 
#' used to characterise the time to invasion of resistant hosts (several values are computed 
#' if several thresholds are given in a vector).
#' \item audpc100S = the audpc in a fully susceptible landscape (used as reference value 
#' for graphics).
#' }
#' @details "epid_outputs" is a character string (or a vector of character strings if several 
#' outputs are to be computed) specifying the type of epidemiological and economic outputs 
#' to generate:  
#' \itemize{
#' \item "audpc" : Area Under Disease Progress Curve (average number of diseased host individuals
#' per time step and square meter) 
#' \item "audpc_rel" : Relative Area Under Disease Progress Curve (average proportion of 
#' diseased host individuals relative to the total number of existing hosts)
#' \item "gla" : Green Leaf Area (average number of healthy host individuals per square meter)
#' \item "gla_rel" : Relative Green Leaf Area (average proportion of healthy host individuals 
#' relative to the total number of existing hosts)
#' \item "eco_yield" : total crop yield (in weight or volume units per ha) 
#' \item "eco_cost" : operational crop costs (in monetary units per ha) 
#' \item "eco_product" : total crop products (in monetary units per ha) 
#' \item "eco_margin" : Margin (products - costs, in monetary units per ha) 
#' \item "contrib": contribution of pathogen genotypes to LIR dynamics
#' \item "HLIR_dynamics", "H_dynamics", "L_dynamics", "IR_dynamics", "HLI_dynamics", etc.: 
#' Epidemic dynamics related to the specified sanitary status (H, L, I or R and all their 
#' combinations). Graphics only, works only if graphic=TRUE.
#' \item "all" : compute all these outputs (default)
#' \item "" : none of these outputs will be generated.
#' }
#' "evol_outputs" is a character string (or a vector of character strings if several outputs 
#' are to be computed) specifying the type of evolutionary outputs to generate :\itemize{
#' \item "evol_patho": Dynamics of pathogen genotype frequencies
#' \item "evol_aggr": Evolution of pathogen aggressiveness
#' \item "durability": Durability of resistance genes
#' \item "all": compute all these outputs (default)
#' \item "": none of these outputs will be generated.
#' }
#' 
#' @return a LandsepiParams object.
#' @seealso \link{loadOutputs}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setOutputs(simul_params, loadOutputs())
#' simul_params@Outputs
#' }
#' @export
setOutputs <- function(params, output_list){
  params@Outputs <- output_list
  
  checkOutputs(params) 
  
  return(params)
}

#' @name checkOutputs
#' @title Check outputs
#' @description Checks outputs validity.
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkOutputs <- function(params) {
  ret <- TRUE
  
  if( !is.character(params@Outputs$epid_outputs) ||
      !is.character(params@Outputs$evol_outputs)) {
    warning("Invalid epidemiological or evolutionary outputs")
    ret <- FALSE
  }
  
  if (!is.na(params@Outputs$audpc100S)){
    if( !is.numeric(params@Outputs$audpc100S) ||
        sum(!is.strict.positive(params@Outputs$audpc100S) > 0) ){
      warning("AUDPC in a fully susceptible landscape must be > 0")
      ret <- FALSE
    }
  }

  if (!is.na(params@Outputs$thres_breakdown)){
    if( !is.wholenumber(params@Outputs$thres_breakdown) || 
        !is.strict.positive(params@Outputs$thres_breakdown) ) {
      warning("Threshold for resistance breakdown must be a whole number > 0")
      ret <- FALSE
    }
  }
  return(ret)
}

###### PRIVATE ######

#' params2CroptypeBDD
#' @description Converts a LandsepiParams object to a value compatible with BDD croptype Table
#' @param params a LandsepiParams object.
#' @return a data.frame BDD compatible
#' @noRd
params2CroptypeBDD <- function(params) {
  if (is.null(params@Croptypes$croptypeName)) {
    croptypes <- params@Croptypes
  } else {
    croptypes <- params@Croptypes[, -which(colnames(params@Croptypes) == "croptypeName")]
  }

  colnames(croptypes) <- 
    c("croptypeID", which(rownames(params@Cultivars) 
                          %in% colnames(croptypes)[-which(.croptypesColNames 
                                                          %in% colnames(croptypes))]))
  res <- data.frame()
  for (i in 1:nrow(croptypes)) {
    for (culti in which(croptypes[i, -1] > 0)) { ## without croptypeID and croptypeName index
      # message(paste0(croptypes[i,1]," ",culti," ",croptypes[i,culti+1]))
      res <- rbind(res, c(croptypes[i, "croptypeID"]
                          , culti - 1, croptypes[i, which(colnames(croptypes) == culti)]))
    }
  }
  res <- cbind(0:(nrow(res) - 1), res)
  colnames(res) <- c("rowid", "croptypeID", "cultivarID", "proportion")
  return(res)
}

#' CroptypeBDD2Params
#' @description Converts BDD table to Croptype LandspeParams object
#' @param inputGPKG a GPKG filename
#' @return a data.frame LandsepiParams@@Croptypes compatible
#' @noRd
CroptypeBDD2Params <- function(inputGPKG) {
  return(getGPKGCroptypes(inputGPKG))
}

#' params2CultivarBDD
#' @description Converts a LandsepiParams object to a value compatible with BDD cultivar Table
#' @param params a LandsepiParam object.
#' @return a data.frame BDD compatible
#' @noRd
params2CultivarBDD <- function(params) {
  cultivars <- data.frame(params@Cultivars, stringsAsFactors = FALSE)
  cultivars <- data.frame("cultivarID" = 0:(nrow(cultivars) - 1), cultivars
                          , stringsAsFactors = FALSE)
  return(cultivars)
}

#' CultivarBDD2Params
#' @description Converts BDD table Cultivar to a Cultivar LandsepiParams object
#' @param inputGPKG a GPKG filename
#' @return a data.frame LandsepiParams@@Cultivars compatible
#' @noRd
CultivarBDD2Params <- function(inputGPKG) {
  return(getGPKGCultivars(inputGPKG))
}

#' params2GeneBDD
#' @description Converts a LandsepiParams object to a value compatible with BDD Gene Table
#' @param params a LandsepiParam object.
#' @return a data.frame BDD compatible
#' @noRd
params2GeneBDD <- function(params) {
  gene <- data.frame(params@Genes, stringsAsFactors = FALSE)
  gene <- data.frame("geneID" = 0:(nrow(gene) - 1), gene, stringsAsFactors = TRUE)
  return(gene)
}

#' GeneBDD2Params
#' @description Converts BDD table Gene to LandsepiParams object
#' @param inputGPKG a LandsepiParams
#' @return a data.frame LandsepiParams@@Genes compatible
#' @noRd
GeneBDD2Params <- function(inputGPKG) {
  return(getGPKGGenes(inputGPKG))
}

#' params2GeneListBDD
#' @description Converts a LandsepiParams object to a value compatible with BDD cultivarsGenes Table
#' @param params a LandsepiParam object.
#' @return a data.frame BDD compatible
#' @noRd
params2GeneListBDD <- function(params) {
  cgList <- params@CultivarsGenes
  res <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:nrow(cgList)) {
    culti <- rownames(cgList)[i]
    for (gene in which(cgList[i, ] > 0)) {
      res <- rbind(res, c(which(rownames(params@Cultivars) == culti) - 1
                          , which(rownames(params@Genes) == colnames(cgList)[gene]) - 1))
    }
  }
  res <- cbind(0:(nrow(res) - 1), res)
  colnames(res) <- c("rowid", "cultivarID", "geneID")
  return(res)
}

#' CultivarGeneBDD2Params
#' @description Converts BDD table to LandsepiParams object
#' @param inputGPKG a GPKG filename
#' @return a data.frame LandsepiParams@@CultivarsGenes compatible
#' @noRd
CultivarGeneBDD2Params <- function(inputGPKG) {
  return(getGPKGCultivarsGenes(inputGPKG))
}
