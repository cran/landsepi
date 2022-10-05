## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- results="hide", message=FALSE-------------------------------------------
library(landsepi)

## ---- results="hide", message="FALSE"-----------------------------------------
simul_params <- createSimulParams(outputDir = getwd())

## -----------------------------------------------------------------------------
simul_params@Seed
simul_params <- setSeed(simul_params, seed = 1)
simul_params@Seed
Nyears = 6 
nTSpY = 120
simul_params <- setTime(simul_params, Nyears = Nyears, nTSpY = nTSpY)

## -----------------------------------------------------------------------------
simul_params <- setInoculum(simul_params, val = 5e-4)

## -----------------------------------------------------------------------------
landscape <- loadLandscape(id = 1)
simul_params <- setLandscape(simul_params, land = landscape)

## -----------------------------------------------------------------------------
basic_patho_param <- loadPathogen(disease = "rust")

## -----------------------------------------------------------------------------
basic_patho_param$repro_sex_prob <- 1  ## at every time step all pathogen individuals reproduces sexually
basic_patho_param
basic_patho_param$repro_sex_prob <- 0.5 ##  at every time step half of the pathogen population
                                        ## reproduce clonally and half sexually
basic_patho_param

## -----------------------------------------------------------------------------
repro_sex_probs <- c(rep(0.0, nTSpY), 1.0)  

## -----------------------------------------------------------------------------
simul_params <- setReproSexProb(simul_params, repro_sex_probs)
simul_params@ReproSexProb 

## -----------------------------------------------------------------------------
basic_patho_param$sex_propagule_release_mean = 1
basic_patho_param$sex_propagule_viability_limit = 5
simul_params <- setPathogen(simul_params, basic_patho_param)  

## -----------------------------------------------------------------------------
basic_patho_param$clonal_propagule_gradual_release = TRUE ## clonal propagules are progressively 
                                                          ## released during the next cropping season
basic_patho_param$clonal_propagule_gradual_release = FALSE ## clonal propagules are released at 
                                                          ## the first day of the next cropping season

## -----------------------------------------------------------------------------
disp_patho <- loadDispersalPathogen(id = 1)

## -----------------------------------------------------------------------------
disp_patho_asex <- disp_patho[[1]]
disp_patho_sex <- disp_patho[[2]]
head(disp_patho_asex)
head(disp_patho_sex)

## -----------------------------------------------------------------------------
disp_patho_asex <- disp_patho[[1]]
disp_patho_sex <- disp_patho[[1]]
head(disp_patho_asex)
head(disp_patho_sex)

## -----------------------------------------------------------------------------
simul_params <- setDispersalPathogen(simul_params, disp_patho_asex, disp_patho_sex)

## -----------------------------------------------------------------------------
# Resistance genes

gene1 <- loadGene(name = "gene 1", type = "majorGene")
gene2 <- loadGene(name = "gene 2", type = "QTL")

#gene2$recombination_sd <- 0.8
gene2$Nlevels_aggressiveness <- 3
genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)


## -----------------------------------------------------------------------------
# Cultivars

cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3)
                        , stringsAsFactors = FALSE)

# Allocating genes to cultivars

simul_params <- setGenes(simul_params, dfGenes = genes)
simul_params <- setCultivars(simul_params, dfCultivars = cultivars)

simul_params <- allocateCultivarGenes(simul_params
                                      , cultivarName = "Resistant1"
                                      , listGenesNames = c("gene 1"))
simul_params <- allocateCultivarGenes(simul_params
                                      , cultivarName = "Resistant2"
                                      , listGenesNames = c("gene 2"))

# Allocating cultivars to croptypes

croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
                                                   , "Resistant crop 1"
                                                   , "Resistant crop 2"))

croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Susceptible crop"
                                       , cultivarsInCroptype = "Susceptible")
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Resistant crop 1"
                                       , cultivarsInCroptype = "Resistant1")
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Resistant crop 2"
                                       , cultivarsInCroptype = "Resistant2")

simul_params <- setCroptypes(simul_params, dfCroptypes = croptypes)

# Allocating croptypes to fields of the landscape

rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
rotation_period <- 0  # number of years before rotation of the landscape
prop <- c(1/3,1/3,1/3)  # proportion (in surface) of each croptype
aggreg <- 0    # level of spatial aggregation
simul_params <- allocateLandscapeCroptypes(simul_params
                                           , rotation_period = rotation_period
                                           , rotation_sequence = rotation_sequence
                                           , prop = prop
                                           , aggreg = aggreg
                                           , graphic = FALSE)

# Define fungicide treatments
# treatment <- list(treatment_reduction_rate = 1,
#                     treatment_efficiency = 0.8,
#                     treatment_timesteps =  vector() ,
#                     #treatment_timesteps =  c(7,21,35) ,
#                     treatment_cultivars =  vector(),
#                     #cultivars= c(0),
#                     treatment_cost = 0) 
# simul_params <- setTreatment(simul_params, treatment)

# Choosing output variables

outputlist <- loadOutputs(epid_outputs = "all", evol_outputs = "all")

simul_params <- setOutputs(simul_params, outputlist)

## ---- eval=FALSE--------------------------------------------------------------
#  checkSimulParams(simul_params)
#  runSimul(simul_params, graphic = TRUE, videoMP4 = FALSE)

## ---- include=FALSE-----------------------------------------------------------
system(paste("rm -rf ", simul_params@OutputDir))

