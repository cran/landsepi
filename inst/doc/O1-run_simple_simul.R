## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----results="hide", message=FALSE--------------------------------------------
library(landsepi)

## ----results="hide", message=FALSE--------------------------------------------
simul_params <- createSimulParams(outputDir = getwd())

## -----------------------------------------------------------------------------
simul_params@Seed
simul_params <- setSeed(simul_params, seed = 1)
simul_params@Seed

## -----------------------------------------------------------------------------
simul_params <- setTime(simul_params, Nyears = 6, nTSpY = 120)
simul_params@TimeParam

## -----------------------------------------------------------------------------
basic_patho_param <- loadPathogen(disease = "rust")
basic_patho_param
basic_patho_param <- loadPathogen(disease = "mildew")
basic_patho_param
basic_patho_param <- loadPathogen(disease = "sigatoka")
basic_patho_param

## -----------------------------------------------------------------------------
basic_patho_param <- loadPathogen("rust")
basic_patho_param$infection_rate <- 0.5
basic_patho_param

## -----------------------------------------------------------------------------
basic_patho_param <- list(infection_rate = 0.4
                          , latent_period_mean = 10
                          , latent_period_var = 9
                          , propagule_prod_rate = 3.125
                          , infectious_period_mean = 24
                          , infectious_period_var = 105
                          , survival_prob = 1e-4
                          , repro_sex_prob = 0
                          , sigmoid_kappa = 5.333
                          , sigmoid_sigma = 3
                          , sigmoid_plateau = 1
                          , sex_propagule_viability_limit = 1
                          , sex_propagule_release_mean = 1
                          , clonal_propagule_gradual_release = 0)

## -----------------------------------------------------------------------------
simul_params <- setPathogen(simul_params, patho_params = basic_patho_param)
simul_params@Pathogen

## ----fig.alt="The first Landscape available"----------------------------------
landscape <- loadLandscape(id = 1)
length(landscape)
plot(landscape, main = "Landscape structure")

## -----------------------------------------------------------------------------
disp_patho_clonal <- loadDispersalPathogen(id = 1)[[1]]
head(disp_patho_clonal)
length(landscape)^2 == length(disp_patho_clonal)

## -----------------------------------------------------------------------------
simul_params <- setLandscape(simul_params, land = landscape)
simul_params <- setDispersalPathogen(simul_params, disp_patho_clonal)

## -----------------------------------------------------------------------------
cultivar1 <- loadCultivar(name = "Susceptible", type = "wheat")
cultivar2 <- loadCultivar(name = "Resistant1", type = "wheat")
cultivar3 <- loadCultivar(name = "Resistant2", type = "banana")
cultivar4 <- loadCultivar(name = "Resistant3", type = "pepper")
cultivar5 <- loadCultivar(name = "Forest", type = "nonCrop")
cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3, cultivar4, cultivar5)
                        , stringsAsFactors = FALSE)
cultivars

## -----------------------------------------------------------------------------
cultivars[cultivars$cultivarName == "Susceptible", "growth_rate"] <- 0.2
cultivars

## -----------------------------------------------------------------------------
cultivars_new <- data.frame(cultivarName = c("Susceptible", "Resistant"),
                            initial_density =   c(0.1, 0.2),
                            max_density =       c(2.0, 3.0),
                            growth_rate =       c(0.1, 0.2),
                            reproduction_rate = c(0.0, 0.0),
                            yield_H =           c(2.5, 2.0),
                            yield_L =           c(0.0, 0.0),
                            yield_I =           c(0.0, 0.0),
                            yield_R =           c(0.0, 0.0),
                            planting_cost =   c(225, 300),
                            market_value =      c(200, 150),
                            stringsAsFactors = FALSE)
cultivars_new

## -----------------------------------------------------------------------------
gene1 <- loadGene(name = "MG 1", type = "majorGene")
gene2 <- loadGene(name = "Lr34", type = "APR")
gene3 <- loadGene(name = "gene 3", type = "QTL")
gene4 <- loadGene(name = "nonhost resistance", type = "immunity")
genes <- data.frame(rbind(gene1, gene2, gene3, gene4), stringsAsFactors = FALSE)
genes

## -----------------------------------------------------------------------------
genes[genes$geneName == "MG 1", "mutation_prob"] <- 1e-3
genes

## -----------------------------------------------------------------------------
genes_new <- data.frame(geneName =               c("MG1", "MG2"),
                        efficiency =             c(1.0  , 0.8  ),
                        age_of_activ_mean =      c(0.0  , 0.0  ),
                        age_of_activ_var =       c(0.0  , 0.0  ),
                        mutation_prob =          c(1E-7 , 1E-4),
                        Nlevels_aggressiveness = c(2    , 2    ),
                        adaptation_cost =        c(0.50 , 0.75 ),
                        relative_advantage =     c(0.50 , 0.75 ),
                        tradeoff_strength =      c(1.0  , 1.0  ),
                        target_trait =           c("IR" , "LAT"),
                        recombination_sd =       c(1.0,1.0),
                        stringsAsFactors = FALSE)
genes_new

## -----------------------------------------------------------------------------
simul_params <- setGenes(simul_params, dfGenes = genes)
simul_params <- setCultivars(simul_params, dfCultivars = cultivars)
simul_params@Genes
simul_params@Cultivars

## -----------------------------------------------------------------------------
simul_params <- allocateCultivarGenes(simul_params
                                      , cultivarName = "Resistant1"
                                      , listGenesNames = c("MG 1"))
simul_params <- allocateCultivarGenes(simul_params
                                      , cultivarName = "Resistant2"
                                      , listGenesNames = c("Lr34", "gene 3"))
simul_params <- allocateCultivarGenes(simul_params
                                      , cultivarName = "Resistant3"
                                      , listGenesNames = c("nonhost resistance"))
simul_params@CultivarsGenes

## -----------------------------------------------------------------------------
croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
                                                   , "Pure resistant crop"
                                                   , "Mixture"
                                                   , "Other"))
croptypes

## -----------------------------------------------------------------------------
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Susceptible crop"
                                       , cultivarsInCroptype = "Susceptible")
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Pure resistant crop"
                                       , cultivarsInCroptype = "Resistant1")
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Mixture"
                                       , cultivarsInCroptype = c("Resistant2","Resistant3")
                                       , prop = c(0.4, 0.6))
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Other"
                                       , cultivarsInCroptype = "Forest")
croptypes

## -----------------------------------------------------------------------------
simul_params <- setCroptypes(simul_params, dfCroptypes = croptypes)
simul_params@Croptypes

## -----------------------------------------------------------------------------
croptypes <- data.frame(croptypeID = c(0, 1, 2, 3) ## must start at 0 and match with values from landscape "croptypeID" layer
                        , croptypeName = c("Susceptible crop"
                                           , "Pure resistant crop"
                                           , "Mixture"
                                           , "Other")
                        , Susceptible = c(1,0,0  ,0)
                        , Resistant1  = c(0,1,0  ,0)
                        , Resistant2  = c(0,0,0.5,0)
                        , Resistant3  = c(0,0,0.5,0)
                        , Forest      = c(0,0,0  ,1)
                        , stringsAsFactors = FALSE)
simul_params <- setCroptypes(simul_params, croptypes)

## -----------------------------------------------------------------------------
# croptypeIDs cultivated in each element of the rotation sequence:
rotation_sequence <- list(c(0,1,3), c(0,2,3))
rotation_period <- 2  # number of years before rotation of the landscape
prop <- list(rep(1/3, 3), rep(1/3, 3)) # proportion (in surface) of each croptype 
aggreg <-1 # level of spatial aggregation
simul_params <- allocateLandscapeCroptypes(simul_params
                                           , rotation_period = rotation_period
                                           , rotation_sequence = rotation_sequence
                                           , prop = prop
                                           , aggreg = aggreg
                                           , graphic = TRUE)
# plot(simul_params@Landscape)

## -----------------------------------------------------------------------------
getMatrixGenePatho(simul_params)
getMatrixCultivarPatho(simul_params)
getMatrixCroptypePatho(simul_params)
getMatrixPolyPatho(simul_params)[1:10,]

## -----------------------------------------------------------------------------
# Option 1: simply use the default parameterisation
simul_params <- setInoculum(simul_params, 5E-4)
 
# Option 2: use loadInoculum()
Npatho <- prod(simul_params@Genes$Nlevels_aggressiveness)
Nhost <- nrow(simul_params@Cultivars)
pI0 <- loadInoculum(simul_params, pI0_all=5E-4, pI0_host=c(1,rep(0, Nhost-1)), pI0_patho=c(1,rep(0, Npatho-1)))
simul_params <- setInoculum(simul_params, pI0)

## ----eval=FALSE---------------------------------------------------------------
# Npatho <- prod(simul_params@Genes$Nlevels_aggressiveness)  ## Nb of pathogen genotypes
# Nhost <- nrow(simul_params@Cultivars)  ## Nb of cultivars
# Npoly <- nrow(simul_params@Landscape)  ## Nb of polygons in the landscape
# Npoly_inoc <- 5  ## number of inoculated polygons
# compatible_poly <- getMatrixPolyPatho(simul_params)[,1]  ## whether the avr pathogen can infect the polygons
# id_poly <- sample(grep(1, compatible_poly), Npoly_inoc)  ## random polygon picked among compatible ones
# pI0_poly <- as.numeric(1:Npoly %in% id_poly)
# pI0 <- loadInoculum(simul_params, pI0_all=5E-4, pI0_host=c(1,rep(0, Nhost-1)), pI0_patho=c(1,rep(0, Npatho-1)),
#                     pI0_poly=pI0_poly)
# simul_params <- setInoculum(simul_params, pI0)

## ----eval=FALSE---------------------------------------------------------------
# ## Example with 4 pathogen genotypes and 2 cultivars
# pI0 <- loadInoculum(simul_params, pI0_patho=c(1E-3,1E-4,1E-4,1E-5), pI0_host=c(1,1))
# simul_params <- setInoculum(simul_params, pI0)

## ----eval=FALSE---------------------------------------------------------------
# Npoly <- nrow(simul_params@Landscape)
# Npoly_inoc <- 5  ## number of inoculated polygons
# id_poly <- sample(1:Npoly, Npoly_inoc)  ## random polygon
# pI0_poly <- as.numeric(1:Npoly %in% id_poly)
# pI0 <- loadInoculum(simul_params, pI0_patho=c(1E-3,1E-4,1E-4,1E-5), pI0_host=c(1,1), pI0_poly=pI0_poly)
# simul_params <- setInoculum(simul_params, pI0)

## ----eval=FALSE---------------------------------------------------------------
# ## example with 2 cultivars, 4 pathogen genotypes and 5 fields
# Nhost=2
# Npatho=4
# Npoly=5
# pI0 <- array(data = 1:40 / 100, dim = c(Nhost, Npatho, Npoly))
# simul_params <- setInoculum(simul_params, pI0)

## ----eval=FALSE---------------------------------------------------------------
# corrected_pI0 <- loadInoculum(simul_params, pI0_mat=pI0)
# simul_params <- setInoculum(simul_params, corrected_pI0)

## -----------------------------------------------------------------------------
inoculumToMatrix(simul_params)[,,1:5]

## -----------------------------------------------------------------------------
Ncroptypes <- nrow(simul_params@Croptypes)
Nyears <- simul_params@TimeParam$Nyears
## Same probability in every croptype:
simul_params <- updateSurvivalProb(simul_params, mat_year=1:Nyears/100)
simul_params@Pathogen
## Same probability every year:
simul_params <- updateSurvivalProb(simul_params, mat_croptype=1:Ncroptypes/10)
simul_params@Pathogen
## specific probability for different croptypes and years:
simul_params <- updateSurvivalProb(simul_params, mat_year=1:Nyears/100, mat_croptype=1:Ncroptypes/10)
simul_params@Pathogen
## One probability per year and per croptype:
simul_params <- updateSurvivalProb(simul_params, mat=matrix(runif(Nyears*Ncroptypes), ncol=Ncroptypes))
simul_params@Pathogen

## -----------------------------------------------------------------------------
survivalProbToMatrix(simul_params)

## -----------------------------------------------------------------------------
## Via function loadTreatment:
treatment <- loadTreatment(disease="mildew")
## Direct implementation:
treatment <- list(treatment_degradation_rate = 0.1,
                  treatment_efficiency = 0.8,
                  treatment_timesteps =  seq(1,120,14) ,
                  treatment_cultivars  = c(0),
                  treatment_cost = 0,
                  treatment_application_threshold = c(0.0))
## Set the parameters
simul_params <- setTreatment(simul_params, treatment)

## -----------------------------------------------------------------------------
outputlist <- loadOutputs(epid_outputs = "all", evol_outputs = "all", disease="rust")
outputlist

## -----------------------------------------------------------------------------
audpc100S <- compute_audpc100S("rust", "wheat", area=1E6)
audpc100S <- compute_audpc100S("mildew", "grapevine", area=1E6)
audpc100S <- compute_audpc100S("sigatoka", "banana", area=1E6, nTSpY=182)

## -----------------------------------------------------------------------------
simul_params <- setOutputs(simul_params, outputlist)

## ----eval=FALSE---------------------------------------------------------------
# checkSimulParams(simul_params)
# simul_params <- saveDeploymentStrategy(simul_params)

## ----eval=FALSE---------------------------------------------------------------
# runSimul(simul_params, graphic = TRUE, videoMP4 = FALSE)

## ----include=FALSE------------------------------------------------------------
system(paste("rm -rf ", simul_params@OutputDir))

