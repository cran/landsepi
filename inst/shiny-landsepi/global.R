#library(shinycssloaders)
library(shiny)
library(shinyBS)
library(DT)
library(shinyjs)
#library(slickR)
library(gridExtra)
library(png)
library(tiff)
library(grid)
library(future)
library(promises)
library(tools)
library(shinyalert)

## Video directory
if(!dir.exists("./www/tmp/")) dir.create("./www/tmp/")
addResourcePath("video", "./www/tmp/")
addResourcePath("graphics", "./www/tmp/")

# Active debug level
# 0 : no print
# 1 : warning
# 2 : error
# 3 : all
ACTIVE_DEBUG <- 0

library("landsepi")
data(package = "landsepi")

source("modules/editableDT.R")

VALUEMAX <- 10000


## del all file and directory of a path
cleanDir <- function(path) {
  files <- dir(path, full.names=TRUE, no..=TRUE)
  lapply(files, FUN = function(file){
      if( dir.exists(file) ) cleanDir(file)
      file.remove(file, recursive=TRUE)
    })
}

ROOT_PATH <- getwd()

if(!dir.exists(paste0(ROOT_PATH,"/www/tmp/"))) dir.create(paste0(ROOT_PATH,"/www/tmp/"))
setwd(paste0(ROOT_PATH,"/www/tmp/"))

cleanDir(paste0(ROOT_PATH,"/www/tmp/"))

## User mode
advanced_mode <- reactiveVal(TRUE)
can_rm_croptypes <- reactiveVal(TRUE)

## Croptypes proportions in landscape
croptypes_proportions <- shiny::reactiveVal(c(1))

## simul params reactive for view
## use to update view 
simul_params_croptypes <- shiny::reactiveVal()
default_croptype <- c()
simul_params_cultivars <- shiny::reactiveVal()
default_croptype <- c()
simul_params_cultivarsgenes <- shiny::reactiveVal()
simul_params_genes <- shiny::reactiveVal()
default_gene <- c()

## Pathogens name reactive value
simul_pathogen <-  shiny::reactiveVal("rust")

## Treatment is active
treatment_is_active <- shiny::reactiveVal(FALSE)


################################
# About Text
################################
aboutText <- paste0("<h1>Landsepi: Landscape Epidemiology and Evolution</h1><img src='landsepi-logo.png'  align='right' alt='' width='120'/>
                          <h3> A stochastic, spatially-explicit, demo-genetic model
                         simulating the spread and evolution of a plant pathogen in a heterogeneous landscape
                         to assess resistance deployment strategies. It is based on a spatial geometry for describing
                         the landscape and allocation of different cultivars, a dispersal kernel for the
                         dissemination of the pathogen, and a SEIR ('Susceptible-Exposed-Infectious-Removed’)
                         structure with a discrete time step. It provides a useful tool to assess the performance
                         of a wide range of deployment options with respect to their epidemiological,
                         evolutionary and economic outcomes.</h3>
                          <h3> Authors:</h3> Rimbaud Loup, Zaffaroni Marta, Rey Jean-François, Papaïx Julien
                          <h3>Package project:</h3><a href='https://CRAN.R-project.org/package=landsepi' target='_blank'> CRAN package</a><br/><a href='https://forge.inrae.fr/landsepi/landsepi' target='_blank'> Source code</a>
                          <br/><a href='https://landsepi.pages-forge.inrae.fr/landsepi' target='_blank'> Package Documentation</a>
                          <br/> License GPL-3
                          <h3> How to cite the package:</h3> <b>Rimbaud L, Zaffaroni M, Rey J, Papaïx J (",sub("-.*","",utils::packageDate("landsepi")),").</b> landsepi: Landscape Epidemiology and Evolution. R package version ",utils::packageVersion("landsepi"),", &lt;URL: https://cran.r-project.org/package=landsepi&gt;.
                          <h3> Funding</h3>
			              This work benefited from ANR project 'ArchiV' (2019–2023, grant n°ANR-18-CE32-0004-01), AFB Ecophyto II-Leviers Territoriaux Project ”Médée” (2020–2023), GRDC grant CSP00192 and the CSIRO/INRA linkage program,
ANR project 'COMBINE' (2022-2026, grant n°ANR-22-CE32-0004), SPE project 'DYNAMO'(2022-2024), INRA Program ASC (2008-2014)
                          <div>
                          <img src='Republique_Francaise_RVB.jpg' alt='RF' style='width:50px; margin-left: 10px;' />
                          <img src='LogoINRAE_Fr_rvb_web.png' alt='INRAE' style='width:50px; margin-left: 10px;' />
                          <img src='logoBIOSP.jpeg' alt='BioSP' style='width:50px; margin-left: 10px;'/>
                          <img src='PATHO_inra_logo.png' alt='Pathologie végétale' style='width:50px; margin-left: 10px;'/>
                          <img src='2022_LogoSAVE_entier_court_bleu.png' alt='SAVE' style='width:50px; margin-left: 10px;'/>
                          <img src='CSIRO_Logo.png' alt='CSIRO' style='width:40px; margin-left: 10px;'/>
                          </div>")

################################
# List Values
################################
listGenes <- c("majorGene", "APR", "QTL" , "immunity")
listCultivarsType <- c("crop","nonCrop")



##################################################################
# Functions
##################################################################

# Print message in terminal
# msg : message to print
# level : level of output (0 no output, 3 all output)
printVerbose <- function(msg, level=ACTIVE_DEBUG) {
  if(ACTIVE_DEBUG < 1 || level < 1) return()
  if(level == 1) warnings(msg)
  if(level == 2) print(msg)
  if(level >= 3) cat(file=stderr(),"### DEBUG ",msg,"\n")
}

## Show message
## id : message id
## selectorafter : id element to place message after
## message : the message
showErrorMessage <- function(id="errorMessage", selectorafter = "#", message = "a message") {
  shiny::insertUI(
    selector = selectorafter,
    where = "afterEnd",
    ui = tags$div(
      id = id,
      class = "alert alert-danger",
      shiny::HTML(message)
      )
    )
}

## Check croptypes Table
# col :
#  1 : ID
#  2 : Name
#  3:ncol : cultivars
#  ncol+1 : Landscape proportions
checkCroptypesTable <- function(data) {
  isok <- TRUE
  ## croptype ID
  # shiny::removeUI(selector = "#croptypeIdError")
  # if( sum(is.na(as.numeric(data[,"croptypeID"]))) != 0 || length(unique(data[,"croptypeID"])) != length(data[,"croptypeID"]) ) {
  #   showErrorMessage(id = "croptypeIdError", selectorafter= "#generateLandscape",
  #                    message = "Each croptype index must be a unique numeric")
  #   isok <- FALSE
  # }
  
  ## croptype name
  shiny::removeUI(selector = "#croptypeNameError")
  if( sum(as.character(data[,"croptypeName"]) == "") != 0 | sum(grepl("^\\s*$",as.character(data[,"croptypeName"]))) != 0) {
    
    showErrorMessage(id = "croptypeNameError", selectorafter= "#generateLandscape",
                     message = "Croptype names must be strings")
    #shinyjs::disable(id = "generateLandscape")
  }
  
  ## croptypes cultivars proportions
    shiny::removeUI(selector = "#croptypeError")
    ## no cultivars
    if( ncol(data) <= 2) {
      showErrorMessage( id = "croptypeError", selectorafter = "#generateLandscape",
                        message = paste0("No cultivar has been defined"))
      isok <- FALSE
    }
    else{
      value <- as.matrix(data[,c(3:ncol(data))], nrow = nrow(data) )
      
      if( sum(as.numeric(value) < 0.0) != 0 | sum(as.numeric(value) > 1.0) != 0 ) {
        showErrorMessage( id = "croptypeError", selectorafter = "#generateLandscape",
                          message = paste0("The proportion of every cultivar within a croptype should be between 0 and 1 (0% and 100%)"))
        isok <- FALSE
      }
      else {
        sum_prop <- sapply(1:nrow(value), FUN = function(i) { !isTRUE(all.equal(sum(as.numeric(value[i,])),1))}) 
        #message(sum_prop)
        if (sum(sum_prop) != 0 ) {
          showErrorMessage(id = "croptypeError", selectorafter = "#generateLandscape",
                           message =  paste0("The sum of the proportions of all cultivars composing a Croptype should be equal to 1 (100%)"))
          isok <- FALSE
        }
      }
    }
    
  return(invisible(isok))
  
}

## Check croptypes Table
# col :
#  1 : Name
#  2:ncol : parameters
checkCultivarsTable <- function(data) {
  isok <- TRUE
  
  shiny::removeUI(selector = "#cultivarsNameError")
  if( sum(as.character(data[,"cultivarName"]) == "") != 0 
      || sum(grepl("^\\s*$",as.character(data[,"cultivarName"]))) != 0) {
    
    showErrorMessage(id = "cultivarsNameError", selectorafter= "#generateLandscape",
                     message = "Cultivar names must be strings without spaces")
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#cultivarsZeroError")
  if( sum(data[,-1] < 0) != 0 ) {
    
    showErrorMessage(id = "cultivarsZeroError", selectorafter= "#generateLandscape",
                     message = "Values in the 'Cultivar' table should be >=0")
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#cultivarsStricZeroError")
  if( sum(data[,"max_density"] <= 0) != 0 || sum(data[,"max_density"] > VALUEMAX) != 0) {
    
    showErrorMessage(id = "cultivarsStricZeroError", selectorafter= "#generateLandscape",
                     message = paste0("Cultivar max_density values should be between 0 and ",VALUEMAX))
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#cultivarsmaxvalueError")
  if( sum(data[,c("growth_rate","reproduction_rate")] > 1) != 0) {
    
    showErrorMessage(id = "cultivarsmaxvalueError", selectorafter= "#generateLandscape",
                     message = paste0("Cultivar 'growth_rate' / 'reproduction_rate' values should be lower than 1"))
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#cultivarsValueError")
  if( sum(data[,- which(c("cultivarName", "growth_rate","reproduction_rate") %in% colnames(data))] < 0) != 0
      || sum(data[,- which(c("cultivarName", "growth_rate","reproduction_rate") %in% colnames(data))] > VALUEMAX) != 0) {
    
    showErrorMessage(id = "cultivarsValueError", selectorafter= "#generateLandscape",
                     message = paste0("Values in the 'Cultivar' table should be lower than ",VALUEMAX))
    isok <- FALSE
  }
  
  return(invisible(isok))
  
}


## Check Cultivars Genes Table
# col :
#  rownames : Cultivars Name
#  1:ncol : Genes names
#  ncol : remove button
# value 0 or 1
checkCultivarsGenesTable <- function(data){
  isok <- TRUE
  if( "delete" %in% colnames(data) ) data_tmp <- data[,-ncol(data),drop=FALSE]
  else data_tmp <- data
  shiny::removeUI(selector = "#cultivarsGenesValueError")
  if( sum(data_tmp != 0) + sum(data_tmp != 1) != nrow(data_tmp)*ncol(data_tmp) ) {
    showErrorMessage(id = "cultivarsGenesValueError", selectorafter= "#generateLandscape",
                     message = paste0("Values in the 'Cultivars and Genes' table should be either 0 or 1"))
    isok <- FALSE
  }
  
  return(invisible(isok))
}

## Check Genes Table
# col :
#  1 : Genes name Name
#  2:ncol : Genes parameters
checkGenesTable <- function(data){
  isok <- TRUE
  #if(nrow(data) == 0 || sum(is.na(data) > 0)) return(invisible(isok))
  
  shiny::removeUI(selector = "#GenesNameError")
  if( sum(as.character(data[,1]) == "") != 0 | sum(grepl("^\\s*$",as.character(data[,1]))) != 0) {
    
    showErrorMessage(id = "GenesNameError", selectorafter= "#generateLandscape",
                     message = "Gene names must be strings")
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#GenesNegatifError")
  if( sum(data[,- which(colnames(data) %in% c("genesName","target_trait"))] < 0) != 0 ) {
    
    showErrorMessage(id = "GenesNegatifError", selectorafter= "#generateLandscape",
                     message = "Values in the 'Genes' table should be >= 0")
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#GenesValueMaxError")
  if( sum(data[, c("age_of_activ_mean", "age_of_activ_var", "Nlevels_aggressiveness", "tradeoff_strength")] > VALUEMAX) != 0 ) {

    showErrorMessage(id = "GenesValueMaxError", selectorafter= "#generateLandscape",
                     message = paste0("Values in the 'Genes' table should be lower than ",VALUEMAX))
    isok <- FALSE
  }

  shiny::removeUI(selector = "#GenesStrictZeroError")
  if( sum(data[, c("tradeoff_strength")] <= 0) != 0 ) {
    
    showErrorMessage(id = "GenesStrictZeroError", selectorafter= "#generateLandscape",
                     message = paste0("Gene 'tradeoff_strength' values should be greater than 0"))
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#GenesUpper1Error")
  if( sum(data[, c("efficiency", "mutation_prob", "adaptation_cost", "relative_advantage", "recombination_sd")] > 1) != 0) {
    showErrorMessage(id = "GenesUpper1Error", selectorafter= "#generateLandscape",
                     message = paste0("Gene 'efficiency' / 'mutation_prob' / 'adaptation_cost' / 'relative_advantage' / 'recombination_sd' values should lower than 1"))
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#GenesZeroError")
  if( sum(data[, c("recombination_sd")] < 0) != 0 ) {
    
    showErrorMessage(id = "GenesZeroError", selectorafter= "#generateLandscape",
                     message = paste0("Gene 'recombination_sd' values should be greater or equal to 0"))
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#GenesTraitError")
  if( sum(data[, c("target_trait")] == "IR") +
      sum(data[, c("target_trait")] == "LAT") +
      sum(data[, c("target_trait")] == "PR") +
      sum(data[, c("target_trait")] == "IP") != nrow(data)) {
    showErrorMessage(id = "GenesTraitError", selectorafter= "#generateLandscape",
                     message = paste0("Gene 'target_trait' values can be either IR, LAT, PR, or IP"))
    isok <- FALSE
  }
  
  return(invisible(isok))
}


### CheckAll Tables
### Will check tables value visible by user
checkAllTables <- function(){
  isok <- TRUE
  
  isok <- isok && checkCroptypesTable(simul_params_croptypes())
  isok <- isok && checkCultivarsTable(simul_params_cultivars())
  isok <- isok && checkCultivarsGenesTable(simul_params_cultivarsgenes())
  isok <- isok && checkGenesTable(simul_params_genes())
  
  return(isok)
}

#
# Select cultivar type depending the pathogen name
# cultivarTypeDisease2type <- function(disease="no pathogen", type="growingHost"){
#     if(disease == "no pathogen"){
#         type <- switch( type,
#                        "growingHost" = "growingHost",
#                        "nongrowingHost" = "nongrowingHost",
#                        "nonCrop" = "nonCrop"
#                        )
#         return(type)
#     }
#     
#     if(disease == "rust"){
#         type <- switch( type,
#                 "growingHost" = "wheat",
#                 "nongrowingHost" = "nongrowingHost",
#                 "nonCrop" = "nonCrop"
#         )
#         return(type)
#     }
#     
#     if(disease == "mildew"){
#         type <- switch( type,
#                 "growingHost" = "grapevine",
#                 "nongrowingHost" = "nongrowingHost",
#                 "nonCrop" = "nonCrop"
#         )
#         return(type)
#     }
#   
#   if(disease == "sigatoka"){
#     type <- switch( type,
#                     "growingHost" = "banana",
#                     "nongrowingHost" = "nongrowingHost",
#                     "nonCrop" = "nonCrop"
#     )
#     return(type)
#   }
#     
#     return("")
# }

disease2CultivarType <- function(disease="rust", type="crop"){
  if (type=="crop"){
    cultivar <- switch(disease,
                       "rust" = "wheat",  
                       "mildew" = "grapevine",
                       "sigatoka" = "banana",
                       "CMV" = "pepper",
                       "")
  }else{ 
    cultivar <- "nonCrop"
  }

  return(cultivar)
}


## Get some parameters specific to a pathogen (pathosystem)
pathosystemParams <- function(disease = "rust") {
  
  options = list(
    nTSpY = 120
      
    )
    
  if(disease == "rust"){
    options$nTSpY = 120
  }
  
  if(disease == "mildew"){
    options$nTSpY = 120
  }
  
  if(disease == "sigatoka"){
    options$nTSpY = 182
  }
  
  printVerbose( paste0("set pathosystem parameters for ",disease), level=3)
  #print(options)
  
  return(options)
}

## Create inoculum vector for local infection
inoculumRandomLocal <- function(params, pI0_base=5E-4, nb_poly_inoc=1){
  
  Npatho <- prod(params@Genes$Nlevels_aggressiveness)
  Nhost <- nrow(params@Cultivars)
  Npoly <- nrow(params@Landscape)
  Npoly_inoc <- nb_poly_inoc  ## number of inoculated polygons
   ## whether the avr pathogen can infect the polygons
  compatible_poly <- getMatrixPolyPatho(params)[,1]
   ## random polygon picked among compatible ones
  id_poly <- sample(grep(1, compatible_poly), Npoly_inoc)
  pI0_poly <- as.numeric(1:Npoly %in% id_poly)  
  pI0 <- loadInoculum(params,
                       pI0_all=pI0_base,
                       pI0_host=c(1,rep(0, Nhost-1)),
                       pI0_patho=c(1,rep(0, Npatho-1)), 
                       pI0_poly=pI0_poly)
  return(pI0)
}

loadDemoMO <- function(params, disease){
  gene1 <- loadGene(name="majorGene_1", type="majorGene")
  gene2 <- loadGene(name="majorGene_2", type="majorGene")
    
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  
  susceptibleName <- disease2CultivarType(simul_pathogen())
  cultivar1 <- loadCultivar(name=susceptibleName, type=disease2CultivarType(disease))
  cultivar2 <- loadCultivar(name="Resistant1", type=disease2CultivarType(disease))
  cultivar3 <- loadCultivar(name="Resistant2", type=disease2CultivarType(disease))
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("majorGene_1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("majorGene_2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible_crop", "Resistant_crop_1", "Resistant_crop_2"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible_crop", susceptibleName)
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant_crop_1", "Resistant1")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant_crop_2", "Resistant2")
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoMI <- function(params,disease){
  gene1 <- loadGene(name="majorGene_1", type="majorGene")
  gene2 <- loadGene(name="majorGene_2", type="majorGene")
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type=disease2CultivarType(disease))
  cultivar2 <- loadCultivar(name="Resistant1", type=disease2CultivarType(disease))
  cultivar3 <- loadCultivar(name="Resistant2", type=disease2CultivarType(disease))
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("majorGene_1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("majorGene_2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible_crop", "Mixture_crop"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible_crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Mixture_crop", c("Resistant1","Resistant2"))
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoRO <- function(params,disease){
  gene1 <- loadGene(name="majorGene_1", type="majorGene")
  gene2 <- loadGene(name="majorGene_2", type="majorGene")
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type=disease2CultivarType(disease))
  cultivar2 <- loadCultivar(name="Resistant1", type=disease2CultivarType(disease))
  cultivar3 <- loadCultivar(name="Resistant2", type=disease2CultivarType(disease))
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("majorGene_1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("majorGene_2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible_crop", "Resistant_crop_1", "Resistant_crop_2"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible_crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant_crop_1", "Resistant1")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant_crop_2", "Resistant2")
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoPY <- function(params, disease){
  gene1 <- loadGene(name="majorGene_1", type="majorGene")
  gene2 <- loadGene(name="majorGene_2", type="majorGene")
  gene1$mutation_prob <- 1E-4
  gene2$mutation_prob <- 1E-4
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type=disease2CultivarType(disease))
  cultivar2 <- loadCultivar(name="Resistant", type=disease2CultivarType(disease))
  cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant", c("majorGene_1", "majorGene_2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible_crop", "Pyramid"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible_crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Pyramid", "Resistant")
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

# Create a numeric input made for positive integer value
IntegerInput <- function(inputId, label, value, max) {
  shiny::numericInput(
    inputId = inputId,
    label = label,
    value = value,
    min = 0,
    max = max,
    step = 1
  )
}

# Create a numeric input made for percentage value
PercentageInput <- function(inputId, label, value) {
  shiny::numericInput(
    inputId = inputId,
    label = label,
    value = value,
    min = 0.0,
    max = 1.0,
    step = 0.02
  )
}


#################################################################
### Tooltip message
#################################################################

SEX_PROPAGULE_VIABILITY_LIMIT <- "Maximum number of cropping seasons up to which a sexual propagule is viable"
SEX_PROPAGULE_RELEASE_MEAN <- "Average number of cropping seasons after which a sexual propagule is released."
SIGMOID_SIGMA <- "Sigma parameter of the sigmoid contamination function (0 to relax density-dependence, 1 for linear dependence)"
SIGMOID_KAPPA <- "Kappa parameter of the sigmoid contamination function"
INFECTIOUS_PERIOD_VAR <- "Variance of the infectious period duration"
INFECTIOUS_PERIOD_MEAN <- "Maximal expected infectious period duration"
LATENT_PERIOD_VAR <- "Variance of the latent period duration"
LATENT_PERIOD_MEAN <- "Minimal expected latent period duration"
PROPAGULE_PROD_RATE <- "Maximal expected effective propagule production rate per timestep and per infectious individual"
INFECTION_RATE <- "Maximal expected infection rate of a propagule on a healthy individual"
SURVIVAL_PROB <- "Off-season survival probability of a propagule"
INOCULUM <- "Initial probability for the first susceptible host (usually indexed by 0) to be infectious (state I) at the beginning of the simulation"
GENERATE_LANDSCAPE <- "Generates a landscape composed of fields where croptypes are allocated with controlled proportions and spatio-temporal aggregation"
RUN_SIMULATION <- "Run the simulation (depending of the parameters it can be long)"
STOP_SIMULATION <- "Force to stop the simulation"
EXPORT_SIMULATION <- "Download a GPKG and txt files containing most of the parameters"
ROTATION_PERIOD <- "Rotation period, in years, between two configurations: (1) croptypes 0 and 1; and (2) croptypes 0 and 2. If Rotation period is 0, there is no rotation"

CULTIVARS_TOOLTIP <- c("Name of the cultivar (DO NOT INCLUDE SPACES)",
                      "Host individuals density (in pure crop) per surface unit at the beginning of the cropping season",
                      "Maximum host individuals density (in pure crop) per surface unit at the end of the cropping season",
                      "Growth rate of the host",
                      "Reproduction rate",
                      "Theoretical yield in pure crop (in weight or volume unit / ha / cropping season) associated with the sanitary status ‘H’",
                      "Theoretical yield in pure crop (in weight or volume unit / ha / cropping season) associated with the sanitary status ‘L’",
                      "Theoretical yield in pure crop (in weight or volume unit / ha / cropping season) associated with the sanitary status ‘I’",
                      "Theoretical yield in pure crop (in weight or volume unit / ha / cropping season) associated with the sanitary status ‘R’",
                      "Planting costs in pure crop (in monetary units / ha / cropping season)",
                      "Market value of the product (in monetary units / weight or volume unit)")
GENES_TOOLTIP <- c("Name of the resistance gene",
                   "Efficiency of the resistance gene (percentage of reduction of the targeted aggressiveness component: IR, 1/LAT, IP or PR)",
                   "Expected delay to resistance activation (for APRs)",
                   "Variance of the delay to resistance activation (for APRs)",
                   "Probability for a pathogenicity gene to mutate",
                   "Number of adaptation levels related to each resistance gene (i.e. 1 + number of required mutations for a pathogenicity gene to fully adapt to the corresponding resistance gene)",
                   "Fitness penalty paid by a pathogen genotype fully adapted to the resistance gene on all hosts",
                   "Fitness advantage of a pathogen genotype fully adapted to the resistance gene on hosts carrying this gene, relative to those that do not carry this gene",
                   "Strength of the trade-off relationship between the level of aggressiveness hosts that do and do not carry the resistance gene",
                   "Aggressiveness component targeted by the resistance gene",
                   "Variance for sexual recombination (only QTL)")
