#library(shinycssloaders)
library(shiny)
library(shinyBS)
library(DT)
library(shinyjs)
#library(slickR)
library(gridExtra)
library(png)
library(grid)
library(future)
library(promises)
library(tools)
library(shinyalert)

# Active debug level
# 0 : no print
# 1 : warning
# 2 : error
# 3 : all
ACTIVE_DEBUG <- 3

library("landsepi")
data(package = "landsepi")

source("modules/editableDT.R")

VALUEMAX <- 1000


## del all file and directory of a path
cleanDir <- function(path) {
  files <- dir(path, full.names=TRUE, no..=TRUE)
  lapply(files, FUN = function(file){
      if( dir.exists(file) ) cleanDir(file)
      file.remove(file)
    })
}

ROOT_PATH <- getwd()

if(!dir.exists(paste0(ROOT_PATH,"/www/tmp/"))) dir.create(paste0(ROOT_PATH,"/www/tmp/"))
setwd(paste0(ROOT_PATH,"/www/tmp/"))

cleanDir(paste0(ROOT_PATH,"/www/tmp/"))

## User mode
advanced_mode <- reactiveVal(FALSE)

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

##################################################################
# Functions
##################################################################

# Print message in terminal
# msg : message to print
# level : level of output (0 no output, 3 all output)
printVerbose <- function(msg, level=ACTIVE_DEBUG) {
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
      paste0(message)
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
                     message = "Cultivar names must be strings")
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
  if( sum(data[,c("growth_rate","reproduction_rate","death_rate")] > 1) != 0) {
    
    showErrorMessage(id = "cultivarsmaxvalueError", selectorafter= "#generateLandscape",
                     message = paste0("Cultivar 'growth_rate' / 'reproduction_rate' / 'death_rate' values should be lower than 1"))
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#cultivarsValueError")
  if( sum(data[,- which(c("cultivarName", "growth_rate","reproduction_rate","death_rate") %in% colnames(data))] < 0) != 0
      || sum(data[,- which(c("cultivarName", "growth_rate","reproduction_rate","death_rate") %in% colnames(data))] > VALUEMAX) != 0) {
    
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
# value 0 or 1
checkCultivarsGenesTable <- function(data){
  isok <- TRUE
  
  shiny::removeUI(selector = "#cultivarsGenesValueError")
  if( sum(data != 0) + sum(data != 1) != nrow(data)*ncol(data) ) {
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
                     message = "Values in the 'Genes' table should be > 0")
    isok <- FALSE
  }
  
  shiny::removeUI(selector = "#GenesValueMaxError")
  if( sum(data[, c("time_to_activ_exp", "time_to_activ_var", "Nlevels_aggressiveness", "tradeoff_strength")] > VALUEMAX) != 0 ) {

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
  if( sum(data[, c("efficiency", "mutation_prob", "fitness_cost")] > 1) != 0) {
    showErrorMessage(id = "GenesUpper1Error", selectorafter= "#generateLandscape",
                     message = paste0("Gene 'efficiency' / 'mutation_prob' / 'fitness_cost' values should lower than 1"))
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


loadDemoMO <- function(params){
  gene1 <- loadGene(name="MG 1", type="majorGene")
  gene2 <- loadGene(name="MG 2", type="majorGene")
    
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type="growingHost")
  cultivar2 <- loadCultivar(name="Resistant1", type="growingHost")
  cultivar3 <- loadCultivar(name="Resistant2", type="growingHost")
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("MG 1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("MG 2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible crop", "Resistant crop 1", "Resistant crop 2"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoMI <- function(params){
  gene1 <- loadGene(name="MG 1", type="majorGene")
  gene2 <- loadGene(name="MG 2", type="majorGene")
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type="growingHost")
  cultivar2 <- loadCultivar(name="Resistant1", type="growingHost")
  cultivar3 <- loadCultivar(name="Resistant2", type="growingHost")
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("MG 1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("MG 2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible crop", "Mixture"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Resistant1","Resistant2"))
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoRO <- function(params){
  gene1 <- loadGene(name="MG 1", type="majorGene")
  gene2 <- loadGene(name="MG 2", type="majorGene")
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type="growingHost")
  cultivar2 <- loadCultivar(name="Resistant1", type="growingHost")
  cultivar3 <- loadCultivar(name="Resistant2", type="growingHost")
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("MG 1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("MG 2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible crop", "Resistant crop 1", "Resistant crop 2"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoPY <- function(params){
  gene1 <- loadGene(name="MG 1", type="majorGene")
  gene2 <- loadGene(name="MG 2", type="majorGene")
  gene1$mutation_prob <- 1E-4
  gene2$mutation_prob <- 1E-4
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type="growingHost")
  cultivar2 <- loadCultivar(name="Resistant", type="growingHost")
  cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant", c("MG 1", "MG 2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible crop", "Pyramid"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
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
SIGMOID_SIGMA <- "Sigma parameter of the sigmoid contamination function (0 to relax density-dependence, 1 for linear dependence)"
SIGMOID_KAPPA <- "Kappa parameter of the sigmoid contamination function"
INFECTIOUS_PERIOD_VAR <- "Variance of the infectious period duration"
INFECTIOUS_PERIOD_EXP <- "Maximal expected infectious period duration"
LATENT_PERIOD_VAR <- "Variance of the latent period duration"
LATENT_PERIOD_EXP <- "Minimal expected latent period duration"
PROPAGULE_PROD_RATE <- "Maximal expected effective propagule production rate per timestep and per infectious individual"
INFECTION_RATE <- "Maximal expected infection rate of a propagule on a healthy individual"
SURVIVAL_PROB <- "Off-season survival probability of a propagule"
INOCULUM <- "Initial probability for the first susceptible host (usually indexed by 0) to be infectious (state I) at the beginning of the simulation"
GENERATE_LANDSCAPE <- "Generates a landscape composed of fields where croptypes are allocated with controlled proportions and spatio-temporal aggregation"
RUN_SIMULATION <- "Run the simulation (depending of the parameters it can be long)"
STOP_SIMULATION <- "Force to stop the simulation"
EXPORT_SIMULATION <- "Download a GPKG containing most of the parameters"
ROTATION_PERIOD <- "Rotation period, in years, between two configurations: (1) croptypes 0 and 1; and (2) croptypes 0 and 2. If Rotation period is 0, there is no rotation"

CULTIVARS_TOOLTIP <- c("Name of the cultivar",
                      "Host individuals density (in pure crop) per surface unit at the beginning of the cropping season",
                      "Maximum host individuals density (in pure crop) per surface unit at the end of the cropping season",
                      "Growth rate",
                      "Reproduction rate",
                      "Death rate",
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
                   "maximal fitness penalty paid by a pathogen genotype fully adapted to the resistance gene on hosts that do not carry this gene",
                   "Strength of the trade-off relationship between the level of aggressiveness hosts that do and do not carry the resistance gene",
                   "Aggressiveness component targeted by the resistance gene")