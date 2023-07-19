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

#' @title Simulation with input parameters as data.frames.
#' @name simul_landsepi
#' @description Stochastic, spatially-explicit, demo-genetic model simulating the spread and evolution
#' of a pathogen in a heterogeneous landscape and generating a wide range of epidemiological, evolutionary
#' and economic outputs.
#' @param seed an integer used as seed value (for random number generator).
#' @param time_param a list of simulation parameters:\itemize{
#' \item Nyears = number cropping seasons,
#' \item nTSpY = number of time-steps per cropping season.
#' }
#' @param croptype_names a vector of croptypes names.
#' @param croptypes_cultivars_prop a dataframe with three columns named 'croptypeID' for croptype index,
#' 'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype.
#' @param cultivars a dataframe of parameters associated with each host genotype (i.e. cultivars)
#' when cultivated in pure crops. Columns of the dataframe are:\itemize{
#' \item cultivarName: cultivar names,
#' \item initial_density: host densities (per square meter) at the beginning of the cropping season 
#' as if cultivated in pure crop,
#' \item max_density: maximum host densities (per square meter) at the end of the cropping season 
#' as if cultivated in pure crop,
#' \item growth_rate: host growth rates,
#' \item reproduction rate: host reproduction rates,
#' \item yield_H: theoretical yield (in weight or volume units / ha / cropping season) associated with 
#' hosts in sanitary status H as if cultivated in pure crop,
#' \item yield_L: theoretical yield (in weight or volume units / ha / cropping season) associated with 
#' hosts in sanitary status L as if cultivated in pure crop,
#' \item yield_I: theoretical yield (in weight or volume units / ha / cropping season) associated with 
#' hosts in sanitary status I as if cultivated in pure crop,
#' \item yield_R: theoretical yield (in weight or volume units / ha / cropping season) associated with 
#' hosts in sanitary status R as if cultivated in pure crop,
#' \item planting_cost = planting costs (in monetary units / ha / cropping season) as if cultivated in pure crop,
#' \item market_value = market values of the production (in monetary units / weight or volume unit).
#' }
#' @param cultivars_genes_list a list containing, for each host genotype, the indices of carried resistance genes.
#' @param genes a data.frame of parameters associated with each resistance gene and with the evolution of
#' each corresponding pathogenicity gene. Columns of the dataframe are:\itemize{
#' \item geneName: names of resistance genes,
#' \item target_trait: aggressiveness components (IR, LAT, IP, or PR) targeted by resistance genes,
#' \item efficiency: resistance gene efficiencies (percentage of reduction of targeted aggressiveness components: 
#' IR, 1/LAT, IP and PR),
#' \item age_of_activ_mean: expected delays to resistance activation (for APRs),
#' \item age_of_activ_var: variances of the delay to resistance activation (for APRs),
#' \item mutation_prob: mutation probabilities for pathogenicity genes (each of them corresponding to a resistance gene),
#' \item Nlevels_aggressiveness: number of adaptation levels related to each resistance gene (i.e. 1 + number
#' of required mutations for a pathogenicity gene to fully adapt to the corresponding resistance gene),
#' \item adaptation_cost: fitness penalties paid by pathogen genotypes fully adapted
#' to the considered resistance genes on host that do not carry the resistance genes,
#' \item tradeoff_strength: strengths of the trade-off relationships between the level of aggressiveness
#' on hosts that do and do not carry the resistance genes.
#' \item recombination_sd: standard deviation of the normal distribution used for recombination of quantitative traits during sexual reproduction (infinitesimal model)
#' }
#' @param landscape a sp object containing the landscape (required only if videoMP4=TRUE).
#' @param area a vector containing polygon areas (must be in square meters).
#' @param rotation a dataframe containing for each field (rows) and year (columns, named "year_1", "year_2", etc.),
#' the index of the cultivated croptype. Importantly, the matrix must contain 1 more column than the real number
#' of simulated years.
#' @param basic_patho_param a list of i. pathogen aggressiveness parameters on a susceptible host
#' for a pathogen genotype not adapted to resistance and ii. sexual reproduction parameters: \itemize{
#' \item infection_rate = maximal expected infection rate of a propagule on a healthy host,
#' \item propagule_prod_rate = maximal expected effective propagule production rate of an infectious host per time step,
#' \item latent_period_mean = minimal expected duration of the latent period,
#' \item latent_period_var = variance of the latent period duration,
#' \item infectious_period_mean = maximal expected duration of the infectious period,
#' \item infectious_period_var = variance of the infectious period duration,
#' \item survival_prob = probability for a propagule to survive the off-season,
#' \item repro_sex_prob = probability for an infectious host to reproduce via sex rather than via cloning,
#' \item sigmoid_kappa = kappa parameter of the sigmoid contamination function,
#' \item sigmoid_sigma = sigma parameter of the sigmoid contamination function,
#' \item sigmoid_plateau = plateau parameter of the sigmoid contamination function,
#' \item sex_propagule_viability_limit = maximum number of cropping seasons up to which a sexual propagule is viable
#' \item sex_propagule_release_mean = average number of seasons after which a sexual propagule is released,
#' \item clonal_propagule_gradual_release = Whether or not clonal propagules surviving the bottleneck are gradually released along the following cropping season.
#' }
#' @param repro_sex_prob a vector of size \code{time_param$nSTpY+1} giving the probability for an infectious host to reproduce
#'  via sex rather than via cloning for every time step,
#' @param disp_patho_clonal a vectorized matrix giving the probability of pathogen dispersal
#' from any field of the landscape to any other field.
#' @param disp_patho_sex a vectorized matrix giving the probability of pathogen dispersal for sexual propagules
#' from any field of the landscape to any other field.
#' @param disp_host a vectorized matrix giving the probability of host dispersal
#' from any field of the landscape to any other field
#' @param treatment list of parameters related to pesticide treatments: \itemize{ 
#' \item treatment_degradation_rate = degradation rate (per time step) of chemical concentration,
#' \item treatment_efficiency = maximal efficiency of chemical treatments (i.e. fractional reduction 
#' of pathogen infection rate at the time of application),
#' \item treatment_timesteps = vector of time-steps corresponding to treatment application dates,
#' \item treatment_cultivars = vector of indices of the cultivars that receive treatments,
#' \item treatment_cost = cost of a single treatment application (monetary units/ha)
#' \item treatment_application_threshold = vector of thresholds (i.e. disease severity, one for each treated cultivar) 
#' above which the treatment is applied in a polygon
#' }
#' @param pI0 initial probability for the first host (whose index is 0) to be infectious (i.e. state I)
#' at the beginning of the simulation. Must be between 0 and 1.
#' @param epid_outputs a character string (or a vector of character strings if several outputs are to be computed)
#' specifying the type of epidemiological and economic outputs to generate (see details):  
#' \itemize{
#' \item "audpc" : Area Under Disease Progress Curve (average number of diseased host individuals
#' per time step and square meter) 
#' \item "audpc_rel" : Relative Area Under Disease Progress Curve (average proportion of diseased host
#' individuals relative to the total number of existing hosts)
#' \item "gla" : Green Leaf Area (average number of healthy host individuals per time step and square meter)
#' \item "gla_rel" : Relative Green Leaf Area (average proportion of healthy host individuals relative to the 
#' total number of existing hosts)
#' \item "eco_yield" : total crop yield (in weight or volume units per ha) 
#' \item "eco_cost" : operational crop costs (in monetary units per ha) 
#' \item "eco_product" : total crop products (in monetary units per ha) 
#' \item "eco_margin" : Margin (products - operational costs, in monetary units per ha) 
#' \item "contrib": contribution of pathogen genotypes to LIR dynamics
#' \item "HLIR_dynamics", "H_dynamics", "L_dynamics", "IR_dynamics", "HLI_dynamics", etc.: Epidemic dynamics
#' related to the specified sanitary status (H, L, I or R and all their combinations). Graphics only,
#' works only if graphic=TRUE.
#' \item "all" : compute all these outputs (default)
#' \item "" : none of these outputs will be generated.
#' }
#' @param evol_outputs a character string (or a vector of character strings if several outputs are to be computed)
#' specifying the type of evolutionary outputs to generate :\itemize{
#' \item "evol_patho": Dynamics of pathogen genotype frequencies
#' \item "evol_aggr": Evolution of pathogen aggressiveness
#' \item "durability": Durability of resistance genes
#' \item "all": compute all these outputs (default)
#' \item "": none of these outputs will be generated.
#' }
#' @param thres_breakdown an integer (or vector of integers) giving the threshold (i.e. number of infections)
#' above which a pathogen genotype is unlikely to go extinct, used to characterise the time to invasion
#' of resistant hosts (several values are computed if several thresholds are given in a vector).
#' @param audpc100S the audpc in a fully susceptible landscape (used as reference value for graphics).
#' @param writeTXT a logical indicating if outputs must be written in text files (TRUE, default) or not (FALSE).
#' @param graphic a logical indicating if graphics must be generated (TRUE, default) or not (FALSE).
#' @param videoMP4 a logical indicating if a video must be generated (TRUE) or not (FALSE, default).
#' Works only if graphic=TRUE and epid_outputs="audpc_rel" (or epid_outputs="all").
#' @param keepRawResults a logical indicating if binary files must be kept after the end of the simulation (default=FALSE).
#' Careful, many files may be generated if keepRawResults=TRUE.
#' @details See ?landsepi for details on the model and assumptions.  
#' Briefly, the model is stochastic, spatially explicit (the basic spatial unit is an individual field), based on a SEIR
#'  (‘susceptible-exposed-infectious-removed’, renamed HLIR for 'healthy-latent-infectious-removed' to avoid confusions
#'  with 'susceptible host') structure with a discrete time step. It simulates the spread and
#'  evolution (via mutation, recombination through sexual reproduction, selection and drift) 
#'  of a pathogen in a heterogeneous cropping landscape, across cropping seasons split by host harvests which impose
#'  potential bottlenecks to the pathogen. A wide array of resistance deployment strategies 
#'  (possibly including chemical treatments) can be simulated and evaluated using several possible 
#'  outputs to assess the epidemiological, evolutionary and economic performance
#'  of deployment strategies (See ?epid_output and ?evol_output for details).
#'
#' @return A list containing all outputs that have been required via "epid_outputs" and "evol_outputs".
#' A set of text files, graphics and a video showing epidemic dynamics can be generated.
#' If keepRawResults=TRUE, a set of binary files is generated for every year of simulation and
#' every compartment: \itemize{
#'  \item H: healthy hosts,
#'  \item Hjuv: juvenile healthy hosts (for host reproduction),
#'  \item L: latently infected hosts,
#'  \item I: infectious hosts,
#'  \item R: removed hosts,
#'  \item P: propagules.}
#' Each file indicates for every time-step the number of individuals in each field, and when appropriate for
#' each host and pathogen genotype. Additionally, a binary file called TFI is 
#' generated and gives the Treatment Frequency Indicator (expressed as the number of treatment applications 
#'  per polygon).
#' @seealso \link{model_landsepi}, \link{epid_output}, \link{evol_output}, \link{video}, \link{runSimul}
#' @examples
#' \dontrun{
#' #### Spatially-implicit simulation with a single 1-km^2 patch 100% cultivated 
#' # with a susceptible cultivar
#' 
#' simul_landsepi()
#' 
#' #### Spatially-implicit simulation with 2 patches (S + R) during 3 years ####
#'
#' ## Simulation parameters
#' time_param <- list(Nyears = 3, nTSpY = 120)
#' area <- c(100000, 100000)
#' rotation <- data.frame(year_1 = c(0, 1), year_2 = c(0, 1), year_3 = c(0, 1), year_4 = c(0, 1))
#' croptype_names <- c("Susceptible crop", "Resistant crop")
#' croptypes_cultivars_prop <- data.frame(
#'   croptypeID = c(0, 1),
#'   cultivarID = c(0, 1),
#'   proportion = c(1, 1)
#' )
#' cultivars <- rbind(
#'   loadCultivar(name = "Susceptible", type = "growingHost"),
#'   loadCultivar(name = "Resistant", type = "growingHost")
#' )
#' genes <- loadGene(name = "MG", type = "majorGene")
#' cultivars_genes_list <- list(numeric(0), 0)
#'
#' ## Run simulation
#' simul_landsepi(
#'   seed = 12345, time_param, croptype_names, croptypes_cultivars_prop, cultivars,
#'   cultivars_genes_list, genes, landscape = NULL, area, rotation,
#'   basic_patho_param = loadPathogen(disease = "rust"),
#'   repro_sex_prob = rep(0, time_param$nTSpY+1),
#'   disp_patho_clonal = c(0.99, 0.01, 0.01, 0.99),
#'   disp_patho_sex = c(0.99, 0.01, 0.01, 0.99),
#'   disp_host = c(1, 0, 0, 1),
#'   pI0 = 5e-4
#' )
#'
#'
#' #### Spatially-explicit simulation with built-in landscape during 10 years ####
#' # Generate a mosaic of four croptypes in balanced proportions
#' # and medium level of spatial aggregation
#'
#' ## Simulation and Landscape parameters
#' Nyears <- 10
#' landscape <- loadLandscape(1)
#' Npoly <- length(landscape)
#' library(sf)
#' area <- st_area(st_as_sf(landscape))
#' rotation <- AgriLand(landscape, Nyears,
#'   rotation_period = 1, rotation_realloc = FALSE,
#'   rotation_sequence = c(0, 1, 2, 3),
#'   prop = rep(1 / 4, 4), aggreg = 0.5, graphic = TRUE, outputDir = getwd()
#' )
#' rotation <- data.frame(rotation)[, 1:(Nyears + 1)]
#' croptype_names <- c("Susceptible crop"
#' , "Resistant crop 1"
#' , "Resistant crop 2"
#' , "Resistant crop 3")
#' croptypes_cultivars_prop <- data.frame(croptypeID = c(0, 1, 2, 3), cultivarID = c(0, 1, 2, 3),
#'                                        proportion = c(1, 1, 1, 1))
#' cultivars <- data.frame(rbind(
#'   loadCultivar(name = "Susceptible", type = "growingHost"),
#'   loadCultivar(name = "Resistant1", type = "growingHost"),
#'   loadCultivar(name = "Resistant2", type = "growingHost"),
#'   loadCultivar(name = "Resistant3", type = "growingHost")
#' ), stringsAsFactors = FALSE)
#' genes <- data.frame(rbind(
#'   loadGene(name = "MG 1", type = "majorGene"),
#'   loadGene(name = "MG 2", type = "majorGene"),
#'   loadGene(name = "MG 3", type = "majorGene")
#' ), stringsAsFactors = FALSE)
#' cultivars_genes_list <- list(numeric(0), 0, 1, 2)
#'
#' ## Run simulation
#' simul_landsepi(
#'   seed = 12345, time_param = list(Nyears = Nyears, nTSpY = 120),
#'   croptype_names, croptypes_cultivars_prop, cultivars,
#'   cultivars_genes_list, genes, landscape, area, rotation,
#'   basic_patho_param = loadPathogen(disease = "rust"),
#'   repro_sex_prob = rep(0, 120+1),
#'   disp_patho_clonal = loadDispersalPathogen(1)[[1]],
#'   disp_patho_sex = as.numeric(diag(Npoly)),
#'   disp_host = as.numeric(diag(Npoly)),
#'   pI0 = 5e-4
#' )
#' }
#'
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018).
#' Assessing the durability andefficiency of landscape-based strategies to deploy 
#' plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#' @include RcppExports.R
#' @export
simul_landsepi <- function(seed = 12345, time_param = list(Nyears = 5, nTSpY = 120)
                           , croptype_names=c("Susceptible crop")
                           , croptypes_cultivars_prop = data.frame(croptypeID = 0, cultivarID = 0, proportion = 1)
                           , cultivars = data.frame(cultivarName = "Susceptible", initial_density = 0.1, max_density = 2.0
                                              , growth_rate = 0.1, reproduction_rate = 0.0 #, death_rate = 0.0
                                              , yield_H = 2.5, yield_L = 0.0, yield_I = 0.0, yield_R = 0.0
                                              , planting_cost = 225, market_value = 200)
                           , cultivars_genes_list=list(numeric(0))
                           , genes = data.frame(geneName = character(0), mutation_prob = numeric(0)
                                             , efficiency = numeric(0) , tradeoff_strength = numeric(0)
                                             , Nlevels_aggressiveness = numeric(0), adaptation_cost = numeric(0)
                                             , age_of_activ_mean = numeric(0) , age_of_activ_var = numeric(0)
                                             , target_trait = character(0)
                                             , recombination_sd = numeric(0))
                           , landscape = NULL, area=1E6
                           , rotation = data.frame(year_1 = c(0), year_2 = c(0), year_3 = c(0), year_4 = c(0), year_5 = c(0), year_6 = c(0))
                           , basic_patho_param = list(name = "rust", survival_prob = 1e-4, repro_sex_prob = 0
                                                      , infection_rate = 0.4, propagule_prod_rate = 3.125
                                                      , latent_period_mean = 10, latent_period_var = 9
                                                      , infectious_period_mean = 24, infectious_period_var = 105
                                                      , sigmoid_kappa = 5.333, sigmoid_sigma = 3, sigmoid_plateau = 1
                                                      , sex_propagule_viability_limit  = 1, sex_propagule_release_mean = 1
                                                      , clonal_propagule_gradual_release = 0)
                           , repro_sex_prob=rep(0, time_param$nTSpY+1)
                           , disp_patho_clonal=c(1), disp_patho_sex=c(1), disp_host=c(1)
                           , treatment=list(treatment_degradation_rate=0.1, 
                                          treatment_efficiency=0, 
                                          treatment_timesteps=logical(0),
                                          treatment_cultivars=logical(0),
                                          treatment_cost=0,
                                          treatment_application_threshold = logical(0))
                           , pI0 = 5e-4
                           , epid_outputs = "all", evol_outputs = "all", thres_breakdown = 50000
                           , audpc100S = 0.76 #8.48 for downy mildew,
                           , writeTXT = TRUE, graphic = TRUE, videoMP4 = FALSE, keepRawResults = FALSE) {

  # Host parameters
  cultivars_param <- list(name = as.character(cultivars$cultivarName),
                          initial_density = as.numeric(cultivars$initial_density),
                          max_density = as.numeric(cultivars$max_density),
                          growth_rate = as.numeric(cultivars$growth_rate),
                          reproduction_rate = as.numeric(cultivars$reproduction_rate),
                          # death_rate = as.numeric(cultivars$death_rate),
                          sigmoid_kappa_host = 0,
                          sigmoid_sigma_host = 1,
                          sigmoid_plateau_host = 1,
                          cultivars_genes_list = cultivars_genes_list,
                          ## Yield relative to H (with a security to avoid 0/0)
                          relative_yield_H = as.numeric(cultivars$yield_H / (cultivars$yield_H + as.numeric(cultivars$yield_H==0))),
                          relative_yield_L = as.numeric(cultivars$yield_L / (cultivars$yield_H + as.numeric(cultivars$yield_H==0))),
                          relative_yield_I = as.numeric(cultivars$yield_I / (cultivars$yield_H + as.numeric(cultivars$yield_H==0))),
                          relative_yield_R = as.numeric(cultivars$yield_R / (cultivars$yield_H + as.numeric(cultivars$yield_H==0)))
                          
  )
  
  # Evolution parameters
  genes_param <- list(name = as.character(genes$geneName),
                      adaptation_cost = as.numeric(genes$adaptation_cost),
                      mutation_prob = as.numeric(genes$mutation_prob),
                      efficiency = as.numeric(genes$efficiency),
                      tradeoff_strength = as.numeric(genes$tradeoff_strength),
                      Nlevels_aggressiveness = as.numeric(genes$Nlevels_aggressiveness),
                      age_of_activ_mean = as.numeric(genes$age_of_activ_mean),
                      age_of_activ_var = as.numeric(genes$age_of_activ_var),
                      target_trait = as.character(genes$target_trait),
                      recombination_sd = as.numeric(genes$recombination_sd)
  )

  # Economic parameters
  eco_param <- list(yield_perHa = cbind(H = as.numeric(cultivars$yield_H),
                                        L = as.numeric(cultivars$yield_L),
                                        I = as.numeric(cultivars$yield_I),
                                        R = as.numeric(cultivars$yield_R)),
                    planting_cost_perHa = as.numeric(cultivars$planting_cost),
                    market_value = as.numeric(cultivars$market_value)
  )

  
  basic_patho_param$repro_sex_prob <- repro_sex_prob
  
  ####
  #### Run the model (C++)
  ####
  
  print("Run the C++ model")

  model_landsepi(
    time_param = time_param,
    area_vector = area,
    rotation_matrix = as.matrix(rotation),
    croptypes_cultivars_prop = as.matrix(croptypes_cultivars_prop),
    dispersal = list(disp_patho_clonal = disp_patho_clonal, disp_patho_sex = disp_patho_sex, disp_host = disp_host),
    inits = list(pI0 = pI0),
    seed = seed,
    cultivars_param = cultivars_param,
    basic_patho_param = basic_patho_param,
    genes_param = genes_param,
    treatment_param = treatment
  )


  ####
  ####  Generate outputs
  ####
  print("Compute model outputs")
  Npoly <- length(area)

  
  ## Evolutionary output
  Npatho <- prod(genes_param$Nlevels_aggressiveness)
  if (Npatho > 1 & evol_outputs[1] != "") {
    evol_res <- evol_output(evol_outputs, time_param, Npoly, cultivars_param, genes_param, thres_breakdown, 
                            writeTXT = writeTXT, graphic = graphic)
  } else {
    evol_res <- NULL
  }

  ## Epidemiological output
  if (epid_outputs[1] != "") {
    epid_res <- epid_output(epid_outputs, time_param, Npatho, area, rotation, croptypes_cultivars_prop,
                            cultivars_param, eco_param, treatment, basic_patho_param, 
                            audpc100S, writeTXT = writeTXT, graphic = graphic)
  } else {
    epid_res <- NULL
  }
  

  ## Video
  if (videoMP4 & !is.null(epid_res[["audpc_rel"]])) {
    video(
      epid_res[["audpc_rel"]], time_param, Npatho, landscape, area, rotation, croptypes_cultivars_prop,
      croptype_names, cultivars_param
      # , audpc100S
      , evol_res$durability
    )
  }

  if (!keepRawResults) {
    print("remove binary files")
    file.remove(list.files(pattern = "\\.bin$")) # delete all the binary files
  }
  print(paste("model outputs stored in :", getwd()))

  return(list(evol_outputs = evol_res, epid_outputs = epid_res))
}
