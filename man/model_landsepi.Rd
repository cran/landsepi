% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{model_landsepi}
\alias{model_landsepi}
\title{Model for Landscape Epidemiology & Evolution}
\usage{
model_landsepi(
  time_param,
  area_vector,
  rotation_matrix,
  croptypes_cultivars_prop,
  dispersal,
  inits,
  seed,
  cultivars_param,
  basic_patho_param,
  genes_param,
  treatment_param
)
}
\arguments{
\item{time_param}{list of simulation parameters:\itemize{
\item Nyears = number cropping seasons,
\item nTSpY = number of time-steps per cropping season.
}}

\item{area_vector}{a vector containing areas of polygons (i.e. fields), in surface units.}

\item{rotation_matrix}{a matrix containing for each field (rows) and year (columns, named "year_1", "year_2", etc.),
the index of the cultivated croptype. Importantly, the matrix must contain 1 more column than the real number
of simulated years.}

\item{croptypes_cultivars_prop}{a matrix with three columns named 'croptypeID' for croptype index,
'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype.}

\item{dispersal}{list of dispersal parameters:\itemize{
\item disp_patho_clonal = vectorised dispersal matrix of the pathogen (clonal propagules),
\item disp_patho_sex = vectorised dispersal matrix of the pathogen (sexual propagules),
\item disp_host = vectorised dispersal matrix of the host.
}}

\item{inits}{list of initial conditions:\itemize{
\item pI0 = vector of length Npoly\emph{Npatho}Nhost giving the probability to be infectious (i.e. state I) at
t=0 pr each polygon, pathogen genotype and host.
}}

\item{seed}{seed (for random number generation).}

\item{cultivars_param}{list of parameters associated with each host genotype (i.e. cultivars)
when cultivated in pure crops:\itemize{
\item initial_density = vector of host densities (per surface unit) at the beginning of the cropping season,
\item max_density = vector of maximum host densities (per surface unit) at the end of the cropping season,
\item growth rate = vector of host growth rates,
\item reproduction rate = vector of host reproduction rates,
\item relative_yield_H = Yield of H individuals relative to H individuals (100\%)
\item relative_yield_L = Yield of L individuals relative to H individuals
\item relative_yield_I = Yield of I individuals relative to H individuals
\item relative_yield_R = Yield of R individuals relative to H individuals
\item sigmoid_kappa_host = kappa parameter for the sigmoid invasion function (for host dispersal),
\item sigmoid_sigma_host = sigma parameter for the sigmoid invasion function (for host dispersal),
\item sigmoid_plateau_host = plateau parameter for the sigmoid invasion function (for host dispersal),
\item cultivars_genes_list = a list containing, for each host genotype, the indices of carried resistance genes,
}}

\item{basic_patho_param}{list of i. pathogen aggressiveness parameters on a susceptible host
for a pathogen genotype not adapted to resistance and ii. sexual reproduction parameters: \itemize{
\item infection_rate = maximal expected infection rate of a propagule on a healthy host,
\item propagule_prod_rate = maximal expected reproduction_rate of an infectious host per timestep,
\item latent_period_mean = minimal expected duration of the latent period,
\item latent_period_var = variance of the latent period duration,
\item infectious_period_mean = maximal expected duration of the infectious period,
\item infectious_period_var = variance of the infectious period duration,
\item survival_prob = matrix giving the probability for a propagule to survive the off-season, for each croptype (rows) and each year (columns)
\item repro_sex_prob = vector of probabilities for an infectious host to reproduce via sex rather than via cloning for each timestep,
\item sigmoid_kappa = kappa parameter of the sigmoid contamination function,
\item sigmoid_sigma = sigma parameter of the sigmoid contamination function,
\item sigmoid_plateau = plateau parameter of the sigmoid contamination function,
\item sex_propagule_viability_limit = maximum number of cropping seasons up to which a sexual propagule is viable
\item sex_propagule_release_mean = average number of cropping seasons after which a sexual propagule is released.
\item clonal_propagule_gradual_release = whether or not clonal propagules surviving the bottleneck are gradually released along the following cropping season.
}}

\item{genes_param}{list of parameters associated with each resistance gene and with the evolution of
each corresponding pathogenicity gene:\itemize{
\item target_trait = vector of aggressiveness components (IR, LAT, IP, or PR) targeted by resistance genes,
\item efficiency = vector of resistance gene efficiencies (percentage of reduction of the targeted
aggressiveness component: IR, 1/LAT, IP and PR),
\item age_of_activ_mean = vector of expected delays to resistance activation (for APRs),
\item age_of_activ_var = vector of variances of the delay to resistance activation (for APRs),
\item mutation_prob = vector of mutation probabilities for pathogenicity genes (each of them corresponding to a resistance gene),
\item Nlevels_aggressiveness = vector of number of adaptation levels related to each resistance gene (i.e. 1 + number
of required mutations for a pathogenicity gene to fully adapt to the corresponding resistance gene),
\item adaptation_cost = vector of adaptation penalties paid by pathogen genotypes fully adapted
to the considered resistance genes on all hosts,
\item relative_advantage = vector of fitness advantages of a pathogen genotype fully adapted to the resistance genes
on hosts carrying these genes, relative to those that do not carry these genes,
\item tradeoff_strength = vector of strengths of the trade-off relationships between the level of aggressiveness
on hosts that do and do not carry the resistance genes.
}}

\item{treatment_param}{list of parameters related to pesticide treatments: \itemize{
\item treatment_degradation_rate = degradation rate (per time step) of chemical concentration,
\item treatment_efficiency = maximal efficiency of chemical treatments (i.e. fractional reduction
of pathogen infection rate at the time of application),
\item treatment_timesteps = vector of time-steps corresponding to treatment application dates,
\item treatment_cultivars = vector of indices of the cultivars that receive treatments,
\item treatment_cost = cost of a single treatment application (monetary units/ha),
\item treatment_application_threshold = vector of thresholds (i.e. disease severity, one for each treated cultivar) above which the treatment is applied
}}
}
\value{
A set of binary files is generated for every year of simulation and every compartment:
\itemize{
\item H: healthy hosts,
\item Hjuv: juvenile healthy hosts (for host reproduction),
\item L: latently infected hosts,
\item I: infectious hosts,
\item R: removed hosts,
\item P: propagules.}
Each file indicates for every time-step the number of individuals in each field, and when
appropriate for each host and pathogen genotypes). Additionally, a binary file called TFI is
generated and gives the Treatment Frequency Indicator (expressed as the number of treatment applications
per polygon).
}
\description{
Stochastic, spatially-explicit, demo-genetic model simulating the spread and evolution of a
plant pathogen in a heterogeneous landscape.
}
\details{
See \code{?landsepi} for details on the model and assumptions.
Briefly, the model is stochastic, spatially explicit (the basic spatial unit is an individual field), based on a SEIR
(‘susceptible-exposed-infectious-removed’, renamed HLIR for 'healthy-latent-infectious-removed' to avoid confusions
with 'susceptible host') structure with a discrete time step. It simulates the spread and
evolution (via mutation, recombination through sexual reproduction, selection and drift)
of a pathogen in a heterogeneous cropping landscape, across cropping seasons split by host harvests which impose
potential bottlenecks to the pathogen. A wide array of resistance deployment strategies
(possibly including chemical treatments) can be simulated.
}
\examples{
\dontrun{
#### Spatially-implicit simulation with 2 patches (S + R) during 3 years ####

## Simulation parameters
time_param <- list(Nyears=3, nTSpY=120)
Npoly=2
Npatho=2
area <- c(100000, 100000)
basic_patho_param <- loadPathogen(disease = "rust")
basic_patho_param$repro_sex_prob <- rep(0, time_param$nTSpY+1)
     cultivars <- as.list(rbind(loadCultivar(name="Susceptible", type="growingHost")
, loadCultivar(name="Resistant", type="growingHost")))
names(cultivars)[names(cultivars)=="cultivarName"] <- "name"
yield0 <- cultivars$yield_H + as.numeric(cultivars$yield_H==0)
cultivars <- c(cultivars, list(relative_yield_H = as.numeric(cultivars$yield_H / yield0)
, relative_yield_L = as.numeric(cultivars$yield_L / yield0)
, relative_yield_I = as.numeric(cultivars$yield_I / yield0)
, relative_yield_R = as.numeric(cultivars$yield_R / yield0)
, sigmoid_kappa_host=0.002, sigmoid_sigma_host=1.001, sigmoid_plateau_host=1
, cultivars_genes_list=list(numeric(0),0)))
rotation <- data.frame(year_1=c(0,1), year_2=c(0,1), year_3=c(0,1), year_4=c(0,1))
croptypes_cultivars_prop <- data.frame(croptypeID=c(0,1), cultivarID=c(0,1), proportion=c(1,1))
genes <- as.list(loadGene(name="MG", type="majorGene"))
treatment=list(treatment_degradation_rate=0.1,
treatment_efficiency=0, 
treatment_timesteps=logical(0),
treatment_cultivars=logical(0),
treatment_cost=0,
treatment_application_threshold = logical(0))
  
## run simulation
model_landsepi(seed=1,
time_param = time_param,
basic_patho_param = basic_patho_param,
inits = list(pI0=c(0.1, rep(0, 7))),
area_vector = area,
dispersal = list(disp_patho_clonal=c(0.99,0.01,0.01,0.99),
disp_patho_sex=c(1,0,0,1),
disp_host=c(1,0,0,1)),
rotation_matrix = as.matrix(rotation),
croptypes_cultivars_prop = as.matrix(croptypes_cultivars_prop),
cultivars_param = cultivars, 
genes_param = genes,
treatment_param = treatment)

## Compute outputs
eco_param <- list(yield_perHa = cbind(H = as.numeric(cultivars$relative_yield_H),
L = as.numeric(cultivars$relative_yield_L),
I = as.numeric(cultivars$relative_yield_I),
R = as.numeric(cultivars$relative_yield_R)),
planting_cost_perHa = as.numeric(cultivars$planting_cost),
market_value = as.numeric(cultivars$market_value))

evol_res <- evol_output(, time_param, Npoly, cultivars, genes)
epid_res <-  epid_output(, time_param, Npatho, area, rotation
, croptypes_cultivars_prop, cultivars, eco_param, treatment, basic_patho_param)



#### 1-year simulation of a rust epidemic in pure susceptible crop in a single 1-km2 patch ####
## Simulation and pathogen parameters
time_param <- list(Nyears=1, nTSpY=120)
area <- c(1E6)
basic_patho_param = loadPathogen(disease = "rust")
basic_patho_param$repro_sex_prob <- rep(0, time_param$nTSpY+1)
## croptypes, cultivars and genes
rotation <- data.frame(year_1=c(0), year_2=c(0))
croptypes_cultivars_prop <- data.frame(croptypeID=c(0), cultivarID=c(0), proportion=c(1))
       cultivars <- as.list(rbind(loadCultivar(name="Susceptible", type="growingHost")))
names(cultivars)[names(cultivars)=="cultivarName"] <- "name"
yield0 <- cultivars$yield_H + as.numeric(cultivars$yield_H==0)
cultivars <- c(cultivars, list(relative_yield_H = as.numeric(cultivars$yield_H / yield0)
, relative_yield_L = as.numeric(cultivars$yield_L / yield0)
    , relative_yield_I = as.numeric(cultivars$yield_I / yield0)
, relative_yield_R = as.numeric(cultivars$yield_R / yield0)
, sigmoid_kappa_host=0.002, sigmoid_sigma_host=1.001, sigmoid_plateau_host=1
, cultivars_genes_list=list(numeric(0))))
genes <-   list(geneName = character(0) , adaptation_cost = numeric(0)
, relative_advantage = numeric(0)
, mutation_prob = numeric(0)
, efficiency = numeric(0) , tradeoff_strength = numeric(0)
, Nlevels_aggressiveness = numeric(0)
, age_of_activ_mean = numeric(0) , age_of_activ_var = numeric(0)
, target_trait = character(0)
, recombination_sd = numeric(0))
treatment=list(treatment_degradation_rate=0.1
                              , treatment_efficiency=0
, treatment_timesteps=logical(0)
, treatment_cultivars=logical(0)
, treatment_cost=0
, treatment_application_threshold = logical(0))

## run simulation
model_landsepi(seed=1, time_param = time_param
, basic_patho_param = basic_patho_param
, inits = list(pI0=5E-4), area_vector = area
, dispersal = list(disp_patho_clonal=c(1), disp_patho_sex=c(1), disp_host=c(1))
, rotation_matrix = as.matrix(rotation)
, treatment_param = treatment
                                , croptypes_cultivars_prop = as.matrix(croptypes_cultivars_prop)
, cultivars_param = cultivars,  genes_param = genes)
 }
}
\references{
Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018).
Assessing the durability andefficiency of landscape-based strategies to deploy plant
resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
}
