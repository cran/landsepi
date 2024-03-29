% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{epid_output}
\alias{epid_output}
\title{Generation of epidemiological and economic model outputs}
\usage{
epid_output(
  types = "all",
  time_param,
  Npatho,
  area,
  rotation,
  croptypes,
  cultivars_param,
  eco_param,
  treatment_param,
  pathogen_param,
  audpc100S = 0.76,
  writeTXT = TRUE,
  graphic = TRUE,
  path = getwd()
)
}
\arguments{
\item{types}{a character string (or a vector of character strings if several outputs are to be computed)
specifying the type of outputs to generate (see details):\itemize{
\item "audpc": Area Under Disease Progress Curve
\item "audpc_rel": Relative Area Under Disease Progress Curve
\item "gla": Green Leaf Area
\item "gla_rel": Relative Green Leaf Area
\item "eco_yield": Total crop yield
\item "eco_cost": Operational crop costs
\item "eco_product": Crop products
\item "eco_margin": Margin (products - operational costs)
\item "contrib": contribution of pathogen genotypes to LIR dynamics
\item "HLIR_dynamics", "H_dynamics", "L_dynamics", "IR_dynamics", "HLI_dynamics", etc.: Epidemic dynamics
related to the specified sanitary status (H, L, I or R and all their combinations).
Graphics only, works only if graphic=TRUE.
\item "all": compute all these outputs (default).
}}

\item{time_param}{list of simulation parameters:\itemize{
\item Nyears = number cropping seasons,
\item nTSpY = number of time-steps per cropping season.
}}

\item{Npatho}{number of pathogen genotypes.}

\item{area}{a vector containing polygon areas (must be in square meters).}

\item{rotation}{a dataframe containing for each field (rows) and year (columns, named "year_1", "year_2", etc.),
the index of the cultivated croptype. Importantly, the matrix must contain 1 more column than the real number
of simulated years.}

\item{croptypes}{a dataframe with three columns named 'croptypeID' for croptype index,
'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype.}

\item{cultivars_param}{list of parameters associated with each host genotype (i.e. cultivars):
\itemize{
\item name = vector of cultivar names,
\item initial_density = vector of host densities (per square meter) at the beginning of the cropping season
as if cultivated in pure crop,
\item max_density = vector of maximum host densities (per square meter) at the end of the cropping season
as if cultivated in pure crop,
\item cultivars_genes_list = a list containing, for each host genotype, the indices of carried resistance genes.
}}

\item{eco_param}{a list of economic parameters for each host genotype as if cultivated in pure crop:\itemize{
\item yield_perHa = a dataframe of 4 columns for the theoretical yield associated with hosts in sanitary status H, L, I and R,
as if cultivated in pure crops, and one row per host genotype
(yields are expressed in weight or volume units / ha / cropping season),
\item planting_cost_perHa = a vector of planting costs (in monetary units / ha / cropping season),
\item market_value = a vector of market values of the production (in monetary units / weight or volume unit).
}}

\item{treatment_param}{list of parameters related to pesticide treatments: \itemize{
\item treatment_degradation_rate = degradation rate (per time step) of chemical concentration,
\item treatment_efficiency = maximal efficiency of chemical treatments (i.e. fractional reduction
of pathogen infection rate at the time of application),
\item treatment_timesteps = vector of time-steps corresponding to treatment application dates,
\item treatment_cultivars = vector of indices of the cultivars that receive treatments,
\item treatment_cost = cost of a single treatment application (monetary units/ha)
\item treatment_application_threshold = vector of thresholds (i.e. disease severity, one for each treated cultivar)
above which the treatment is applied in a polygon
}}

\item{pathogen_param}{a list of i. pathogen aggressiveness parameters on a susceptible host
for a pathogen genotype not adapted to resistance and ii. sexual reproduction parameters: \itemize{
\item infection_rate = maximal expected infection rate of a propagule on a healthy host,
\item propagule_prod_rate = maximal expected effective propagule production rate of an infectious host per time step,
\item latent_period_mean = minimal expected duration of the latent period,
\item latent_period_var = variance of the latent period duration,
\item infectious_period_mean = maximal expected duration of the infectious period,
\item infectious_period_var = variance of the infectious period duration,
\item survival_prob = probability for a propagule to survive the off-season,
\item repro_sex_prob = probability for an infectious host to reproduce via sex rather than via cloning,
\item sigmoid_kappa = kappa parameter of the sigmoid contamination function,
\item sigmoid_sigma = sigma parameter of the sigmoid contamination function,
\item sigmoid_plateau = plateau parameter of the sigmoid contamination function,
\item sex_propagule_viability_limit = maximum number of cropping seasons up to which a sexual propagule is viable
\item sex_propagule_release_mean = average number of seasons after which a sexual propagule is released,
\item clonal_propagule_gradual_release = whether or not clonal propagules surviving the bottleneck are gradually released along the following cropping season.
}}

\item{audpc100S}{the audpc in a fully susceptible landscape (used as reference value for graphics).}

\item{writeTXT}{a logical indicating if the output is written in a text file (TRUE) or not (FALSE).}

\item{graphic}{a logical indicating if a tiff graphic of the output is generated (only if more than one year is simulated).}

\item{path}{path of text file (if writeTXT = TRUE) and tiff graphic (if graphic = TRUE) to be generated.}
}
\value{
A list containing, for each required type of output, a matrix summarising the output for each year and cultivar
(as well as the whole landscape).
Each matrix can be written in a txt file (if writeTXT=TRUE), and illustrated in a graphic (if graphic=TRUE).
}
\description{
Generates epidemiological and economic outputs from model simulations.
}
\details{
Outputs are computed every year for every cultivar as well as for the whole landscape. \describe{
\item{\strong{Epidemiological outputs.}}{
The epidemiological impact of pathogen spread can be evaluated by different measures: \enumerate{
\item Area Under Disease Progress Curve (AUDPC): average number of diseased host individuals (status I + R)
per time step and square meter.
\item Relative Area Under Disease Progress Curve (AUDPCr): average proportion of diseased host individuals
(status I + R) relative to the total number of existing hosts (H+L+I+R).
\item Green Leaf Area (GLA): average number of healthy host individuals (status H) per time step and per square meter.
\item Relative Green Leaf Area (GLAr): average proportion of healthy host individuals (status H) relative to the total number
of existing hosts (H+L+I+R).
\item Contribution of pathogen genotypes: for every year and every host (as well as for the whole landscape and the whole
simulation duration), fraction of cumulative LIR infections attributed to each pathogen genotype.
}
}
\item{\strong{Economic outputs.}}{
The economic outcome of a simulation can be evaluated using: \enumerate{
\item Crop yield: yearly crop yield (e.g. grains, fruits, wine) in weight (or volume) units
per hectare (depends on the number of productive hosts and associated theoretical yield).
\item Crop products: yearly product generated from sales, in monetary units per hectare
(depends on crop yield and market value). Note that when disease = "mildew" a price reduction
between 0\% and 5\% is applied to the market value depending on disease severity.
\item Operational crop costs: yearly costs associated with crop planting (depends on initial
host density and planting cost) and pesticide treatments (depends on the number of applications and
the cost of a single application) in monetary units per hectare.
\item Crop margin, i.e. products - operational costs, in monetary units per hectare.
}
}
}
}
\examples{
\dontrun{
demo_landsepi()
}
}
\references{
Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018). Assessing the durability and efficiency of
landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
}
\seealso{
\link{evol_output}
}
