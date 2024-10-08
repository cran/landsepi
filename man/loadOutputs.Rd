% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods-LandsepiParams.R
\name{loadOutputs}
\alias{loadOutputs}
\title{Load outputs}
\usage{
loadOutputs(epid_outputs = "all", evol_outputs = "all", disease = "rust")
}
\arguments{
\item{epid_outputs}{a character string (or a vector of character strings if several outputs
are to be computed) specifying the type of epidemiological and economic outputs to generate
(see details):\itemize{
\item "audpc" : Area Under Disease Progress Curve (average number of diseased host individuals
per time step and square meter)
\item "audpc_rel" : Relative Area Under Disease Progress Curve (average proportion of
diseased host individuals relative to the total number of existing hosts)
\item "gla" : Green Leaf Area (average number of healthy host individuals per time step and
square meter)
\item "gla_rel" : Relative Green Leaf Area (average proportion of healthy host individuals
relative to the total number of existing hosts)
\item "eco_yield" : total crop yield (in weight or volume units per ha)
\item "eco_cost" : operational crop costs (in monetary units per ha)
\item "eco_product" : total crop products (in monetary units per ha)
\item "eco_margin" : Margin (products - operational costs, in monetary units per ha)
\item "contrib": contribution of pathogen genotypes to LIR dynamics
\item "HLIR_dynamics", "H_dynamics", "L_dynamics", "IR_dynamics", "HLI_dynamics", etc.:
Epidemic dynamics related to the specified sanitary status (H, L, I or R and all their
combinations). Graphics only, works only if graphic=TRUE.
\item "all" : compute all these outputs (default)
\item "" : none of these outputs will be generated.
}}

\item{evol_outputs}{a character string (or a vector of character strings if several
outputs are to be computed) specifying the type of evolutionary outputs to generate :\itemize{
\item "evol_patho": Dynamics of pathogen genotype frequencies
\item "evol_aggr": Evolution of pathogen aggressiveness
\item "durability": Durability of resistance genes
\item "all": compute all these outputs (default)
\item "": none of these outputs will be generated.
}}

\item{disease}{a disease name, among "rust" (default), "mildew", "sigatoka" and "no pathogen"}
}
\value{
a list of outputs and parameters for output generation
}
\description{
Creates an output list
}
\examples{
outputList <- loadOutputs(epid_outputs = "audpc", evol_outputs = "durability")
outputList
}
\seealso{
\link{setOutputs}, \link{compute_audpc100S}
}
