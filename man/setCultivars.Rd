% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods-LandsepiParams.R
\name{setCultivars}
\alias{setCultivars}
\title{Set cultivars}
\usage{
setCultivars(params, dfCultivars)
}
\arguments{
\item{params}{a landsepiParams object.}

\item{dfCultivars}{a data.frame defining the cultivars (see details). It can be generated
manually or, alternatively, via \code{\link{loadCultivar}}.}
}
\value{
a LandsepiParams object
}
\description{
Updates a LandsepiParams object with cultivars parameters
}
\details{
dfCultivars is a dataframe of parameters associated with each host genotype
(i.e. cultivars, lines) when cultivated in pure crops. Columns of the dataframe are:\itemize{
\item cultivarName: cultivar names (cannot accept space),
\item initial_density: host densities (per square meter) at the beginning of the cropping season
as if cultivated in pure crop,
\item max_density: maximum host densities (per square meter) at the end of the cropping season
as if cultivated in pure crop,
\item growth rate: host growth rates,
\item reproduction rate: host reproduction rates,
\item yield_H: theoretical yield (in weight or volume units / ha / cropping season)
associated with hosts in sanitary status H as if cultivated in pure crop,
\item yield_L: theoretical yield (in weight or volume units / ha / cropping season)
associated with hosts in sanitary status L as if cultivated in pure crop,
\item yield_I: theoretical yield (in weight or volume units / ha / cropping season)
associated with hosts in sanitary status I as if cultivated in pure crop,
\item yield_R: theoretical yield (in weight or volume units / ha / cropping season)
associated with hosts in sanitary status R as if cultivated in pure crop,
\item planting_cost = planting costs (in monetary units / ha / cropping season) as if
cultivated in pure crop,
\item market_value = market values of the production (in monetary units / weight or volume unit).
}

The data.frame must be defined as follow (example):\tabular{lllllllllll}{
   cultivarName \tab initial_density \tab max_density \tab growth_rate \tab reproduction_rate \tab yield_H \tab yield_L \tab yield_I \tab yield_R \tab planting_cost \tab market_value \cr
   Susceptible \tab 0.1 \tab 2.0 \tab 0.1 \tab 0.0 \tab 2.5 \tab 0.0 \tab 0.0 \tab 0.0 \tab 225 \tab 200 \cr
   Resistant1 \tab 0.1 \tab 2.0 \tab 0.1 \tab 0.0 \tab 2.5 \tab 0.0 \tab 0.0 \tab 0.0 \tab 225 \tab 200 \cr
   Resistant2 \tab 0.1 \tab 2.0 \tab 0.1 \tab 0.0 \tab 2.5 \tab 0.0 \tab 0.0 \tab 0.0 \tab 225 \tab 200 \cr
}
}
\examples{
\dontrun{
simul_params <- createSimulParams()
cultivar1 <- loadCultivar(name = "Susceptible", type = "wheat")
cultivar2 <- loadCultivar(name = "Resistant", type = "wheat")
cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
simul_params <- setCultivars(simul_params, cultivars)
simul_params@Cultivars
}
}
\seealso{
\link{loadCultivar}
}
