# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
#                    Julien Papaix <julien.papaix@inrae.fr>
#                    Jean-François Rey <jean-francois.rey@inrae.fr>
#                    Jean-Loup Gaussen <jean-loup-thomas.gaussen@inrae.fr>
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
#

#' @encoding UTF-8
#' @title Landscape Epidemiology and Evolution
#' @description A stochastic, spatially-explicit, demo-genetic model simulating 
#' the spread and evolution of a plant pathogen in a heterogeneous landscape to assess 
#' resistance deployment strategies.
#' @aliases landsepi-package landsepi
#'
#' @author Loup Rimbaud \email{loup.rimbaud@@inrae.fr}
#' @author Marta Zaffaroni \email{marta.zaffaroni@@inrae.fr}
#' @author Jean-Francois Rey \email{jean-francois.rey@@inrae.fr}
#' @author Julien Papaix \email{julien.papaix@@inrae.fr}
#' @author Jean-Loup Gaussen \email{jean-loup-thomas.gaussen@@inrae.fr}
#' @author Manon Couty \email{manon.couty@@insa-lyon.fr}
#'
#' Maintainer: Jean-Francois Rey \email{jean-francois.rey@@inrae.fr}
#' @docType package
#' @name landsepi-package
#' @details \tabular{ll}{
#'          Package: \tab landsepi\cr
#'          Type: \tab Package\cr
#'          Version: \tab 1.5.2\cr
#'          Date: \tab 2025-07-30\cr
#'          License: \tab GPL (>=2)\cr
#'          }
#'
#' The landsepi package implements a spatially explicit stochastic model able to assess the epidemiological,
#' evolutionary and economic outcomes of strategies to deploy plant resistance to pathogens. 
#' It also helps investigate the effect of landscape organisation, the considered pathosystem and
#' the epidemio-evolutionary context on the performance of a given strategy.
#' 
#' It is based on a spatial geometry for describing the landscape and allocation of different cultivars,
#' a dispersal kernel for the dissemination of the pathogen,
#' and a SEIR (‘susceptible-exposed-infectious-removed’, renamed HLIR for 'healthy-latent-infectious-removed'
#' to avoid confusions with 'susceptible host') structure with a discrete time step. It simulates the spread and
#' evolution (via mutation, recombination through sexual reproduction, selection and drift) 
#' of a pathogen in a heterogeneous cropping landscape, across cropping seasons split by host harvests which impose
#' potential bottlenecks to the pathogen.
#'
#' The lansdcape is represented by a set of polygons where the pathogen can disperse
#' (the basic spatial unit is an individual polygon; an agricultural field may be composed of a single 
#' or several polygons). *landsepi* includes built-in simulated landscapes 
#' (and associated dispersal matrices for rust pathogens, see below), but is it possible 
#' to use your own landscape (in shapefile format) and dispersal matrix.  
#' 
#' A wide array of resistance deployment strategies can be simulated in landsepi: fields of the 
#' landscape are cultivated with different croptypes that can rotate through time; each croptype is
#' composed of either a pure cultivar or a mixture; and each cultivar may carry one or several resistance 
#' genes. Thus, all combinations of rotations, mosaics, mixtures and pyramiding strategies are 
#' possible. Resistance genes affect several possible pathogen aggressiveness components: 
#' infection rate, durations of the latent period and the infectious period, and propagule
#' production rate. Resistance may be complete (i.e. complete inhibition of the targeted aggressiveness component) or partial
#' (i.e. the targeted aggressiveness component is only softened), and expressed from the beginning of the season, or later
#' (to simulate Adult Plant Resistance (APR), also called Mature Plant Resistance). 
#' Cultivar allocation can be realised via an algorithm (\code{allocateCroptypeCultivars()}) 
#' but it is possible to use your own cultivar allocation if it is included in the shapefile 
#' containing the landsape. 
#' Additionally, any cultivar may be treated with contact pesticides, which reduce the pathogen infection rate 
#' with an efficiency gradually decreasing with time and host growth.
#'
#' To each resistance gene in the host (whether it may be a major gene or a QTL for quantitative resistance) 
#' is associated a pathogenicity gene in the pathogen.
#' Through mutation of pathogenicity genes, the pathogen can restore its aggressiveness on resistance hosts and thus
#' adapt to resistance (leading to sudden breakdown or gradual erosion of resistance genes).
#' Pathogenicity genes may also be reassorted via sexual reproduction or gene recombination.
#' Increased aggressiveness on a resistant host (i.e. adaptation to the corresponding resistance genes)
#' can be penalised by a fitness cost, either on all hosts, or only on susceptible hosts (in the latter case, 
#' pathogen genotypes adapted to a resistance gene have a reduced aggressiveness on hosts that do not carry this gene, 
#' and a 'relative advantage' on host that do carry such gene). 
#' The relation between pathogen aggressiveness on susceptible and resistant hosts
#' is defined by a trade-off relationship whose shape depends on the strength of the trade-off. 
#' Strong trade-off means that the gain in fitness on resistant hosts is smaller than the cost on susceptible hosts.
#' 
#' The package includes five examples of landscape structures and a default parameterisation to represent
#' plant pathogens as typified by rusts of cereal crops (genus \emph{Puccinia},
#' e.g. stripe rust, stem rust and leaf rust of wheat and barley). A parameterisation to
#' downy mildew of grapevine (\emph{Plasmopara viticola}) and black sigatoka of banana 
#' (\emph{Pseudocercospora fijiensis}) are also available.
#' The main function of the package is \code{runSimul()}.
#' It can be parameterised to simulate various resistance deployment strategies using either the provided 
#' landscapes and parameters for cereal rusts, or landscapes and parameters set by the user. 
#' See \code{demo_landsepi()} for a demonstration, and our tutorials (\code{browseVignettes("landsepi")}) 
#' for details on how to use landsepi.
#'
#' \describe{
#' \item{\strong{Assumptions} (in bold those that can be relaxed with appropriate parameterization): }{
#' \enumerate{
#' 
#' 
#'  \item The spatial unit is a polygon, i.e. a piece of land delimited by boundaries and possibly 
#'  cultivated with a crop. Such crop may be host or non-host, and the polygon is considered a homogeneous 
#'  mixture of host individuals (i.e. there is no intra-polygon structuration). 
#'  An agricultural field may be composed of a single or several polygons.
#'  \item A host ‘individual’ is an infection unit (i.e. it can be infected by one and only one 
#'  pathogen propagule, there is no co-infection) and may correspond to **a given amount of plant tissue 
#'  (where a local infection may develop, e.g. fungal lesion) or a whole plant (e.g. systemic viral infection). 
#'  In the first case, plant growth increases the amount of available plant tissue 
#'  (hence the number of individuals) during the cropping season.** Plant growth is deterministic (logistic growth) 
#'  and **only healthy individuals (state H) contribute to plant growth (castrating pathogen)**.
#'  \item Host individuals are in one of these four categories: H (healthy), E (exposed and latent, 
#'  i.e. infected but not infectious nor symptomatic), I (infectious and symptomatic), or R 
#'  (removed, i.e. epidemiologically inactive).
#'  \item **The decreasing availability of healthy host tissues (as epidemics spread) makes pathogen 
#'  infection less likely (i.e. density-dependence due to plant architecture).**
#'  \item **Hosts are cultivated (i.e. sown/planted and harvested), thus there is no host reproduction, 
#'  dispersal and natural death.**
#'  \item Environmental and climate conditions are constant, and host individuals of a given genotype 
#'  are equally susceptible to disease from the first to the last day of every cropping season.
#'  \item Crop yield depends on the average amount of producing host individuals during the cropping 
#'  season and does not depend on the time of epidemic peak. **Only healthy individuals (state H) 
#'  contribute to crop yield.**
#'  \item Cultivars may be treated with chemicals which reduce the pathogen infection rate (contact treatment). 
#'  Treatment efficiency decreases with host growth (i.e. new biomass is not protected by treatments) 
#'  **and time (i.e. pesticide degradation)**. Cultivars to be treated and dates of chemical applications 
#'  are fixed prior to simulations but only polygons where disease severity exceeds a given threshold (possibly 0) are treated. 
#'  \item Components of a mixture are independent each other (i.e. there is neither plant-plant 
#'  interaction nor competition for space, and harvests are segregated). If one component is treated 
#'  with a chemical, it does not affect other components. 
#'  \item The pathogen is haploid.
#'  \item **Initially, the pathogen is not adapted to any source of resistance, and is only present on 
#'  susceptible hosts (at state I).**
#'  \item **Pathogen dispersal is isotropic (i.e. equally probable in every direction).**
#'  \item **Boundaries of the landscape are reflective: propagules stay in the system as if it was closed.** 
#'  \item Pathogen reproduction can be purely clonal, purely sexual, or mixed (alternation of clonal 
#'  and sexual reproduction).
#'  \item If there is sexual reproduction (or gene recombination), it occurs only between parental 
#'  infections located in the same polygon and the same host genotype (i.e. cultivar). 
#'  At that scale, the pathogen population is panmictic (i.e. all pairs of parents have the 
#'  same probability to occur). 
#'  The propagule production rate of a parental pair is the sum of the propagule production rates of 
#'  the parents. For a given parental pair, the genotype of each propagule is issued from random loci 
#'  segregation of parental qualitative resistance genes. For each quantitative resistance gene, the 
#'  value of each propagule trait is issued from a normal distribution around the average of the 
#'  parental traits, following the infinitesimal model (Fisher 1919).
#'  \item All types of propagules (i.e. clonal and sexual) share the same pathogenicity parameters 
#'  (e.g. infection rate, latent period duration, etc.) but each of them has their own dispersal and survival 
#'  abilities (see after). 
#'  \item At the end of each cropping season, pathogens experience a bottleneck representing the 
#'  off-season and then propagules are produced (either via clonal or sexual reproduction). 
#'  **The probability of survival is the same every year and in every polygon.** 
#'  Clonal propagules are released during the following season only, either altogether at the first day of 
#'  the season, or progressively (in that case the day of release of each propagule is sampled from 
#'  a uniform distribution). Sexual propagules are gradually released during several of the following 
#'  seasons (between-season release). The season of release of each propagule is sampled from an 
#'  exponential distribution, truncated by a maximum viability limit. Then, the day of release in a 
#'  given season is sampled from a uniform distribution (within-season release).
#'  \item Pathogenicity genes mutate independently from each other.
#'  \item **Pathogen adaptation to a given resistance gene consists in restoring the same aggressiveness 
#'  component as the one targeted by the resistance gene.**
#'  \item If a fitness cost penalises pathogen adaptation to a given resistance gene, this cost is paid 
#'  on all hosts with possibly a relative advantage on hosts carrying the resistance gene. 
#'  It consists in a reduction in the same aggressiveness 
#'  component as the one targeted by the resistance gene.
#'  \item When there is a delay for activation of a given resistance gene (APR), the age of activation 
#'  is the same for all hosts carrying this gene and located in the same polygon.
#'  \item Variances of the durations of the latent and the infectious periods of the pathogen are 
#'  not affected by plant resistance. 
#'   }
#'  }

#' \item{\strong{Epidemiological outputs}}{
#' The epidemiological outcome of a deployment strategy is evaluated using: \enumerate{
#' \item the area under the disease progress curve (AUDPC) to measure disease severity
#' (i.e. the average number of diseased plant tissue -status I and R- per time step and square meter),
#' \item the relative area under the disease progress curve (AUDPCr) to measure the average proportion 
#' of diseased tissue (status I and R) relative to the total number of existing host individuals (H+L+I+R).
#' \item the Green Leaf Area (GLA) to measure the average amount of healthy plant tissue (status H) per time step and square meter,
#' \item the relative Green Leaf Area (GLAr) to measure the average proportion of healthy tissue (status H)
#' relative to the total number of existing host individuals (H+L+I+R).
#' \item the yearly contribution of pathogen genotypes to LIR dynamics on every host as well as the whole landscape.
#'   }
#' A set of graphics and a video showing epidemic dynamics can also be generated.
#'  }

#' \item{\strong{Evolutionary outputs}}{
#' The evolutionary outcome is assessed by measuring: \enumerate{
#' \item the dynamics of pathogen genotype frequencies,
#' \item the evolution of pathogen aggressiveness,
#' \item the durability of resistance genes. Durability can be estimated using the time until the pathogen reaches the three
#' steps to adapt to plant resistance: (1) first appearance of adapted mutants,
#' (2) initial migration to resistant hosts and infection, and
#' (3) broader establishment in the resistant host population (i.e. the point at which extinction becomes unlikely).
#'   }
#'  }

#' \item{\strong{Economic outputs}}{
#' The economic outcome of a simulation can be evaluated using: \enumerate{
#' \item the crop yield: yearly crop production (e.g. grains, fruits, wine) in weight (or volume) units
#' per hectare (depends on the number of productive hosts and associated theoretical yield),
#' \item the crop products: yearly products generated from sales, in monetary units per hectare
#' (depends on crop yield and market value),
#' \item the crop operational costs: yearly costs associated with crop planting (depends on initial 
#' host density and planting cost) and pesticide treatments (depends on the number of applications and 
#' the cost of a single application) in monetary units per hectare.
#' \item the margin, i.e. products - operational costs, in monetary units per hectare.
#'   }
#'  }
#' }
#'

#' \strong{Future versions:}
#'
#' Future versions of the package will include in particular:\itemize{
#' \item Sets of pathogen parameters to simulate other pathosystems (e.g. Cucumber mosaic virus on pepper, potato virus Y on pepper).
#' \item An updated version of the shiny interface.
#' }
#' \strong{Dependencies:}
#'
#' The package for compiling needs:\itemize{
#' \item g++
#' \item libgsl2
#' \item libgsl-dev}
#' and the following R packages:\itemize{
#' \item Rcpp
#' \item sp
#' \item stats
#' \item Matrix
#' \item mvtnorm
#' \item fields
#' \item splancs
#' \item sf
#' \item DBI
#' \item RSQLite
#' \item foreach
#' \item parallel
#' \item doParallel
#' \item deSolve}
#' In addition, to generate videos the package will need ffmpeg.
#' @keywords model spatial demo-genetic deployment resistance durability stochastic SEIR
#' @references
#' ## When referencing the simulation model, please cite the following article:
#'
#' Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018). Assessing the durability and efficiency of
#' landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#'
#' ## When referencing the R package, please cite the following package:
#'
#' Rimbaud L., Papaïx J. and Rey J.-F. (2018). landsepi: Landscape Epidemiology and Evolution. \emph{R package},
#' url: https://cran.r-project.org/package=landsepi.
#' @examples
#' \dontrun{
#' library("landsepi")
#'
#' ## Run demonstrations (in 10-year simulations) for different deployment strategies:
#' demo_landsepi(strat = "MO") ## for a mosaic of cultivars
#' demo_landsepi(strat = "MI") ## for a mixture of cultivars
#' demo_landsepi(strat = "RO") ## for a rotation of cultivars
#' demo_landsepi(strat = "PY") ## for a pyramid of resistance genes
#' }
#' @useDynLib landsepi, .registration = TRUE
#' @import methods
#' @import graphics
#' @import stats
#' @import Rcpp
#' @importFrom Matrix nearPD
#' @importFrom mvtnorm rmvnorm
#' @import fields
#' @import sp
#' @import sf
#' @importFrom utils data
#' @import RSQLite
#' @import DBI
"_PACKAGE"

# @title Print package information
# @name getInfo
# @description Displays some information about the package
# @importFrom utils packageVersion
getInfo <- function() {
  packageStartupMessage("Package: landsepi | Landscape Epidemiology and Evolution")
  packageStartupMessage("Version: ", appendLF = FALSE)
  packageStartupMessage(utils::packageVersion("landsepi"))
  packageStartupMessage("License: GPL (>= 2)")
}

# @title Things to do at package attach
# @name .onAttach
# @param libname a character string giving the library directory where
#  the package defining the namespace was found.
# @param pkgname a character string giving the name of the package.
# @description Print package information and check dependencies
.onAttach <- function(libname, pkgname) {
  getInfo()
}
# .onLoad <- function(libname, pkgname) {
#   getInfo()
# }
