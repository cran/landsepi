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
#' @author Julien Papaix \email{julien.papaix@@inrae.fr}
#' @author Jean-Francois Rey \email{jean-francois.rey@@inrae.fr}
#' @author Jean-Loup Gaussen \email{jean-loup-thomas.gaussen@@inrae.fr}
#'
#' Maintainer: Jean-Francois Rey \email{jean-francois.rey@@inrae.fr}
#' @docType package
#' @name landsepi-package
#' @details \tabular{ll}{
#'          Package: \tab landsepi\cr
#'          Type: \tab Package\cr
#'          Version: \tab 1.0.2\cr
#'          Date: \tab 2020-11-26\cr
#'          License: \tab GPL (>=2)\cr
#'          }
#'
#' The landsepi package implements a spatially explicit stochastic model able to assess the epidemiological,
#' evolutionary and economic outcomes of strategies to deploy plant resistance to pathogens.
#' It is based on a spatial geometry for describing the landscape and allocation of different cultivars,
#' a dispersal kernel for the dissemination of the pathogen,
#' and a SEIR (‘susceptible-exposed-infectious-removed’, renamed HLIR for 'healthy-latent-infectious-removed'
#' to avoid confusions with 'susceptible host') structure with a discrete time step. It simulates the spread and
#' evolution of a pathogen in a heterogeneous cropping landscape, across cropping seasons split by host harvests which impose
#' potential bottlenecks to the pathogen.
#'
#' The lansdcape is represented by a set of polygons where the pathogen can disperse
#' (the basic spatial unit is an individual field). *landsepi* includes built-in simulated landscapes 
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
#'
#' To each resistance gene in the host (whether it may be a major gene or a QTL for quantitative resistance) 
#' is associated a pathogenicity gene in the pathogen.
#' Through mutation of pathogenicity genes, the pathogen can restore its aggressiveness on resistance hosts and thus
#' adapt to resistance (leading to sudden breakdown or gradual erosion of resistance genes).
#' Pathogenicity genes may also be reassorted via sexual reproduction or gene recombination.
#' Increased in aggressiveness on a resistant host (i.e. adaptation to the corresponding resistance genes)
#' can be penalised by a fitness cost on susceptible hosts, i.e. pathogen genotypes adapted to a resistance gene have
#' a reduced aggressiveness on hosts that do not carry this gene. 
#' The relation between pathogen aggressiveness on susceptible and resistant hosts
#' is defined by a trade-off relationship whose shape depends on the strength of the trade-off. 
#' Strong trade-off means that the gain in fitness on resistant hosts is smaller than the cost on susceptible hosts.
#'
#' This model provides a useful tool to assess the performance of a wide range of deployment options
#' via epidemiological, evolutionary and economic outputs.
#' It also helps investigate the effect of landscape organisation, the considered pathosystem and
#' the epidemio-evolutionary context on the performance of a given strategy.
#'
#' The package includes five examples of landscape structures and a default parameterisation to represent
#' plant pathogens as typified by rusts of cereal crops (genus \emph{Puccinia},
#' e.g. stripe rust, stem rust and leaf rust of wheat and barley).
#' The main function of the package is \code{runSimul()}.
#' It can be parameterised to simulate various resistance deployment strategies using either the provided 
#' landscapes and parameters for cereal rusts, or landscapes and parameters set by the user. 
#' See \code{demo_landsepi()} for a demonstration, and our tutorials (\code{browseVignettes("landsepi")}) 
#' for details on how to use landsepi.
#'
#' \describe{
#' \item{\strong{Assumptions} (in bold those that can be relaxed with appropriate parameterization): }{
#' \enumerate{
#'  \item The spatial unit is the field, i.e. a piece of land delimited by boundaries and possibly cultivated with a crop.
#'  Such crop may be host or non-host, and the field is considered a homogeneous mixture of individuals (i.e. there is no
#'  intra-field structuration).
#'  \item Host individuals are in one of these four categories: H (healthy), 
#'  E (latent, i.e. infected but not infectious nor symptomatic), I (infectious and symptomatic), 
#'  or R (removed, i.e. epidemiologically inactive). 
#'  \item **A host `individual' is an infection unit and may correspond to a given amount of plant tissue
#'  (where a local infection may develop, e.g. fungal lesion) or a whole plant (e.g. systemic viral infection). 
#'  In the first case, plant growth increases the amount of available plant tissue (hence the number of individuals) 
#'  during the cropping season.** Plant growth is deterministic (logistic growth) and 
#'  only healthy hosts (state H) contribute to plant growth (castrating pathogen).
#'  \item **The decreasing availability of healthy host tissues (as epidemics spread) makes pathogen infection less likely
#'  (i.e. density-dependence due to plant architecture).**
#'  \item **Host are cultivated, thus there is no host reproduction, dispersal and natural death.**
#'  \item Environmental and climate conditions are constant, and host individuals of a given genotype are equally
#'  susceptible to disease from the first to the last day of every cropping season.
#'  \item Crop yield depends on the average amount of producing host individuals during the cropping season 
#'  and does not depend on the time of epidemic peak. **Only healthy individuals (state H) contribute to crop yield.**
#'  \item Initially, the pathogen is not adapted to any source of resistance, and is only present on
#'  susceptible hosts (at state I).
#'  \item **Pathogen dispersal is isotropic (i.e. equally probable in every direction).**
#'  \item **Pathogen reproduction is clonal.**
#'  \item Pathogenicity genes mutate independently from each other.
#'  \item **Pathogen adaptation to a given resistance gene consists in restoring the same aggressiveness component
#'  as the one targeted by the resistance gene.**
#'  \item If a fitness cost penalises pathogen adaptation to a given resistance gene, this cost is paid on 
#'  hosts that do not carry this gene, and consists in a reduction in the same aggressiveness component as 
#'  the one targeted by the resistance gene.
#'  \item When there is a delay for activation of a given resistance gene (APR), the time to activation is the same for
#'  all hosts carrying this gene and located in the same field.
#'  \item Variances of the durations of the latent and the infectious periods of the pathogen are not affected by plant resistance.
#'  \item If there is sexual reproduction (or gene recombination), it occurs only between parental infections located in the same field
#'  and the same host genotype. The propagule production rate of a couple is the sum of the propagule production rates of the parents.
#'  The genotype of each daughter propagule is issued from random loci segregation between parental loci.
#'   }
#'  }

#' \item{\strong{Epidemiological outputs}}{
#' The epidemiological outcome of a deployment strategy is evaluated using: \enumerate{
#' \item the area under the disease progress curve (AUDPC) to measure disease severity
#' (i.e. the average proportion of diseased hosts -status I and R- relative to the carrying capacity),
#' \item the absolute Green Leaf Area (GLAa) to measure the average amount of healthy tissue (status H),
#' \item the relative Green Leaf Area (GLAr) to measure the average proportion of healthy tissue (status H)
#' relative to the total number of existing hosts (H+L+I+R).
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
#' \item the crop production: yearly crop production (e.g. grains, fruits, wine) in weight (or volume) units
#' per hectare (depends on the number of productive hosts and associated yield),
#' \item the crop benefits: yearly benefits generated from product sales, in monetary units per hectare
#' (depends on crop production and market value of the product),
#' \item the crop costs: yearly costs associated with crop production (including planting, amortisation, labour, ...)
#'  in monetary units per hectare (depends on initial host density and production cost),
#' \item the gross margin, i.e. benefits - costs, in monetary units per hectare.
#'   }
#'  }
#' }
#'

#' \strong{Future versions:}
#'
#' Future versions of the package will include in particular:\itemize{
#' \item Sets of pathogen parameters to simulate other pathosystems (e.g. canola blackleg, grapevine downy mildew, potato virus Y on pepper).
#' \item More flexible initial conditions (e.g. size, location and composition of pathogen inoculum at the beginning of the simulation).
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
#' \item maptools
#' \item fields
#' \item splancs
#' \item sf
#' \item DBI
#' \item RSQLite
#' \item foreach
#' \item parallel
#' \item doParallel}
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
#' ## Run demonstrations (in 20-year simulations) for different deployment strategies:
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
