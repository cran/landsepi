# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@csiro.au>
#                    Julien Papaix <julien.papaix@csiro.au>
#                    Jean-François Rey <jean-francois.rey@inra.fr>
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
#' @description A spatio-temporal stochastic model to assess resistance deployment strategies against plant pathogens.
#' The model is based on stochastic geometry for describing the landscape and the resistant hosts,
#' a dispersal kernel for the dissemination of the pathogen, and a SEIR (Susceptible-Exposed-Infectious-Removed) 
#' architecture to simulate plant response to disease.
#' @aliases landsepi-package landsepi
#' 
#' @author Loup Rimbaud \email{loup.rimbaud@@csiro.au}
#' @author Julien Papaix \email{julien.papaix@@inra.fr}
#' @author Jean-Francois Rey \email{jean-francois.rey@@inra.fr}
#' 
#' Maintainer: Jean-Francois Rey \email{jean-francois.rey@@inra.fr}
#' @docType package
#' @name landsepi-package
#' @details \tabular{ll}{
#'          Package: \tab lansepi\cr
#'          Type: \tab Package\cr
#'          Version: \tab 0.0.2\cr
#'          Date: \tab 2018-03-19\cr
#'          License: \tab GPL (>=2)\cr
#'          }
#'
#' The landsepi package implements a spatially explicit stochastic model able to assess the epidemiological and evolutionary outcomes of four major strategies
#' to deploy plant resistance to pathogens. These strategies include the combination of several resistance sources across time (crop rotations) or space.
#' The spatial scale of deployment can vary from multiple resistance sources occurring in a single cultivar (pyramiding), 
#' in different cultivars within the same field (cultivar mixtures) or in different fields (mosaics). The simulated sources of resistance can 
#' consist of qualitative resistance (i.e. major genes) or quantitative resistance traits against several components of pathogen aggressiveness: 
#' infection rate, latent period duration, propagule production rate, and infectious period duration. This model provides a useful tool to assess 
#' the performance of a wide range of deployment options, and helps investigate the effect of landscape, epidemiological and evolutionary parameters 
#' on the performance of a given strategy. 
#' 
#' The simulation model is based on a SEIR (Susceptible-Exposed-Infectious-Removed) architecture to describe host response to disease. The lansdcape is 
#' represented by a set of polygons where the pathogen can disperse. Initially, the pathogen is not adapted to any source of resistance, and is only 
#' present on susceptible hosts. However, through mutation, it can evolve and may acquire infectivity genes (which leads to breakdown of major resistance genes) 
#' or increase aggressiveness (which leads to the erosion of the relevant quantitative resistance traits). 
#' However, evolution of a pathogen toward infectivity or increased aggressiveness on a resistant host may be penalised 
#' by a fitness cost on susceptible hosts. Consequently, pathogens carrying infectivity genes may have reduced infection rate (cost of infectivity) on susceptible 
#' hosts relative to pathogens that do not carry these genes. Similarly, a gain in pathogen aggressiveness on quantitatively resistant hosts is penalised by a 
#' decreased aggressiveness on susceptible hosts, leading to a trade-off.
#' 
#' The evolutionary outcome of a deployment strategy is assessed by measuring the time until the pathogen reaches the three steps to adapt to plant resistance: \itemize{
#' \item (d1) first appearance of adapted mutants, 
#' \item (d2) initial migration to resistant hosts and infection, and 
#' \item (d3) broader establishment in the resistant host population (i.e. the point at which extinction becomes unlikely). }
#' Epidemiological outcomes are evaluated using: \itemize{
#' \item (e1) the Green Leaf Area (GLA) as a proxy for yield, and 
#' \item (e2) the area under the disease progress curve (AUDPC) to measure disease severity.}
#' 
#' The package includes five examples of landscape structures. The demonstration function is parameterise to roughly represent biotrophic foliar fungi of cereal crops, as
#'  typified by rusts of wheat (genus \emph{Puccinia}).

#' @keywords model spatial demo-genetic deployment resistance durability stochastic SEIR
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (in press). Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology}.
#' @examples \dontrun{
#' library("landsepi")
#' ## Run a demonstration
#' demo_landsepi() 
#' ## Run a simulation
#' simul_landsepi()
#' }
#' @useDynLib landsepi
#' @import methods
#' @import graphics
#' @import stats
#' @import Rcpp
#' @import rgdal
#' @import rgeos
#' @importFrom Matrix nearPD
#' @import MASS
#' @import maptools
#' @import fields
#' @import sp
#' @import splancs
#' @importFrom utils data
"_PACKAGE"
