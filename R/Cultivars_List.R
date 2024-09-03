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

#' Cultivars Type list
#'
#' A set of configurated cultivars types
#'
#' @format A list of list indexed by type name
#' \itemize{
#' \item cultivarName: cultivar names (cannot accept space),
#' \item initial_density: host densities (per square meter) at the beginning of the cropping season 
#' as if cultivated in pure crop,
#' \item max_density: maximum host densities (per square meter) at the end of the cropping season 
#' as if cultivated in pure crop,
#' \item growth rate: host growth rates,
#' \item reproduction rate: host reproduction rates,
#' \item yield_H: theoretical yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status H as if cultivated in pure crop,
#' \item yield_L: theoretical yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status L as if cultivated in pure crop,
#' \item yield_I: theoretical yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status I as if cultivated in pure crop,
#' \item yield_R: theoretical yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status R as if cultivated in pure crop,
#' \item planting_cost = planting costs (in monetary units / ha / cropping season) as if cultivated in pure crop,
#' \item market_value = market values of the production (in monetary units / weight or volume unit).
#' }
# "Cultivars_list"
Cultivars_list <- list(
  # "growingHost" = list(
  #   "cultivarName" = "",
  #   "initial_density" = 0.1,
  #   "max_density" = 2.0,
  #   "growth_rate" = 0.1,
  #   "reproduction_rate" = 0.0,
  #   "yield_H" = 2.5,
  #   "yield_L" = 0.0,
  #   "yield_I" = 0.0,
  #   "yield_R" = 0.0,
  #   "planting_cost" = 225,
  #   "market_value" = 200
  # ),
  "wheat" = list(
    "cultivarName" = "",
    "initial_density" = 0.1,
    "max_density" = 2.0,
    "growth_rate" = 0.1,
    "reproduction_rate" = 0.0,
    "yield_H" = 2.5,
    "yield_L" = 0.0,
    "yield_I" = 0.0,
    "yield_R" = 0.0,
    "planting_cost" = 225,
    "market_value" = 200
  ),
  "grapevine" = list(
    "cultivarName" = "",
    "initial_density" = 1.0,
    "max_density" = 20.0,
    "growth_rate" = 0.1,
    "reproduction_rate" = 0.0,
    "yield_H" = 6.7,
    "yield_L" = 6.7,
    "yield_I" = 6.7,
    "yield_R" = 0.0,
    "planting_cost" = 5481,
    "market_value" = 600
  ),
  "banana" = list(
    "cultivarName" = "Cavendish",
    "initial_density" = 0.9,
    "max_density" = 1.8,
    "growth_rate" = 0.02,
    "reproduction_rate" = 0.0,
    "yield_H" = 46.8,
    "yield_L" = 46.8,
    "yield_I" = 0.0,
    "yield_R" = 0.0,
    "planting_cost" = 0,
    "market_value" = 0
  ),
  "pepper" = list(
    "cultivarName" = "Gorria",
    "initial_density" = 1.75,
    "max_density" = 1.75,
    "growth_rate" = 0.0,
    "reproduction_rate" = 0.0,
    "yield_H" = 0.5,
    "yield_L" = 0.5,
    "yield_I" = 0.0,
    "yield_R" = 0.0,
    "planting_cost" = 0,
    "market_value" = 0
  ),
  "nonCrop" = list(
    "cultivarName" = "",
    "initial_density" = 0.0,
    "max_density" = 2.0,
    "growth_rate" = 0.0,
    "reproduction_rate" = 0.0,
    "yield_H" = 0.0,
    "yield_L" = 0.0,
    "yield_I" = 0.0,
    "yield_R" = 0.0,
    "planting_cost" = 0,
    "market_value" = 0
  )
)
