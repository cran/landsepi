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


#' @title Periodic covariance function
#' @name periodic_cov
#' @description Periodic function used to compute the variance-covariance matrix of the fields of the landscape. 
#' @param d a numeric object containing pairwise distances between the centroids of the fields
#' @param range range (half-period of oscillations)
#' @param phi amplitude of the oscillations
#' @details The periodic covariance is defined by exp(-2 * sin(d*pi/(2*range))^2 / phi^2). It is used to generate highly fragmented or highly aggregated landscapes.
#' @return An object of the same type as d. 
#' @examples
#' periodic_cov(10, range=5)
#' @export
periodic_cov <- function(d,range,phi=1) {  return( exp(-2 * sin(d*pi/(2*range))^2 / phi^2) ) }
