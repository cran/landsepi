# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inra.fr>
#                    Julien Papaix <julien.papaix@inra.fr>
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


#' @title Inverse logit function
#' @name invlogit
#' @description Given a numeric object return the invlogit of the values. Missing values (NAs) are allowed. 
#' @param x a numeric object
#' @details The invlogit is defined by exp(x) / (1+exp(x)). Values in x of -Inf or Inf return invlogits of 0 or 1 respectively. Any NAs in the input will also be NAs in the output. 
#' @return An object of the same type as x containing the invlogits of the input values.
#' @examples
#' invlogit(10)
#' @export
invlogit <- function(x) exp(x) / (1+exp(x))
