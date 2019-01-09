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


#' @title Logit function
#' @name logit
#' @description Given a numeric object return the logit of the values. Missing values (NAs) are allowed. 
#' @param x a numeric object containing values between 0 and 1
#' @details The logit is defined by log(x/(1-x)). Values in x of 0 or 1 return logits of -Inf or Inf respectively. Any NAs in the input will also be NAs in the output. 
#' @return An object of the same type as x containing the logits of the input values.
#' @examples
#' logit(0.5)
#' @export
logit <- function(x) log(x/(1-x))
