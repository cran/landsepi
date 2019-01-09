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


#' Run demo landsepi
#' @title Package demonstration
#' @name demo_landsepi
#' @description Run a demonstration of the package.
#' @param seed an integer used as seed value (for random number generator)
#' @details Run a 30-year simulated example of mosaic deployment strategy of two resistant cultivars in balanced proportions and 
#' high level of spatial aggregation. The generated model outputs are text files, graphics and a video.
#' @include RcppExports.R AgriLand.R graphLand.R  multiN.R  periodic_cov.R 
#' @importFrom utils data
#' @export
demo_landsepi <- function(seed=12345){
    pathRES <- getwd()
    simul_landsepi(seed, nYears=30, idLan=1, propSR=2/3, isolSR=3, propRR=1/2, isolRR=3, strat="MO", Nhost=3
                               , resistance1=c(1,0,0,0,0,0,0,0), resistance2=c(0,1,0,0,0,0,0,0), taumut=1e-7
                               , graphic=TRUE, video=TRUE)
}
