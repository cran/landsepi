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


#' @title Allocation of cultivars
#' @name multiN
#' @description Algorithm based on latent Gaussian fields to allocate two different types of crops across the simulated landscapes. 
#' @param d a matrix of the pairwise distances between the centroids of the fields of the landscape.
#' @param area vector containing the areas of the fields.
#' @param range range of spatial autocorrelation between fields (<0 for fragmented landscapes, >0 for aggregated landscapes, =0 for random allocation).
#' @param prop relative proportion of the second type of crop. 
#' @details This algorithm allows the control of the proportions of each type of crop in terms of surface coverage, and their level of spatial aggregation. 
#' A random vector of values is drawn from a multivariate normal distribution with expectation 0 and a variance-covariance matrix 
#' which depends on the pairwise distances between the centroids of the fields. The variance-covariance matrix is computed from a periodic function for highly 
#' fragmented or highly aggregated types of crops, an exponential function for moderately aggregated types of crops, and from a normal distribution for a random allocation of 
#' the types of crops. Next, the two types of crops are allocated to different fields depending on whether the each value drawn from the multivariate normal distribution 
#' is above or below a threshold. The proportion of each type of crop in the landscape is controlled by the value of this threshold (parameter prop). 
#' @return A dataframe containing the index of each field (column 1) and the index (0 or 1) of the type of crop grown on these fields (column 2).
#' @examples
#' \dontrun{
#' d <- matrix(rpois(100,100), nrow=10)
#' area <- data.frame(num=1:10, area=10)
#' ## Fragmented landscape
#' multiN(d, area, range=-2, prop=0.5)
#' }
#' @export
multiN <- function(d, area, range, prop) {
     nPoly <- nrow(area)

     ## Multivariate normal distribution
     if (range != 0) { 
          if (range > 0)  ## aggregation
               covMat <- Exponential(d, range)  ## function from 'fields' package
          if (range < 0) {           ## repulsion
               covMat_tmp <- periodic_cov(d, range)
               covMatPD <- nearPD(covMat_tmp, keepDiag=TRUE, ensureSymmetry=TRUE, maxit=1000)    ## conversion in a positive-definite matrix
               if (covMatPD$converge==TRUE) {
                    print("OK: Periodic covariance could be converted in a positive-definite matrix")
                    covMat <- covMatPD$mat
               } else {
                    print(paste("WARNING: covariance cannot be converted in a positive-definite matrix -- range =", range))
                    print("Exponential kernel used instead")
                    range2 <- 1*(range > -200) + 2000*(range < -200)
                    covMat <- Exponential(d, abs(range2))
               }
          }  ## if range < 0

          s <- mvrnorm(1,mu=rep(0,nPoly),Sigma=covMat)  ## function from 'MASS' package  

     } else {
          s <- rnorm(nPoly,mean=0,sd=1)  ## random allocation
     }
     
     area$s <- s
     ## Ordered gaussian beam
     ord_S <- order(s)
     area.ord_S <- area[ord_S,]
     area.ord_S$cumsum <- cumsum(area.ord_S$area)
     
     ## Threshold for allocation of the second type of crop
     prop.areaTot <- prop * area.ord_S$cumsum[nPoly]   
     area.ord_S$dif <- abs(area.ord_S$cumsum - prop.areaTot)
     th_index <- which.min(area.ord_S$dif)
     area.ord_S$cultivar <- as.numeric( (1:nPoly) < th_index )
     
     ## Final landscape
     habitat <- area.ord_S[order(area.ord_S$num), c("num", "cultivar")]
     
     return(habitat)
     # plot(centroid, col=grey((s-min(s))/diff(range(s))), pch=16, cex=2)
     # plot(centroid, col=(s<0)+1, pch=16, cex=2)
}


