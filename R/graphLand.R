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



#' @title Plotting the landscape
#' @name plotland
#' @description Plot a landscape with colors or hatched lines to represent different types of fields
#' @param landscape a spatialpolygon object containing field coordinates 
#' @param COL vector containing the color of each field
#' @param DENS vector containing the density of hatched lines for each field
#' @param ANGLE vector containing the angle of hatched lines for each field
#' @param COL.LEG vector containing the colors in the first legend
#' @param DENS.LEG vector containing the density of hatched lines in the second legend
#' @param ANGLE.LEG vector containing the angle of hatched lines in the second legend
#' @param TITLE title of the graphic
#' @param SUBTITLE subtitle of the graphic
#' @param LEGEND1 labels in the first legend (colors)
#' @param LEGEND2 labels in the second legend (hatched lines)
#' @param TITLE.LEG2 title for the second legend
#' @param XMAX maximal coordinate on horizontal axis
#' @param YMAX maximal coordinate on vertical axis
#' @examples
#' \dontrun{
#' ## Draw a landscape with various colours
#' landscapeTEST1
#' plotland(landscapeTEST1, COL=1:length(landscapeTEST1),
#' DENS=rep(0,length(landscapeTEST1)), ANGLE=rep(30,length(landscapeTEST1)))
#' }
#' @include RcppExports.R logit.R invlogit.R
# @S3method plot land
#' @export
plotland <- function(landscape, COL=rep(0,length(landscape)), DENS=rep(0,length(landscape)), ANGLE=rep(30,length(landscape))
                      , COL.LEG=unique(COL), DENS.LEG=unique(DENS), ANGLE.LEG=unique(ANGLE)
                      , TITLE="", SUBTITLE="", LEGEND1=rep("", length(COL.LEG)), LEGEND2=rep("", length(COL.LEG)), TITLE.LEG2="", XMAX=2000, YMAX=2000) {
     par(cex=2, xpd=NA, bg="white", mar=c(5,4,4,2))
     nPoly <- length(landscape)
     plot(0,0, xlim=c(0,XMAX), ylim=c(0,YMAX), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", type="n", main=TITLE)    # Empty graph
     mtext(SUBTITLE[1], side=3, line=0, padj=-.5, las=1, cex=1.4)
     if (length(SUBTITLE)>1)
          mtext(SUBTITLE[2], side=3, line=0, padj=1.5, las=1, cex=1.4)
     for (i in 1:nPoly) {
          polymap(landscape@polygons[[i]]@Polygons[[1]]@coords, add=TRUE, col=COL[i], border="black", lwd=2.5)
          polymap(landscape@polygons[[i]]@Polygons[[1]]@coords, add=TRUE, col="black", density=DENS[i], angle=ANGLE[i], border=NA)
     }
     if (LEGEND1[1]!=""){
         if (TITLE.LEG2=="")
             legend(XMAX/2.66, -YMAX/40, legend=LEGEND1, fill=COL.LEG, bty="n")
         else {
             legend(XMAX/2.66, -YMAX/40, legend=LEGEND1, col="black", density=2*DENS.LEG, angle=ANGLE.LEG, bty="n")
             legend(-XMAX/5, YMAX, legend=LEGEND2, fill=COL.LEG, bty="n", title=TITLE.LEG2)
         }
     }
}
