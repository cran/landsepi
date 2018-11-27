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


#' @title Plotting pathotype frequencies
#' @name plotevolQR
#' @description Plot the dynamic of pathotype frequencies in a tiff file. 
#' @param pathRES the path of the repository where the graphic will be generated
#' @param nIncr number of pathotypes
#' @param trait index of the evolving trait (1 to 4 for major genes, 5 to 8 for quantitative resistance traits against infection rate (5), latent period (6), reproduction rate (7) and infectious period duration (8)) 
#' @param I_aggrProp a matrix giving the frequency of every pathotype (row) for every time-step (columns)
#' @param D durability of the trait
#' @param nTS number of simulated time-steps
#' @param nYears number of simulated years
#' @param nTSpY number of time-steps per year
#' @examples
#' \dontrun{
#' freqMatrix <- matrix(0, nrow=2, ncol=100)
#' freqMatrix[2,26:100] <- (26:100)/100
#' freqMatrix[1,] <- 1-freqMatrix[2,]
#' plotevolQR(getwd(), nIncr=2, trait=1, freqMatrix, D=25, nTS=100,nYears=10,nTSpY=10)
#' }
#' @export
plotevolQR <- function(pathRES, nIncr,trait, I_aggrProp, D, nTS,nYears,nTSpY) {
    COL.grey <- gray(0:150/150)
    COL.grey <- COL.grey[length(COL.grey):1]
    TITLE <- "Aggressiveness trait"
    LABELS <- c("\nS specialist",rep(NA,nIncr-2),"\nR specialist")
    if (trait < 5) {
        TITLE <- "Infectivity gene"
        LABELS <- c("\nNot infective","\nInfective")
    }
    
    tiff(filename=paste(pathRES,"/EVOLpatho_", names(D), ".tiff",sep=""),width=100,height=100,units='mm',compression='lzw',res=300)
    par(xpd=F, mar=c(4,4,2,2))
    image(x=1:nIncr, y=1:nTS, z=I_aggrProp, col=COL.grey, ylab="", xlab="Phenotype", main=paste(TITLE,names(D)), axes=F, zlim=c(0,1), ylim=c(1,nTS+1))
    box()
    axis(side=1, at=1:nIncr, labels=LABELS)
    if (nYears==1) {
        axis(2, at=round(seq(1,nTS,length.out=8)), las=1)
        title(ylab="Evolutionnary time (days)")
    } else {
        axis(2, at=seq(1,nTS+1,nTSpY*((nYears-1)%/%8+1)), labels=seq(0,nYears,((nYears-1)%/%8+1)), las=1)
        title(ylab="Evolutionnary time (years)")
    }
    if (!is.na(D) & D<=nTS)
        abline(h=D, col="blue", lty=trait+1, lwd=2.5)  ## durability of the trait
    dev.off()     
}





