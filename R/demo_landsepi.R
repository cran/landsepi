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


#' Run demo landsepi
#' @title Package demonstration
#' @name demo_landsepi
#' @description Run a demonstration of the package.
#' @param seed an integer used as seed value (for random number generator)
#' @details Run a 30-year simulated example of mosaic deployment strategy of two resistant cultivars in balanced proportions and 
#' high level of spatial aggregation. The generated model outputs are graphics and text files.
#' @include RcppExports.R AgriLand.R graphLand.R  multiN.R  periodic_cov.R 
#' @importFrom utils data
#' @export
demo_landsepi <- function(seed=12345){
    
    pathRES <- getwd()
    graphOn <- 1
    
    #landscape parameters
    propSR <- 2/3    ## proportion of cultivars >0
    isolSR <- 3      ## Class of aggregation between S and R (0:3 with 1=low, 2=medium, 3=high, 4=random aggregation)
    propRR <- 1/2    ## relative proportion of cultivars >1
    isolRR <- 3      ## Class of aggregation between R1 and R2 (0:3 with 1=low, 2=medium, 3=high, 4=random aggregation)
    strat  <- "MO"
    Nhote  <- 3 
    nYears <- 30
    Cmax0  <- 2
    Cmax1  <- 2
    
    ## hack for cran check / do not use it
    landscapeTEST1 <- get("landscapeTEST1")
    ## use landscapeTEST1 directly
    landscape <- AgriLand(landscapeTEST1,filename="landscapeTEST1",propSR,isolSR,propRR,isolRR,strat,Nhote,nYears,Cmax0,Cmax1,seed,graphOn)

    #time parameters
    nTSpY = 120
    #timesStep=1:120#times step to be saved
    paramT <- list(nYears=nYears,nTSpY=nTSpY)#,timesStep=timesStep)

    #dispersal    
    ## Hack for cran check 
    dispP_1 <- get("dispP_1")
    dispH <- get("dispH")

    #host parameters
    croisH0 <- 0.10
    reproH0 <- 0.0
    croisH1 <- 0.10
    reproH1 <- 0.0
    deathH <- 0.0
    RESISTANCE0 <- c(0,0,0,0,0,0,0,0)
    RESISTANCE1 <- c(1,0,0,0,0,0,0,0)
    RESISTANCE2 <- c(0,1,0,0,0,0,0,0)
    RESISTANCE <- as.vector(cbind(RESISTANCE0,RESISTANCE1,RESISTANCE2))
    khost <- 0.002
    sighost <- 1.001
    shost <- 1.0
    paramH <- list(Nhote=Nhote,croisH0=croisH0,reproH0=reproH0,croisH1=croisH1,reproH1=reproH1,deathH=deathH,resistance=as.integer(RESISTANCE),khost=khost,sighost=sighost,shost=shost)

    #epidemiology
    PSURV <- 1e-4
    EFF <- 0.4
    REPROP <- 3.125
    TLATEXP <- 10
    TLATVAR <- 9
    TSPOEXP <- 24
    TSPOVAR <- 105
    kpatho <- 5.333
    sigpatho <- 3
    spatho <- 1
    paramepi <- list(psurv=PSURV,eff=EFF,reproP=REPROP,TlatEXP=TLATEXP,TlatVAR=TLATVAR,TspoEXP=TSPOEXP,TspoVAR=TSPOVAR,kpatho=kpatho,sigpatho=sigpatho,spatho=spatho)
    
    #evolution
    COSTINFECT <- 0.75
    COSTAGGR <- 0.5
    TAUMUT <- 1e-7
    MGEFF <- 1.0
    QREFF <- 0.5
    BETA <- 1.0
    NAGGR <- 6
    ADAPTATION <- c(1,1,0,0,0,0,0,0)
    paramevol <- list(costinfect=COSTINFECT,costaggr=COSTAGGR,taumut=TAUMUT,MGeff=MGEFF,QReff=QREFF,beta=BETA,Naggr=NAGGR,adaptation=as.integer(ADAPTATION))
    
    #initial conditions
    C_0 <- 0.1
    PI0 <- 5e-4
    
    #run the model!
    modelLandsEPI(paramT,
                   landscape,
                   dispersal=list(dispP=dispP_1,dispH=dispH),
                   inits=list(C_0=C_0, PI0=PI0),
                   val_seed=seed,
                   hostP=paramH,
                   epiP=paramepi,
                   evolP=paramevol)

    ## generate the output
    HLIRdynamics(pathRES, graphOn, paramT, landscape, paramH, epiP=paramepi, evolP=paramevol,nMapPY=0)
}
