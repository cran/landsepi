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


#' @title Simulation with provided data
#' @name simul_landsepi
#' @description Simulation of the deployment of plant resistance, using landscape structures provided with the package
#' and a parameterisation of the model to represent pathogens as typified by rusts of cereals (e.g. stripe rust, stem rust
#' , and leaf rust of wheat and barley). All parameters are optional. See details for explanations.
#' @param seed an integer used as seed value (for random number generator).
#' @param idLan an integer giving the index of landscape structure (1 to 5).
#' @param propSR proportion of fields where resistance is deployed: (RC)/(SC+RC) or (RC1+RC2)/(SC+RC1+RC2). Must be between 0 and 1.
#' @param isolSR an integer giving the spatial aggregation of fields where resistance is deployed (1=highly fragmented, 2=balanced, 3=highly aggregated).
#' @param propRR when applicable (mixtures and mosaics only), relative proportion of the second resistant cultivar: (RC2)/(RC1+RC2). Must be between 0 and 1.
#' @param isolRR when applicable, an integer specifying the spatial (for mosaics: 1=highly fragmented, 2=balanced, 3=highly aggregated) or 
#' temporal (for rotations: 1=every year, 2=every two years, 3=every three years) aggregation of fields cultivated with the second resistant cultivar.
#' @param strat a character string specifying the deployment strategy ("MO"=mosaic, "MI"=mixture, "RO"=rotations, "PY"=pyramiding).
#' @param nHost an integer giving the number of cultivars (1, 2 or 3).
#' @param nYears an integer giving the number of simulated years.
#' @param pI0 initial probability of infection of the susceptible cultivar. Must be between 0 and 1.
#' @param resistance1 a logical vector of size 8 giving the resistance formula of the 2nd cultivar (see details)
#' @param resistance2 when applicable, a logical vector of size 8 giving the resistance formula of the 3rd cultivar (see details)
#' @param costInfect cost of infectivity paid by infective pathogens (i.e. adapted to plant cultivars carrying a major gene) on susceptible hosts. Must be between 0 and 1.
#' @param costAggr cost of aggressiveness paid by fully adapted pathogens (relative to plant cultivars carrying a quantitative resistance trait) on fully susceptible hosts. Must be between 0 and 1.
#' @param taumut mutation probability: probability for a propagule to change its infectivity or its aggressiveness on a resistant cultivar 
#' carrying a major gene or a quantitative resistance trait. Must be above 0. If equal to 0, then the pathogen cannot evolve.
#' @param MGeff efficiency of major-gene resistance on the infection rate of non-adapted pathogens. Must be between 0 and 1.
#' @param QReff efficiency of quantitative resistance on the target aggressiveness trait (infection rate, latent period duration, 
#' sporulation rate, or sporulation duration) of non-adapted pathogens. Must be between 0 and 1.
#' @param beta trade-off strength for pathogen adaptation to quantitative resistance (<1 for weak, =1 for linear, >1 for strong). Must be above 0.
#' @param nAggr an integer specifying the number of increments to completely adapt to quantitative resistance. Must be greater or equal 2.
#' @param graphOn a logical indicating if graphics must be generated (1) or not (0).
#' @details   \describe{
#' 
#' \item{Landscape structure}{The landscape structure is the physical structure of the area, defined as the spatial arrangement of fields.}
#' 
#' \item{Deployment strategies}{Deployment strategies include the deployment of a susceptible cultivar (SC)
#' and one (RC) or two (RC1 and RC2) resistant cultivars carrying up to four major resistance genes or up to four
#' quantitative resistance traits (against infection rate, latent period, sporulation rate and sporulation duration 
#' of the pathogen). In addition, the different resistance sources can be combined in time (crop rotation: recurrent
#' succession of cultivars in the same field), or space within a single cultivar (pyramiding), in different cultivars
#' of the same field (mixtures) or in different fields (mosaics).}
#' 
#' \item{Resistance formulas}{The genetic resistance carried by a plant cultivar is specified by a vector of size 8: the
#' four first elements indicate whether the cultivar carries major resistance genes #1, #2, #3 and #4, respectively.
#' The following four elements indicate whether the cultivar carried a quantitative resistance trait against the
#'  infection rate, the latent period duration, the sporulation rate, or the sporulation duration of the pathogen, respectively.
#'  For example, the formula c(1,0,0,0,0,1,0,0) indicates the presence of major gene #1 and a quantitative resistance which 
#'  increases the duration of the latent period of the pathogen.}
#' 
#' \item{Model outputs}{ 
#' 
#' \describe{
#' \item{\strong{Evolutionary outputs.}}{ \describe{
#' \item{\emph{Durability of qualitative resistance:}}{
#' For a given major gene, several computations are 
#' performed: \itemize{
#' \item (d1) time to first appearance of a pathogen mutant;
#' \item (d2) time to first true infection of a resistant host by such mutants; and 
#' \item (d3) time when the number of infections of resistant hosts by these mutants reaches a threshold above which mutant pathogens are unlikely to go extinct.
#' } 
#' }
#' \item{\emph{Erosion of quantitative resistance:}}{
#' pathogen adaptation to quantitative resistance is gradual, so the three measures described above are computed for every step 
#' towards complete erosion of resistance (i.e. nAgw-1 levels).} 
#' \item{\emph{Durability of a deployment strategy:}}{ a simulation run is divided into three periods: \enumerate{
#' \item the initial short-term period when all resistance sources are at their highest potential; 
#' \item a transitory period during which a given deployment strategy is only partially effective; and 
#' \item a longer-term period when all the resistances have been overcome or completely eroded. 
#' }
#' To assess the end of the short-term period, the time to establishment (durability measure (d3)) is computed for every major gene, and every quantitative trait 
#' at the first level of erosion (agw(p)=2). The minimal value of these measures, denoted by D1, delimitates short-term and transitory periods. 
#' Similarly, the time to establishment is computed for every major gene, and for every quantitative trait at the highest level of erosion (agw(p)=nAgw). 
#' The maximal value of these measures, termed D2, delimits transitory and long-term periods.
#' }
#' }
#' }
#' \item{\strong{Epidemiological outputs.}}{ 
#' The epidemiological impact of pathogen spread is evaluated by two different measures: \enumerate{
#' \item Green Leaf Area (GLA): The GLA represents the average number of productive hosts per time step and per surface unit.
#' \item Area Under Disease Progress Curve (AUDPC): The AUDPC is the average proportion of diseased hosts relative to the carrying capacity and represents disease severity.
#' }
#' \describe{
#' \item{\emph{Global epidemiological control:}}{
#' The GLA and AUDPC of every cultivar as well as the whole landscape are averaged across the whole simulation run, 
#'  to measure the global epidemiological performance of a deployment strategy.} 
#'  \item{\emph{Short-term epidemiological control:}}{
#'   The average GLA and AUDPC of the susceptible 
#'  cultivar is computed on whole cropping seasons from the beginning of the simulation until the end of 
#'  the season preceding year before D1.}
#'  \item{\emph{Epidemiological control during the transitory period:}}{
#'  The average GLA and AUDPC of the susceptible cultivar is computed on whole seasons from the 
#'  beginning of the season following year after D1 to the end of the season year before preceding D2.} 
#'  \item{\emph{Long-term epidemiological control:}}{
#'  The average GLA and AUDPC of the whole landscape is computed on whole seasons 
#'  from the beginning of the year after D2 to the end of the simulation.}
#'  }
#'  }
#'  }
#' }
#' }
#' @return A set of binary files is generated for every year of simulation and every compartment: \itemize{
#' \item H: healthy hosts,
#' \item Hjuv: juvenile healthy hosts,
#' \item L: latently infected hosts,
#' \item I: infectious hosts,
#' \item R: removed hosts,
#' \item S: propagules.
#' } 
#' Each file indicates for every time-step the number of individuals in each field, and when appropriate for each cultivar and pathotype)
#' These binary files are used to generate a set of text files containing all outputs of the simulations (see details).
#' A set of graphics and epidemic maps can also be generated.
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (in press). Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology}.
#' @examples \donttest{
#' ## Default parameterisation
#' simul_landsepi()
#' 
#' ## Mosaic of two major genes
#' simul_landsepi(seed=1, idLan=1, propSR=2/3, isolSR=3, propRR=1/2, isolRR=3, strat="MO", nHost=3
#' , nYears=50, resistance1=c(1,0,0,0,0,0,0,0), resistance2=c(0,1,0,0,0,0,0,0)
#' , costInfect=0.5, taumut=1e-7)
#' 
#' ## Mixture of two major genes
#' simul_landsepi(seed=1, idLan=1, propSR=2/3, isolSR=3, propRR=1/2, strat="MI", nHost=3
#' , nYears=50, resistance1=c(1,0,0,0,0,0,0,0), resistance2=c(0,1,0,0,0,0,0,0)
#' , costInfect=0.5, taumut=1e-7)
#' 
#' ## Rotations of two major genes
#' simul_landsepi(seed=1, idLan=1, propSR=2/3, isolSR=3, isolRR=1, strat="RO", nHost=3
#' , nYears=50, resistance1=c(1,0,0,0,0,0,0,0), resistance2=c(0,1,0,0,0,0,0,0)
#' , costInfect=0.5, taumut=1e-7)
#' 
#' ## Pyramiding of two major genes
#' simul_landsepi(seed=1, idLan=1, propSR=2/3, isolSR=3, strat="PY", nHost=2
#' , nYears=50, resistance1=c(1,1,0,0,0,0,0,0), costInfect=0.5, taumut=1e-7)
#' 
#' ## Combination of a major gene with a quantitative resistance against the latent period    
#' simul_landsepi(seed=1, idLan=1, propSR=0.8, isolSR=1, strat="PY", nHost=2
#' , nYears=50, resistance1=c(1,0,0,0,0,1,0,0)
#' , costInfect=0.5, costAggr=0.5, taumut=1e-7, MGeff=1.0, QReff=0.5, beta=1.0, nAggr=6)
#' }
#' @include RcppExports.R AgriLand.R graphLand.R  multiN.R  periodic_cov.R 
#' @importFrom utils data
#' @export
simul_landsepi <- function(seed=12345, idLan=1, propSR=2/3, isolSR=3, propRR=1/2, isolRR=3, strat="MO", nHost=3, nYears=5, pI0=5e-4
                          , resistance1=c(1,0,0,0,0,0,0,0), resistance2=c(1,1,0,0,0,0,0,0)
                          , costInfect=0.75, costAggr=0.5, taumut=1e-7, MGeff=1.0, QReff=0.5, beta=1.0, nAggr=6
                          , graphOn=1){

    pathRES <- getwd()

    ## Time parameters
    nTSpY <- 120
    paramT <- list(nYears=nYears,nTSpY=nTSpY)
    
    ## Landscape and dispersal 
    dispH <- get("dispH")    # Hack for cran check 
    
    landscape <- NULL
    dispP <- NULL
    if (idLan==1){
        landscapeTEST1 <- get("landscapeTEST1")
        dispP_1 <- get("dispP_1")   # Hack for cran check 
        landscape <- landscapeTEST1
        dispP <- dispP_1
    }else if (idLan==2){
        landscapeTEST2 <- get("landscapeTEST2")
        dispP_2 <- get("dispP_2")   # Hack for cran check 
        landscape <- landscapeTEST2
        dispP <- dispP_2
    }else if (idLan==3){
        landscapeTEST3 <- get("landscapeTEST3")
        dispP_3 <- get("dispP_3")   # Hack for cran check 
        landscape <- landscapeTEST3
        dispP <- dispP_3
    }else if (idLan==4){
        landscapeTEST4 <- get("landscapeTEST4")
        dispP_4 <- get("dispP_4")   # Hack for cran check 
        landscape <- landscapeTEST4
        dispP <- dispP_4
    }else if (idLan==5){
        landscapeTEST5 <- get("landscapeTEST5")
        dispP_5 <- get("dispP_5")   # Hack for cran check 
        landscape <- landscapeTEST5
        dispP <- dispP_5 }


    ## Host parameters
    C_0 <- 0.1
    Cmax0 <- 2
    Cmax1 <- 2
    croisH0 <- 0.10
    reproH0 <- 0.0
    croisH1 <- 0.10
    reproH1 <- 0.0
    deathH <- 0.0
    resistance0 <- c(0,0,0,0,0,0,0,0)
    resistanceMatrix <- cbind(resistance0,resistance1,resistance2)
    resistance <- as.vector(resistanceMatrix)
    if (nHost>1){
        adaptation <- apply(resistanceMatrix[,1:nHost],1,sum)
    }else{
        adaptation <- resistance0
    }
    adaptation[adaptation>0] <- 1
    khost <- 0.002
    sighost <- 1.001
    shost <- 1.0
    paramH <- list(Nhote=nHost,croisH0=croisH0,reproH0=reproH0,croisH1=croisH1,reproH1=reproH1,deathH=deathH,resistance=as.integer(resistance),khost=khost,sighost=sighost,shost=shost)

    ## Pathogen parameters
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
    
    ## Evolution parameters
    paramevol <- list(costinfect=costInfect,costaggr=costAggr,taumut=taumut,MGeff=MGeff,QReff=QReff,beta=beta,Naggr=nAggr,adaptation=as.integer(adaptation))
    
    ## Landscape
    # hack for cran check 
    landscape <- AgriLand(landscape,filename="landscape",propSR,isolSR,propRR,isolRR,strat,nHost,nYears,Cmax0,Cmax1,seed,graphOn)
    
    #run the model!
    modelLandsEPI(paramT,
                   landscape,
                   dispersal=list(dispP=dispP,dispH=dispH),
                   inits=list(C_0=C_0, PI0=pI0),
                   val_seed=seed,
                   hostP=paramH,
                   epiP=paramepi,
                   evolP=paramevol)

    ## generate the output
    HLIRdynamics(pathRES, graphOn, paramT, landscape, paramH, epiP=paramepi, evolP=paramevol,nMapPY=0)
}
