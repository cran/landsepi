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


#' @title Generation of model outputs
#' @name HLIRdynamics
#' @description Generate epidemiological and evolutionary outputs from model simulations.
#' @param pathRES a character string indicating the path of the repository where outputs will be generated.
#' @param timeP a list of simulation parameters (number of years, number of time-steps per year).
#' @param landscape a shapefile containing the agricultural landscape (can be generated through function AgriLand).
#' @param hostP a list of host parameters (number of cultivars, initial planting density, maximal carrying capacity, growth rate, reproduction rate, 
#' death rate, number of possible resistance sources (8), resistance formula, parameters of the sigmoid invasion function: kappa, sigma and s).
#' @param evolP a list of evolution parameters (cost of infectivity, cost of aggressiveness, mutation rate, efficiency of major 
#' resistance genes, efficiency of quantitative resistance, trade-off strength, number of increments of quantitative 
#' resistance erosion, adaptation formula).
#' @param ecoP a list of economic parameters (yield, purchase price, selling price).
#' @param graphic a logical indicating if graphics of the outputs must be generated (TRUE) or not (FALSE).
#' @param video a logical indicating if a video must be generated (TRUE) or not (FALSE). Works only if graphic is TRUE as well. Note that ffmpeg is required to generate videos.
#' @param th_break an integer giving the threshold (number of infections) above which mutant pathogen are unlikely to go extinct, used to 
#' characterise resistance breakdown.
#' @param nMapPY an integer specifying the number of epidemic maps per year in the video (video must be TRUE).
#' @details \describe{
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
#' @return A set of text files containing all outputs of the simulations (see details). 
#' A set of graphics and a video with epidemic maps can also be generated.
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018). Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#' @examples \donttest{
#' demo_landsepi()
#' }
#' @include RcppExports.R logit.R invlogit.R graphLand.R
#' @importFrom sf st_read
#' @importFrom grDevices colorRampPalette dev.off graphics.off gray png tiff
#' @importFrom utils write.table
#' @export
HLIRdynamics <- function(pathRES, timeP, landscape, hostP, evolP, ecoP, graphic, video=FALSE, nMapPY=5, th_break=50000){

    nYears <- timeP$nYears
    nTSpY <- timeP$nTSpY
    
    Nhost <- hostP$Nhost
    resistance <- matrix(hostP$resistance,ncol=3)
    
    landscape_tmp <- st_read(landscape$shapefilename,layer=landscape$layername_hab)
    propRR <- landscape$propRR
    strat <- landscape$strat
    
    Naggr <- evolP$Naggr    
    adaptation <- evolP$adaptation  
    
    C0 <- hostP$C0
    Cmax <- hostP$Cmax
    yield <- ecoP$yield
    sellPrice <- ecoP$sellPrice
    purchPrice <- ecoP$purchPrice

#### LANDSCAPE ####

    area <- landscape_tmp$area
    hab1 <- landscape_tmp$habitat0
    hab2 <- landscape_tmp$habitat1
    habitat <- cbind(hab1, hab2)+1  ## +1 due to 0 as first index in C
    rotation <- landscape$rotation+1
    Npoly <- length(area)

    trans.landscape_tmp <- as(landscape_tmp,"Spatial")
#    trans.landscape <- list()
#    ipoly <- 1
#    for (i in 1:length(trans.landscape_tmp)){
#        toto <- trans.landscape_tmp@polygons[[i]]
#        for (j in 1:length(toto)){
#            titi <- toto@Polygons[[j]]@coords
#            trans.landscape[[ipoly]] <- as.poly(titi)
#            ipoly <- ipoly+1
#        }
#    }
#    class(trans.landscape) <- "listpoly"
    
K <- matrix(nrow=Npoly, ncol=Nhost)
for (host in 1:Nhost)
     K[,host] <- floor(area * Cmax[host])   ## Carrying capacity of each host in each paddock
 

# plot(landscape, color=F, border="darkgrey", lwd=1.5)
# title("Landscape structure")
# mtext("m", side=2, line=2.5, padj=0, cex=1, las=1)
# mtext("m", side=1, line=0, padj=3.5, cex=1, las=1)
# summary(unlist(landscape))   ## Bounds of the landscape

#### IMPORTATION OF THE SIMULATION OUTPUT ####
nTS <- nYears * nTSpY          ## number of time-steps
adapt_vect <- adaptation
Npatho <- 2^sum(adapt_vect[1:4]) * Naggr^sum(adapt_vect[5:8])
res_vect <- t(resistance)
Rhosts <- (1:nrow(res_vect))[apply(res_vect,1,sum)>0]
Rhosts <- Rhosts[Rhosts<=Nhost]

H <- as.list(1:nTS)
# Hjuv <- as.list(1:nTS)
S <- as.list(1:nTS)
L <- as.list(1:nTS)
I <- as.list(1:nTS)
R <- as.list(1:nTS)
index <- 0

print("Reading binary files to compute model outputs")
for (year in 1:nYears) {
     print(paste("year", year, "/", nYears))
     binfileH <- file(paste(pathRES,sprintf("/H-%02d", year), ".bin",sep=""), "rb")
      H.tmp <- readBin(con=binfileH, what="int", n=Npoly*Nhost*nTSpY, size = 4, signed=T, endian="little")
     # binfileHjuv = file(paste(pathRES, sprintf("/Hjuv-%02d", year), ".bin",sep=""), "rb")
     #  Hjuv.tmp <- readBin(con=binfileHjuv, what="int", n=Npoly*Nhost*nTSpY, size = 4, signed=T,endian="little")
     binfileS = file(paste(pathRES, sprintf("/S-%02d", year), ".bin",sep=""), "rb")
      S.tmp <- readBin(con=binfileS, what="int", n=Npoly*Npatho*nTSpY, size = 4, signed=T,endian="little")
     binfileI = file(paste(pathRES, sprintf("/I-%02d", year), ".bin",sep=""), "rb")
      I.tmp <- readBin(con=binfileI, what="int", n=Npoly*Npatho*Nhost*nTSpY, size = 4, signed=T,endian="little")
     binfileL = file(paste(pathRES, sprintf("/L-%02d", year), ".bin",sep=""), "rb")
      L.tmp <- readBin(con=binfileL, what="int", n=Npoly*Npatho*Nhost*nTSpY, size = 4, signed=T,endian="little")
     binfileR = file(paste(pathRES, sprintf("/R-%02d", year), ".bin",sep=""), "rb")
      R.tmp <- readBin(con=binfileR, what="int", n=Npoly*Npatho*Nhost*nTSpY, size = 4, signed=T,endian="little")
     
     for (t in 1:nTSpY) {
          H[[t + index]] <- matrix(H.tmp[((Nhost*Npoly)*(t-1)+1):(t*(Nhost*Npoly))],ncol=Nhost,byrow=T)
          # Hjuv[[t + index]] <- matrix(Hjuv.tmp[((Nhost*Npoly)*(t-1)+1):(t*(Nhost*Npoly))],ncol=Nhost,byrow=T)
          S[[t + index]] <- matrix(S.tmp[((Npatho*Npoly)*(t-1)+1):(t*(Npatho*Npoly))],ncol=Npatho,byrow=T)
          L[[t + index]] <- array(data=L.tmp[((Npatho*Npoly*Nhost)*(t-1)+1):(t*(Npatho*Npoly*Nhost))],dim=c(Nhost,Npatho,Npoly))
          I[[t + index]] <- array(data=I.tmp[((Npatho*Npoly*Nhost)*(t-1)+1):(t*(Npatho*Npoly*Nhost))],dim=c(Nhost,Npatho,Npoly))
          R[[t + index]] <- array(data=R.tmp[((Npatho*Npoly*Nhost)*(t-1)+1):(t*(Npatho*Npoly*Nhost))],dim=c(Nhost,Npatho,Npoly))
     }
     index <- index + nTSpY
     close(binfileH)
     # close(binfileHjuv)
     close(binfileS)
     close(binfileL)
     close(binfileI)
     close(binfileR)
}

####  HLIR DYNAMIC  ####
H_host <- NULL
L_host <- NULL
I_host <- NULL
R_host <- NULL
S_patho <- NULL
L_patho <- NULL
I_patho <- NULL

IL_patho_Rhost <- NULL ## only on resistant hosts
for (t in 1:nTS){
    H_host <- cbind(H_host, apply(H[[t]],2,sum))
    L_host <- cbind(L_host, apply(L[[t]],1,sum))
    I_host <- cbind(I_host, apply(I[[t]],1,sum))
    R_host <- cbind(R_host, apply(R[[t]],1,sum))  
    N_host <- H_host + L_host + I_host + R_host
    S_patho <- cbind(S_patho, apply(S[[t]],2,sum))
    L_patho <- cbind(L_patho, apply(L[[t]],2,sum))
    I_patho <- cbind(I_patho, apply(I[[t]],2,sum))
    if (length(Rhosts)>0 & Npatho>1)
        IL_patho_Rhost <- cbind(IL_patho_Rhost, apply(L[[t]][Rhosts,,]+I[[t]][Rhosts,,],1 + (length(Rhosts)>1),sum))  ## length(Rhosts)=2 --> need to use the 2nd dimension
}

I_pathoProp <- matrix(NA, nrow=Npatho, ncol=nTS)
Itot <- apply(matrix(I_patho, ncol=nTS),2,sum)  ## matrix to avoid pb if only 1 row
for (p in 1:Npatho){
    for (t in 1:nTS){
        if (Itot[t] > 0)
            I_pathoProp[p,t] <- I_patho[p,t] / Itot[t]
    }
}


## Calculation of the carrying capacity and planting quantity
K_host <- matrix(nrow=Nhost, ncol=nTS)
C_host <- matrix(nrow=Nhost, ncol=nYears)
for (y in 1:nYears){
     hab <- habitat[,rotation[y]]
     ts_year <- ((y-1)*nTSpY+1) : (y*nTSpY)
     for (host in 1:Nhost) {
         C_host[host, y] <- C0[host] * sum(area*(hab==host))
          K_host[host,ts_year] <- sum(K[,host]*(hab==host)) 
          if (strat=="MI" & host>1){
              C_host[host, y] <- C_host[host, y] + C0[host] * sum(area*(habitat[,2]==host))
               K_host[host,ts_year] <- K_host[host,ts_year] + sum(K[,host]*(habitat[,2]==host))
          } ## if strat=="MI"
     } ## for host
}  ## for y

## Write computed results
# write(as.vector(H_host),file=paste(pathRES,"/Hhost.txt",sep=""),sep=",")
# write(as.vector(L_host),file=paste(pathRES,"/Lhost.txt",sep=""),sep=",")
# write(as.vector(I_host),file=paste(pathRES,"/Ihost.txt",sep=""),sep=",")
# write(as.vector(R_host),file=paste(pathRES,"/Rhost.txt",sep=""),sep=",")
# write(as.vector(I_pathoProp),file=paste(pathRES,"/IpathoProp.txt",sep=""),sep=",")

#############################################
#####               OUTPUTS              ####
#############################################

#### Extinction of the pathogen if not present any more at the last time-step (L=I=0)
#### No breakdown of resistance if no mutant is present at the last time-step
EXT <- c(ext=NA, noMut=NA)
EXT["ext"] <- as.numeric(sum(L_host[,nTS], I_host[,nTS]) == 0)
if (Npatho>1)
    EXT["noMut"] <- as.numeric(sum(L_patho[2:Npatho,nTS],I_patho[2:Npatho,nTS]) == 0)
EXT

####      A)    AUDPC  &   GLA       ####
## as.numeric to allow memory allocation of long integer

## _I=only infectious sites
## _IR=infectious + removed sites
audpc_IR <- matrix(NA, nrow=Nhost, ncol=nYears)  ## Average proportion of diseased host relative to carrying capacity (per cultivar, per year)
audpc_year <- matrix(NA, ncol=nYears, nrow=1)   ## --- for whole landscape (per year)
# audpc_I <- matrix(NA, nrow=Nhost, ncol=nYears)
# peak_I <- matrix(NA, nrow=Nhost, ncol=nYears)   ## elevation of the peak
# peak_IR <- matrix(NA, nrow=Nhost, ncol=nYears)
# timeToPeak_I <- matrix(NA, nrow=Nhost, ncol=nYears)   ## time to peak
# timeToPeak_IR <- matrix(NA, nrow=Nhost, ncol=nYears)

GLAabs <- matrix(NA, ncol=nYears, nrow=Nhost)  ## Average absolute number of healthy host per timestep and square meter (per cultivar, per year)
GLArel <- matrix(NA, ncol=nYears, nrow=Nhost)  ## Average proportion of healthy hosts relative to the number of hosts (H+L+I+R) (per cultivar, per year)

product <- matrix(NA, ncol=nYears, nrow=Nhost)    ## Crop production (weight unit of fruits/grains/etc per host and per year)

for (host in 1:Nhost){
     for (y in 1:nYears){
          ts_year <- ((y-1)*nTSpY+1):(y*nTSpY)
          GLAabs[host,y] <- sum(as.numeric(H_host[host,ts_year])) / (nTSpY * sum(area))
          audpc_year[,y] <- sum(as.numeric(I_host[,ts_year]),as.numeric(R_host[,ts_year])) / sum(as.numeric(K_host[,ts_year]))

          product[host,y] <- sum(yield[host,"H"] * as.numeric(H_host[host,ts_year]), yield[host,"L"] * as.numeric(L_host[host,ts_year]), yield[host,"I"] * as.numeric(I_host[host,ts_year]), yield[host,"R"] * as.numeric(R_host[host,ts_year]))
          
          if (sum(as.numeric(K_host[host,ts_year])) > 0){   ## to avoid pb with rotations
               # propI <- I_host[host,ts_year] / K_host[host,ts_year]
               # propIR <- (I_host[host,ts_year]+R_host[host,ts_year]) / K_host[host,ts_year]
               # audpc_I[host,y] <- sum(as.numeric(I_host[host,ts_year])) / sum(as.numeric(K_host[host,ts_year]))
               audpc_IR[host,y] <- sum(as.numeric(I_host[host,ts_year]),as.numeric(R_host[host,ts_year])) / sum(as.numeric(K_host[host,ts_year]))
               # peak_I[host,y] <- max(propI)
               # peak_IR[host,y] <- max(propIR)
               # timeToPeak_I[host,y] <- order(propI, decreasing=TRUE)[1]
               # timeToPeak_IR[host,y] <- order(propIR, decreasing=TRUE)[1]
          } ## if K>0
          
          if ( sum(as.numeric(N_host[host,ts_year])) > 0 ) ## (i.e. if the cultivar is cultivated this year)
              GLArel[host,y] <- sum(as.numeric(H_host[host,ts_year])) / sum(as.numeric(N_host[host,ts_year]))
     } ## for y
} ## for host

audpc <- audpc_IR
audpc_host <- c(S=NA, R1=NA, R2=NA)
audpc_host[1:Nhost] <- apply(audpc, 1, mean, na.rm=TRUE)  ## Average AUDPC for whole simulation (per cultivar)
audpc_host
product

## Write AUDPC & GLA
# write.table(GLAabs,file=paste(pathRES,"/GLAabs_detail.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
# write.table(GLArel,file=paste(pathRES,"/GLArel_detail.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
# write.table(audpc,file=paste(pathRES,"/audpc_detail.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
# write.table(audpc_year,file=paste(pathRES,"/audpc_year.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(matrix(audpc_host, nrow=1),file=paste(pathRES,"/audpc_host.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)

####       B) Durability of resistant sources       ####
## (=NA if resistant gene not present ; =nYears + 1 timestep if longer than the simulated period)
traits <- c("mg1", "mg2", "mg3", "mg4", "infection", "latency", "repro", "infPeriod")

## pathoToAggr matrix
nIGcum <- cumsum(adapt_vect[1:4])
nAGcum <- cumsum(adapt_vect[5:8])
pathoToAggr <- matrix(0, nrow=Npatho, ncol=8)
colnames(pathoToAggr) <- traits
for (p in 1:Naggr^(adapt_vect[8])-1) {
    for (o in 1:Naggr^(adapt_vect[7])-1) {
        for (n in 1:Naggr^(adapt_vect[6])-1) {
            for (m in 1:Naggr^(adapt_vect[5])-1) {
                for (l in 1:2^(adapt_vect[4])-1) {
                    for (k in 1:2^(adapt_vect[3])-1) { 
                        for (j in 1:2^(adapt_vect[2])-1) {
                            for (i in 1:2^(adapt_vect[1])-1) {
                                line = 1 + i + j*2^nIGcum[1] + k*2^nIGcum[2] + l*2^nIGcum[3] + m*(2^nIGcum[4])
                                line = line + n*(2^nIGcum[4])*(Naggr^nAGcum[1]) + o*(2^nIGcum[4])*(Naggr^nAGcum[2]) + p*(2^nIGcum[4])*(Naggr^nAGcum[3])
                                pathoToAggr[line,1] = i+1
                                pathoToAggr[line,2] = j+1
                                pathoToAggr[line,3] = k+1
                                pathoToAggr[line,4] = l+1
                                pathoToAggr[line,5] = m+1
                                pathoToAggr[line,6] = n+1
                                pathoToAggr[line,7] = o+1
                                pathoToAggr[line,8] = p+1
                            }
                        }
                    }
                }
            }
        }
    }
}

## _1mut: time to 1st appearance of an adapted mutant
## _1inf: time to 1st infection of an adapted mutant on R hosts
## _esta: time to establishment (50000 infection of R hosts by adapted mutants)
## col 1:4 for major genes
## col 5:8 for QR traits
## col 9 for the superpathogen
## lines: levels of resistance erodion
nLevErosion <- Naggr-1
D_1mut <- matrix(NA, nrow=nLevErosion, ncol=9)
colnames(D_1mut) <- c(traits, "super")
rownames(D_1mut) <- 1:nLevErosion
D_1inf <- D_1mut
D_esta <- D_1inf

if (length(Rhosts)>0 & Npatho>1) {
    for (i in 1:nLevErosion) {
        for (j in 1:9){
            if (j<=4) ## qualitative resistance
                indexPatho <- pathoToAggr[,j] > 1
            if (j>4 & j<9)  ## quantitative resistance
                indexPatho <- pathoToAggr[,j] > i
            if (j==9) ## superpathogen
                indexPatho <- Npatho
            
            if (sum(indexPatho>0)) {  ## i.e. the trait is evolving
                S_indexPatho <- apply(matrix(S_patho[indexPatho,], ncol=nTS), 2, sum)  ## conversion in matrix to avoid problem if only 1 row
                IL_indexPatho <- apply(matrix(IL_patho_Rhost[indexPatho,], ncol=nTS), 2, sum)
                D_1mut[i,j] <- min(which(S_indexPatho>0), nTS+1, na.rm=TRUE)
                D_1inf[i,j] <- min(which(IL_indexPatho>0), nTS+1, na.rm=TRUE)
                D_esta[i,j] <- min(which(IL_indexPatho>th_break), nTS+1, na.rm=TRUE)
            }
        } ## for j
    } ## for i
} ## if Rhosts
D_1mut
D_1inf
D_esta

# rotation_host <- numeric()
# for (y in 1:nYears)
#      rotation_host[y] <- max(habitat[,rotation[y]])


####       C) Epidemiological control       ####
## ST=short-term control (on S): integration from the beginning of the simulation until the first bound (e.g. 1st breakdown)
## (=NA if 1st breakdown before 1 year)
## TP=transitory period (on S): integration from first to 2nd bound
## (=NA if no complete seasons between bounds)
## LT=long-term control (on whole landscape): integration from 2nd bound (e.g. last breakdown) until the end
## (=NA if last bound higher than nYears-1 ~complete durability of at least one source)
## TOT=global (on whole landscape): integration from 0 to the end

GLA <- rep(NA,4)
names(GLA) <- c("glaST","glaTP","glaLT","glaTOT")
AUDPC <- rep(NA,4)
names(AUDPC) <- c("audpcST","audpcTP","audpcLT","audpcTOT")
indexHost <- list(1,1,1:Nhost,1:Nhost)
## Criterion to define the bounds
# durab <- D_esta[nLevErosion,1:8]
bounds <- range(D_esta[,1:8], na.rm=TRUE)   ## time to start breakdown/erosion -- time to complete breakdown/erosion
indexYears <- data.frame(rbind(ST=c(1, floor(bounds[1]/nTSpY))
                               , TP=c(ceiling(bounds[1]/nTSpY)+1, floor(bounds[2]/nTSpY))
                               , LT=c(ceiling(bounds[2]/nTSpY)+1, nYears)))
colnames(indexYears) <- c("y0","yf")
indexYears$y0[!is.finite(indexYears$y0)] <- nYears+1  ## avoid calculation if infinite bound
indexYears$yf[!is.finite(indexYears$yf)] <- 0
for (i in 1:3){
    if (indexYears$y0[i] <= indexYears$yf[i]) {
        GLA[i] <- sum(GLAabs[indexHost[[i]],indexYears$y0[i]:indexYears$yf[i]]) / (indexYears$yf[i]-indexYears$y0[i]+1)
        AUDPC[i] <- sum(audpc[indexHost[[i]],indexYears$y0[i]:indexYears$yf[i]]) / (indexYears$yf[i]-indexYears$y0[i]+1)
        if (i==3)
            AUDPC[3] <- mean(audpc_year[indexYears$y0[i]:indexYears$yf[i]]) 
    } 
} ## for i
AUDPC["audpcTOT"] <- mean(audpc_year)  ## Average AUDPC of the simulation of the whole landscape
GLA["glaTOT"] <- sum(GLAabs)/nYears
AUDPC
GLA



####       D) Economic performance       ####
benefit <- sellPrice  %*% product     ## yearly benefit from product sales
cost    <- purchPrice %*% C_host   ## yearly cost from crop planting
balance <- benefit - cost   ## yearly economic balance



####  Writing the output  ####
write.table(matrix(EXT,nrow=1),file=paste(pathRES,"/extinction.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(matrix(D_1mut,nrow=1),file=paste(pathRES,"/dur_1mut.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(matrix(D_1inf,nrow=1),file=paste(pathRES,"/dur_1inf.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(matrix(D_esta,nrow=1),file=paste(pathRES,"/dur_esta.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(matrix(GLA,nrow=1),file=paste(pathRES,"/GLA.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(matrix(AUDPC,nrow=1),file=paste(pathRES,"/AUDPC.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(matrix(balance,nrow=1),file=paste(pathRES,"/grossMargin.txt",sep=""),sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)


#############################################
#####               GRAPHICS             ####
#############################################
if (graphic) {
    graphics.off()
    nCol <- 21
    grad_red <- colorRampPalette(c("white","#FF5555"))
    RED <- grad_red(nCol)
    grad_blue <- colorRampPalette(c("white","#4F94CD"))
    BLUE <- grad_blue(nCol)
    grad_green <- colorRampPalette(c("white","darkolivegreen4"))
    GREEN <- grad_green(nCol)
    grad_grey <- colorRampPalette(c("white", "gray30"))
    GRAY <- grad_grey(Nhost)
    GRAY[1] <- "black"
    
    grad_whiteYellowRed <- colorRampPalette(c("white","#FFFF99","#990000"))
    grad_greenYellowRed <- colorRampPalette(c("#CCFF99","#FFFF99","#990000"))
    COL1 <- grad_whiteYellowRed(nCol)
    COL2 <- grad_greenYellowRed(nCol)
    
    
    legend.hab <- c("Susceptible","Resistant")
    if ((strat=="MO" | strat=="RO") & Nhost>2){
        # legend.hab <- c("Susceptible","Resistant 1","Resistant 2")
        legend.hab <- "Susceptible"
        for (i in 1:length(Rhosts))
            legend.hab <- c(legend.hab, paste("Resistant",i))
    }
    if (strat=="MI")
        # legend.hab <- c("Susceptible","Resistant 1 + Resistant 2")
        legend.hab <- c("Susceptible","Resistant 1","Resistant 2")
    if (strat=="PY")
        legend.hab <- c("Susceptible","Resistant 1+2")
    if (Nhost==1)
        legend.hab <- c("Susceptible")
    
    
     ####
     #### A) Host - Pathogen dynamic
     ####
     ####  1) Host dynamic  ####
     keyDates <- data.frame(start=D_esta[1,1:8], full=c(rep(NA,4),D_esta[nLevErosion,5:8]))  ## dataframe of durabilities
     keyDates$start[keyDates$start==nTS+1] <- NA
     keyDates$full[keyDates$full==nTS+1] <- NA
     keyDates$res0 <- res_vect[1,]
     keyDates$res1 <- res_vect[2,]
     keyDates$res2 <- res_vect[3,]
     keyDates$COL <- rep(c(BLUE[nCol],RED[nCol]), each=4)
     keyDates$LTY <- 1*keyDates$res0 + 2*keyDates$res1 + 3*keyDates$res2
     COL.tot <- "gray50"
     COL.periods <- c(GREEN[nCol], "gray50", RED[nCol])
     COL.periods2 <- c(GREEN[nCol%/%3], "gray80", RED[nCol%/%3])
     
     if (nYears==1) {
         ## Dynamic of the host (H, L, I, R)
         tiff(filename=paste(pathRES,"/DYNhost.tiff",sep=""),width=180,height=110,units='mm',compression='lzw',res=300)
         par(xpd=NA, mar=c(4,4,0,9))
         plot(0,0, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(1,nTS), ylim=c(0,1), ylab="", xlab="")
          hab <- habitat[,rotation[1]]
          for (host in 1:Nhost){
               lines(1:nTS, H_host[host,]/K_host[host,], col=GREEN[nCol], lty=host)
               lines(1:nTS, L_host[host,]/K_host[host,], col=BLUE[nCol], lty=host)
               lines(1:nTS, I_host[host,]/K_host[host,], col=RED[nCol], lty=host)
               lines(1:nTS, R_host[host,]/K_host[host,], col="black", lty=host)
          }
          axis(1, at=round(seq(1,nTS,length.out=8)))
          title(xlab="Time (days)", ylab="Proportion of hosts relative to carrying capacity")
          legend(nTS*1.05, .8, bty="n", title="Status:", legend=c("H","L","I","R"), border=NA, fill=c(GREEN[nCol],BLUE[nCol],RED[nCol],"black"))
          legend(nTS*1.05, .3, title.adj=-.01, bty="n", title="Cultivar:", legend=legend.hab, lty=1:Nhost, lwd=2, col="black")
          axis(2, at=seq(0,1,.2), las=1)
          dev.off()
     } else {
         ## Dynamic of Healthy hosts (H)
         tiff(filename=paste(pathRES,"/DYNhost_H.tiff",sep=""),width=180,height=110,units='mm',compression='lzw',res=300)
         par(xpd=NA, mar=c(4,4,0,9))
         plot(0,0, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(1,nTS), ylim=c(0,1), ylab="", xlab="")
         # for (i in 1:3){
         #     if (indexYears$y0[i] <= indexYears$yf[i])
         #         polygon(rep(c((indexYears$y0[i]-1)*nTSpY,indexYears$yf[i]*nTSpY),each=2),par("usr")[c(3,4,4,3)],col=COL.periods2[i], border=NA)
         # }
          for (host in 1:Nhost)
               lines(1:nTS, H_host[host,]/K_host[host,], col=GRAY[host], lty=host)
         for (trait in 1:nrow(keyDates)) {
             abline(v=keyDates$start[trait], col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
             abline(v=keyDates$full[trait], col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
          }
          axis(1, at=seq(1,nTS+1,nTSpY*((nYears-1)%/%8+1)), labels=seq(0,nYears,((nYears-1)%/%8+1)))
          title(xlab="Time (years)", ylab="Proportion of healthy hosts: H/K")
          legend(nTS*1.05, .5, title.adj=-.01, bty="n", title="Cultivar:", legend=legend.hab, lty=1:Nhost, lwd=2, col=GRAY)
          axis(2, at=seq(0,1,.2), las=1)
          for (i in 1:3){
              if (indexYears$y0[i] <= indexYears$yf[i])
                  segments((indexYears$y0[i]-1)*nTSpY, par("usr")[3], indexYears$yf[i]*nTSpY, par("usr")[3], col=COL.periods[i], lwd=12, xpd=FALSE)
          }
          dev.off()
          
          ## Dynamic of Diseased hosts (I+R)
          tiff(filename=paste(pathRES,"/DYNhost_IR.tiff",sep=""),width=180,height=110,units='mm',compression='lzw',res=300)
          par(xpd=NA, mar=c(4,4,0,9))
          plot(0,0, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(1,nTS), ylim=c(0,1), ylab="", xlab="")
          # for (i in 1:3){
          #     if (indexYears$y0[i] <= indexYears$yf[i])
          #         polygon(rep(c((indexYears$y0[i]-1)*nTSpY,indexYears$yf[i]*nTSpY),each=2),par("usr")[c(3,4,4,3)],col=COL.periods2[i], border=NA)
          # }
          for (host in 1:Nhost)
              lines(1:nTS, (I_host[host,]+R_host[host,])/K_host[host,], col=GRAY[host], lty=host)
          for (trait in 1:nrow(keyDates)) {
              abline(v=keyDates$start[trait], col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
              abline(v=keyDates$full[trait], col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
          }
          axis(1, at=seq(1,nTS+1,nTSpY*((nYears-1)%/%8+1)), labels=seq(0,nYears,((nYears-1)%/%8+1)))
          title(xlab="Time (years)", ylab="Proportion of diseased hosts: (I+R)/K")
          for (i in 1:3){
              if (indexYears$y0[i] <= indexYears$yf[i])
                  segments((indexYears$y0[i]-1)*nTSpY, par("usr")[3], indexYears$yf[i]*nTSpY, par("usr")[3], col=COL.periods[i], lwd=12, xpd=FALSE)
          }
          legend(nTS*1.05, .5, title.adj=-.01, bty="n", title="Cultivar:", legend=legend.hab, lty=1:Nhost, lwd=2, col=GRAY)
          axis(2, at=seq(0,1,.2), las=1)
          dev.off()
     }


     ####  2) Pathogen dynamic and evolution  ####
     ## Graphic with grey scale
     
     
     D_freq <- rep(NA,9)  ## Durability relative to pathogen frequencies
     names(D_freq) <- c(traits,"super")
     I_aggrProp <- list()
     for (trait in (1:8)[as.logical(adapt_vect)]){
          nIncr <- 2*(trait < 5) + Naggr*(trait>=5)  ## number of increments (2 for GFG, Naggr for QR)
          I_aggrProp[[trait]] <- matrix(0, nrow=nIncr, ncol=nTS)
          for (t in 1:nTS) {
               for (p in 1:Npatho) {
                    aggr <- pathoToAggr[p,trait]
                    I_aggrProp[[trait]][aggr, t] <- I_aggrProp[[trait]][aggr, t] + I_pathoProp[p, t]
               }
          }
          D_freq[trait] <- which(I_aggrProp[[trait]][1,]<0.95)[1] - 1   ## last time-step the 1st pathogen is above 0.95
          if (graphic)
               plotevolQR(pathRES,nIncr,trait,I_aggrProp[[trait]],D_freq[trait],nTS,nYears,nTSpY)   
     }
     D_freq[9] <- which(I_pathoProp[Npatho,]>0.05)[1] - 1   ## durability relative to the superpathogen
     
     ## Graphic with curves (for qualitative resistance)
     if (sum(adapt_vect[1:4])>0) {
         COL <- c("#FF5555","#4F94CD","darkolivegreen4","#CD950C","black")  ## colors: red, blue, green, orange, black
         tiff(filename=paste(pathRES,"/EVOLpatho_curve.tiff",sep=""),width=180,height=110,units='mm',compression='lzw',res=300)
         par(xpd=FALSE, mar=c(4,4,0,9))
         idIG <- (1:4)[as.logical(adapt_vect[1:4])]
         plot(0,0, type="n", xlim=c(1, nTS), bty="n", las=1, ylim=c(0,1), xaxt="n", xlab="", ylab="Frequency of infective genotypes")
         for (trait in idIG) {
             lines(I_aggrProp[[trait]][2,], col=COL[trait], lwd=1.5)
             abline(v=D_freq[trait], col=COL[trait], lty=trait+1, lwd=1)  ## durability of the trait
         }
         if (length(idIG)>1) {
             lines(I_pathoProp[Npatho,], lty=2, col=COL[5], lwd=1.5)
             abline(v=D_freq[9], col=COL[5], lty=3, lwd=1)  ## durability relative to the superpathogen
         }
         if (nYears==1) {
             axis(1, at=round(seq(1,nTS,length.out=8)), las=1)
             title(xlab="Evolutionnary time (days)")
         } else {
             axis(1, at=seq(1,nTS+1,nTSpY*((nYears-1)%/%8+1)), labels=seq(0,nYears,((nYears-1)%/%8+1)), las=1)
             title(xlab="Evolutionnary time (years)")
         }
         par(xpd=TRUE)
         if (length(idIG)>1) {
             legend(nTS*1.05, .5, title.adj=-.01, bty="n", title="Infectivity gene:", legend=c(idIG,paste(idIG, collapse="+")), lty=c(rep(1,length(idIG)),2), lwd=2, col=COL[c(idIG,5)])
         } else {
             legend(nTS*1.05, .5, title.adj=-.01, bty="n", title="Infectivity gene:", legend=c(idIG), lty=c(rep(1,length(idIG))), lwd=2, col=COL[c(idIG)])
         }
         dev.off()
     } ## if there are some MG

     
     ####
     ####  B) GLA & AUDPC  ####
     ####
     GLAnoDis <- 1.48   ## GLA in absence of disease
     audpc100S <- 0.38  ## AUDPC in a 100% susceptible landscape
     GLAabs_year <- apply(GLAabs, 2, sum)    ## Total GLA of the landscape (per year)
     GLAabs_host <- apply(GLAabs, 1, mean)   ## Average GLA of the simulation (per cultivar)
     
     PCH <- 15:(15+Nhost-1)
     PCH.tot <- 0
     LTY.tot <- 1
     COL.tot <- RED[nCol] #"gray50"
     legend.GLA <- c(legend.hab,"Whole landscape")
     if (strat=="MI")  ## re-adjust legend of maps for mixtures
         legend.hab <- c("Susceptible","Resistant 1 + Resistant 2")
     
     if (nYears>1) {
         ## GLA
         tiff(filename=paste(pathRES, "/GLAabs.tiff",sep=""), width=180, height=110, units='mm', compression='lzw', res=300)
         # m <- matrix(c(rep(1,5),2),6,1)    # matrix(c(rep(1,5),3,rep(2,5),3),6,2)
         # layout(m)
         # par(xpd=F, cex=.9,mar=c(5,4.1,4.1,2.1))
         par(xpd=NA, cex=.9, mar=c(4,4,0,9))
         plot(0,0, type="n", bty="n", xlim=c(1, nYears), ylim=c(0, GLAnoDis), xaxt="n", yaxt="n"
              # , main="Absolute green leaf area (GLAa)"
              , xlab="Years", ylab="Nb of healthy hosts / day / m^2") 
         axis(1, at=seq(1,nYears,by=(nYears-1)%/%8+1))
         axis(2, at=round(seq(0, GLAnoDis,length.out=5),2), las=1)
         for (host in 1:Nhost) {
             lines(1:nYears, GLAabs[host,], lwd=2, lty=host, col=GRAY[host])
             points(1:nYears, GLAabs[host,], lwd=1, pch=PCH[host], col=GRAY[host])  ## need to add the points for rotations
             points(par("usr")[1], GLAabs_host[host], pch=PCH[host], col=GRAY[host], xpd=TRUE)
         }
         lines(1:nYears, GLAabs_year, lwd=2, lty=LTY.tot, col=COL.tot)
         points(1:nYears, GLAabs_year, lwd=1, pch=PCH.tot, col=COL.tot)
         points(par("usr")[1], GLA["glaTOT"], pch=PCH.tot, col=COL.tot, cex=1, xpd=TRUE)
         for (trait in 1:nrow(keyDates)) {
             abline(v=ceiling(keyDates$start[trait]/nTSpY), col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
             abline(v=ceiling(keyDates$full[trait]/nTSpY), col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
         }
         for (i in 1:3){
             if (indexYears$y0[i] <= indexYears$yf[i])
                 segments(indexYears$y0[i]-1/4, par("usr")[3], indexYears$yf[i]+1/4, par("usr")[3], col=COL.periods[i], lwd=12, xpd=FALSE)
         }
          # legend(nYears*.5,max(GLAabs)*.9, title.adj=0.1, bty="n", title="Cultivar:", legend=c("Resistant","Susceptible"), lty=1:2, lwd=2, col="black")
         legend(nYears*1.05, .5*GLAnoDis, title.adj=-.01, bty="n", title="Cultivar:", legend=legend.GLA, cex=0.9
                , lty=c(1:Nhost,LTY.tot), lwd=2, pch=c(PCH,PCH.tot), pt.cex=c(rep(1,Nhost),1), col=c(GRAY, COL.tot), seg.len=2.5)
         
         # par(xpd=NA, mar=c(0,0,0,0))
         # plot(0,0, type="n", axes=F) 
         # legend("top", bty="n", title="Cultivar:", legend=legend.GLA
         #        , lty=c(1:Nhost,4), lwd=2, pch=PCH[c(1:Nhost,4)], pt.cex=c(rep(1,Nhost),1.5), col=c(rep("black",Nhost), COL.tot), seg.len=2.5)
         # 
         dev.off()

         tiff(filename=paste(pathRES, "/GLArel.tiff",sep=""), width=180, height=110, units='mm', compression='lzw', res=300)
         par(xpd=NA, cex=.9,mar=c(4,4,0,9))
         plot(0, 0, type="l", bty="n", lwd=2, xlim=c(1,nYears), ylim=c(0,1), xaxt="n", yaxt="n"
              # , main="Relative green leaf area (GLAr)"
              , xlab="Years", ylab="Proportion of healthy hosts: H/N") 
         for (host in 1:Nhost) {
             lines(1:nYears, GLArel[host,], lwd=2, lty=host, col=GRAY[host])
             points(1:nYears, GLArel[host,], lwd=1, pch=PCH[host], col=GRAY[host])
         }
         for (trait in 1:nrow(keyDates)) {
             abline(v=ceiling(keyDates$start[trait]/nTSpY), col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
             abline(v=ceiling(keyDates$full[trait]/nTSpY), col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
         }
         axis(1, at=seq(1,nYears,by=(nYears-1)%/%8+1))
         axis(2, at=seq(0,1,.2), las=1)
         legend(nYears*1.05, .5, title.adj=-.01, bty="n", title="Cultivar:", legend=legend.GLA[1:(length(legend.GLA)-1)], cex=0.9
                , lty=1:Nhost, lwd=2, pch=PCH, pt.cex=rep(1,Nhost), col=GRAY, seg.len=2.5)
         dev.off()
          
          ## AUDPC
          tiff(filename=paste(pathRES, "/AUDPC.tiff",sep=""), width=180, height=110, units='mm', compression='lzw', res=300)
          par(xpd=NA, cex=.9,mar=c(4,4,0,9))
          plot(0,0, type="n", bty="n", xlim=c(1, nYears), ylim=c(0, audpc100S), xaxt="n", yaxt="n"
               # , main="AUDPC"
               , xlab="Years", ylab="Proportion of diseased hosts: (I+R)/K") 
          axis(1, at=seq(1,nYears,by=(nYears-1)%/%8+1))
          axis(2, at=round(seq(0, audpc100S,length.out=5),2), las=1)
          for (host in 1:Nhost) {
              lines(1:nYears, audpc[host,], lwd=2, lty=host, col=GRAY[host])
              points(1:nYears, audpc[host,], lwd=1, pch=PCH[host], col=GRAY[host])  ## need to add the points for rotations
              points(par("usr")[1], audpc_host[host], pch=PCH[host], col=GRAY[host], xpd=TRUE)
          }
          for (trait in 1:nrow(keyDates)) {
              abline(v=ceiling(keyDates$start[trait]/nTSpY), col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
              abline(v=ceiling(keyDates$full[trait]/nTSpY), col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=2.6, xpd=FALSE)
          }
          lines(1:nYears, audpc_year, lwd=2, lty=LTY.tot, col=COL.tot)
          points(1:nYears, audpc_year, lwd=1, pch=PCH.tot, col=COL.tot)
          points(par("usr")[1], AUDPC["audpcTOT"], pch=PCH.tot, col=COL.tot, cex=1, xpd=TRUE)
          for (i in 1:3){
              if (indexYears$y0[i] <= indexYears$yf[i])
                  segments(indexYears$y0[i]-1/4, par("usr")[3], indexYears$yf[i]+1/4, par("usr")[3], col=COL.periods[i], lwd=12, xpd=FALSE)
          }
          legend(nYears*1.05, .5*audpc100S, title.adj=-.01, bty="n", title="Cultivar:", legend=legend.GLA, cex=0.9
                 , lty=c(1:Nhost,LTY.tot), lwd=2, pch=c(PCH,PCH.tot), pt.cex=c(rep(1,Nhost),1), col=c(GRAY, COL.tot), seg.len=2.5)
          dev.off()
          
          ## Production
          productTot <- apply(product,2,sum)
          maxProd <- max(productTot)
          tiff(filename=paste(pathRES, "/production.tiff",sep=""), width=180, height=110, units='mm', compression='lzw', res=300)
          par(xpd=NA, cex=.9,mar=c(4,4,0,9))
          plot(0,0, type="n", bty="n", xlim=c(1, nYears), ylim=c(0, maxProd), xaxt="n", yaxt="n"
               , xlab="Years", ylab="Production (weight units)") 
          axis(1, at=seq(1,nYears,by=(nYears-1)%/%8+1))
          axis(2, at=round(seq(0, maxProd,length.out=5)), las=1)
          for (host in 1:Nhost) {
              lines(1:nYears, product[host,], lwd=2, lty=host, col=GRAY[host])
              points(1:nYears, product[host,], lwd=1, pch=PCH[host], col=GRAY[host])  ## need to add the points for rotations
              points(par("usr")[1], mean(product[host,]), pch=PCH[host], col=GRAY[host], xpd=TRUE)
          }
          lines(1:nYears, productTot, lwd=2, lty=LTY.tot, col=COL.tot)
          points(1:nYears, productTot, lwd=1, pch=PCH.tot, col=COL.tot)
          points(par("usr")[1], mean(productTot), pch=PCH.tot, col=COL.tot, cex=1, xpd=TRUE)
          legend(nYears*1.05, .5*maxProd, title.adj=-.01, bty="n", title="Cultivar:", legend=legend.GLA, cex=0.9
                 , lty=c(1:Nhost,LTY.tot), lwd=2, pch=c(PCH,PCH.tot), pt.cex=c(rep(1,Nhost),1), col=c(GRAY, COL.tot), seg.len=2.5)
          dev.off()
          
          ## Economic performance
          maxEco <- max(benefit, cost, balance)
          minEco <- min(0, benefit, cost, balance)
          tiff(filename=paste(pathRES, "/grossMargin.tiff",sep=""), width=180, height=110, units='mm', compression='lzw', res=300)
          par(xpd=NA, cex=.9,mar=c(4,4,0,9))
          plot(0,0, type="n", bty="n", xlim=c(1, nYears), ylim=c(minEco, maxEco), xaxt="n", yaxt="n"
               , xlab="Years", ylab="Monetary units") 
          axis(1, at=seq(1,nYears,by=(nYears-1)%/%8+1))
          axis(2, at=round(seq(minEco, maxEco,length.out=5)), las=1)
          if (minEco<0)
              abline(h=0, lwd=1, lty=2, col="gray", xpd=FALSE)
          
          lines(1:nYears, benefit, lwd=2, lty=1, col="black")
          points(1:nYears, benefit, lwd=1, pch=PCH[1], col="black")  ## need to add the points for rotations
          points(par("usr")[1], mean(benefit), pch=PCH[1], col="black", xpd=TRUE)
          lines(1:nYears, cost, lwd=2, lty=2, col="gray50")
          points(1:nYears, cost, lwd=1, pch=PCH[1]+1, col="gray50")  ## need to add the points for rotations
          points(par("usr")[1], mean(cost), pch=PCH[1]+1, col="gray50", xpd=TRUE)
          
          lines(1:nYears, balance, lwd=2, lty=LTY.tot, col=COL.tot)
          points(1:nYears, balance, lwd=1, pch=PCH.tot, col=COL.tot)
          points(par("usr")[1], mean(balance), pch=PCH.tot, col=COL.tot, cex=1, xpd=TRUE)
          
          legend(nYears*1.05, .5*maxEco, title.adj=-.01, bty="n", legend=c("Benefit","Cost","Gross Margin"), cex=0.9
                 , lty=c(1:2,LTY.tot), lwd=2, pch=c(PCH[1],PCH[1]+1,PCH.tot), pt.cex=rep(1,3), col=c("black", "gray50", COL.tot), seg.len=2.5)
          dev.off()
          
     }
     
     ####
     ####  C) VIDEO (MAPS of HLIR dynamic)  ####
     ####
     if (video){
         
         if( system("ffmpeg -version", ignore.stdout=TRUE) == 127) 
             stop("You need to install ffmpeg before generating videos")
         
         
     graphics.off()
     # library("sf")
     # library(animation)
     dir.create("maps")
     
     # title.H <- "Proportion of healty hosts"
     title.IR <- "Proportion of diseased hosts"
     invlogitbound <- 3
     maxColorScale <- 1
     # intvls <- maxColorScale * sort(c(0, invlogit(seq(-invlogitbound,invlogitbound,l=nCol-2)), 1), decreasing=FALSE)   ## intervals to define coloration
     intvls <- maxColorScale * sort(seq(0, 1, l=nCol), decreasing=FALSE)
     # colMap <- heat.colors(nCol, alpha=.3)    ## color for each interval (alpha = transparency)
     # plot(1:nCol~rep(0,nCol),col=colMap, pch=15, cex=3)
     dens.H <- c(0,8,15,12,18)[1:Nhost]
     angle.H <- c(0,30,120,75,165)[1:Nhost]
     legend.H <- sprintf("%.3f",intvls)
     legend.H <- paste(intvls*100,"%", sep="")
 
     for (y in 1:nYears) {
         print(paste("Year",y,"/",nYears, " -  habitat",rotation[y]))
         hab <- habitat[,rotation[y]]#habitat[[rotation[y]]][,y]
         K_poly <- rep(0,Npoly)
         for (host in 1:Nhost) {
             K_poly <- K_poly + (K[,host]*(hab==host))
             if (strat=="MI" & host>1)
                 K_poly <- K_poly + (K[,host]*(habitat[,2]==host))
         }
         for (d in round(seq(1,nTSpY,length.out=nMapPY))) {
             subtitle.H <- paste("Year =", y, "   Day =", d)
             ## Proportion of each type of host relative to carrying capacity
             ts <- (y-1)*nTSpY + d
             # N_poly <- apply(H[[ts]],1,sum) + apply(L[[ts]],3,sum) + apply(I[[ts]],3,sum) + apply(R[[ts]],3,sum)
             
             propH <- apply(H[[ts]],1,sum) / K_poly
             # propL <- apply(L[[ts]],3,sum) / K_poly
             propI <- apply(I[[ts]],3,sum) / K_poly
             propR <- apply(R[[ts]],3,sum) / K_poly
             intvlsH <- findInterval(propH, intvls)
             
             # propIR <- (apply(I[[ts]],3,sum) + apply(R[[ts]],3,sum)) / N_poly
             propIR <- (propI + propR)
             intvlsIR <- findInterval(propIR, intvls)

             png(filename=paste(pathRES, "/maps/HLIR_",sprintf("%02d",y),"-",sprintf("%03d",d),".png",sep=""),width=700,height=350,units="mm",res=72)
             par(mfrow=c(1,2), cex=2, xpd=NA, mar=c(9.5,5,4,2))
             ## Map dynamic of healthy hosts
             # plot.land(landscape, GREEN[intvlsH], dens.H[hab], angle.H[hab], GREEN, dens.H, angle.H, title.H, subtitle.H, legend.hab, legend.H, "H/K")
             
             ## moving AUDPC
             plot(0,0, type="n", bty="n", xlim=c(1, nYears), ylim=c(0,1), xaxt="n", yaxt="n"
                  , xlab="", ylab="", main="Disease severity relative to a fully susceptible landscape") 
             axis(1, at=seq(1,nYears,by=(nYears-1)%/%8+1))
             ticksMarks <- round(seq(0,1,length.out=5),2)
             axis(2, at=ticksMarks, labels=paste(100*ticksMarks, "%"), las=1)
             title(xlab="Years", mgp=c(2,1,0))
             for (host in 1:Nhost) {
                 lines(1:y, audpc[host,1:y]/audpc100S, lwd=2, lty=host, col=GRAY[host])
                 points(1:y, audpc[host,1:y]/audpc100S, lwd=1, pch=PCH[host], col=GRAY[host])  ## need to add the points for rotations
             }
             for (trait in 1:nrow(keyDates)) {
                 date1 <- ceiling(keyDates$start[trait]/nTSpY)
                 date2 <- ceiling(keyDates$full[trait]/nTSpY)
                 if (!is.na(date1) & date1 <= y)
                     abline(v=date1, col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=4, xpd=FALSE)
                 if (!is.na(date2) & date2 <= y)
                     abline(v=date2, col=keyDates$COL[trait], lty=keyDates$LTY[trait], lwd=4, xpd=FALSE)
             }
             lines(1:y, audpc_year[1:y]/audpc100S, lwd=2, lty=LTY.tot, col=COL.tot)
             points(1:y, audpc_year[1:y]/audpc100S, lwd=1, pch=PCH.tot, col=COL.tot)
             legend(nYears/3, -1/5, bty="n", legend=legend.GLA, cex=1
                    , lty=c(1:Nhost,LTY.tot), lwd=2, pch=c(PCH,PCH.tot), pt.cex=c(rep(1,Nhost),1), col=c(GRAY, COL.tot), seg.len=2.5)
             
             ## Map dynamic of diseased hosts
             plotland(trans.landscape_tmp, COL1[intvlsIR], dens.H[hab], angle.H[hab], COL1, dens.H, angle.H, title.IR, subtitle.H, legend.hab, legend.H, " ")
             dev.off()
         } ## for d
     }  ## for y
     
     ## Convert png files to mp4 video
     setwd("maps")
     # title_video <- paste("simul_",strat,Nhost, "_prop", round(propSR,2), "_agreg", isolSR,"_rep",seed,".mp4", sep="")
     title_video <- paste("video_",strat,Nhost,".mp4", sep="")
     command_line <- paste("ffmpeg -f image2 -framerate ", nMapPY/2, " -pattern_type glob -i \"*.png\" ", title_video,sep="")
     system(command_line)

     system(paste("mv ", title_video, " ", "../."))
     setwd("../")
     print("remove directory maps")
     system("rm -rf maps")
     
     ### Not used  (gif too big)
     # convert the .png files to one .gif file using ImageMagick. 
     # The system() function executes the command as if it was done
     # in the terminal. the -delay flag sets the time between showing
     # the frames, i.e. the speed of the animation.
     
     #system("convert -delay 10 *.png video_simul.gif")
     #saveGIF({xxxxx},movie.name="test.gif",interval=0.2) 
     #system("ffmpeg -framerate 1 -i *.png -r 10 -c:v libx264 -crf 15 -pix_fmt yuv420p out.mp4")  
 
     
     } ## if video
     
} ## if graphic


print("remove binary files")
file.remove(list.files(pattern="\\.bin$"))   # delete all the binary files
print(paste("model outputs stored in :", pathRES))

}

