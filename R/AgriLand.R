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


#' @title Landscape generation
#' @name AgriLand
#' @description Generate a landscape composed of fields where a susceptible (SC) and one (RC) or two (RC1 and RC2) resistant cultivars 
#' are allocated with controlled proportions and spatio-temporal aggregation.
#' @param landscape a spatialpolygon object containing field coordinates.
#' @param filename a character string specifying the output layer name.
#' @param propSR proportion of fields where resistance is deployed: (RC)/(SC+RC) or (RC1+RC2)/(SC+RC1+RC2). Must be between 0 and 1.
#' @param isolSR an integer giving the spatial aggregation of fields where resistance is deployed (1=highly fragmented, 2=balanced, 3=highly aggregated).
#' @param propRR when applicable (mixtures and mosaics only), relative proportion of the second resistant cultivar: (RC2)/(RC1+RC2). Must be between 0 and 1.
#' @param isolRR when applicable, an integer giving the spatial (for mosaics) or temporal (for rotations) aggregation of fields 
#' cultivated with the second resistant cultivar (1=highly fragmented, 2=balanced, 3=highly aggregated).
#' @param strat a character string specifying the deployment strategy ("MO"=mosaic, "MI"=mixture, "RO"=rotations, "PY"=pyramiding).
#' @param Nhote an integer giving the number of cultivars (1, 2 or 3).
#' @param nYears an integer giving the number of simulated years.
#' @param Cmax0 carrying capacity of the susceptible cultivar in number of hosts per meter square.
#' @param Cmax1 carrying capacity of the resistant cultivars in number of hosts per meter square.
#' @param seed an integer giving the seed value (for random number generator).
#' @param graphOn a logical indicating if a graph of the landscape must be generated (1) or not (0).
#' @details An algorithm based on latent Gaussian fields is used to allocate two different crop cultivars across the simulated landscapes 
#' (e.g. a susceptible and a resistant cultivar, denoted as SC and RC, respectively). This algorithm allows the control of the proportions 
#' of each cultivar in terms of surface coverage, and their level of spatial aggregation. A random vector of values is drawn from a 
#' multivariate normal distribution with expectation 0 and a variance-covariance matrix which depends on the pairwise distances between 
#' the centroids of the fields. Next, the crop cultivars are allocated to different fields depending on whether each value drawn from the 
#' multivariate normal distribution is above or below a threshold. The proportion of each cultivar in the landscape is controlled by the value 
#' of this threshold. The sequential use of this algorithm allows the allocation of more than two crop cultivars (e.g. SC, RC1 and RC2). 
#' Therefore, deployment strategies involving two sources of resistance is simulated by: 
#' \enumerate{
#' \item running the allocation algorithm once to segregate the fields where the susceptible cultivar is grown, and 
#' \item applying one of the following deployment strategies to the remaining candidate fields:
#' \itemize{
#' \item Mosaics: two resistant cultivars (RC1 and RC2, carrying the first and the second resistance sources, respectively) are assigned to candidate 
#' fields by re-running the allocation algorithm;
#' \item Mixtures: both RC1 and RC2 are allocated to all candidate fields;
#' \item Rotations: RC1 and RC2 are alternatively cultivated in candidate fields, depending on the number of cropping 
#' seasons over which a given cultivar is grown before being rotated;
#' \item Pyramiding: all candidate fields are cultivated with RC12, a resistant cultivar carrying both resistance sources.
#' }
#' }
#' @return a shapefile containing the landscape structure (i.e. coordinates of field boundaries) and composition (i.e. cultivars).
#' @importFrom sf st_as_sf
#' @importFrom sf st_write
#' @importFrom grDevices dev.off graphics.off png tiff
#' @include multiN.R periodic_cov.R graphLand.R 
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (in press). Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology}.
#' @examples 
#' ## Generate a landscape consisting in a mosaic of fields cultivated with a susceptible cultivar
#' ## and two resistant cultivars in balanced proportions and high level of spatial aggregation
#' \dontrun{ 
#' landscapeTEST1
#' AgriLand(landscapeTEST1,filename="landscapeTEST1",propSR=2/3,isolSR=3,
#' propRR=1/2,isolRR=3,strat="MO",Nhote=3,nYears=30,Cmax0=2,Cmax1=2,seed=12345,graphOn=1)
#' }
#' @export
AgriLand <- function(landscape,filename="landscapeTEST1",propSR,isolSR,propRR,isolRR,strat,Nhote,nYears,Cmax0,Cmax1,seed,graphOn){
    set.seed(seed)
nPoly.tmp <- length(landscape)
prop <- c(propSR, propRR)
isol <- c(isolSR, isolRR)

## Parameters of cultivar allocation
nAlloc <- (Nhote>1)    ## basic number of allocations to perform
if (strat=="MO" & Nhote>2){nAlloc <- 2}

## Isolation/aggregation parameter
aggreg <- c(-180, 200, -2000, 0)
## selection of the appropriate aggregation parameter with isol
## isol = 1 --> low
## isol = 2 --> moderate
## isol = 3 --> high
## isol = 4 --> random

## Centroid of the paddocks
centroid <- NULL
area <- NULL
for (i in 1:nPoly.tmp){
    for (j in 1:length(landscape@polygons[[i]])){
        centroid <- rbind(centroid,apply(landscape@polygons[[i]]@Polygons[[j]]@coords,2,mean))
        area <- c(area,landscape@polygons[[i]]@Polygons[[j]]@area)
    }
}
nPoly <- nrow(centroid)
d <- as.matrix(dist(centroid))   ## 2-by-2 distance between centroid of each paddock
neigh <- (d <= 399)               ## matrix of neighborood

## Multivariate distribution
area.df <- data.frame(num=1:nPoly, area)
habitat <- data.frame(num=1:nPoly, cultivar=rep(0,nPoly))
i=0
num_index=1

while( sum(num_index)>0 & i<nAlloc ) {
     ## update area, d, qnd i for allocation of the cultivar
     area.tmp <- area.df[habitat$cultivar==i,]
     d.tmp <- d[habitat$cultivar==i, habitat$cultivar==i]
     i <- i+1
     ## Compute the multivariate normal distribution
     habitat.tmp <- multiN(d.tmp, area.tmp, aggreg[isol[i]], prop[i])
     ## update habitat by allocating cultivar i
     num_index <- habitat.tmp$cultivar==1
     habitat[habitat.tmp$num[num_index],"cultivar"] <- i
}

habitat1 <- habitat$cultivar
habitat2 <- habitat1

if (strat=="RO" | strat=="TO"| strat=="MI"){habitat2[habitat2==1] <- 2}

## Crop rotations & Turn-over: calculation of time series
rotation <- data.frame(y=1:(nYears+1), habitat=rep(0,nYears+1))  ## need to have nYears+1 because of C algorithm for host plantation
if (strat=="RO" | strat=="TO") {
    ## isolRR gives the number of years for each habitat
    rotation$habitat <- rep(c(1,0), each=isol[2], length=nYears+1)
}

#shape file for the landscape
landscapeINIT <- SpatialPolygonsDataFrame(landscape, data.frame(habitat1=habitat1,habitat2=habitat2,area=area), match.ID=T)
habitat <- st_as_sf(landscapeINIT)
results <- st_as_sf(landscape)
st_write(habitat, paste0(filename,".gpkg"), "habitat",layer_options="OVERWRITE=yes")
st_write(results, paste0(filename,".gpkg"), "results", update = TRUE,layer_options="OVERWRITE=yes")
  
## Graphic representing the landscape
if (graphOn) {
    title.hab <- "Simulated landscape"
    isol.name <- c("low", "medium", "high", "random")
    alloc.name <- c("R/(S+R)","R2/(R1+R2)")
    subtitle.hab <- paste("Cropping ratio",alloc.name,"=",round(prop,2), "   Aggregation",alloc.name,"=",isol.name[isol])
    if (strat!="MO" | Nhote==2)
        subtitle.hab <- subtitle.hab[1]
    col.hab <- c("white", "gray80","gray45")  ## habitat 0 = S ; habitat 1 = R1 ; habitat 2 = R2
    dens.hab <- c(0,0,0)
    angle.hab <- c(0,0,0)
    legend.hab <- c("Susceptible","Resistant")
    if ((strat=="MO" | strat=="RO" | strat=="TO") & Nhote>2)
        legend.hab <- c("Susceptible","Resistant 1","Resistant 2")
    if (strat=="MI")
        legend.hab <- c("Susceptible","Resistant 1 + Resistant 2")
    if (strat=="PY")
        legend.hab <- c("Susceptible","Resistant 1+2")
    png(filename="landscape.png",width=1000,height=1000)
    plotland(landscape, col.hab[habitat1+1], dens.hab[habitat1+1], angle.hab[habitat1+1], col.hab, dens.hab, angle.hab, title.hab, subtitle.hab, legend.hab)
    dev.off()
}
return(list(shapefilename=paste0(filename,".gpkg"),layername_hab="habitat",layername_res="results",rotation=as.integer(rotation$habitat),Cmax0=Cmax0,Cmax1=Cmax1,propRR=propRR,strat=strat))
}



