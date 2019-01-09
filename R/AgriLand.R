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
#' @param Nhost an integer giving the number of cultivars (1, 2 or 3).
#' @param nYears an integer giving the number of simulated years.
#' @param seed an integer giving the seed value (for random number generator).
#' @param graphic a logical indicating if a graph of the landscape must be generated (TRUE) or not (FALSE).
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
#' @return a shapefile containing the landscape structure (i.e. coordinates of field boundaries) and composition (i.e. cultivars) in time (i.e. each year) and space (i.e. each field).
#' @importFrom sf st_as_sf
#' @importFrom sf st_write
#' @importFrom grDevices dev.off graphics.off png tiff
#' @include multiN.R periodic_cov.R graphLand.R 
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018). Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#' @examples 
#' ## Generate a landscape consisting in a mosaic of fields cultivated with a susceptible cultivar
#' ## and two resistant cultivars in balanced proportions and high level of spatial aggregation
#' \dontrun{ 
#' landscapeTEST1
#' AgriLand(landscapeTEST1,filename="landscapeTEST1",propSR=2/3,isolSR=3,
#' propRR=1/2,isolRR=3,strat="MO",Nhost=3,nYears=30,seed=12345,graphic=TRUE)
#' }
#' @export

AgriLand <- function(landscape,filename="landscapeTEST1",propSR,isolSR,propRR,isolRR,strat,Nhost,nYears,seed,graphic=FALSE){
     set.seed(seed)
     nPoly.tmp <- length(landscape)
     prop <- c(propSR, propRR)
     aggreg <- c(isolSR, isolRR)
     
     ## Parameters of cultivar allocation
     nAlloc <- (Nhost>1)    ## basic number of allocations to perform
     if (strat=="MO" & Nhost>2){nAlloc <- Nhost - 1}
     
     ## aggregation parameter related to the range of the covariance matrix
     range <- c(-180, 200, -2000, 0)
     ## selection of the appropriate range parameter with aggreg
     ## aggreg = 1 --> low
     ## aggreg = 2 --> moderate
     ## aggreg = 3 --> high
     ## aggreg = 4 --> random
     
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
          ## update area, d, and i for allocation of the cultivar
          area.tmp <- area.df[habitat$cultivar==i,]
          d.tmp <- d[habitat$cultivar==i, habitat$cultivar==i]
          i <- i+1
          ## Compute the multivariate normal distribution
          habitat.tmp <- multiN(d.tmp, area.tmp, range[aggreg[i]], prop[i])
          ## update habitat by allocating cultivar i
          num_index <- habitat.tmp$cultivar==1
          habitat[habitat.tmp$num[num_index],"cultivar"] <- i
     }
     
     ## Writing the habitat files for C function
     habitat0 <- habitat$cultivar
     habitat1 <- habitat0
     if (strat=="MI" | strat=="RO") { habitat1[habitat1==1] <- 2 }
     
     ## Crop rotations: calculation of time series
     rotation <- data.frame(y=1:(nYears+1), habitat=rep(0,nYears+1))  ## need to have nYears+1 because of C algorithm for host plantation
     if (strat=="RO") {
          ## aggreg[2] gives the number of years for each habitat
          rotation$habitat <- rep(sample(c(1,0)), each=aggreg[2], length=nYears+1)
     }
     
     ## shape file for the landscape
     landscapeINIT <- SpatialPolygonsDataFrame(landscape, data.frame(habitat0=habitat0,habitat1=habitat1,area=area), match.ID=T)
     habitat <- st_as_sf(landscapeINIT)
     results <- st_as_sf(landscape)
     st_write(habitat, paste0(filename,".gpkg"), "habitat",layer_options="OVERWRITE=yes", driver= "GPKG")
     st_write(results, paste0(filename,".gpkg"), "results", update = TRUE,layer_options="OVERWRITE=yes",driver= "GPKG")
     
     ## Graphic representing the landscape
     if (graphic) {
          title.hab <- "Simulated landscape"
          aggreg.name <- c("low", "medium", "high", "random")
          # alloc.name <- c("R/(S+R)","R2/(R1+R2)")
          # subtitle.hab <- paste("Cropping ratio",alloc.name,"=",round(prop,2), "   Aggregation",alloc.name,"=",aggreg.name[aggreg])
          subtitle.hab <- paste("Cropping ratios = (", paste(round(prop[1:nAlloc],2),collapse="; "), ")"
                                , "   Aggregations = (", paste(aggreg.name[aggreg[1:nAlloc]],collapse="; "), ")", sep="")
          
          # if (strat!="MO" | Nhost==2)
          #      subtitle.hab <- subtitle.hab[1]
          # col.hab <- c("white", "gray80","gray45")  ## habitat 0 = S ; habitat 1 = R1 ; habitat 2 = R2
          # dens.hab <- c(0,0,0)
          # angle.hab <- c(0,0,0)
          colfunc <- colorRampPalette(c("white", "gray30"))
          col.hab <- colfunc(nAlloc+1+as.numeric(strat=="RO"))
          dens.hab <- rep(0,nAlloc+1+as.numeric(strat=="RO"))
          angle.hab <- rep(0,nAlloc+1+as.numeric(strat=="RO"))
          legend.hab <- c("Susceptible","Resistant")
          if ((strat=="MO" | strat=="RO") & Nhost>2){
               # legend.hab <- c("Susceptible","Resistant 1","Resistant 2")
               legend.hab <- "Susceptible"
               for (i in 1:(nAlloc+as.numeric(strat=="RO")))
                    legend.hab <- c(legend.hab, paste("Resistant",i))
          }
          if (strat=="MI")
               legend.hab <- c("Susceptible","Resistant 1 + Resistant 2")
          if (strat=="PY")
               legend.hab <- c("Susceptible","Resistant 1+2")
          png(filename="landscape.png", width=1000, height=1000)
          plotland(landscape, col.hab[habitat0+1], dens.hab[habitat0+1], angle.hab[habitat0+1], col.hab, dens.hab, angle.hab, title.hab, subtitle.hab, legend.hab)
          dev.off()
     }
     return(list(shapefilename=paste0(filename,".gpkg"),layername_hab="habitat",layername_res="results",rotation=as.integer(rotation$habitat),propRR=propRR,strat=strat))
}



