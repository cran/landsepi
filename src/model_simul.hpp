/* 
 * Part of the landsepi R package.
 * Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inra.fr>
 *                    Julien Papaix <julien.papaix@inra.fr>
 *                    Jean-François Rey <jean-francois.rey@inra.fr>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,i
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#ifndef __MODEL_SIMUL__
#define __MODEL_SIMUL__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gdal.h>
#include <ogrsf_frmts.h>

#ifndef GDALV2
#if GDAL_VERSION_MAJOR >= 2
# define GDALV2 1
#endif
#endif

#include <Rcpp.h>
#include "memory.hpp"
#include "functions.hpp"            /* functions used in print.c + sigmoid function + trade-off function */
#include "printReadWrite.hpp"       /* Printing, reading, writing functions                              */
#include "initialisation.hpp"

#define NLOCI 8          /* number of resistance/adaptation loci */
#define ID_E 4           /* index of QR on infection efficiency    */
#define ID_L 5           /* index of QR on latent period duration  */
#define ID_R 6           /* index of QR on reproduction rate       */
#define ID_S 7           /* index of QR on sporulation duration    */

using namespace Rcpp;
using namespace std;
//' @title Model Landscape Epidemiology & Evolution
//' @name model_landsepi
//' @description Stochastic, spatially-explicit, demo-genetic model simulating the spread and evolution of a pathogen in a heterogeneous landscape.
//' @param timeP list of simulation parameters (number of years, number of time-steps per year)
//' @param landscape landscape generated through AgriLand
//' @param dispersal list of dispersal parameters (vectorised dispersal matrix of the pathogen, vectorised dispersal matrix of the host)
//' @param inits list initial conditions (initial probability of infection by the pathogen)
//' @param val_seed seed (for random number generation)
//' @param hostP list of host parameters (number of cultivars, initial planting density, maximal carrying capacity, growth rate, reproduction rate, death rate, resistance formula, 
//' parameters of the sigmoid invasion function: kappa, sigma and s)
//' @param pathoP list of pathogen parameters (probability to survive the off-season, probability to reproduce via sex rather than via cloning, 
//' infection rate, reproduction rate, average latent period duration, variance of the latent period, average infectious period duration
//' , variance of the infectious period duration, parameters of the sigmoid contamination function: kappa, sigma, s)
//' @param evolP list of evolution parameters (cost of infectivity, cost of aggressiveness, mutation rate, efficiency of major 
//' resistance genes, efficiency of quantitative resistance, trade-off strength, number of increments of quantitative 
//' resistance erosion, average time to expression of quantitative resistance, Variance of the time to expression of quantitative resistance,
//' adaptation formula)
//' @details \itemize{
//' \item The model is stochastic, spatially explicit (the basic spatial unit is an individual field), based on a SEIR 
//' (‘susceptible-exposed-infectious-removed’) structure with a discrete time step. It simulates the spread and evolution 
//' of a pathogen in an agricultural landscape, across cropping seasons split by host harvests which impose potential bottlenecks 
//' to the pathogen. 
//' \item A wide array of deployment strategies can be simulated: mosaics, mixtures, rotations and pyramiding of multiple
//'  major resistance genes which affect pathogen infectivity, and up to four quantitative resistance traits. 
//'  These traits target different aggressiveness components of the pathogen, i.e. the infection rate, the duration of the latent 
//'  period and the infectious period, and the propagule production rate. Quantitative resistance may be expressed from the time of planting, 
//'  or later in the cropping season (Adult Plant Resistance or Mature Plant Resistance). 
//'  \item The genotype of cultivated plant cultivars is specified using 
//'  the "resistance formulas", i.e. a vector of size 8. the four first elements indicate whether the cultivar carries major resistance 
//'  genes #1, #2, #3 and #4, respectively. The following four elements indicate whether the cultivar carried a quantitative resistance 
//'  trait against the infection rate, the latent period duration, the sporulation rate, or the sporulation duration of the pathogen, respectively. 
//'  For example, the formula c(1,0,0,0,0,1,0,0) indicates the presence of major gene #1 and a quantitative resistance 
//'  which increases the duration of the latent period of the pathogen.
//'  \item Initially, the pathogen is not adapted to any source of 
//'  resistance, and is only present on susceptible hosts. However, through mutation, it can evolve and may acquire infectivity 
//'  genes (which leads to breakdown of major resistance genes) or increase aggressiveness (which leads to the erosion of the 
//'  relevant quantitative resistance traits). These genes may also be reassorted via sexual reproduction.
//'  \item Evolution of a pathogen toward infectivity or increased aggressiveness on 
//'  a resistant host is often penalised by a fitness cost on susceptible hosts. Consequently, in the present model, pathogens 
//'  carrying infectivity genes may have reduced infection rate (cost of infectivity) on susceptible hosts relative 
//'  to pathogens that do not carry these genes. Similarly, a gain in pathogen aggressiveness on quantitatively resistant hosts 
//'  is penalised by a decreased aggressiveness on susceptible hosts, leading to a trade-off.
//'  }
//' @return A set of binary files is generated for every year of simulation and every compartment: \itemize{ 
//' \item H: healthy hosts,
//' \item Hjuv: juvenile healthy hosts,
//' \item L: latently infected hosts,
//' \item I: infectious hosts,
//' \item R: removed hosts,
//' \item S: propagules.
//' } 
//' Each file indicates for every time-step the number of individuals in each field, and when appropriate for each cultivar and pathotype)
//' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018). Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
//' @export
// [[Rcpp::export]]
void model_landsepi(Rcpp::List timeP, Rcpp::List landscape, Rcpp::List dispersal, Rcpp::List inits ,int val_seed, Rcpp::List hostP, Rcpp::List pathoP, Rcpp::List evolP); 

void split_IclonalIsex(const gsl_rng *gen, int Npatho, int Nhost, int **I, int **Iclonal_poly, int **Isex_poly, double probSex);
void reproClonal(const gsl_rng *gen, int t, int Nhost, int Npatho, int *S, int **I, double repro, int **resistance, int **pathoToAggr, double **aggr, int activeQR);
void reproSex(const gsl_rng *gen, int t, int Nhost, int Npatho, int *S, int **I, double repro, int **resistance, int **pathoToAggr, int ********aggrToPatho, double **aggr, int activeQR);
void mutation_locus(const gsl_rng *gen, int patho, int Npatho, int Naggr, int trait_mut, double **mutkernel, int **pathoToAggr, int ********aggrToPatho, int **SpathoMut);
void mutation(const gsl_rng *gen, int Npatho, int Naggr, int *S, int *adaptation, double **mutkernelMG,  double **mutkernelQR, int **pathoToAggr, int ********aggrToPatho);
void dispersal(const gsl_rng *gen, int Npoly, int Nhost, int Npatho, int **H, int **Hjuv, int **S, double *reproH, double **disphost, double **disppatho);
void bottleneck(const gsl_rng *gen, int Npoly, int t, int Nhost, int Npatho, double pSurv, double *Tspo, int ***L, int ***I, double **aggr, int **resistance, int **pathoToAggr, int ***eqIsurv, int *activeQR);
void host_dynamic(const gsl_rng *gen, int Nhost, int Npatho, int area, double *Cmax, int *H, int *Hjuv, int **L, int **I, int **R, int *N, double *mortH, double *delta, double khost, double sighost, double shost);
void contamination(const gsl_rng *gen, int Nhost, int Npatho, int *H, int *S, int *N, int **Hcontaminated, double kpatho, double sigpatho, double spatho);
void infection(const gsl_rng *gen, int t, int nTSpY, int Nhost, int Npatho, int *H, int **Hcontaminated, int **L, int **I, int **R, int ***L2I, int ***I2R, double eff, double *Tlat, double *Tspo, int **resistance, int **pathoToAggr, double **infect, double **aggr, int activeQR);
void dynepi(const gsl_rng *gen, int nYears, int nTSpY, int Npoly, int Nhost, int Npatho, int *area, int **habitat, int *rotation, double *C0, double pI0, double pSurv, double probSex, double *Cmax, double **mort, double *delta, double *reproH, double **disphost, double **disppatho, double eff, double *Tlat, double repro, double *Tspo, double **mutkernelMG,  double **mutkernelQR, int **pathoToAggr, int ********aggrToPatho, int Naggr, int **activeQR, double **infect, double **aggr, char *strat, int **resistance, int *adaptation, int printOn, OGRLayer *poLayer, double khost, double sighost, double shost, double kpatho, double sigpatho, double spatho);

/*void sortie_layer(OGRLayer *poLayer,int **H,int Npoly, int Nhost, int timeStep, int year);*/

                
#endif
