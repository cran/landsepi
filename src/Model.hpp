/*
 * Part of the landsepi R package.
 * Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
 *                    Julien Papaix <julien.papaix@inrae.fr>
 *                    Jean-François Rey <jean-francois.rey@inrae.fr>
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

#include <Rcpp.h>
using namespace Rcpp;
#include <gsl/gsl_randist.h> // gsl_ran_binomial, gsl_ran_poisson...
#include <gsl/gsl_rng.h> // gsl_rng_mt19937, gsl_rng*...
#include <gsl/gsl_cdf.h> // for cumulativee distribution
#include <gsl/gsl_linalg.h>
#include <math.h> // pow...
#include <stdio.h> // FILE*...
#include <array>
#include <fstream> // ofstream
#include <string>
#include <vector>
#include <chrono> 
using namespace std::chrono;
#include <iostream>
#include "Basic_patho.hpp"
#include "Croptype.hpp"
#include "Cultivar.hpp"
#include "Gene.hpp"
#include "Treatment.hpp"
#include "functions.hpp" // find_paramGamma, sigmoid, sample...
#include "printReadWrite.hpp" // write_HHjuvPLIR

//#include <Eigen/Eigen>
//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include "eigenmvn.h"

//#include <armadillo>


class Model {
  public:
    const int Nyears;               // Duration of the simulation (number of cropping seasons)
    const int time_steps_per_year;  // Number of timesteps per cropping season
    const int Npoly;                // Number of polygons (i.e. pieces of land where crops may be cultivated)
    const int Nhost;                // Number of host genotypes
    const int Npatho;               // Number of pathogen genotypes
    const int Ngene;                // Number of genes
    const std::vector<double> area;        // Vector of areas of polygons
    const Vector2D<int> rotation;          // Array of croptypes for each polygon and each year
    const gsl_rng* random_generator;       // Randomness generator
    const std::vector<Cultivar> cultivars; // Array of available cultivars
    const std::vector<Gene> genes;         // Array of available genes
    const Basic_patho basic_patho;         // Aggressiveness of a wild-type pathogen on a S cultivar and parameters relative to release of sexual propagules
    const Treatment treatment;             // Effect of chemical treatments on the pathogen
    const std::map<int,Croptype> croptypes; // List of pairs {cultivar, relative proportion} for each croptype
    const double sigmoid_kappa_host;       // Kappa parameter for the sigmoid invasion function (for host dispersal)
    const double sigmoid_sigma_host;       // Sigma parameter for the sigmoid invasion function (for host dispersal)
    const double sigmoid_plateau_host;     // Plateau parameter for the sigmoid invasion function (for host dispersal)
    const Vector3D<double> pI0;         // Initial inoculum as an array of dimensions (Npoly, Npatho, Nhost) with the probability to be I at t=0
    const Vector2D<double> disp_patho_clonal;     // Pathogen dispersal matrix (clonal propagules)
    const Vector2D<double> disp_patho_sex;     // Pathogen dispersal matrix (sexual propagules)
    const Vector2D<double> disp_host;      // Host dispersal matrix

    Model(const int& Nyears, const int& time_steps_per_year, const int& Npoly, const int& Nhost, const int& Npatho,
          const int& Ngene, const std::vector<double>& area, const Vector2D<int>& rotation,
          const gsl_rng* random_generator, const std::vector<Cultivar>& cultivars, const std::vector<Gene>& genes,
          const Basic_patho& basic_patho, const Treatment& treatment, const std::map<int,Croptype>& croptypes, const double& sigmoid_kappa_host,
          const double& sigmoid_sigma_host, const double& sigmoid_plateau_host, const Vector3D<double>& pI0,
          const Vector2D<double>& disp_patho_clonal, const Vector2D<double>& disp_patho_sex, const Vector2D<double>& disp_host, const int& seed);

    void dynepi();
    
    void infection(const int& t, std::vector<int>& H, const Vector2D<int>& Hcontaminated, Vector2D<int>& L,
                   Vector2D<int>& I, Vector2D<int>& R, Vector3D<int>& L2I, Vector3D<int>& I2R,
                   const std::vector<int>& activeR, const std::vector<int>& N,  const std::vector<int>& Nspray,
                   const std::vector<int>& t_lastspray);
    Vector2D<int> contamination(const std::vector<int>& H, const std::vector<int>& P, const std::vector<int>& N);
    void host_dynamic(const int& poly, const int& year, const int& t, std::vector<int>& H, std::vector<int>& Hjuv,
                             Vector2D<int>& L, Vector2D<int>& I, Vector2D<int>& R, Vector3D<int>& L2I, Vector3D<int>& I2R, 
                             std::vector<int>& N, std::vector<int>& Nspray, std::vector<int>& t_lastspray, std::vector<int>& TFI);
    Vector3D<int> bottleneck(const int& t, const Vector3D<int>& L, const Vector3D<int>& I,
                             const Vector2D<int>& activeQR, const int& year);
    //void dormancy(std::vector<int>& P, const int& year, const int& poly, Vector2D<int>& P_rest, std::vector<int>& P_germ);
    void dispersal(Vector2D<int>& Propagules, const Vector2D<double>& disp_matrix, const int& Ngeno);
    void mutation(std::vector<int>& P);
    void mutation_locus(const int& patho, const int& trait_mut, Vector2D<int>& PpathoMut);
    void reproSex(const int& t, std::vector<int>& P, const Vector2D<int>& I, const std::vector<int>& activeQR, 
                  const std::vector<int>& Nlevels_aggressiveness, const int& Nquali_gene);
    void between_season_pr_inoc(std::vector<int>& P_sex_primary_tmp, Vector2D<int>& P_stock, int& year);
    void in_season_pr_inoc(std::vector<int>& P_stock_release, Vector2D<int>& P_sex, const bool& distr);
    void reproClonal(const int& t, std::vector<int>& P, const Vector2D<int>& I, const std::vector<int>& activeR);
    std::array<Vector2D<int>, 2> split_IclonalIsex(const int& t, const Vector2D<int>& I);
    bool get_resistance(const int& index_gene, const int& host, const int& t, const int& activeR);
    double get_treat_effect(const int& Nt, const int& Nspray, const int& t, const int& t_lastspray);
    std::vector<int> get_P_stock_release(Vector2D<int>& P_stock_poly, const int& year);
    void get_P_daily(Vector2D<int>& P_daily, Vector3D<int>& P_primary, const int& t);
    Vector2D<int> get_sum_Vector2D(Vector2D<int>& M1, Vector2D<int>& M2);
        
    /* Init functions */
    int switch_aggr_to_patho(const std::vector<int>& aggr);
    std::vector<int> switch_patho_to_aggr(const int& index_patho);
    std::vector<double> switch_aggr_to_trait(const std::vector<int>& aggr, const std::vector<bool>& activeR);
    std::vector<int> switch_trait_to_aggr(const std::vector<double>& trait, const std::vector<bool>& activeR);
    Vector2D<int> init_activeR();
    void init_HjuvLIR(Vector2D<int>& Hjuv, Vector3D<int>& L, Vector3D<int>& I, Vector3D<int>& R);
   // void init_PgermPrest(Vector2D<int>& P_germ, Vector3D<int>& P_rest);
    void init_P(Vector2D<int>& P, Vector2D<int>& P_sex_secondary, Vector2D<int>& P_clonal_secondary, 
                 Vector3D<int>& P_sex_primary, Vector2D<int>& P_sex_primary_tmp, Vector3D<int>& P_clonal_primary, 
                 Vector2D<int>& P_clonal_primary_tmp, Vector2D<int>& P_sex_daily, Vector2D<int>& P_clonal_daily,
                 Vector3D<int>& P_stock);
    void init_Nspray_t_lastspray(Vector2D<int>& Nspray, Vector2D<int>& t_lastspray);
    void init_TFI(Vector3D<int>& TFI);
    void init_Nlevels_aggressiveness(std::vector<int>& Nlevels_aggressiveness);
    void init_L2I2R(Vector4D<int>& L2I, Vector4D<int>& I2R);
    Vector2D<int> intro_H(const int& year);
    void intro_I(Vector2D<int>& H, Vector3D<int>& I, Vector4D<int>& I2R, const Vector2D<int>& activeR);

    /* Random functions */
    double rng_uniform();
    int ran_poisson(const double& mu);
    double ran_gamma(const double& a, const double& b);
    int ran_binomial(const double& p, const int& n);
    std::vector<int> ran_multinomial(const int& N, const std::vector<double>& p);
    double ran_gaussian(const double& mu, const double& sigma);
    double ran_exponential(const double& mu);
    double ran_exponential_trunc(const double& mu, const double& limit_sup);
    std::vector<double>ran_multivariate_gaussian(const std::vector<double>& mu, const std::vector<std::vector<double>>& A);
    std::vector<std::vector<double>> ran_multisample_multivariate_gaussian(unsigned n, const std::vector<double>& mu, const std::vector<std::vector<double>>& A);

    /* Print parameters in an output .txt file */
    void print_param(const int& seed, const std::vector<double>& mutation_prob, const std::vector<double>& efficiency,
                     const std::vector<double>& adaptation_cost, const std::vector<double>& relative_advantage, 
                     const std::vector<double>& tradeoff_strength);

    /* Write model output in .txt files and print output on screen */
    void write_HHjuvPLIR(const Vector2D<int>& H, const Vector2D<int>& Hjuv, const Vector2D<int>& P,
                         const Vector3D<int>& L, const Vector3D<int>& I, const Vector3D<int>& R, FILE* fH, FILE* fHjuv,
                         FILE* fP, FILE* fL, FILE* fI, FILE* fR);

    /* Write model output in .txt files and print output on screen (Pbefinter ONLY) */
    void write_Pbefinter(const Vector3D<int>& eqIsurv, FILE* feqIsurv, const Vector2D<int>& Pbefinter, FILE* fPbefinter);
    
    /* Write model output in .txt files and print output on screen */
    void write_TFI(const Vector2D<int>& TFI, FILE* fTFI);

};


/********************
 * Inline functions
 ********************/

inline double Model::rng_uniform() {
    return gsl_rng_uniform(this->random_generator);
}

inline int Model::ran_poisson(const double& mu) {
    return gsl_ran_poisson(this->random_generator, mu);
}

inline double Model::ran_gamma(const double& a, const double& b) {
    return gsl_ran_gamma(this->random_generator, a, b);
}

inline int Model::ran_binomial(const double& p, const int& n) {
    return gsl_ran_binomial(this->random_generator, p, n);
}

inline std::vector<int> Model::ran_multinomial(const int& N, const std::vector<double>& p) {
    const unsigned int K = p.size();
    std::vector<int> res(K);
    gsl_ran_multinomial(this->random_generator, K, N, p.data(), (unsigned int*)res.data());
    return res;
}

inline double Model::ran_gaussian(const double& mu, const double& sigma) {
    /* This function returns a Gaussian random variate, with mean mu and standard deviation sigma */
    return mu+gsl_ran_gaussian(this->random_generator, sigma);
}

inline double Model::ran_exponential(const double& mu) {
    return gsl_ran_exponential(this->random_generator, mu);
}

inline double Model::ran_exponential_trunc(const double& mu, const double& limit_sup) {
    /* This function returns a truncated exponential random variate, with mean mu, 
    lower limit = 0 and upper limit = limit_sup 
    (following approach https://www.r-bloggers.com/2020/08/generating-data-from-a-truncated-distribution/ )*/
    const double& limit_inf = 0.0; 
    double F_inf = gsl_cdf_exponential_P(limit_inf, mu);
    double F_sup = gsl_cdf_exponential_P(limit_sup, mu);
    //double u = gsl_ran_flat(this->random_generator, F_inf, F_sup);
    double u = (F_sup - F_inf)*rng_uniform()+F_inf;
    double f_u = gsl_cdf_exponential_Pinv(u, mu);
    return f_u;
}


inline std::vector<double> Model::ran_multivariate_gaussian(const std::vector<double>& mu, const std::vector<std::vector<double>>& A){
    /* This function generates a random vector satisfying the k-dimensional multivariate
     Gaussian distribution with mean \mu and variance-covariance matrix A.*/
    
    //std::vector<std::vector<int>> L(K,std::vector<int>(K));
    gsl_vector * mup = gsl_vector_alloc(mu.size());
    gsl_vector * res = gsl_vector_alloc(mu.size());
    for(size_t i = 0; i<mu.size(); i++)
        gsl_vector_set(mup, i, mu[i]);
    gsl_matrix * L = gsl_matrix_alloc(A.size(), A[0].size());
    for(size_t i = 0; i<A.size(); i++){
        for(size_t j=0; j<A[0].size(); j++){
            gsl_matrix_set(L, i, j, A[i][j]);
        }
    }
    gsl_linalg_cholesky_decomp(L);
    //std::vector<double> res(K);
    
    gsl_ran_multivariate_gaussian(this->random_generator, mup, L, res);
    //std::vector<double> resp(res->data, &res->data[res->size-1]);
    std::vector<double> resp(mu.size());
    for(unsigned int i=0; i<res->size; i++ ){
        resp[i]=res->data[i];
    }
    
    gsl_vector_free(res);
    gsl_vector_free(mup);
    gsl_matrix_free(L);
    return resp;
}

inline std::vector<std::vector<double>> Model::ran_multisample_multivariate_gaussian(unsigned n, const std::vector<double>& mu, const std::vector<std::vector<double>>& A){
    /* This function generates a random vector satisfying the k-dimensional multivariate
     Gaussian distribution with mean \mu and variance-covariance matrix A.*/
    
    //std::vector<std::vector<int>> L(K,std::vector<int>(K));
    std::vector<std::vector<double>> resp(n, std::vector<double>(mu.size()));
    gsl_vector * mup = gsl_vector_alloc(mu.size());
    gsl_vector * res = gsl_vector_alloc(mu.size());
    for(size_t i = 0; i<mu.size(); i++)
        gsl_vector_set(mup, i, mu[i]);
    gsl_matrix * L = gsl_matrix_alloc(A.size(), A[0].size());
    for(size_t i = 0; i<A.size(); i++){
        for(size_t j=0; j<A[0].size(); j++){
            gsl_matrix_set(L, i, j, A[i][j]);
        }
    }
    gsl_linalg_cholesky_decomp(L);
    //std::vector<double> res(K);
    
    for(size_t i = 0; i<n ; i++){
        gsl_ran_multivariate_gaussian(this ->random_generator, mup, L, res);
        std::copy(res->data, &(res->data[res->size]), resp[i].begin());
        /*for(unsigned int i=0; i<res->size; i++ ){
         resp[i]=res->data[i];
        }*/
        
    }
    
    gsl_vector_free(res);
    gsl_vector_free(mup);
    gsl_matrix_free(L);
    return resp;
}


//' @title Model for Landscape Epidemiology & Evolution
//' @name model_landsepi
//' @description Stochastic, spatially-explicit, demo-genetic model simulating the spread and evolution of a 
//' plant pathogen in a heterogeneous landscape.
//' 
//' @param time_param list of simulation parameters:\itemize{ 
//' \item Nyears = number cropping seasons, 
//' \item nTSpY = number of time-steps per cropping season.
//' }
//' @param area_vector a vector containing areas of polygons (i.e. fields), in surface units.
//' @param rotation_matrix a matrix containing for each field (rows) and year (columns, named "year_1", "year_2", etc.), 
//' the index of the cultivated croptype. Importantly, the matrix must contain 1 more column than the real number 
//' of simulated years. 
//' @param croptypes_cultivars_prop a matrix with three columns named 'croptypeID' for croptype index, 
//' 'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype. 
//' @param dispersal list of dispersal parameters:\itemize{ 
//' \item disp_patho_clonal = vectorised dispersal matrix of the pathogen (clonal propagules), 
//' \item disp_patho_sex = vectorised dispersal matrix of the pathogen (sexual propagules), 
//' \item disp_host = vectorised dispersal matrix of the host.
//' }
//' @param inits list of initial conditions:\itemize{
//' \item pI0 = vector of length Npoly*Npatho*Nhost giving the probability to be infectious (i.e. state I) at 
//' t=0 pr each polygon, pathogen genotype and host.
//' }
//' @param seed seed (for random number generation).
//' @param cultivars_param list of parameters associated with each host genotype (i.e. cultivars) 
//' when cultivated in pure crops:\itemize{   
//' \item initial_density = vector of host densities (per surface unit) at the beginning of the cropping season,  
//' \item max_density = vector of maximum host densities (per surface unit) at the end of the cropping season, 
//' \item growth rate = vector of host growth rates, 
//' \item reproduction rate = vector of host reproduction rates, 
//' \item relative_yield_H = Yield of H individuals relative to H individuals (100%)
//' \item relative_yield_L = Yield of L individuals relative to H individuals
//' \item relative_yield_I = Yield of I individuals relative to H individuals
//' \item relative_yield_R = Yield of R individuals relative to H individuals
//' \item sigmoid_kappa_host = kappa parameter for the sigmoid invasion function (for host dispersal),
//' \item sigmoid_sigma_host = sigma parameter for the sigmoid invasion function (for host dispersal),
//' \item sigmoid_plateau_host = plateau parameter for the sigmoid invasion function (for host dispersal),
//' \item cultivars_genes_list = a list containing, for each host genotype, the indices of carried resistance genes,
//' } 
//' @param basic_patho_param list of i. pathogen aggressiveness parameters on a susceptible host 
//' for a pathogen genotype not adapted to resistance and ii. sexual reproduction parameters: \itemize{
//' \item infection_rate = maximal expected infection rate of a propagule on a healthy host, 
//' \item propagule_prod_rate = maximal expected reproduction_rate of an infectious host per timestep, 
//' \item latent_period_mean = minimal expected duration of the latent period, 
//' \item latent_period_var = variance of the latent period duration, 
//' \item infectious_period_mean = maximal expected duration of the infectious period, 
//' \item infectious_period_var = variance of the infectious period duration,
//' \item survival_prob = matrix giving the probability for a propagule to survive the off-season, for each croptype (rows) and each year (columns) 
//' \item repro_sex_prob = vector of probabilities for an infectious host to reproduce via sex rather than via cloning for each timestep, 
//' \item sigmoid_kappa = kappa parameter of the sigmoid contamination function, 
//' \item sigmoid_sigma = sigma parameter of the sigmoid contamination function, 
//' \item sigmoid_plateau = plateau parameter of the sigmoid contamination function,
//' \item sex_propagule_viability_limit = maximum number of cropping seasons up to which a sexual propagule is viable
//' \item sex_propagule_release_mean = average number of cropping seasons after which a sexual propagule is released.
//' \item clonal_propagule_gradual_release = whether or not clonal propagules surviving the bottleneck are gradually released along the following cropping season.
//' }
//' @param genes_param list of parameters associated with each resistance gene and with the evolution of 
//' each corresponding pathogenicity gene:\itemize{ 
//' \item target_trait = vector of aggressiveness components (IR, LAT, IP, or PR) targeted by resistance genes, 
//' \item efficiency = vector of resistance gene efficiencies (percentage of reduction of the targeted 
//' aggressiveness component: IR, 1/LAT, IP and PR), 
//' \item age_of_activ_mean = vector of expected delays to resistance activation (for APRs), 
//' \item age_of_activ_var = vector of variances of the delay to resistance activation (for APRs),  
//' \item mutation_prob = vector of mutation probabilities for pathogenicity genes (each of them corresponding to a resistance gene), 
//' \item Nlevels_aggressiveness = vector of number of adaptation levels related to each resistance gene (i.e. 1 + number 
//' of required mutations for a pathogenicity gene to fully adapt to the corresponding resistance gene), 
//' \item adaptation_cost = vector of adaptation penalties paid by pathogen genotypes fully adapted 
//' to the considered resistance genes on all hosts, 
//' \item relative_advantage = vector of fitness advantages of a pathogen genotype fully adapted to the resistance genes 
//' on hosts carrying these genes, relative to those that do not carry these genes,
//' \item tradeoff_strength = vector of strengths of the trade-off relationships between the level of aggressiveness 
//' on hosts that do and do not carry the resistance genes.
//' }
//' @param treatment_param list of parameters related to pesticide treatments: \itemize{ 
//' \item treatment_degradation_rate = degradation rate (per time step) of chemical concentration,
//' \item treatment_efficiency = maximal efficiency of chemical treatments (i.e. fractional reduction 
//' of pathogen infection rate at the time of application),
//' \item treatment_timesteps = vector of time-steps corresponding to treatment application dates,
//' \item treatment_cultivars = vector of indices of the cultivars that receive treatments,
//' \item treatment_cost = cost of a single treatment application (monetary units/ha),
//' \item treatment_application_threshold = vector of thresholds (i.e. disease severity, one for each treated cultivar) above which the treatment is applied
//' }
//' 
//' @details See \code{?landsepi} for details on the model and assumptions. 
//' Briefly, the model is stochastic, spatially explicit (the basic spatial unit is an individual field), based on a SEIR
//' (‘susceptible-exposed-infectious-removed’, renamed HLIR for 'healthy-latent-infectious-removed' to avoid confusions 
//' with 'susceptible host') structure with a discrete time step. It simulates the spread and 
//' evolution (via mutation, recombination through sexual reproduction, selection and drift) 
//' of a pathogen in a heterogeneous cropping landscape, across cropping seasons split by host harvests which impose
//' potential bottlenecks to the pathogen. A wide array of resistance deployment strategies 
//' (possibly including chemical treatments) can be simulated.
//'  
//' @return A set of binary files is generated for every year of simulation and every compartment: 
//' \itemize{
//'  \item H: healthy hosts,
//'  \item Hjuv: juvenile healthy hosts (for host reproduction),
//'  \item L: latently infected hosts,
//'  \item I: infectious hosts,
//'  \item R: removed hosts,
//'  \item P: propagules.}
//' Each file indicates for every time-step the number of individuals in each field, and when 
//' appropriate for each host and pathogen genotypes). Additionally, a binary file called TFI is 
//' generated and gives the Treatment Frequency Indicator (expressed as the number of treatment applications 
//'  per polygon).
//' 
//' @examples
//' \dontrun{
//' #### Spatially-implicit simulation with 2 patches (S + R) during 3 years ####
//' 
//' ## Simulation parameters
//' time_param <- list(Nyears=3, nTSpY=120)
//' Npoly=2
//' Npatho=2
//' area <- c(100000, 100000)
//' basic_patho_param <- loadPathogen(disease = "rust")
//' basic_patho_param$repro_sex_prob <- rep(0, time_param$nTSpY+1)
//'      cultivars <- as.list(rbind(loadCultivar(name="Susceptible", type="growingHost")
//' , loadCultivar(name="Resistant", type="growingHost")))
//' names(cultivars)[names(cultivars)=="cultivarName"] <- "name"
//' yield0 <- cultivars$yield_H + as.numeric(cultivars$yield_H==0)
//' cultivars <- c(cultivars, list(relative_yield_H = as.numeric(cultivars$yield_H / yield0)
//' , relative_yield_L = as.numeric(cultivars$yield_L / yield0)
//' , relative_yield_I = as.numeric(cultivars$yield_I / yield0)
//' , relative_yield_R = as.numeric(cultivars$yield_R / yield0)
//' , sigmoid_kappa_host=0.002, sigmoid_sigma_host=1.001, sigmoid_plateau_host=1
//' , cultivars_genes_list=list(numeric(0),0)))
//' rotation <- data.frame(year_1=c(0,1), year_2=c(0,1), year_3=c(0,1), year_4=c(0,1))
//' croptypes_cultivars_prop <- data.frame(croptypeID=c(0,1), cultivarID=c(0,1), proportion=c(1,1))
//' genes <- as.list(loadGene(name="MG", type="majorGene"))
//' treatment=list(treatment_degradation_rate=0.1,
//' treatment_efficiency=0, 
//' treatment_timesteps=logical(0),
//' treatment_cultivars=logical(0),
//' treatment_cost=0,
//' treatment_application_threshold = logical(0))
//'   
//' ## run simulation
//' model_landsepi(seed=1,
//' time_param = time_param,
//' basic_patho_param = basic_patho_param,
//' inits = list(pI0=c(0.1, rep(0, 7))),
//' area_vector = area,
//' dispersal = list(disp_patho_clonal=c(0.99,0.01,0.01,0.99),
//' disp_patho_sex=c(1,0,0,1),
//' disp_host=c(1,0,0,1)),
//' rotation_matrix = as.matrix(rotation),
//' croptypes_cultivars_prop = as.matrix(croptypes_cultivars_prop),
//' cultivars_param = cultivars, 
//' genes_param = genes,
//' treatment_param = treatment)
//'
//' ## Compute outputs
//' eco_param <- list(yield_perHa = cbind(H = as.numeric(cultivars$relative_yield_H),
//' L = as.numeric(cultivars$relative_yield_L),
//' I = as.numeric(cultivars$relative_yield_I),
//' R = as.numeric(cultivars$relative_yield_R)),
//' planting_cost_perHa = as.numeric(cultivars$planting_cost),
//' market_value = as.numeric(cultivars$market_value))
//'
//' evol_res <- evol_output(, time_param, Npoly, cultivars, genes)
//' epid_res <-  epid_output(, time_param, Npatho, area, rotation
//' , croptypes_cultivars_prop, cultivars, eco_param, treatment, basic_patho_param)
//'
//'
//'
//' #### 1-year simulation of a rust epidemic in pure susceptible crop in a single 1-km2 patch ####
//' ## Simulation and pathogen parameters
//' time_param <- list(Nyears=1, nTSpY=120)
//' area <- c(1E6)
//' basic_patho_param = loadPathogen(disease = "rust")
//' basic_patho_param$repro_sex_prob <- rep(0, time_param$nTSpY+1)
//' ## croptypes, cultivars and genes
//' rotation <- data.frame(year_1=c(0), year_2=c(0))
//' croptypes_cultivars_prop <- data.frame(croptypeID=c(0), cultivarID=c(0), proportion=c(1))
//'        cultivars <- as.list(rbind(loadCultivar(name="Susceptible", type="growingHost")))
//' names(cultivars)[names(cultivars)=="cultivarName"] <- "name"
//' yield0 <- cultivars$yield_H + as.numeric(cultivars$yield_H==0)
//' cultivars <- c(cultivars, list(relative_yield_H = as.numeric(cultivars$yield_H / yield0)
//' , relative_yield_L = as.numeric(cultivars$yield_L / yield0)
//'     , relative_yield_I = as.numeric(cultivars$yield_I / yield0)
//' , relative_yield_R = as.numeric(cultivars$yield_R / yield0)
//' , sigmoid_kappa_host=0.002, sigmoid_sigma_host=1.001, sigmoid_plateau_host=1
//' , cultivars_genes_list=list(numeric(0))))
//' genes <-   list(geneName = character(0) , adaptation_cost = numeric(0)
//' , relative_advantage = numeric(0)
//' , mutation_prob = numeric(0)
//' , efficiency = numeric(0) , tradeoff_strength = numeric(0)
//' , Nlevels_aggressiveness = numeric(0)
//' , age_of_activ_mean = numeric(0) , age_of_activ_var = numeric(0)
//' , target_trait = character(0)
//' , recombination_sd = numeric(0))
//' treatment=list(treatment_degradation_rate=0.1
//'                               , treatment_efficiency=0
//' , treatment_timesteps=logical(0)
//' , treatment_cultivars=logical(0)
//' , treatment_cost=0
//' , treatment_application_threshold = logical(0))
//' 
//' ## run simulation
//' model_landsepi(seed=1, time_param = time_param
//' , basic_patho_param = basic_patho_param
//' , inits = list(pI0=5E-4), area_vector = area
//' , dispersal = list(disp_patho_clonal=c(1), disp_patho_sex=c(1), disp_host=c(1))
//' , rotation_matrix = as.matrix(rotation)
//' , treatment_param = treatment
//'                                 , croptypes_cultivars_prop = as.matrix(croptypes_cultivars_prop)
//' , cultivars_param = cultivars,  genes_param = genes)
//'  }
//' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018).
//' Assessing the durability andefficiency of landscape-based strategies to deploy plant 
//' resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
//' 
//' @export
// [[Rcpp::export]]
void model_landsepi(Rcpp::List time_param
                        , Rcpp::NumericVector area_vector
                        , Rcpp::NumericMatrix rotation_matrix
                        , Rcpp::NumericMatrix croptypes_cultivars_prop
                        , Rcpp::List dispersal
                        , Rcpp::List inits
                        , int seed
                        , Rcpp::List cultivars_param
                        , Rcpp::List basic_patho_param
                        , Rcpp::List genes_param
                        , Rcpp::List treatment_param
                        );
#endif
