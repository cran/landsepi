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

#include "functions.hpp"
#include "Model.hpp"


// int i not equal to 0
bool all_is_zero(int i)
{
  return i==0;
}

/****************************************************************/
/*                         functions.c                          */
/****************************************************************/
/* Useful functions to manipulate C objects */

/* Convert an integer to its binary equivalent using a recursive procedure */
int as_binary(const int& num) {
    if(num == 0) {
        return 0;
    } else {
        return (num % 2) + 10 * as_binary(num / 2);
    }
}

/*-----------------------------------*/
/*    Gamma & beta distributions     */
/*-----------------------------------*/
/* Find the values of alpha1 and alpha2 from expectation and variance of a gamma distribution */
std::array<double, 2> find_paramGamma(const double& exp, const double& var) {
    std::array<double, 2> alpha;
    alpha[0] = exp * exp / var; /* shape */
    alpha[1] = var / exp; /* scale */
    return alpha;
}

/* Find the values of alpha1 and alpha2 from expectation and variance of a beta distribution */
std::array<double, 2> find_paramBeta(const double& exp, const double& var) {
    std::array<double, 2> alpha;
    alpha[0] = (exp * exp * (1 - exp) / var) - exp; /* shape 1 */
    alpha[1] = alpha[0] * (1 - exp) / exp; /* shape 2 */
    return alpha;
}

/*-------------------------*/
/*    Sigmoid function     */
/*-------------------------*/
/* Sigmoid function used for spore contamination and host dispersal */
double sigmoid(const double& plateau, const double& kappa, const double& sigma, const double& x) {
    if(x < 1) {
        return plateau * (1 - ((pow(M_E, -kappa * (pow(x, sigma))) - pow(M_E, -kappa)) / (1 - pow(M_E, -kappa))));
    } else {
        return plateau;
    }
}

/*--------------------------*/
/*    Trade-off function    */
/*--------------------------*/
/* Trade-off function on aggressiveness traits (cf Debarre et al. JEB 2010) */
std::vector<double> tradeoff(const std::vector<double>& x, const double& strength) {
    std::vector<double> y(x.size());
    for(unsigned int i = 0; i < x.size(); i++) {
        y[i] = 1 - pow(1 - pow(x[i], 1 / strength), strength);
    }
    return y;
}

double cdf_gaussian_P(const double& x, const double& sigma) {
    return gsl_cdf_gaussian_P(x,sigma);
}

/*--------------------------*/
/*     Sample functions     */
/*--------------------------*/
/* Single random draw using a vector of probabilities   */
/* cumProb must be ordered with last element equal to 1 */
int sample_multinomial_once(const gsl_rng* gen, const std::vector<double>& cumProb) {
    int j = 0;
    const double randNum = gsl_rng_uniform(gen);
    while(randNum > cumProb[j]) {
        j++;
    }
    return j;
}

/* Samples an integer array without replacement until entire array has been sampled once */
/* then continues by sampling with replacement */
/* The output sequence is also randomised (thus can be used to randomise a vector) */
std::vector<int> sample(const gsl_rng* gen, const std::vector<int>& inArray) {
    const unsigned int length = inArray.size();
    std::vector<int> outArray(length);
    double inIndicesLength = static_cast<double>(length); /* Temporary length of indices array */
    double outIndicesLength = static_cast<double>(length); /* Temporary length of indices array */
    std::vector<int> inArrayCopy(length); /* To store copy of inArray */
    std::vector<int> outIndices(length); /* Array of indices for empty slots in outArray */

    for(unsigned int i = 0; i < length; i++) {
        inArrayCopy[i] = inArray[i];
        outIndices[i] = i;
    }
    for(unsigned int i = 0; i < length; i++) {
        /* NOTE: (--) is done AFTER affectation */
        const int l = static_cast<int>(outIndicesLength-- * gsl_rng_uniform(gen));
        const int j = outIndices[l]; /* Choose an Index for output */
        /* Replace index value with final index in array */
        outIndices[l] = outIndices[static_cast<int>(outIndicesLength)];
        if(i < length) {
            const int k = static_cast<int>(inIndicesLength-- * gsl_rng_uniform(gen)); /* Choose an Index for input */
            outArray[j] = inArrayCopy[k]; /* Copy input value to output */
            /* Replace index value with final index in array */
            inArrayCopy[k] = inArrayCopy[static_cast<int>(inIndicesLength)];
        } else {
            const int k = static_cast<int>(length * gsl_rng_uniform(gen)); /* Choose an Index for input */
            outArray[j] = inArray[k]; /* Copy input value to output */
        }
    }
    return outArray;
}



std::vector<double> mu_transformation(const std::vector<double>& mu_old, const int& n_propagule) {
  
  std::vector<double> mu_trans;
  for (unsigned int i=0; i<mu_old.size(); i++ ){
      mu_trans.insert(mu_trans.end(), n_propagule, mu_old[i]);
  }
  return mu_trans;
}

std::vector<std::vector<double>> cov_transformation(std::vector<std::vector<double>>& cov_old, const int& n_propagule) {
  std::vector<std::vector<double>> cov_trans(n_propagule*cov_old.size(),std::vector<double>(n_propagule*cov_old.size(),0.0));
  for(unsigned int r=0; r < cov_old.size(); r++){
    for(unsigned int c=0; c < cov_old.size(); c++){
      for(int p=0; p < n_propagule; p++){
        cov_trans[n_propagule*r+p][n_propagule*c+p] = cov_old[r][c];
      }
    }
  }
  return cov_trans;
}
/*
std::vector<std::vector<double>> reshape(const std::vector<double>& flat_vec, std::size_t ncols )
{ //reshape a 1d vector into a matrix with ncols  
  // sanity check
  if( ncols == 0 || flat_vec.size()%ncols != 0 ) throw std::domain_error( "bad #cols" ) ;
  
  const auto nrows = flat_vec.size() / ncols ;
  
  std::vector<std::vector<double>> mtx ;
  const auto begin = std::begin(flat_vec) ;
  for( std::size_t row = 0 ; row < nrows ; ++row ){
    for(std::size_t col = 0 ; col < ncols ; ++col ){
    //mtx.push_back( { begin + row*ncols, begin + (row+1)*ncols } ) ;
    mtx.push_back( flat_vec[begin + (row+nrows)*col]) ;
    }
  }
  return mtx ;
}
*/

/*--------------------------*/
/*     erfinv function      */
/*--------------------------*/
/* erfinv(X) computes the inverse error function of X.  */

  double erfinv(double x) {
    return( gsl_cdf_ugaussian_Pinv((1+x)/2) / sqrt(2) );
  }


