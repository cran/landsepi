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

#ifndef __BASIC_PATHO__
#define __BASIC_PATHO__

#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include "functions.hpp"

// Pathogen aggressiveness components on a susceptible host for a pathogen genotype not adapted to resistance
struct Basic_patho {
    const double infection_rate;        // Maximal expected infection rate of a propagule on a healthy host
    const double propagule_prod_rate;   // Maximal expected propagule production rate per timestep and per infectious host
    const double latent_period_mean;     // Minimal expected latent period duration
    const double latent_period_var;     // Variance of the latent period duration
    const double infectious_period_mean; // Maximal expected infectious period duration
    const double infectious_period_var; // Variance of the infectious period duration
    const Vector2D<double> survival_prob;         // Off-season survival probability of a propagule for each croptype and each year
    const std::vector<double> repro_sex_prob;        // Probability for an infectious host to reproduce via sex
    const double sigmoid_kappa;         // Kappa parameter of the sigmoid contamination function
    const double sigmoid_sigma;         // Sigma parameter of the sigmoid contamination function
    const double sigmoid_plateau;       // Plateau parameter of the sigmoid contamination function
    const int sex_propagule_viability_limit;   // Maximum number of cropping seasons up to which a sexual propagule is viable
    const double sex_propagule_release_mean;       // average number of cropping seasons after which a sexual propagule is released.
    const bool clonal_propagule_gradual_release;        // whether or not clonal propagules surviving the bottleneck are gradually released along the following cropping season. 

    Basic_patho();
    Basic_patho(const double& infection_rate, const double& propagule_prod_rate, const double& latent_period_mean,
                const double& latent_period_var, const double& infectious_period_mean,
                const double& infectious_period_var, const Vector2D<double>& survival_prob, const std::vector<double>& repro_sex_prob,
                const double& sigmoid_kappa, const double& sigmoid_sigma, const double& sigmoid_plateau,
                const int& sex_propagule_viability_limit, const double& sex_propagule_release_mean,
                const bool& clonal_propagule_gradual_release);
    std::string to_string() const;
};

#endif
