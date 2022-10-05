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

#include "Basic_patho.hpp"

// Default constructor
Basic_patho::Basic_patho()
    : infection_rate(0.0),
      propagule_prod_rate(0.0),
      latent_period_mean(0.0),
      latent_period_var(0.0),
      infectious_period_mean(0.0),
      infectious_period_var(0.0),
      survival_prob(0.0),
      repro_sex_prob(0.0),
      sigmoid_kappa(0.0),
      sigmoid_sigma(0.0),
      sigmoid_plateau(0.0),
      sex_propagule_viability_limit(0),
      sex_propagule_release_mean(0.0),
      clonal_propagule_gradual_release(0){}

// Constructor
Basic_patho::Basic_patho(const double& infection_rate, const double& propagule_prod_rate,
                         const double& latent_period_mean, const double& latent_period_var,
                         const double& infectious_period_mean, const double& infectious_period_var,
                         const double& survival_prob, const std::vector<double>& repro_sex_prob, const double& sigmoid_kappa,
                         const double& sigmoid_sigma, const double& sigmoid_plateau,
                         const int& sex_propagule_viability_limit, const double& sex_propagule_release_mean,
                         const bool& clonal_propagule_gradual_release
                         )
    : infection_rate(infection_rate),
      propagule_prod_rate(propagule_prod_rate),
      latent_period_mean(latent_period_mean),
      latent_period_var(latent_period_var),
      infectious_period_mean(infectious_period_mean),
      infectious_period_var(infectious_period_var),
      survival_prob(survival_prob),
      repro_sex_prob(repro_sex_prob),
      sigmoid_kappa(sigmoid_kappa),
      sigmoid_sigma(sigmoid_sigma),
      sigmoid_plateau(sigmoid_plateau),
      sex_propagule_viability_limit(sex_propagule_viability_limit),
      sex_propagule_release_mean(sex_propagule_release_mean),
      clonal_propagule_gradual_release(clonal_propagule_gradual_release){}

// Transform Basic_patho attributs to string
std::string Basic_patho::to_string() const {
  std::ostringstream oss;
  std::copy(this->repro_sex_prob.begin(), this->repro_sex_prob.end(), std::ostream_iterator<int>(oss, ", "));
  
    std::string str("");
    str += "  infection_rate:        " + std::to_string(this->infection_rate) + "\n";
    str += "  propagule_prod_rate:   " + std::to_string(this->propagule_prod_rate) + "\n";
    str += "  latent_period_mean:     " + std::to_string(this->latent_period_mean) + "\n";
    str += "  latent_period_var:     " + std::to_string(this->latent_period_var) + "\n";
    str += "  infectious_period_mean: " + std::to_string(this->infectious_period_mean) + "\n";
    str += "  infectious_period_var: " + std::to_string(this->infectious_period_var) + "\n";
    str += "  survival_prob:         " + std::to_string(this->survival_prob) + "\n";
    str += "  repro_sex_prob:        " + oss.str() + "\n";
    str += "  sigmoid_kappa:         " + std::to_string(this->sigmoid_kappa) + "\n";
    str += "  sigmoid_sigma:         " + std::to_string(this->sigmoid_sigma) + "\n";
    str += "  sigmoid_plateau:       " + std::to_string(this->sigmoid_plateau) + "\n";
    str += "  sex_propagule_viability_limit:         " + std::to_string(this->sex_propagule_viability_limit) + "\n";
    str += "  sex_propagule_release_mean:       " + std::to_string(this->sex_propagule_release_mean) + "\n";
    str += "  clonal_propagule_gradual_release:       " + std::to_string(this->clonal_propagule_gradual_release) + "\n";
        return str;
}
