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

#include "Gene.hpp"

/* Initialisation of mutation matrix (only called by Gene constructor) 
 values outside the bounds 0-1 are associated to class 0 and Nlevelaggressiveness-1, respectively
 mutation_prob is the probability to mutate at once from the first to the last class of aggressiveness
 */

Vector2D<double> Gene::init_mutkernel(const double& mutation_prob) {
  //Mutation matrix for pathogen genes 
  Vector2D<double> mutkernel(Nlevels_aggressiveness, std::vector<double>(Nlevels_aggressiveness));
  std::vector<double> level_aggress_extremes(Nlevels_aggressiveness+1);
  std::vector<double> mean_level_aggress(Nlevels_aggressiveness);
  level_aggress_extremes[0] = 0.0;
  for(int i = 1; i < Nlevels_aggressiveness+1; i++) {
    level_aggress_extremes[i] = level_aggress_extremes[i-1] + 1.0/Nlevels_aggressiveness;
  }
  
  for(int i = 0; i < Nlevels_aggressiveness; i++) {
    mean_level_aggress[i] = (level_aggress_extremes[i] + level_aggress_extremes[i+1])/2.0;
  }
  
  double sigma = abs(1/sqrt(2)*(level_aggress_extremes[Nlevels_aggressiveness-1] - mean_level_aggress[0])/(erfinv(2*(1-mutation_prob)-1)));
  
  for(int r = 0; r < Nlevels_aggressiveness; r++){
    double mu = mean_level_aggress[r];
    for(int c = 0; c < Nlevels_aggressiveness; c++){
      double z_inf = level_aggress_extremes[c];
      double z_sup = level_aggress_extremes[c+1];
      const double x_inf = z_inf - mu;
      const double x_sup = z_sup - mu;
      if(c == 0){
        mutkernel[r][c] = cdf_gaussian_P(x_sup,sigma);
      } else if (c == Nlevels_aggressiveness-1){
        mutkernel[r][c] = 1-cdf_gaussian_P(x_inf,sigma);
      } else {
        mutkernel[r][c] = cdf_gaussian_P(x_sup,sigma) - cdf_gaussian_P(x_inf,sigma);
      }
      
    }
  }
  return mutkernel;
}


/* Initialisation of aggressiveness matrix (only called by Gene constructor) */
Vector2D<double> Gene::init_aggressiveness_matrix(const double& efficiency, const double& adaptation_cost,
                                                  const double& relative_advantage,
                                                  const double& tradeoff_strength) {
  // Leonard KJ (1977) Selection Pressures and Plant Pathogens. 
  // Annals of the New York Academy of Sciences, 287, 207-222. https://doi.org/10.1111/j.1749-6632.1977.tb34240.x
  
    /* Aggressiveness matrix */
    Vector2D<double> aggressiveness_matrix(Nlevels_aggressiveness, std::vector<double>(2));
    const double aggressiveness_0 = 1 - efficiency;

    /* (Nlevels_aggressiveness-1): intervals between Nlevels_aggressiveness values */
    const double step = (Nlevels_aggressiveness > 1) ? 1 / static_cast<double>(Nlevels_aggressiveness - 1) : 0;

    std::vector<double> gain(Nlevels_aggressiveness);
    for(int i = 0; i < Nlevels_aggressiveness; i++) {
        gain[i] = i * step;
    }

    const std::vector<double> cost = tradeoff(gain, tradeoff_strength);
    for(int i = 0; i < Nlevels_aggressiveness; i++) {
        aggressiveness_matrix[i][1] = aggressiveness_0 + gain[i] * (efficiency - adaptation_cost + relative_advantage);
        aggressiveness_matrix[i][0] = 1 - cost[i] * adaptation_cost;
    }
    return aggressiveness_matrix;
}

// Default constructor
Gene::Gene()
    : age_of_activ_mean(0.0),
      age_of_activ_var(0.0),
      Nlevels_aggressiveness(0),
      target_trait(""),
      mutkernel(Vector2D<double>()),
      aggressiveness_matrix(Vector2D<double>()),
      recombination_sd(0.0){}

// Constructor
Gene::Gene(const double& age_of_activ_mean, const double& age_of_activ_var, const int& Nlevels_aggressiveness,
           const std::string& target_trait, const double& mutation_prob, const double& efficiency,
           const double& adaptation_cost, const double& relative_advantage, const double& tradeoff_strength, 
           const double& recombination_sd)
    : age_of_activ_mean(age_of_activ_mean),
      age_of_activ_var(age_of_activ_var),
      Nlevels_aggressiveness(Nlevels_aggressiveness),
      target_trait(target_trait),
      mutkernel(init_mutkernel(mutation_prob)),
      aggressiveness_matrix(init_aggressiveness_matrix(efficiency, adaptation_cost, relative_advantage, tradeoff_strength)),
      recombination_sd(recombination_sd){}

// Transform Gene attributs to string
std::string Gene::to_string() const {
    std::string str("");
    str += "  age_of_activ_mean:      " + std::to_string(this->age_of_activ_mean) + "\n";
    str += "  age_of_activ_var:       " + std::to_string(this->age_of_activ_var) + "\n";
    str += "  Nlevels_aggressiveness: " + std::to_string(this->Nlevels_aggressiveness) + "\n";
    str += "  target_trait:           " + this->target_trait + "\n";
    str += "  recombination_sd:       " + std::to_string(this->recombination_sd) + "\n";
    str += "  mutkernel:\n";
    for(int i = 0; i < this->Nlevels_aggressiveness; i++) {
        str += "    " + std::to_string(i) + ": ";
        for(int j = 0; j < this->Nlevels_aggressiveness; j++) {
            str += std::to_string(mutkernel[i][j]) + " ";
        }
        str += "\n";
    }
    str += "  aggressiveness_matrix:\n";
    for(int i = 0; i < this->Nlevels_aggressiveness; i++) {
        str += "    " + std::to_string(i) + ": ";
        str += std::to_string(aggressiveness_matrix[i][0]) + " ";
        str += std::to_string(aggressiveness_matrix[i][1]) + "\n";
    }
    return str;
}
