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

#include "Treatment.hpp"

// Default constructor
Treatment::Treatment()
  : treatment_degradation_rate(0),
    treatment_efficiency(0),
    treatment_timesteps({}),
    treatment_cultivars({}),
    treatment_cost(0),
    treatment_application_threshold({}){}

// Constructor
Treatment::Treatment(const double& treatment_degradation_rate, const double& treatment_efficiency, const std::vector<int>& treatment_timesteps,
          const std::vector<int>& treatment_cultivars, const double& treatment_cost, const std::vector<double>& treatment_application_threshold)
  : treatment_degradation_rate(treatment_degradation_rate),
    treatment_efficiency(treatment_efficiency),
    treatment_timesteps(treatment_timesteps),
    treatment_cultivars(treatment_cultivars),
    treatment_cost(treatment_cost),
    treatment_application_threshold(treatment_application_threshold){}

// Transform treatment attributs to string
std::string Treatment::to_string() const {
  std::ostringstream oss, oss2, oss3;
  std::copy(this->treatment_timesteps.begin(), this->treatment_timesteps.end(), std::ostream_iterator<int>(oss, ", "));
  std::copy(this->treatment_cultivars.begin(), this->treatment_cultivars.end(), std::ostream_iterator<int>(oss2, ", "));
  std::copy(this->treatment_application_threshold.begin(), this->treatment_application_threshold.end(), std::ostream_iterator<double>(oss3, ", "));
   
  std::string str("");
  str += "  treatment degradation rate:        " + std::to_string(this->treatment_degradation_rate) + "\n";
  str += "  treatment efficiency:         " + std::to_string(this->treatment_efficiency) + "\n";
  str += "  treatment timesteps:         " + oss.str() + "\n";
  str += "  treatment cultivars:         " + oss2.str() + "\n";
  str += "  treatment cost:    " + std::to_string(this->treatment_cost) + "\n";
  str += "  treatment application threshold:    " + oss3.str() + "\n";
  return str;
}
