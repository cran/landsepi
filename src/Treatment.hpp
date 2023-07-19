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

#ifndef __TREATMENT__
#define __TREATMENT__

#include <string>
#include <vector>
#include <sstream>
#include <iterator>

// List of parameters related to pesticide treatments:
struct Treatment {
  const double treatment_degradation_rate; // degradation per time step of treatment concentration
  const double treatment_efficiency; // maximal efficiency of chemical treatments (i.e. fractional reduction of pathogen infection rate at the application date)
  const std::vector<int> treatment_timesteps; // Vector of time-steps corresponding to treatment application dates
  const std::vector<int> treatment_cultivars; // Vector of cultivar indices that receive treatments
  const double treatment_cost; // cost a single treatment application (monetary units/ha)
  const std::vector<double> treatment_application_threshold; // Vector of thresholds (i.e. disease severity, one for each treated cultivar) above which the treatment is applied in a polygon

  Treatment();
  Treatment(const double& treatment_degradation_rate, const double& treatment_efficiency, const std::vector<int>& treatment_timesteps,
              const std::vector<int>& cultivars, const double& treatment_cost, const std::vector<double>& treatment_application_threshold);
  std::string to_string() const;
};

#endif
