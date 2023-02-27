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

#ifndef __CULTIVAR__
#define __CULTIVAR__

#include "Gene.hpp"

// Cultivar parameters
struct Cultivar {
    const double initial_density;    // Host density per surface unit at the beginning of the cropping season (in pure crop)
    const double max_density;        // Maximum host density per surface unit at the end of the cropping season (in pure crop)
    const double growth_rate;        // Growth rate
    const double reproduction_rate;  // Reproduction rate
//    const double death_rate;         // Death rate
    const double relative_yield_H;   // Yield of H individuals relative to H individuals (100%)
    const double relative_yield_L;   // Yield of L individuals relative to H individuals
    const double relative_yield_I;   // Yield of I individuals relative to H individuals
    const double relative_yield_R;   // Yield of R individuals relative to H individuals
    const std::vector<int> genes_id; // Indices of carried resistance genes

    Cultivar();
    Cultivar(const double& initial_density
                 , const double& max_density
                 , const double& growth_rate
                 , const double& reproduction_rate
             //    , const double& death_rate
                 , const double& relative_yield_H
                 , const double& relative_yield_L
                 , const double& relative_yield_I
                 , const double& relative_yield_R
                 , const std::vector<int>& genes_id);
    std::string to_string() const;
};

#endif
