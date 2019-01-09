/* 
 * Part of the landsepi R package.
 * Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inra.fr>
 *                    Julien Papaix <julien.papaix@inra.fr>
 *                    Jean-Frnaçois Rey <jean-francois.rey@inra.fr>
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


#ifndef __INITIALISATION__
#define __INITIALISATION__

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "functions.hpp"
#include "memory.hpp"
#include "printReadWrite.hpp"

/****************************************************************/
/*                    initialisation.c                          */
/****************************************************************/
/* Initialisation of the conversion matrices between pathotype index and aggressiveness formula */
void init_aggrFormula(int Npatho, int Naggr, int *adaptation, int ********aggrToPatho, int **pathoToAggr); 

/* Initialisation of mutation matrices */
void init_mutkernel(double **mutkernelMG, double **mutkernelQR, double taumut, int Naggr);

/* Initialisation of infectivity and aggressiveness matrices */
void init_infectAggr(int Naggr, double MGeff, double QReff, double costInfect, double costAggr, double beta, double **infect, double **aggr);
/* Initialisation of H, L, I and R at 0 */
void init_HHjuvSLIR(int Npoly, int Nhost, int Npatho, int **H, int **Hjuv, int **S, int ***L, int ***I, int ***R); 
/* Initialise L2I and I2R with 0 */
void init_L2I2R(int Npoly, int Npatho, int Nhost, int nTSpY, int ****L2I, int ****I2R);

/* Initialisation of S at 0 */
void init_S(int Npoly, int Npatho, int **S); 
    
/* Initialise SpathoMut at 0 */
void init_SpathoMut(int Npatho, int **SpathoMut);
/* Plantation of H in the beginning of a season */
void intro_H(int Npoly, int Nhost, int **H, int *area, int **habitat, double *C0, char *strat, int id_rotation);
/* Pathogen introduction : infectious sites from pathotype 0 in cultivar 0 */
void intro_I(const gsl_rng *gen, int nTSpY, int Npoly, int Nhost, int Npatho, int **H, int ***I, int ****I2R, double pI0, double *Tspo, int **resistance, double **infect, double **aggr, int ID_E, int ID_S);

/* activeQR: matrix indicating for each year and field, at which time step QR starts to be active */
void init_activeQR(const gsl_rng *gen, int nYears, int Npoly, double *timeToQR, double *timeToQR_alpha, int **activeQR);
#endif
