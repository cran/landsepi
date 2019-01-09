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


#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <gdal.h>
//#include <gdal/org_api.h>
#include <ogrsf_frmts.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

#ifndef GDALV2
#if GDAL_VERSION_MAJOR >= 2
# define GDALV2 1
#endif
#endif


/****************************************************************/
/*                         functions.c                          */
/****************************************************************/
/* Useful functions to manipulate C objects */

/* Convert an integer to its binary equivalent using a recursive procedure */
int as_binary(int num);

/* Perform the sum of a table of integer from pointer type in 3 dimensions on its 1st dimension */
void sum1_i3(int z, int l, int c, int ***t,int **tsum1);

/* Perform the sum of a table of doubles from pointer type in 3 dimensions on its 1st dimension */
void sum1_d3(int z, int l, int c, double ***t, double **tsum1);
/* Perform the sum of a table of integer from pointer type in 2 dimensions on its 2 dimensions */
int sum2_i2(int l, int c, int **t);
/* Perform the sum of a table of float from pointer type in 2 dimensions on its 2 dimensions */
double sum2_d2(int l, int c, double **t);
/* sort a table of integers */
void sort_i(int *a, int *n);
/* sort a table of doubles */
void sort_d(double *a, int *n);
/* Split a big integer into a vector of figures */
int split_i1(int n, int bigNber, int *vector);
/* function returning the max of the elements of a vector */
double max_d(double *a, int *n);
/*-----------------------------------*/
/*    Gamma & beta distributions     */
/*-----------------------------------*/
/* Find the values of alpha1 and alpha2 from expectation and variance of a gamma distribution */
void find_paramGamma(double exp, double var, double *alpha);
/* Find the values of alpha1 and alpha2 from expectation and variance of a beta distribution */
void find_paramBeta(double exp, double var, double *alpha); 
/*-------------------------*/
/*    Sigmoid function     */
/*-------------------------*/
void sigmoid(double s, double k, double sig, double x, double *y);

/*--------------------------*/
/*    Trade-off function    */
/*--------------------------*/
/* Trade-off function on aggressiveness traits (cf Debarre et al. JEB 2010) */
void tradeoff(int n, double *x, double *y, double beta);

/*--------------------------*/
/*     Sample functions     */
/*--------------------------*/
/* Single random draw using a vector of probabilities   */
int sample_multinomial_once(const gsl_rng *gen, double *cumProb);
/* Samples an integer array without replacement until entire array has been sampled once */
void sample(const gsl_rng *gen, int *inArray, int inLength, int *outArray, int outLength);


void addField(OGRLayer * poLayer, const char * fieldname, const double * values);

#endif
