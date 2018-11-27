/* 
 * Part of the landsepi R package.
 * Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@csiro.au>
 *                    Julien Papaix <julien.papaix@csiro.au>
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


#ifndef __PRINT_READ_WRITE__
#define __PRINT_READ_WRITE__

#include <stdio.h>
#include <stdlib.h>
#include "memory.hpp"
#include "functions.hpp"
#include <string.h>

/* ************************************************************************* */
/*                         printReadWrite.c                                  */
/* ************************************************************************* */

/* Print a vector of integer */
void print_i1(FILE *f, int n, int *t, char *title);

/* Print a vector of float */
void print_d1(FILE *f, int n, double *t, char *title); 

/* Print a table of integer */
void print_i2(FILE *f, int l, int c, int **t, const char *title);

/* Print a 3-dimension table of integer */
void print_i3(FILE *f, int z, int l, int c, int ***t, char *title);
/* Print a table of float */
void print_d2(FILE *f, int l, int c, double **t, char *title);
/* Print a 3-dimension table of float */
void print_d3(FILE *f, int z, int l, int c, double ***t, char *title);
/* Print the sum of the 1st dimension a table of integer of dimension 3 */
void print_i3sum1(FILE *f, int z, int l, int c, int ***t, char *title); 
/* Print the sum of the 1st dimension a table of float of dimension 3 */
void print_d3sum1(FILE *f, int z, int l, int c, double ***t, char *title);
/* Print the sum of the 1st dimension a table of integer of dimension 3 */
void print_i2sum2(FILE *f, int l, int c, int **t, char *title);
/* Print the sum of the 1st dimension a table of float of dimension 3 */
void print_d2sum2(FILE *f, int l, int c, double **t, char *title); 
/* Print the parameters in the console (f=stdout) or in an output .txt file */
void  print_param(FILE *f, int seed, char *nomfile_disppatho, char *nomfile_disphote, char *nomfile_habitat1, char *nomfile_habitat2, int nYears, int nTSpY, double pSurv, double eff, double *Tlat, double repro, double *Tspo, int Npoly, int NpolyTot, int Nhote, int Npatho, int Naggr, double *C0, double pI0, double *Cmax, double *delta, double *reproH, int *area, int **habitat, int *rotation, double **mort, char *strat, double propRR, int **resistance, int *adaptation, double MGeff, double QReff, double taumut, double costInfect, double costAggr, double beta);

/*   Function to read tables from .txt   */
//int lire(int MAX_LINE, int N_elements, char *nomfile, char *pdelim, double *y);
/* Write model output in .txt files and print output on screen */
void write_HHjuvSLIR(int Npoly, int Npatho, int Nhote, int t, int **H, int **Hjuv, int **S, int ***L, int ***I, int ***R, FILE *fH, FILE *fHjuv, FILE *fS, FILE *fL, FILE *fI, FILE *fR, int printOn);

#endif

