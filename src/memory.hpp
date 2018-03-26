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


#ifndef __MEMORY__
#define __MEMORY__

#include <stdlib.h>
  /*---------------------------------------*/
  /*    Functions for memory allocation    */
  /*---------------------------------------*/

double *init_d1(int r1);
double **init_d2(int r1, int r2);
void free_d2(double **x, int r1);
double ***init_d3(int r1, int r2, int r3);
void free_d3(double ***x, int r1, int r2);
double ****init_d4(int r1, int r2, int r3, int r4);
void free_d4(double ****x, int r1, int r2, int r3);
int *init_i1(int r1);
int **init_i2(int r1, int r2);
void free_i2(int **x, int r1);
int ***init_i3(int r1, int r2, int r3);
void free_i3(int ***x, int r1, int r2);
int ****init_i4(int r1, int r2, int r3, int r4);
void free_i4(int ****x, int r1, int r2, int r3);
int *****init_i5(int r1, int r2, int r3, int r4, int r5);
void free_i5(int *****x, int r1, int r2, int r3, int r4);
int ******init_i6(int r1, int r2, int r3, int r4, int r5, int r6);
void free_i6(int ******x, int r1, int r2, int r3, int r4, int r5);
int *******init_i7(int r1, int r2, int r3, int r4, int r5, int r6, int r7);
void free_i7(int *******x, int r1, int r2, int r3, int r4, int r5, int r6);
int ********init_i8(int r1, int r2, int r3, int r4, int r5, int r6, int r7, int r8);
void free_i8(int ********x, int r1, int r2, int r3, int r4, int r5, int r6, int r7);
#endif
