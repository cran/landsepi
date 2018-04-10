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


#include "memory.hpp"

  /*---------------------------------------*/
  /*    Functions for memory allocation    */
  /*---------------------------------------*/

double *init_d1(int r1){
  double *x;
  
  x= (double *) malloc(r1*sizeof(double));
  
  return x;
}

double **init_d2(int r1, int r2){
  double **x;
  int i;
  
  x= (double **) malloc(r1*sizeof(double *));
  for (i=0;i<r1;i++)
    x[i]=init_d1(r2);
  
  return x;
}

void free_d2(double **x, int r1){
  int i;
  
  for (i=0;i<r1;i++)
    {free(x[i]);
     x[i]=NULL;}
    
  free(x);
}

double ***init_d3(int r1, int r2, int r3){
  double ***x;
  int i;
  
  x= (double ***) malloc(r1*sizeof(double **));
  for (i=0;i<r1;i++)
    x[i]=init_d2(r2,r3);
  
  return x;
}

void free_d3(double ***x, int r1, int r2){
  int i;
  
  for (i=0;i<r1;i++)
    free_d2(x[i],r2);
  free(x);
}

double ****init_d4(int r1, int r2, int r3, int r4){
  double ****x;
  int i;
  
  x= (double ****) malloc(r1*sizeof(double ***));
  for (i=0;i<r1;i++)
    x[i]=init_d3(r2,r3,r4);
  
  return x;
}

void free_d4(double ****x, int r1, int r2, int r3){
  int i;
  
  for (i=0;i<r1;i++)
    free_d3(x[i],r2,r3);
  free(x);
}

int *init_i1(int r1){
  int *x;
  
  x= (int *) malloc(r1*sizeof(int));
  
  return x;
}

int **init_i2(int r1, int r2){
  int **x;
  int i;
  
  x= (int **) malloc(r1*sizeof(int *));
  for (i=0;i<r1;i++)
    x[i]=init_i1(r2);
  
  return x;
}

void free_i2(int **x, int r1){
  int i;
  
  for (i=0;i<r1;i++)
    {free(x[i]);
     x[i]=NULL;}
    
  free(x);
}

int ***init_i3(int r1, int r2, int r3){
  int ***x;
  int i;
  
  x= (int ***) malloc(r1*sizeof(int **));
  for (i=0;i<r1;i++)
    x[i]=init_i2(r2,r3);
  
  return x;
}

void free_i3(int ***x, int r1, int r2){
  int i;
  
  for (i=0;i<r1;i++)
    free_i2(x[i],r2);
  free(x);
}

int ****init_i4(int r1, int r2, int r3, int r4){
  int ****x;
  int i;
  
  x= (int ****) malloc(r1*sizeof(int ***));
  for (i=0;i<r1;i++)
    x[i]=init_i3(r2,r3,r4);
  
  return x;
}

void free_i4(int ****x, int r1, int r2, int r3){
  int i;
  
  for (i=0;i<r1;i++)
    free_i3(x[i],r2,r3);
  free(x);
}

int *****init_i5(int r1, int r2, int r3, int r4, int r5){
  int *****x;
  int i;
  
  x= (int *****) malloc(r1*sizeof(int ****));
  for (i=0;i<r1;i++)
    x[i]=init_i4(r2,r3,r4,r5);
  
  return x;
}

void free_i5(int *****x, int r1, int r2, int r3, int r4){
  int i;
  
  for (i=0;i<r1;i++)
    free_i4(x[i],r2,r3,r4);
  free(x);
}

int ******init_i6(int r1, int r2, int r3, int r4, int r5, int r6){
  int ******x;
  int i;
  
  x= (int ******) malloc(r1*sizeof(int *****));
  for (i=0;i<r1;i++)
    x[i]=init_i5(r2,r3,r4,r5,r6);
  
  return x;
}

void free_i6(int ******x, int r1, int r2, int r3, int r4, int r5){
  int i;
  
  for (i=0;i<r1;i++)
    free_i5(x[i],r2,r3,r4,r5);
  free(x);
}

int *******init_i7(int r1, int r2, int r3, int r4, int r5, int r6, int r7){
  int *******x;
  int i;
  
  x= (int *******) malloc(r1*sizeof(int ******));
  for (i=0;i<r1;i++)
    x[i]=init_i6(r2,r3,r4,r5,r6,r7);
  
  return x;
}

void free_i7(int *******x, int r1, int r2, int r3, int r4, int r5, int r6){
  int i;
  
  for (i=0;i<r1;i++)
    free_i6(x[i],r2,r3,r4,r5,r6);
  free(x);
}

int ********init_i8(int r1, int r2, int r3, int r4, int r5, int r6, int r7, int r8){
  int ********x;
  int i;
  
  x= (int ********) malloc(r1*sizeof(int *******));
  for (i=0;i<r1;i++)
    x[i]=init_i7(r2,r3,r4,r5,r6,r7,r8);
  
  return x;
}

void free_i8(int ********x, int r1, int r2, int r3, int r4, int r5, int r6, int r7){
  int i;
  
  for (i=0;i<r1;i++)
    free_i7(x[i],r2,r3,r4,r5,r6,r7);
  free(x);
}
