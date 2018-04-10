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

#include "functions.hpp"

/****************************************************************/
/*                         functions.c                          */
/****************************************************************/
/* Useful functions to manipulate C objects */

/* Perform the sum of a table of integer from pointer type in 3 dimensions on its 1st dimension */
void sum1_i3(int z, int l, int c, int ***t,int **tsum1) {
	int i,j,k;
	for (i=0;i<l;i++){
		for (j=0;j<c;j++){
			tsum1[i][j]=0;
			for (k=0;k<z;k++)
				tsum1[i][j] += t[k][i][j];
		} /* for j */
	} /* for i */
}

/* Perform the sum of a table of doubles from pointer type in 3 dimensions on its 1st dimension */
void sum1_d3(int z, int l, int c, double ***t, double **tsum1) {
	int i,j,k;
	for (i=0;i<l;i++){
		for (j=0;j<c;j++){
			tsum1[i][j]=0.0;
			for (k=0;k<z;k++)
				tsum1[i][j] += t[k][i][j];
		} /* for j */
	} /* for i */
}

/* Perform the sum of a table of integer from pointer type in 2 dimensions on its 2 dimensions */
int sum2_i2(int l, int c, int **t) {
	int i,j;
	int tsum=0;
	for (i=0;i<l;i++){
		for (j=0;j<c;j++)
			tsum += t[i][j];
	}
	return tsum;
}

/* Perform the sum of a table of float from pointer type in 2 dimensions on its 2 dimensions */
double sum2_d2(int l, int c, double **t) {
	int i,j;
	double tsum=0.0;
	for (i=0;i<l;i++){
		for (j=0;j<c;j++)
			tsum += t[i][j];
	}
	return tsum;
}

/* sort a table of integers */
void sort_i(int *a, int *n) {
     int i, j, temp;
     for(i=0;i<*n;i++) {   
          for(j=i+1;j<*n;j++) {     
               if(a[i]>a[j]) {	    
                    temp=a[j];
                    a[j]=a[i];
                    a[i]=temp;	    
               }	
          }
     }
     /*printf("table sorted :\n");
      for(i=0;i<*n;i++)
      printf("%d | ", a[i]);*/
}

/* sort a table of doubles */
void sort_d(double *a, int *n) {
     int i, j;
     double temp;
     for(i=0;i<*n;i++) {   
          for(j=i+1;j<*n;j++) {     
               if(a[i]>a[j]) {         
                    temp=a[j];
                    a[j]=a[i];
                    a[i]=temp;	    
               }	
          }
     }
     /*printf("table sorted :\n");
      for(i=0;i<*n;i++)
      printf("%.2f | ", a[i]);*/
}

/* Split a big integer into a vector of figures */
int split_i1(int n, int bigNber, int *vector) {
     int i;
     for (i=0 ; i<n ; i++){
          vector[n-1-i] = bigNber%10;
          bigNber/=10;
     }
  return 0;
}

/* function returning the max of the elements of a vector */
double max_d(double *a, int *n) {
     int i;
     double b[*n];  
     for (i=0;i<*n;i++)
          b[i]=a[i];       /* storage in b avoid to sort a */
     sort_d(b, n);
     return b[*n-1]; 
}


/*-----------------------------------*/
/*    Gamma & beta distributions     */
/*-----------------------------------*/
/* Find the values of alpha1 and alpha2 from expectation and variance of a gamma distribution */
void find_paramGamma(double exp, double var, double *alpha) {
     alpha[0] = exp * exp / var;  /* shape */
     alpha[1] = var / exp;        /* scale */
}

/* Find the values of alpha1 and alpha2 from expectation and variance of a beta distribution */
void find_paramBeta(double exp, double var, double *alpha) {
     alpha[0] = (exp * exp * (1-exp) / var) - exp;  /* shape 1 */
     alpha[1] = alpha[0] * (1-exp) / exp;           /* shape 2 */
}

/*-------------------------*/
/*    Sigmoid function     */
/*-------------------------*/
void sigmoid(double s, double k, double sig, double x, double *y){
     if (x < 1)
          *y = s * (1 - ((pow(M_E,-k * (pow(x,sig))) - pow(M_E,-k)) / (1 - pow(M_E,-k))));
     else 
          *y = s;
}

/*--------------------------*/
/*    Trade-off function    */
/*--------------------------*/
/* Trade-off function on aggressiveness traits (cf Debarre et al. JEB 2010) */
void tradeoff(int n, double *x, double *y, double beta){
     int i;
     for (i=0;i<n;i++)
          y[i] = 1 - pow(1-pow(x[i],1/beta),beta);
}




void addField(OGRLayer * poLayer, const char * fieldname, const double * values) {
  
  OGRFieldDefn oField( fieldname, OFTReal );
  oField.SetPrecision(32);
  if( poLayer->CreateField( &oField ) != OGRERR_NONE )
  {
    Rcpp::Rcerr<< "Creating Name field failed.\n" ;
  }

  OGRFeature *poFeature;

  int ind = 0;
  poLayer->ResetReading();
  while( (poFeature = poLayer->GetNextFeature()) != NULL )
  {
    //cerr<<ind<<" "<< values[ind]<<endl;
    poFeature->SetField(fieldname,values[ind++]);
    if( poLayer->SetFeature(poFeature) != OGRERR_NONE ) Rcpp::Rcerr<<"Error SetFeature layer\n";
    OGRFeature::DestroyFeature( poFeature );
  }

}
