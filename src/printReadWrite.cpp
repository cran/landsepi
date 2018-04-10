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


/* ************************************************************************* */
/*                         printReadWrite.c                                  */
/* ************************************************************************* */
#include "printReadWrite.hpp"

/* Print a vector of integer */
void print_i1(FILE *f, int n, int *t, char *title) {
	int i;
	if (title!="") 
		fprintf(f, "%s : \n", title);
	for (i=0; i<n; i++) 
		fprintf (f, "%7d", t[i]);
	fprintf(f, "\n");
}

/* Print a vector of float */
void print_d1(FILE *f, int n, double *t, char *title) {
	int i;
	if (title!="") 
		fprintf(f, "%s : \n", title);
	for (i=0; i<n; i++)
		fprintf (f, " %.3f ", t[i]);
	fprintf(f, "\n");
}

/* Print a table of integer */
void print_i2(FILE *f, int l, int c, int **t, const char *title) {
	int i ;
	if (title!="")
		fprintf(f, "%s : \n", title);
	for (i=0; i<l; i++)
		print_i1(f,c,t[i],"");
	fprintf(f, "\n");
}

/* Print a 3-dimension table of integer */
void print_i3(FILE *f, int z, int l, int c, int ***t, char *title) {
	int k ;
	if (title!="")
		fprintf(f, "%s : \n", title);
	for (k=0; k<z; k++)
		print_i2(f,l,c,t[k],"");
	fprintf(f, "\n");
}

/* Print a table of float */
void print_d2(FILE *f, int l, int c, double **t, char *title) {
	int i ;
	if (title!="")
		fprintf(f, "%s : \n", title);
	for (i=0; i<l; i++)
		print_d1(f,c,t[i],"");
	fprintf(f, "\n");
}

/* Print a 3-dimension table of float */
void print_d3(FILE *f, int z, int l, int c, double ***t, char *title) {
	int k ;
	if (title!="")
		fprintf(f, "%s : \n", title);
	for (k=0; k<z; k++)
		print_d2(f,l,c,t[k],"");
	fprintf(f, "\n");
}

/* Print the sum of the 1st dimension a table of integer of dimension 3 */
void print_i3sum1(FILE *f, int z, int l, int c, int ***t, char *title) {
	int **tsum;
	tsum = init_i2(l,c);
	sum1_i3(z,l,c,t,tsum); 
	print_i2(f,l,c,tsum,title);
	free_i2(tsum,l);
}

/* Print the sum of the 1st dimension a table of float of dimension 3 */
void print_d3sum1(FILE *f, int z, int l, int c, double ***t, char *title) {
	double **tsum;
	tsum = init_d2(l,c);
	sum1_d3(z,l,c,t,tsum); 
	print_d2(f,l,c,tsum,title);
	free_d2(tsum,l);
}  

/* Print the sum of the 1st dimension a table of integer of dimension 3 */
void print_i2sum2(FILE *f, int l, int c, int **t, char *title) {
	int tsum;
	tsum = sum2_i2(l,c,t); 
	fprintf(f, "%s = %5d\n",title,tsum);
}

/* Print the sum of the 1st dimension a table of float of dimension 3 */
void print_d2sum2(FILE *f, int l, int c, double **t, char *title) {
	double tsum;
	tsum = sum2_d2(l,c,t); 
	fprintf(f, "%s = %.3f\n",title,tsum);
}

/* Print the parameters in the console (f=stdout) or in an output .txt file */
void  print_param(FILE *f, int seed, char *nomfile_disppatho, char *nomfile_disphote, char *nomfile_habitat1, char *nomfile_habitat2, int nYears, int nTSpY, double pSurv, double eff, double *Tlat, double repro, double *Tspo, int Npoly, int NpolyTot, int Nhote, int Npatho, int Naggr, double *C0, double pI0, double *Cmax, double *delta, double *reproH, int *area, int **habitat, int *rotation, double **mort, char *strat, double propRR, int **resistance, int *adaptation, double MGeff, double QReff, double taumut, double costInfect, double costAggr, double beta) {

  int hote;
  
  fprintf(f, "###     MODEL PARAMETERS     ###\n");  
  fprintf(f, "seed:                      %d\n", seed);
  fprintf(f, "pathogen dispersal kernel: %s\n", nomfile_disppatho);
  fprintf(f, "host dispersal kernel:     %s\n", nomfile_disphote);
  fprintf(f, "landscape1:                %s\n", nomfile_habitat1);
  fprintf(f, "landscape2:                %s\n", nomfile_habitat2);
  
  fprintf(f, "\n*****             Seasonality                  *****\n");
  fprintf(f, "nYears:     %d\n", nYears); 
  fprintf(f, "nDays/year: %d\n", nTSpY);
  fprintf(f, "pI0:        %.4f\n", pI0);
  fprintf(f, "Relative bottleneck size: %.0e\n", pSurv);
  
  fprintf(f, "\n*****     Landscape, host parameters and initial conditions     *****\n"); 
  fprintf(f, "Npoly:    %d    simulated (%d available)\n", Npoly, NpolyTot);
  fprintf(f, "Nhote:    %d       Npatho:    %d       Naggr:    %d\n", Nhote, Npatho, Naggr);
  print_d1(f, 2, C0,  "C0 [2]");
  print_d1(f, Nhote, Cmax,"Cmax [hote]");
  print_d1(f, Nhote, delta, "delta [hote]");
  print_d1(f, Nhote, reproH, "reproH [hote]");
  if (Npoly<11){
  	print_i1(f, Npoly, area, "area [poly]");
  	print_i2(f, 2, Npoly, habitat, "habitat [1:2][poly]");	
	print_d2(f, Npoly, Nhote, mort, "mortH [poly][hote]");
  }    
  if (nYears<11)
       print_i1(f, nYears+1, rotation, "rotation [nYears+1]");
  fprintf(f, "*****  Epidemiological and evolution parameters *****\n"); 
  fprintf(f, "eff:    %.2f     Tlat:   %.2f     repro:    %.2f      Tspo:    %.2f\n", eff, Tlat[0], repro, Tspo[0]);
  fprintf(f, "taumut: %.0e\n", taumut);
  fprintf(f, "Tlat_var:    %.2f       Tspo_var:    %.2f\n", Tlat[1], Tspo[1]);
  
  fprintf(f, "\n*****  Deployment and resistance parameters *****\n"); 
  fprintf(f, "deployment strategy:    %s\n", strat);
  fprintf(f, "Cropping ratio 2 (R2/(R1+R2)):    %.2f\n", propRR);
  print_i2(f, Nhote, 8, resistance, "resistance [hote][mg1, mg2, mg3, mg4, qr_E, qr_L, qr_R, qr_S]");
  print_i1(f, 8, adaptation, "adaptation [mg1, mg2, mg3, mg4, qr_E, qr_L, qr_R, qr_S]");
  fprintf(f, "# Major-gene resistance\n");
  fprintf(f, "MGeff:         %.2f\n", MGeff);
  fprintf(f, "costInfect:  %.2f\n", costInfect);
  fprintf(f, "# Quantitative resistance\n");
  fprintf(f, "QReff:       %.2f\n", QReff);
  fprintf(f, "costAggr  :  %.2f\n", costAggr);
  fprintf(f, "beta:         %.2f\n", beta);
  fprintf(f, "\n");
}


/*   Function to read tables from .txt   */
/*int lire(int MAX_LINE, int N_elements, char *nomfile, char *pdelim, double *y) {
     FILE *fpi = NULL;
     char lu[MAX_LINE];
     double x[N_elements]; /* tableau des valeurs */
/*     char *p;
     int i,j;

     fpi = fopen (nomfile, "r");
     if (!fpi) {
          //printf("Peut pas ouvrir le fichier: %s \n",nomfile);
          return(0);
     }

     /* Lecture de la 1iere ligne */
/*     i=0;
     fgets (lu, MAX_LINE, fpi);
     
     /* décodage des valeurs */
     /* La 1iere valeur de la ligne */
/*     p = strtok (lu, pdelim);
     x[i] = atof (p);
     i++;
     
     /* Les valeurs suivantes: 1ier argument de strok=NULL */
/*     while (feof(fpi) == 0) {
          while ((p = strtok (NULL, pdelim)) != NULL)
          {
               x[i] = atof (p);
               i++;
          }
          
          /* Lecture de la ligne suivante */
//          fgets (lu, MAX_LINE, fpi);
          /* Décodage de la 1iere valeur */
/*          p = strtok (lu, pdelim);
          x[i] = atof (p);
          i++;
     } /* while */
     
/*     for (j=0; j<(i-1); j++){
          y[j] = x[j];
     }
     fclose(fpi);
     return(0);
}

/* Write model output in .txt files and print output on screen */
void write_HHjuvSLIR(int Npoly, int Npatho, int Nhote, int t, int **H, int **Hjuv, int **S, int ***L, int ***I, int ***R, FILE *fH, FILE *fHjuv, FILE *fS, FILE *fL, FILE *fI, FILE *fR, int printOn) {
     
     int poly,patho,hote;  
     for (poly = 0; poly < Npoly; poly++) {
          for  (patho = 0 ; patho < Npatho ; patho++){
               fwrite(&S[poly][patho],sizeof(int),1,fS);
               for (hote = 0 ; hote < Nhote ; hote++){
                    fwrite(&L[poly][patho][hote],sizeof(int),1,fL);
                    fwrite(&I[poly][patho][hote],sizeof(int),1,fI);			
                    fwrite(&R[poly][patho][hote],sizeof(int),1,fR); 			
               } /* for hote */
          } /* for patho */
          for (hote = 0 ; hote < Nhote ; hote++){
               fwrite(&Hjuv[poly][hote],sizeof(int),1,fHjuv);
               fwrite(&H[poly][hote],sizeof(int),1,fH);
          } /* for hote */  
     } /* for poly */
     
     /* Print model output on screen */
     /*if (printOn) {
          printf("### t = %d\n", t);
          print_i2(stdout, Npoly,Npatho,S,"S [poly][patho]");
          print_i2(stdout, Npoly,Nhote,H,"H [poly][hote]");
          print_i2(stdout, Npoly,Nhote,Hjuv,"Hjuv [poly][hote]");
          print_i3sum1(stdout, Npoly,Npatho,Nhote,L,"L (sum poly)[patho][hote]");
          print_i3sum1(stdout, Npoly,Npatho,Nhote,I,"I (sum poly)[patho][hote]");
          print_i3sum1(stdout, Npoly,Npatho,Nhote,R,"R (sum poly)[patho][hote]");  
     }*/
} 



