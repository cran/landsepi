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
void print_param(FILE *f, unsigned long int seed, int nYears, int nTSpY, int Npoly, int Nhost, char *strat, double *C0, double *Cmax, double *growthH, double *reproH, double pI0, double pSurv, double kpatho, double spatho, double sigpatho, double eff, double repro, double *Tlat, double *Tspo, int Npatho, int Naggr, int **resistance, int *adaptation, double MGeff, double QReff, double *timeToQR, double taumut, double probSex, double costInfect, double costAggr, double beta) {

  fprintf(f, "###     MODEL PARAMETERS     ###\n");  
  fprintf(f, "seed:                      %d\n", seed);

  fprintf(f, "\n*****             Seasonality                  *****\n");
  fprintf(f, "nYears:     %d\n", nYears); 
  fprintf(f, "nDays/year: %d\n", nTSpY);
  fprintf(f, "pI0:        %.4f\n", pI0);
  fprintf(f, "Relative bottleneck size: %.0e\n", pSurv);
  
  fprintf(f, "\n*****     Landscape & deployment strategy     *****\n"); 
  fprintf(f, "Npoly:    %d\n", Npoly);
  fprintf(f, "Nhost:    %d\n", Nhost);
  fprintf(f, "Deployment strategy:    %s\n", strat);
  
  fprintf(f, "\n*****        Host        *****\n"); 
  print_d1(f, Nhost, C0,  "C0 [host]");
  print_d1(f, Nhost, Cmax,"Cmax [host]");
  print_d1(f, Nhost, growthH, "growthH [host]");
  print_d1(f, Nhost, reproH, "reproH [host]");

  fprintf(f, "\n*****        Pathogen      *****\n"); 
  fprintf(f, "eff:    %.2f     Tlat:   %.2f     repro:    %.2f      Tspo:    %.2f\n", eff, Tlat[0], repro, Tspo[0]);
  fprintf(f, "Tlat_var:    %.2f       Tspo_var:    %.2f\n", Tlat[1], Tspo[1]);
  fprintf(f, "k_patho:    %.2f       s_patho:    %.2f        sig_patho:     %.2f\n", kpatho, spatho, sigpatho);
  
  fprintf(f, "\n*****     Evolution       *****\n"); 
  fprintf(f, "Npatho:       %d\n", Npatho);
  fprintf(f, "probSex: %.2f    taumut: %.0e\n", probSex, taumut);
  print_i2(f, Nhost, NLOCI, resistance, "resistance [host][mg1, mg2, mg3, mg4, qr_E, qr_L, qr_R, qr_S]");
  print_i1(f, NLOCI, adaptation, "adaptation [mg1, mg2, mg3, mg4, qr_E, qr_L, qr_R, qr_S]");
  fprintf(f, "# Major-gene resistance\n");
  fprintf(f, "MGeff:         %.2f\n", MGeff);
  fprintf(f, "costInfect:  %.2f\n", costInfect);
  fprintf(f, "# Quantitative resistance\n");
  fprintf(f, "QReff:       %.2f\n", QReff);
  fprintf(f, "costAggr  :  %.2f\n", costAggr);
  fprintf(f, "Naggr:       %d\n", Naggr);
  fprintf(f, "beta:         %.2f\n", beta);
  fprintf(f, "timeToQR_exp: %.2f\n", timeToQR[0]);
  fprintf(f, "timeToQR_var: %.2f\n", timeToQR[1]);
  fprintf(f, "\n");
}


/* Write model output in .txt files and print output on screen */
void write_HHjuvSLIR(int Npoly, int Npatho, int Nhost, int t, int **H, int **Hjuv, int **S, int ***L, int ***I, int ***R, FILE *fH, FILE *fHjuv, FILE *fS, FILE *fL, FILE *fI, FILE *fR /*, int printOn*/) {
     
     int poly,patho,host;  
     for (poly = 0; poly < Npoly; poly++) {
          for  (patho = 0 ; patho < Npatho ; patho++){
               fwrite(&S[poly][patho],sizeof(int),1,fS);
               for (host = 0 ; host < Nhost ; host++){
                    fwrite(&L[poly][patho][host],sizeof(int),1,fL);
                    fwrite(&I[poly][patho][host],sizeof(int),1,fI);			
                    fwrite(&R[poly][patho][host],sizeof(int),1,fR); 			
               } /* for host */
          } /* for patho */
          for (host = 0 ; host < Nhost ; host++){
               fwrite(&Hjuv[poly][host],sizeof(int),1,fHjuv);
               fwrite(&H[poly][host],sizeof(int),1,fH);
          } /* for host */  
     } /* for poly */
     
     /* Print model output on screen */
     /*if (printOn) {
          printf("### t = %d\n", t);
          print_i2(stdout, Npoly,Npatho,S,"S [poly][patho]");
          print_i2(stdout, Npoly,Nhost,H,"H [poly][host]");
          print_i2(stdout, Npoly,Nhost,Hjuv,"Hjuv [poly][host]");
          print_i3sum1(stdout, Npoly,Npatho,Nhost,L,"L (sum poly)[patho][host]");
          print_i3sum1(stdout, Npoly,Npatho,Nhost,I,"I (sum poly)[patho][host]");
          print_i3sum1(stdout, Npoly,Npatho,Nhost,R,"R (sum poly)[patho][host]");  
     }*/
} 



