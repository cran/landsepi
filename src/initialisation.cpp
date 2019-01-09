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


/****************************************************************/
/*                    initialisation.c                          */
/****************************************************************/
#include "initialisation.hpp"

/* Initialisation of the conversion matrices between pathotype index and aggressiveness formula */
void init_aggrFormula(int Npatho, int Naggr, int *adaptation, int ********aggrToPatho, int **pathoToAggr) {
  
  int i,j,k,l,m,n,o,p;
  int line;
  int ag=0;    /* index of aggressiveness */
  int nIG[4], nAG[4];
  nIG[0]=adaptation[0];
  nAG[0]=adaptation[4];
  for (i=1; i<4; i++) {
       nIG[i] = nIG[i-1] + adaptation[i];
       nAG[i] = nAG[i-1] + adaptation[i+4];
  }

  /* Builging pathoToAggr: give the aggressiveness formula from pathotype index */
  /* Fill the columns corresponding to the genes that can evolve (given by the adaptation formula) */
  /* with a complete experimental design of the possible aggressiveness indices */
  /* and fill the other columns with 0 (meaning absence of infectivity or aggressiveness genes) */
  for (p = 0; p < pow(Naggr,adaptation[7]); p++) {
    for (o = 0; o < pow(Naggr,adaptation[6]); o++) {
      for (n = 0; n < pow(Naggr,adaptation[5]); n++) {
        for (m = 0; m < pow(Naggr,adaptation[4]); m++) {
          for (l = 0; l < pow(2,adaptation[3]); l++) {
            for (k = 0; k < pow(2,adaptation[2]); k++) { 
              for (j = 0; j < pow(2,adaptation[1]); j++) {
                for (i = 0; i < pow(2,adaptation[0]); i++) {
				line = i;
				line += j*pow(2,nIG[0]);
				line += k*pow(2,nIG[1]);
				line += l*pow(2,nIG[2]);
				line += m*pow(2,nIG[3]);
				line += n*pow(2,nIG[3])*pow(Naggr,nAG[0]);
				line += o*pow(2,nIG[3])*pow(Naggr,nAG[1]);
				line += p*pow(2,nIG[3])*pow(Naggr,nAG[2]);
				pathoToAggr[line][0] = i;
				pathoToAggr[line][1] = j;
				pathoToAggr[line][2] = k;
				pathoToAggr[line][3] = l;
				pathoToAggr[line][4] = m;
				pathoToAggr[line][5] = n;
				pathoToAggr[line][6] = o;
				pathoToAggr[line][7] = p;
				/* building aggrToPatho: give the pathotype index from aggressiveness formula */
				aggrToPatho[i][j][k][l][m][n][o][p] = ag++;
	  		 }
  	         }
  	       }
  	     }
        }
      }
    }
  }
  
  //if (Npatho < 101)
    //print_i2(stdout, Npatho, 8, pathoToAggr, (const char *)"pathoToAggr");
}


/* Initialisation of mutation matrices */
void init_mutkernel(double **mutkernelMG, double **mutkernelQR, double taumut, int Naggr) {
  int i,j;
  
  /* Mutation matrix for major genes */
  for (i=0;i<2;i++) {
  	mutkernelMG[i][i] = 1 - taumut;
  	mutkernelMG[i][(int) pow((i-1), 2)] = taumut;
  }
  //print_d2(stdout, 2,2,mutkernelMG,"mutkernelMG");
  
  /* Mutation matrix for quantitative resistance */
  for (i=0 ; i<Naggr ; i++){
    	for (j=0 ; j<Naggr ; j++)
          mutkernelQR[i][j] = 0.0;
    	mutkernelQR[i][i] = 1 - taumut;  /* Diagonal */
  } 
  mutkernelQR[0][1] = taumut;    /* because only one possible direction for mutation */
  mutkernelQR[Naggr-1][Naggr-2] = taumut;
  if (Naggr>2) {
	for (i=1;i<Naggr-1;i++) {   /* two possible directions */
	  mutkernelQR[i][i+1] = taumut/2;  
	  mutkernelQR[i][i-1] = taumut/2;
	}
  }  
  //print_d2(stdout, Naggr,Naggr,mutkernelQR,"mutkernelQR");
}


/* Initialisation of infectivity and aggressiveness matrices */
void init_infectAggr(int Naggr, double MGeff, double QReff, double costInfect, double costAggr, double beta, double **infect, double **aggr) {
  int i,j;
  /* Infectivity matrix */
  infect[0][0] = 1;
  infect[1][0] = 1 - costInfect;
  infect[0][1] = 1 - MGeff;
  infect[1][1] = 1;
  //print_d2(stdout,2,2,infect,"infectivity matrix [2][2]");
  
  /* Aggressiveness matrix */
  double aggr0=1-QReff;
  double step=1/(double) (Naggr-1);  /* (Naggr-1) intervals between Naggr values */
  double *gain; 
  double *cost;
  gain = init_d1(Naggr);
  cost = init_d1(Naggr);
  for (i=0;i<Naggr;i++)
	gain[i] = i*step;	
  tradeoff(Naggr,gain,cost,beta);
  
  for (i=0;i<Naggr;i++){
       aggr[i][1] = aggr0 + gain[i]*QReff;
       aggr[i][0] = 1 - cost[i]*costAggr;
  }
  //print_d2(stdout, Naggr,2,aggr,"aggressiveness matrix [aggr][2]");
  
  free(cost);
  free(gain);
}
 
 
 /* Initialisation of S at 0 */
void init_S(int Npoly, int Npatho, int **S) {
 int poly, patho;
 for (poly=0; poly<Npoly; poly++) {
      for (patho=0; patho<Npatho; patho++)
           S[poly][patho] = 0;
 } /* for poly */
}

/* Initialisation of H, L, I and R at 0 */
void init_HHjuvSLIR(int Npoly, int Nhost, int Npatho, int **H, int **Hjuv, int **S, int ***L, int ***I, int ***R) {
  int poly,host,patho;
  for (poly=0; poly<Npoly; poly++) {
   	for (host=0; host<Nhost; host++) {
  		H[poly][host] = 0; 
   	     Hjuv[poly][host] = 0;
		for (patho=0; patho<Npatho; patho++) {
  			L[poly][patho][host] = 0;
  			I[poly][patho][host] = 0;
  			R[poly][patho][host] = 0;
  		} /* for patho */
   	} /* for host */
     for (patho=0; patho<Npatho; patho++)
          S[poly][patho] = 0;
  } /* for poly */
}


/* Initialise L2I and I2R with 0 */
void init_L2I2R(int Npoly, int Npatho, int Nhost, int nTSpY, int ****L2I, int ****I2R){
     int poly, patho, host, t;
     for (poly=0; poly<Npoly; poly++) {
          for (patho=0; patho<Npatho; patho++) {
               for (host=0; host<Nhost; host++) {
                    for (t=0; t<nTSpY; t++) {
                         L2I[poly][patho][host][t] = 0;
                         I2R[poly][patho][host][t] = 0;
                    } /* for t */
               } /* for host */
          } /* for patho */
     } /* for poly */
}


/* Initialise SpathoMut at 0 */
void init_SpathoMut(int Npatho, int **SpathoMut) {
     int patho_old, patho_mut;
     for (patho_old = 0; patho_old < Npatho; patho_old++) {
          for (patho_mut = 0; patho_mut < Npatho; patho_mut++)
               SpathoMut[patho_old][patho_mut] = 0;
     }
}


/* Plantation of H in the beginning of a season */
void intro_H(int Npoly, int Nhost, int **H, int *area, int **habitat, double *C0, char *strat, int id_rotation) {
     int poly,host,patho;
     for (poly=0; poly<Npoly; poly++) {
          for (host=0; host<Nhost; host++) {
               H[poly][host] = area[poly] * C0[host] * (habitat[id_rotation][poly] == host);
               if (strcmp(strat,"MI") == 0 && host>0)
                    H[poly][host] += area[poly] * C0[host] * (habitat[1][poly] == host);
          } /* for host */
     } /* for poly */
}


/* Pathogen introduction : infectious sites from pathotype 0 in cultivar 0 */
void intro_I(const gsl_rng *gen, int nTSpY, int Npoly, int Nhost, int Npatho, int **H, int ***I, int ****I2R, double pI0, double *Tspo, int **resistance, double **infect, double **aggr, int ID_E, int ID_S) {
  int poly, i;
  int id_MG, ig, ag_E, ag_S, mg, qr_E, qr_S;
  double eff_exp, Tspo_exp;
  double Tspo_alpha[2];
  int lag;
  
  eff_exp = pI0;
  ig = 0;                               /* introduced pathogen has no infectivity gene */
  for (id_MG=0; id_MG<4; id_MG++){
     mg = resistance[0][id_MG];         /* indicate if the cultivar has a major gene */
     eff_exp *= infect[ig][mg];
  }
  ag_E = 0;                            /* introduced pathogen has no aggressiveness gene on infection rate */
  qr_E = resistance[0][ID_E];          /* indicate if the cultivar has a QR on infection rate */
  eff_exp *= aggr[ag_E][qr_E];

  ag_S = 0; 
  Tspo_exp = Tspo[0];
  qr_S = resistance[0][ID_S];      /* indicate if the cultivar has a QR on sportulation duration */
  Tspo_exp *= aggr[ag_S][qr_S];
  Tspo_exp += 0.01*(Tspo_exp == 0);   /* Security to avoid problem in alpha calculation */				
  find_paramGamma(Tspo_exp, Tspo[1], Tspo_alpha);
  
  for (poly=0; poly<Npoly; poly++) {
     I[poly][0][0] = gsl_ran_binomial(gen, eff_exp, H[poly][0]);
     H[poly][0] -= I[poly][0][0];
     
     /* Calculation of latent and infectious period for infected sites at t=0 */
     for (i=0; i<I[poly][0][0]; i++) {	      		
          lag = (int) gsl_ran_gamma(gen, Tspo_alpha[0], Tspo_alpha[1]);
          lag += 1*(lag == 0);            /* Security to avoid infectious duration of 0 day */
          if (lag < nTSpY)
               I2R[poly][0][0][lag]++;
     } /* for i */
  } /* for poly */
}


/* activeQR: matrix indicating for each year and field, at which time step QR starts to be active */
void init_activeQR(const gsl_rng *gen, int nYears, int Npoly, double *timeToQR, double *timeToQR_alpha, int **activeQR){
     int year, poly;
     if (timeToQR[0]>0 && timeToQR[1]>0){
          find_paramGamma(timeToQR[0], timeToQR[1], timeToQR_alpha);
          for (year=0; year<(nYears+1); year++){   /* nYears + 1 to allow last time step of the simulation */
               for (poly=0; poly<Npoly; poly++)
                    activeQR[year][poly] = (int) gsl_ran_gamma(gen, timeToQR_alpha[0], timeToQR_alpha[1]);
          }
      }else{
          for (year=0; year<(nYears+1); year++){
               for (poly=0; poly<Npoly; poly++)
                    activeQR[year][poly] = 0;
          }
     }
}
