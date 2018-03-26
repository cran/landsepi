/* 
 * Part of the landsepi R package.
 * Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@csiro.au>
 *                    Julien Papaix <julien.papaix@csiro.au>
 *                    Jean-François Rey <jean-francois.rey@inra.fr>
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





/* REMINDER                              */
/*  | is bitwise OR,  || is logical OR   */
/*  & is bitwise AND, && is logical AND  */
/*  (int) a|b is the euclidian quotient  */
/*  (int) a%b is the euclidian remainder */
/*  ++variable == variable++ == variable=variable+1 */
/*  ++*pointer == (*pointer)++ == *pointer=*pointer+1    CAREFULL   *pointer++   DOES NOT WORK */

/* CULTIVAR GENOTYPE : mg are major resistant genes, qr are traits of quantitative resistance  */
/* PATHOGEN GENOTYPE : ig are infectivity genes,     ag are aggressiveness genes               */
/* Cultivar 0 is susceptible (no mg, no qr), Cultivars > 0 have one or more resistance sources */
/* Pathogen 0 is avirulent   (no ig, no ag), Pathogen > 0 have some ig or ag                   */

#include "model_simul.hpp"



/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/


/* ---------------------------------- */
/*         MUTATION OF THE SPORES     */
/* ---------------------------------- */  
/* Update SpathoMut with mutation on "trait_mut" through multinomial draw */
void mutation(const gsl_rng *gen, int patho, int Npatho, int Naggr, int trait_mut, double **mutkernel, int **pathoToAggr, int ********aggrToPatho, int **SpathoMut) {
     
     int patho_old, aggr_mut;//patho_mut
     int ig0, ig1, ig2, ig3, agE, agL, agR, agS;  /* index of the infectivity genes (ig) or aggressiveness genes (ag) */
     int id_aggr_old, id_patho_mut;	
     int **SaggrMut;  	
     SaggrMut = init_i2(Npatho, Naggr);

     /* mutation of trait_mut */
     for (patho_old = 0; patho_old < Npatho; patho_old++){
          id_aggr_old = pathoToAggr[patho_old][trait_mut];/* Aggressiveness index before mutation */
          gsl_ran_multinomial (gen, Naggr, SpathoMut[patho][patho_old], mutkernel[id_aggr_old], (unsigned int*) SaggrMut[patho_old]);
          SpathoMut[patho][patho_old] = 0;  /* re-initialisation of SpathoMut */
     } /* for patho_old */

     
//     for (patho_mut = 0; patho_mut < Npatho; patho_mut++)
//          SpathoMut[patho][patho_mut] = 0;
     
     /* update SpathoMut with mutants */  
     for (patho_old = 0; patho_old < Npatho; patho_old++) {
          /* Aggressiveness index relative to the different traits */
          ig0 = pathoToAggr[patho_old][0];	
          ig1 = pathoToAggr[patho_old][1];
          ig2 = pathoToAggr[patho_old][2];
          ig3 = pathoToAggr[patho_old][3];
          agE = pathoToAggr[patho_old][4];
          agL = pathoToAggr[patho_old][5];
          agR = pathoToAggr[patho_old][6];
          agS = pathoToAggr[patho_old][7];
          for (aggr_mut = 0; aggr_mut < Naggr; aggr_mut++) {
               /* Pathotype index after mutation */
               if (trait_mut == 0)
                    id_patho_mut = aggrToPatho[aggr_mut][ig1][ig2][ig3][agE][agL][agR][agS];  
               else if (trait_mut == 1)
                    id_patho_mut = aggrToPatho[ig0][aggr_mut][ig2][ig3][agE][agL][agR][agS];	
               else if (trait_mut == 2)
                    id_patho_mut = aggrToPatho[ig0][ig1][aggr_mut][ig3][agE][agL][agR][agS];  
               else if (trait_mut == 3)
                    id_patho_mut = aggrToPatho[ig0][ig1][ig2][aggr_mut][agE][agL][agR][agS];
               else if (trait_mut == 4)
                    id_patho_mut = aggrToPatho[ig0][ig1][ig2][ig3][aggr_mut][agL][agR][agS];  
               else if (trait_mut == 5)
                    id_patho_mut = aggrToPatho[ig0][ig1][ig2][ig3][agE][aggr_mut][agR][agS];	
               else if (trait_mut == 6)
                    id_patho_mut = aggrToPatho[ig0][ig1][ig2][ig3][agE][agL][aggr_mut][agS];  
               else if (trait_mut == 7)
                    id_patho_mut = aggrToPatho[ig0][ig1][ig2][ig3][agE][agL][agR][aggr_mut];		  
               
               SpathoMut[patho][id_patho_mut] += SaggrMut[patho_old][aggr_mut];  /* Add in index of new pathotype */
          } /* for aggr_mut */
     } /* for patho_old */	
          
     //	print_i2(stdout, Npatho, Npatho, SpathoMut, "SpathoMut");	
     free_i2(SaggrMut, Npatho); 

/* Old code for mutation */
/* --------------------- */
/* mutation of the 1st trait */
//	id_aggr1 = pathoToAggr[patho][trait1];  /* Aggressiveness index of the 1st adaptative trait */
//	gsl_ran_multinomial (gen, Naggr, Sprod,mutkernelQR[id_aggr1],SaggrMut1);

/* mutation of the 2nd trait */
//	id_aggr2 = pathoToAggr[patho][trait2];  /* Aggressiveness index of the 2nd adaptative trait */
//	for (aggr1=0;aggr1<Naggr;aggr1++) {
//		gsl_ran_multinomial (gen, Naggr, SaggrMut1[aggr1],mutkernelQR[id_aggr2],SaggrMut2[aggr1]);
/* Find the pathotype associated with the new traits 1 and 2 */
//		for (aggr2=0;aggr2<Naggr;aggr2++) {
//			id_pathomut = aggrToPatho[aggr1][aggr2];
//			SpathoMut[patho][id_pathomut] = SaggrMut2[aggr1][aggr2];
//		} /* for aggr2 */
//	} /* for aggr1 */

//	printf("*** poly=%d , patho=%d, Sprod_tmp=%d, Sprod=%d, id_aggr1=%d, id_aggr2=%d\n", poly, patho, Sprod_tmp, Sprod, id_aggr1, id_aggr2);
//	print_i1(stdout, Naggr,SaggrMut1,"SaggrMut1");
//	print_i2(stdout, Naggr,Naggr,SaggrMut2,"SaggrMut2");
}


/* -------------------------------------- */
/*        PRODUCTION AND DISPERSAL        */
/* -------------------------------------- */
/* Production, mutation and dispersal of spores */
/* Production and dispersal of new hosts */
void reproDisp(const gsl_rng *gen, int Npoly, int Nhote, int Npatho, int Naggr, int **H, int **Hjuv, int **S, int ***I, double *reproH, double repro, double **disphote, double **disppatho, int **resistance, int *adaptation, double **mutkernelMG,  double **mutkernelQR, int **pathoToAggr, int ********aggrToPatho, double **aggr){
     /* H, Hjuv, S and I are the numbers of individuals in a given poly */ 
     int poly, polyE, patho, patho_mut, hote, id_IG, id_AG;
     int ag_R, qr_R;
     
     double Sprod_tmp;
     double Sprod_exp;
     int Sprod;
     int **SpathoMut;
     int **Smut;
     int ***Sdisp;
     int ***Hjuvtmp;
     SpathoMut = init_i2(Npatho,Npatho);  
     Smut = init_i2(Npoly,Npatho);
     Sdisp = init_i3(Npatho,Npoly,Npoly);      
     Hjuvtmp = init_i3(Nhote,Npoly,Npoly); 
     
     /* Production and mutation of spores */
     for (poly = 0 ; poly < Npoly ; poly++){
          /* Initialisation SpathoMut at 0 */		
          init_SpathoMut(Npatho, SpathoMut);
          
          for (patho = 0 ; patho < Npatho ; patho++){
               /* Poisson draw of the number of produced spores, depending on the effect of aggressiveness on repro */	
               Sprod_exp = 0.0;
               for (hote = 0 ; hote < Nhote ; hote++) {
                    Sprod_tmp = repro * I[poly][patho][hote];
                    ag_R = pathoToAggr[patho][ID_R];        /* indicate if the pathogen has a aggressiveness gene on reproduction rate */
                    qr_R = resistance[hote][ID_R];          /* indicate if the cultivar has a QR on reproduction rate */
               
                    Sprod_tmp *= aggr[ag_R][qr_R];
                    Sprod_exp += Sprod_tmp;           /* expectation of the number of produced spores */
               }	
               Sprod = gsl_ran_poisson(gen, Sprod_exp);	
               SpathoMut[patho][patho] = Sprod;
               
               /* Mutation of the different traits listed in the adaptation formula */
               for (id_IG = 0; id_IG < 4; id_IG++) {
                    if (adaptation[id_IG])
                         mutation(gen, patho, Npatho, 2, id_IG, mutkernelMG, pathoToAggr, aggrToPatho, SpathoMut);
               }
               for (id_AG = 4; id_AG < 8; id_AG++) {
                    if (adaptation[id_AG])
                         mutation(gen, patho, Npatho, Naggr, id_AG, mutkernelQR, pathoToAggr, aggrToPatho, SpathoMut); 
               }
          } /* for patho */
               
          for (patho_mut = 0; patho_mut < Npatho; patho_mut++) {
               Smut[poly][patho_mut] = 0;
               for (patho = 0; patho<Npatho; patho++)
                    Smut[poly][patho_mut] += SpathoMut[patho][patho_mut];
          } /* for patho_mut */     
     } /* for poly */
               
     /* Production and dispersal of the host */
     /* Dispersal of spores */
     for (poly = 0 ; poly < Npoly ; poly++){
          /* Spore dispersal */
          for (patho = 0 ; patho < Npatho ; patho++)    
               gsl_ran_multinomial (gen, Npoly, Smut[poly][patho], disppatho[poly], (unsigned int*)Sdisp[patho][poly]);
          /* host reproduction: production and dispersal of Hjuv */
          for (hote = 0 ; hote < Nhote ; hote++)
               gsl_ran_multinomial (gen, Npoly, (int) (reproH[hote] * H[poly][hote]), disphote[poly], (unsigned int*)Hjuvtmp[hote][poly]);
     } /* for poly */

     for (poly = 0 ; poly < Npoly ; poly++){
          /* Number of spores (S) and Hjuv landing in each field */
          for (patho = 0 ; patho < Npatho ; patho++){
               S[poly][patho] = 0;
               for (polyE = 0 ; polyE < Npoly ; polyE++)
                    S[poly][patho] += Sdisp[patho][polyE][poly];	
          } /* for patho */
          
          for (hote = 0 ; hote < Nhote ; hote++){
               Hjuv[poly][hote] = 0;
               for (polyE = 0 ; polyE < Npoly ; polyE++)
                    Hjuv[poly][hote] += Hjuvtmp[hote][polyE][poly];
          } /* for hote */
     } /* for poly */
                    
//     printf("# Mutation and dispersal\n");
//     print_i2(stdout, Npoly,Npatho,Smut,"Smut (after mutation, before dispersal)");
//	print_i3(stdout, Npatho,Npoly,Npoly,Sdisp,"Sdisp");

     free_i3(Hjuvtmp,Nhote,Npoly);  
     free_i2(SpathoMut,Npatho);  
     free_i2(Smut,Npoly);
     free_i3(Sdisp,Npatho,Npoly);
}          


/* --------------------------------------------- */
/*       BOTTLENECK AT THE END OF THE SEASON     */
/* --------------------------------------------- */
void bottleneck(const gsl_rng *gen, int Npoly, int Nhote, int Npatho, double pSurv, double *Tspo, int ***L, int ***I, double **aggr, int **resistance, int **pathoToAggr, int ***eqIsurv, int printOn){
     int poly, patho, hote;
     int Tspo_exp;
     int ag_S, qr_S;
     
     for (patho = 0 ; patho < Npatho ; patho++){
          for (hote = 0 ; hote < Nhote ; hote++){
               /* calculate the mean sporulation period */
               ag_S = pathoToAggr[patho][ID_S];        /* indicate if the pathogen has a aggressiveness gene on sporulation duration */
               qr_S = resistance[hote][ID_S];          /* indicate if the cultivar has a QR on sporulation duration */
               Tspo_exp = (int) (Tspo[0] * aggr[ag_S][qr_S]); 
               
               for (poly = 0; poly < Npoly; poly++){
                    /* Reduce the number of infected hosts (bottleneck) */
                    eqIsurv[poly][patho][hote] = gsl_ran_binomial(gen, pSurv, L[poly][patho][hote] + I[poly][patho][hote]); 
                    /* Calculate the equivalent number of infectious hosts */
                    eqIsurv[poly][patho][hote] *= Tspo_exp;
               } /* for poly */
          } /* for hote */
     } /* for patho */
     //if (printOn)
     //     print_i3sum1(stdout, Npoly,Npatho,Nhote,eqIsurv,"# Bottleneck: eqIsurv (sum poly)[patho][hote]");
}


/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/* The following functions are nested in a { for poly } */

/* -------------------------- */
/*         HOST DYNAMIC       */
/* -------------------------- */
/* Compute host reproduction, death and growth and updtate the number of H in the concerned poly */
void host_dynamic(const gsl_rng *gen, int Nhote, int Npatho, int area, double *Cmax, int *H, int *Hjuv, int **L, int **I, int **R, int *N, double *mortH, double *delta, double khote, double sighote, double shote){ 
  /* H, Hjuv, L, I and R are the number of host in a given poly */    
  int patho, hote;
  double f1hote;
  int availSites, siteaccess;
  int K[Nhote];
  int L_hote[Nhote];
  int I_hote[Nhote];
  int R_hote[Nhote];
  int H2M, H2Mjuv, Hnewind, H2H;

  for (hote = 0 ; hote < Nhote ; hote++){
       
       /* Calculation of totals for L, I, R and K */
       K[hote] = (int) (area * Cmax[hote]);    /* carrying capacity of the cultivar in the concerned paddock */
       L_hote[hote] = 0;
       I_hote[hote] = 0;
       R_hote[hote] = 0;
       for (patho = 0 ; patho < Npatho ; patho++){
            L_hote[hote] += L[patho][hote];
            I_hote[hote] += I[patho][hote];
            R_hote[hote] += R[patho][hote];
       }
       N[hote] = H[hote] + L_hote[hote] + I_hote[hote] + R_hote[hote];        /* Number of occupied sites */
       
       /* HOST MORTALITY: H2M */
       /* ------------------- */
       H2M = gsl_ran_binomial(gen, (double) mortH[hote], H[hote]);
       H2Mjuv = gsl_ran_binomial(gen, (double) mortH[hote], Hjuv[hote]);    
       Hjuv[hote] -= H2Mjuv;   /* update Hjuv */
       
       /* HOST REPRODUCTION: Hnewind */
       /* -------------------------- */
       /* Hjuv settlement in the field */
       availSites = K[hote] - N[hote] - H2M;                      /* Number of available sites */
       if (availSites < 0)  /* security */
            availSites ==0;
       sigmoid((double) shote, (double) khote, (double) sighote, availSites / (double) K[hote], &f1hote);
       siteaccess = gsl_ran_binomial(gen, f1hote, availSites);
       if (siteaccess < Hjuv[hote])
            Hnewind = siteaccess;
       else
            Hnewind = Hjuv[hote];
       
       /* HOST GROWTH: H2H */
       /* ---------------- */
       H2H = delta[hote] * (H[hote] - H2M) * (1 - ((N[hote] - H2M + Hnewind)/(double) K[hote]));
       if (H2H < 0) {   /* security */
            Rprintf("CAREFUL !!!!!  H2H < 0   one of the areas may be 0: check if Npoly, NpolyTot and idLAN are correct\n");
            H2H = 0;
       } else if ((N[hote] - H2M + Hnewind + H2H) > K[hote]) {
            Rprintf("CAREFUL !!!!!  H2H too big\n");
            H2H = K[hote] - (N[hote] - H2M + Hnewind);
       }
       
       /* UPDATE NUMBER OF HOSTS */
       /* ---------------------- */
       H[hote] = H[hote] - H2M + Hnewind + H2H;
       N[hote] = N[hote] - H2M + Hnewind + H2H;
       
       /*  printf("# Host dynamic: death, reproduction and growth of host %d\n", hote);
        printf("death : H2M = %d   H2Mjuv = %d\n", H2M, H2Mjuv);
        printf("repro : Hnewind = %d\n", Hnewind);
        printf("growth: H2H = %d\n", H2H);
        printf("final H: H = %d\n", H[hote]);
        */       
  } /* for hote */
}


/* -------------------------------------------------------------- */
/*         CONTAMINATION : spore deposition on healthy sites      */
/* -------------------------------------------------------------- */
/* Calculation of the number of contaminated sites */
void contamination(const gsl_rng *gen, int Nhote, int Npatho, int *H, int *S, int *N, int **Hcontaminated, double kpatho, double sigpatho, double spatho){ 
     /* H, S are the number of individuals in a given poly */ 
     int patho, hote;
     double f1patho;
     int Stot=0;
     int Htot=0;
     int S_hote[Nhote+1];
     int *S_hote_patho;
     double probaH[Nhote+1];
     double *probaP;
     double probaHtot=0;
     double probaPtot=0;
     int Hcontaminable;        /* Contaminable site, where a spore may deposit */
     int *Hcontaminated_tmp;
     probaP = init_d1(Npatho+1);
     S_hote_patho = init_i1(Npatho+1);
     Hcontaminated_tmp = init_i1(Npatho+1);

 /* Calculation of total for H and S */
 for (hote = 0 ; hote < Nhote ; hote++)
      Htot += H[hote];
 for (patho = 0 ; patho < Npatho ; patho++)
      Stot += S[patho];

  /* Probability for H to belong to each host */
  if (Htot == 0){
       for (hote=0; hote<Nhote; hote++)
            probaH[hote] = 0;
  }
  else {
       for (hote=0; hote<Nhote; hote++)
            probaH[hote] = (double) H[hote] /(double) Htot;
  }
  for (hote=0 ; hote<Nhote; hote++)
       probaHtot += probaH[hote];
  probaH[Nhote] = 1 - probaHtot;   
  
  /* Probability to belong to the different pathotypes */  
  if (Stot == 0){  //si pas de spores on met tout à 0
       for (patho=0; patho<Npatho; patho++)
            probaP[patho] = 0;
  }
  else {          //sinon on met la proportion de chacun
       for (patho=0; patho<Npatho; patho++)
            probaP[patho] = (double) S[patho] /(double) Stot;
  }  
  for (patho=0 ; patho<Npatho; patho++)
       probaPtot += probaP[patho];    
  probaP[Npatho] = 1 - probaPtot;   
  
  /* distribution of the spores among the different cultivars */
  gsl_ran_multinomial (gen, Nhote+1, Stot, probaH, (unsigned int*)S_hote); 
  
  for (hote=0; hote<Nhote; hote++){
     /* distribution of the spores among the different pathotypes */ 
     gsl_ran_multinomial (gen, Npatho+1, S_hote[hote], probaP, (unsigned int*)S_hote_patho); 
       
     /* Calculation of the number of contaminable sites */
     if (N[hote] > 0)
          sigmoid((double) spatho, (double) kpatho, (double) sigpatho, (H[hote]/(double) N[hote]), &f1patho);
     else
          f1patho = 0.0;
     Hcontaminable = gsl_ran_binomial (gen, f1patho, H[hote]);
     
     /* distribution of the contaminable sites among the different pathotypes */ 
     gsl_ran_multinomial (gen, Npatho+1, Hcontaminable, probaP, (unsigned int*)Hcontaminated_tmp);
     
     /* The true number of contaminated sites is the minimum between sites and spores */
     for (patho=0; patho<Npatho; patho++){
          if (Hcontaminated_tmp[patho] < S_hote_patho[patho])
               Hcontaminated[patho][hote] = Hcontaminated_tmp[patho];
          else
               Hcontaminated[patho][hote] = S_hote_patho[patho];
     } /* for patho */     
  } /* for hote */
     
  //printf("# Contamination and infection of healthy hosts\n");
  //print_i2(stdout, Npatho, Nhote, Hcontaminated,"Hcontaminated [hote][patho]");

 free(probaP);
 free(S_hote_patho);
 free(Hcontaminated_tmp);

}  

       
/* ----------------------------------------------------------- */
/*         INFECTIOUS CYCLE : transitions H -> L -> I -> R     */
/* ----------------------------------------------------------- */    
/* Calculate the number of contaminated sites that become infected, infectious or removed and update H, L, I, R */
void infection(const gsl_rng *gen, int t, int nTSpY, int Nhote, int Npatho, int *H, int **Hcontaminated, int **L, int **I, int **R, int ***L2I, int ***I2R, double eff, double *Tlat, double *Tspo, int **resistance, int **pathoToAggr, double **infect, double **aggr){ 
  /* H, Hcontaminated, L, I, R, L2I and I2R are the number of individuals in a given poly */   
  int patho, hote, h2l, id_MG;
  int ig, mg, ag_E, ag_L, ag_S, qr_E, qr_L, qr_S;
  double eff_exp, Tlat_exp, Tspo_exp;
  int H2L;
  double Tlat_alpha[2];
  double Tspo_alpha[2];
  int lag1, lag2;
  
  for (patho = 0 ; patho < Npatho ; patho++){
    for (hote = 0 ; hote < Nhote ; hote++){
         
      /* infection of healthy sites: H2L */
      /* ------------------------------- */
      eff_exp = eff;
      /* Interaction major genes - Infectivity genes */
      for (id_MG=0; id_MG<4; id_MG++) {
        ig = pathoToAggr[patho][id_MG];        /* indicate if the pathogen has an infectivity gene */
        mg = resistance[hote][id_MG];          /* indicate if the cultivar has an major gene */
        eff_exp *= infect[ig][mg];
      } /* for id_MG */
      /* Interaction Quantitative resistance - Aggressiveness gene */
      ag_E = pathoToAggr[patho][ID_E];        /* indicate if the pathogen has a aggressiveness gene on infection efficiency */
      qr_E = resistance[hote][ID_E];          /* indicate if the cultivar has a QR on infection efficiency */
      eff_exp *= aggr[ag_E][qr_E];
      H2L = gsl_ran_binomial(gen, eff_exp, Hcontaminated[patho][hote]);
     
      /* Latent and infectious periods */
      /* ----------------------------- */
      /* find parameters of gamma distributions from mean and variance */     	    	
      Tlat_exp = Tlat[0];
      ag_L = pathoToAggr[patho][ID_L];        /* indicate if the pathogen has a aggressiveness gene on latent period */
      qr_L = resistance[hote][ID_L];          /* indicate if the cultivar has a QR on reproduction rate */
      Tlat_exp /= (aggr[ag_L][qr_L] + 0.001*(aggr[ag_L][qr_L] == 0));  /* security to avoid division by 0 */
      Tlat_exp += 0.001*(Tlat_exp == 0);       /* Security to avoid problem in alpha calculation */
      find_paramGamma(Tlat_exp, Tlat[1], Tlat_alpha);
      
      Tspo_exp = Tspo[0];
      ag_S = pathoToAggr[patho][ID_S];        /* indicate if the pathogen has a aggressiveness gene on sporulation duration */
      qr_S = resistance[hote][ID_S];          /* indicate if the cultivar has a QR on sporulation duration */
      Tspo_exp *= aggr[ag_S][qr_S];
      Tspo_exp += 0.001*(Tspo_exp == 0);       /* Security to avoid problem in alpha calculation */				
      find_paramGamma(Tspo_exp, Tspo[1], Tspo_alpha);
      
      for (h2l=0; h2l<H2L; h2l++) { /* recently infected hosts */
          /* Latent period */		
          lag1 = (int) gsl_ran_gamma(gen, Tlat_alpha[0], Tlat_alpha[1]);
           //printf("lat = %d (mean=%.2f)  ", lag1, Tlat[0] * (2 - aggr[ag_L][qr_L]));
           if ((t+lag1) < nTSpY)     	
                L2I[patho][hote][t+lag1]++;
           
           /* Sporulation period */		       		
           lag2 = lag1 + (int) gsl_ran_gamma(gen, Tspo_alpha[0], Tspo_alpha[1]);	
           //printf("spo = %d (mean=%.2f)  ", lag2, Tspo[0] * aggr[ag_S][qr_S]);
           if ((t+lag2) < nTSpY)
                I2R[patho][hote][t+lag2]++;
      } /* for h2l */
        
      /* Update H, L, I, R */
      /* ----------------- */ 
      H[hote]        -= H2L;
      L[patho][hote] += H2L;
      L[patho][hote] -= L2I[patho][hote][t];
      I[patho][hote] += L2I[patho][hote][t];
      I[patho][hote] -= I2R[patho][hote][t];
      R[patho][hote] += I2R[patho][hote][t];  

    } /* for patho */
  } /* for hote */

//print_i3(stdout, Npatho, Nhote, nTSpY, L2I, "L2I(poly) [Npatho][Nhote][t]");
//print_i3(stdout, Npatho, Nhote, nTSpY, I2R, "I2R(poly) [Npatho][Nhote][t]"); 
}


 
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
 
/* --------------------------- */
/*   DYNAMIC OF THE EPIDEMIC   */
/* --------------------------- */

void dynepi(const gsl_rng *gen, int nYears, int nTSpY, int Npoly, int Nhote, int Npatho, int *area, int **habitat, int *rotation, double *C0, double pI0, double pSurv, double *Cmax, double **mort, double *delta, double *reproH, double **disphote, double **disppatho, double eff, double *Tlat, double repro, double *Tspo, double **mutkernelMG,  double **mutkernelQR, int **pathoToAggr, int ********aggrToPatho, int Naggr, double **infect, double **aggr, char *strat, int **resistance, int *adaptation, int printOn, OGRLayer *poLayer, double khote, double sighote, double shote, double kpatho, double sigpatho, double spatho){//, NumericVector timesStep,

  int year, t, poly;
  char name_fH[15], name_fHjuv[15], name_fS[15], name_fL[15], name_fI[15], name_fR[15];
  int **S;
  int **H;
  int **Hjuv;
  int **Hcontaminated;      /* Contaminated sites (where a spore is deposited) */
  int ***L;
  int ***I;
  int ***R;
  int N[Nhote];
  int ****L2I;
  int ****I2R;
  int ***eqIsurv;    /* equivalent number of I that survive and produce spores for the next season */ 

  S = init_i2(Npoly,Npatho); 
  H = init_i2(Npoly,Nhote);
  Hjuv = init_i2(Npoly,Nhote);
  Hcontaminated = init_i2(Npatho,Nhote);
  L = init_i3(Npoly,Npatho,Nhote);
  I = init_i3(Npoly,Npatho,Nhote);
  R = init_i3(Npoly,Npatho,Nhote);
  L2I = init_i4(Npoly,Npatho,Nhote,nTSpY);
  I2R = init_i4(Npoly,Npatho,Nhote,nTSpY); 
  eqIsurv = init_i3(Npoly, Npatho, Nhote);
  /* Initialisation (t=0) */ 
  init_HHjuvSLIR(Npoly, Nhote, Npatho, H, Hjuv, S, L, I, R);
  init_L2I2R(Npoly, Npatho, Nhote, nTSpY, L2I, I2R);
  intro_H(Npoly, Nhote, H, area, habitat, C0, strat, rotation[0]); 
  intro_I(gen, nTSpY, Npoly, Nhote, Npatho, H, I, I2R, pI0, Tspo, resistance, infect, aggr, ID_E, ID_S);
       
  for (year = 1; year < (nYears+1); year++) {
  	Rprintf("----------------------------- YEAR %d -----------------------------\n",year);
	/* Create the files to write the output */
	FILE *fH;
	FILE *fHjuv;
	FILE *fL;
	FILE *fI;
	FILE *fS;
	FILE *fR;

	sprintf(name_fH, "H-%02d.bin", year);
	sprintf(name_fHjuv, "Hjuv-%02d.bin", year);
	sprintf(name_fS, "S-%02d.bin", year);
	sprintf(name_fL, "L-%02d.bin", year);
	sprintf(name_fI, "I-%02d.bin", year);
	sprintf(name_fR, "R-%02d.bin", year);

	fH = fopen(name_fH,"wb");
	fHjuv = fopen(name_fHjuv,"wb");
	fL = fopen(name_fL,"wb");
	fI = fopen(name_fI,"wb");
	fS = fopen(name_fS,"wb");
	fR = fopen(name_fR,"wb"); 

	//int INDtime=0;
	/* Loop for all the timesteps of the cropping season */
	for (t = 0 ; t < nTSpY ; t++){
	     
	     /* Writing model output for timestep t */
	     //if(t==timesStep[INDtime]){
	       //sortie_layer(poLayer,H,Npoly,Nhote,t,year);
	       write_HHjuvSLIR(Npoly, Npatho, Nhote, t, H, Hjuv, S, L, I, R, fH, fHjuv, fS, fL, fI, fR, printOn);
	       //INDtime=INDtime+1;
	       //}	     
	     reproDisp(gen, Npoly, Nhote, Npatho, Naggr, H, Hjuv, S, I, reproH, repro, disphote, disppatho, resistance, adaptation, mutkernelMG, mutkernelQR, pathoToAggr, aggrToPatho, aggr);
	     for (poly = 0 ; poly < Npoly ; poly++){     
	          //      printf("poly  %d\n", poly);
	       host_dynamic(gen, Nhote, Npatho, area[poly], Cmax, H[poly], Hjuv[poly], L[poly], I[poly], R[poly], N, mort[poly], delta, khote, sighote, shote);
	       contamination(gen, Nhote, Npatho, H[poly], S[poly], N, Hcontaminated, kpatho, sigpatho, spatho);
	          infection(gen, t, nTSpY, Nhote, Npatho, H[poly], Hcontaminated, L[poly], I[poly], R[poly], L2I[poly], I2R[poly], eff, Tlat, Tspo, resistance, pathoToAggr, infect, aggr);
	     } /* for poly */

	}  /* for t */

     /* last time-step of the season: bottleneck before starting a new season */
     //sortie_layer(poLayer,H,Npoly,Nhote,nTSpY,year);	
     write_HHjuvSLIR(Npoly, Npatho, Nhote, nTSpY, H, Hjuv, S, L, I, R, fH, fHjuv, fS, fL, fI, fR, printOn);  /* Writing model output for last timestep */
     bottleneck(gen, Npoly, Nhote, Npatho, pSurv, Tspo, L, I, aggr, resistance, pathoToAggr, eqIsurv, printOn);     /* Calculation of the equivalent number of I */
     init_HHjuvSLIR(Npoly, Nhote, Npatho, H, Hjuv, S, L, I, R);                                              /* Re-initialisation at 0 */
     init_L2I2R(Npoly, Npatho, Nhote, nTSpY, L2I, I2R);                                                      /* Re-initialisation at 0 */
    
     /* Generate S issued from eqIsurv = (remaining L+I) * Tspo */
     reproDisp(gen, Npoly, Nhote, Npatho, Naggr, H, Hjuv, S, eqIsurv, reproH, repro, disphote, disppatho, resistance, adaptation, mutkernelMG, mutkernelQR, pathoToAggr, aggrToPatho, aggr);
     /* re-plantation --> regenerate H */
     intro_H(Npoly, Nhote, H, area, habitat, C0, strat, rotation[year]);
     /* infection of newly planted hosts to generate the inoculum of the next season */
     for (poly = 0 ; poly < Npoly ; poly++){ 
          contamination(gen, Nhote, Npatho, H[poly], S[poly], H[poly], Hcontaminated, kpatho, sigpatho, spatho);  /* N = H[poly] in beginning of next season */
          infection(gen, 0, nTSpY, Nhote, Npatho, H[poly], Hcontaminated, L[poly], I[poly], R[poly], L2I[poly], I2R[poly], eff, Tlat, Tspo, resistance, pathoToAggr, infect, aggr);
     }
      
	fclose(fH);
	fclose(fHjuv);
	fclose(fL);
	fclose(fI);
	fclose(fS);
	fclose(fR);	    
  } /* for year */ 

  free_i2(Hcontaminated,Npatho);  
  free_i4(L2I, Npoly, Npatho, Nhote);
  free_i4(I2R, Npoly, Npatho, Nhote);	    	    
  free_i2(S,Npoly);
  free_i2(H,Npoly);
  free_i2(Hjuv,Npoly); 
  free_i3(L,Npoly,Npatho);
  free_i3(I,Npoly,Npatho); 
  free_i3(R,Npoly,Npatho);
  free_i3(eqIsurv,Npoly,Npatho);
  
}

/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/

void sortie_layer(OGRLayer *poLayer,int **H,int Npoly, int Nhote, int pastemps, int year){

  char nomH[256];
  double *resH;
  resH=init_d1(Npoly);
     
  for (int i=0; i<Nhote; i++){
    for (int j=0; j<Npoly; j++){
      resH[j] = H[j][i];
    }
    sprintf(nomH, "H%i_%i-%i", i, year, pastemps);
    addField(poLayer, nomH, resH);
  }
  
  free(resH);

}



/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/

void modelLandsEPI(Rcpp::List times, Rcpp::List landscape, Rcpp::List dispersal, Rcpp::List inits ,int val_seed, Rcpp::List hostP, Rcpp::List epiP, Rcpp::List evolP) {


  int nYears = Rcpp::as<int>(times["nYears"]);
  int nTSpY = Rcpp::as<int>(times["nTSpY"]);
  //NumericVector timesStep=Rcpp::as<NumericVector>(times["timesStep"]);

  CharacterVector shapefilename = Rcpp::as<CharacterVector>(landscape["shapefilename"]);
  CharacterVector layername_hab = Rcpp::as<CharacterVector>(landscape["layername_hab"]);
  CharacterVector layername_res = Rcpp::as<CharacterVector>(landscape["layername_res"]);
  char * strat = (char*)Rcpp::as<std::string>(Rcpp::as<CharacterVector>(landscape["strat"])).c_str();
  IntegerVector rotation_tmp = Rcpp::as<IntegerVector>(landscape["rotation"]);
  double propRR = Rcpp::as<double>(landscape["propRR"]);
  double CMAX0 = Rcpp::as<double>(landscape["Cmax0"]);
  double CMAX1 = Rcpp::as<double>(landscape["Cmax1"]);

  //fprintf(stderr,"\ntype rotation: %s\n",typeid(rotation[0]).name());

  int * rotation=init_i1(nYears+1);
  IntegerVector::iterator it;
  int roti;
  for (roti=0, it = rotation_tmp.begin() ; it != rotation_tmp.end(); ++it, roti++){
    rotation[roti]=*it;
  }

  //CharacterVector nomfile_disppatho = Rcpp::as<CharacterVector>(dispersal["dispP"]);
  //CharacterVector nomfile_disphote = Rcpp::as<CharacterVector>(dispersal["dispH"]);
  NumericVector dispP_tmp=Rcpp::as<NumericVector>(dispersal["dispP"]);
  NumericVector dispH_tmp=Rcpp::as<NumericVector>(dispersal["dispH"]);

  double C_0 = Rcpp::as<double>(inits["C_0"]);
  double PI0 = Rcpp::as<double>(inits["PI0"]);

  int Nhote = Rcpp::as<int>(hostP["Nhote"]);
  double CROISH0 = Rcpp::as<double>(hostP["croisH0"]);
  double CROISH1 = Rcpp::as<double>(hostP["croisH1"]);
  double REPROH0 = Rcpp::as<double>(hostP["reproH0"]);
  double REPROH1 = Rcpp::as<double>(hostP["reproH1"]);
  double MORTH = Rcpp::as<double>(hostP["deathH"]);
  double khote = Rcpp::as<double>(hostP["khost"]);
  double sighote = Rcpp::as<double>(hostP["sighost"]);
  double shote = Rcpp::as<double>(hostP["shost"]);
  IntegerVector resistance_tmp = Rcpp::as<IntegerVector>(hostP["resistance"]);
  
  int **resistance;  
  resistance = init_i2(Nhote,NTRAIT);
  for (int ihote=0; ihote<Nhote; ihote++){
    for (int itrait=0; itrait<NTRAIT; itrait++){
      resistance[ihote][itrait] = resistance_tmp[itrait+ihote*NTRAIT];
    }
  }

  double PSURV = Rcpp::as<double>(epiP["psurv"]);
  double EFF = Rcpp::as<double>(epiP["eff"]);
  double REPROP = Rcpp::as<double>(epiP["reproP"]);
  double TLATEXP = Rcpp::as<double>(epiP["TlatEXP"]);
  double TLATVAR = Rcpp::as<double>(epiP["TlatVAR"]);
  double TSPOEXP = Rcpp::as<double>(epiP["TspoEXP"]);
  double TSPOVAR = Rcpp::as<double>(epiP["TspoVAR"]);
  double kpatho = Rcpp::as<double>(epiP["kpatho"]);
  double sigpatho = Rcpp::as<double>(epiP["sigpatho"]);
  double spatho = Rcpp::as<double>(epiP["spatho"]);

  double COSTINFECT = Rcpp::as<double>(evolP["costinfect"]);
  double COSTAGGR = Rcpp::as<double>(evolP["costaggr"]);
  double TAUMUT = Rcpp::as<double>(evolP["taumut"]);
  double MGEFF = Rcpp::as<double>(evolP["MGeff"]);
  double QREFF = Rcpp::as<double>(evolP["QReff"]);
  double BETA = Rcpp::as<double>(evolP["beta"]);
  int NAGGR = Rcpp::as<int>(evolP["Naggr"]);
  IntegerVector adaptation_tmp = Rcpp::as<IntegerVector>(evolP["adaptation"]);
  int adaptation[NTRAIT];  
  for (int itrait=0; itrait<NTRAIT; itrait++){
    adaptation[itrait] = adaptation_tmp[itrait];
  }
  
  int NPoly = 0;

  GDALAllRegister();

  #ifdef GDALV2
    GDALDataset *poDS;
  #else
    OGRDataSource *poDS;
  #endif
  // Open a directory with shapefile each shapefile is a layer
  #ifdef GDALV2
    poDS = (GDALDataset *) GDALOpenEx( Rcpp::as<std::string>(shapefilename).c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE, NULL, NULL, NULL);
  #else
    poDS = OGRSFDriverRegistrar::Open( Rcpp::as<std::string>(shapefilename).c_str(), FALSE);
  #endif
        
  if( poDS == NULL ){
    Rcpp::stop("Open : %s failed.\n",shapefilename);
  }

  OGRLayer *poLayer;
  poLayer = poDS->GetLayerByName( Rcpp::as<std::string>(layername_hab).c_str());

  //GIntBig NPoly = poDS->GetFeatureCount(); 
  poLayer->ResetReading();
  OGRFeature *poFeature;
  while( (poFeature = poLayer->GetNextFeature()) != NULL )
  {
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    if( poGeometry != NULL && (wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon || wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon)) NPoly++;
    else Rcerr<<"WARNING unknow geometry type : "<<poGeometry->getGeometryType()<<endl;
  }

  int ipoly=0;
  int * area = init_i1(NPoly);
  int **habitat;
  habitat = init_i2(2,NPoly);
  poLayer->ResetReading();
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    habitat[0][ipoly] = poFeature->GetFieldAsInteger(0);
    habitat[1][ipoly] = poFeature->GetFieldAsInteger(1);
    area[ipoly] = poFeature->GetFieldAsInteger(2);
    ipoly = ipoly+1;
  }

  
  /*------------------*/
  /* Seed generation  */
  /*------------------*/
  const gsl_rng_type *gen_type = gsl_rng_mt19937;
  /* Makoto Matsumoto and Takuji Nishimura, “Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator”. ACM Transactions on Modeling and Computer Simulation, Vol. 8, No. 1 (Jan. 1998), Pages 3–30 */
  gsl_rng *gen = gsl_rng_alloc(gen_type);
  /* This function returns a pointer to a newly-created instance of a random number generator of type gen_type */
  unsigned long int seed = (int) val_seed;
  gsl_rng_set(gen, seed);
  /* This function initializes (or `seeds') the random number generator. */



  /*-----------------------*/
  /* Load input parameters */
  /*-----------------------*/

  int    i,j,poly,hote,year;
  char   *pdelim = ",";
  int    printOn = (int) PRINTON;
  //int    nYears  = (int) NYEARS;
  //int    nTSpY   = (int) NTSPY;
  //int    Nhote   = (int) NHOTE;
  
  /* Evolution and deployment parameters */
  //char *strat       = "STRAT";
  //double propRR     = (double) PROPRR;
  double costInfect = (double) COSTINFECT;
  double costAggr   = (double) COSTAGGR;
  double taumut     = (double) TAUMUT;
  double MGeff      = (double) MGEFF;
  double QReff      = (double) QREFF;
  double beta       = (double) BETA;
  
  /* Host parameters */
  double C0[2];
  double * Cmax = init_d1(Nhote);
  double * delta = init_d1(Nhote);
  double * reproH = init_d1(Nhote);
  double **mort; 
  
  /* Cultivar V0 */
  C0[0]     = (double) C_0;
  C0[1]     = (double) C_0;  /* only useful for mixtures */
  Cmax[0]   = (double) CMAX0;	
  delta[0]  = (double) CROISH0;
  reproH[0] = (double) REPROH0;  
  /* Other cultivars */
  for (hote=1; hote<Nhote; hote++) {
       Cmax[hote]   = (double) CMAX1;
       delta[hote]  = (double) CROISH1;
       reproH[hote] = (double) REPROH1;	
  }
  if (strcmp(strat,"MI") == 0) {    /* Adjustment of C0 and Cmax for mixtures */
       C0[0]   = C0[0] * (1-propRR);
       C0[1]   = C0[1] * propRR;
       Cmax[1] = CMAX1 * (1-propRR);
       Cmax[2] = CMAX1 * propRR;
  }
  
  mort = init_d2(NPoly,Nhote);
  for (poly = 0 ; poly < NPoly ; poly++){
       for (hote = 0 ; hote < Nhote ; hote++)
            mort[poly][hote] = (double) MORTH;
  }    
  
  /* Epidemiological parameters */
  double pI0     = (double) PI0;
  double pSurv   = (double) PSURV;
  double eff     = (double) EFF;
  double repro   = (double) REPROP;
  double Tlat[2], Tspo[2]; 
  Tlat[0]        = (double) TLATEXP;
  Tlat[1]        = (double) TLATVAR;
  Tspo[0]        = (double) TSPOEXP;
  Tspo[1]        = (double) TSPOVAR;  

  /* Resistance formula for each host (resistance genes that are present) */
  //int **resistance;
  //int resist_tmp[3];
  
  //resistance = init_i2(Nhote,8);
  //resist_tmp[0]=(int) atof(RESISTANCE0);
  //resist_tmp[1]=(int) atof(RESISTANCE1);
  //resist_tmp[2]=(int) atof(RESISTANCE2);
  //for (hote=0; hote<Nhote; hote++)
  //     split_i1(8, resist_tmp[hote], resistance[hote]);
  
  /* Adaptation formula (genes that can evolve) */
  //int adaptation[8];
  //int adapt_tmp=(int) atof(ADAPTATION);
  //split_i1(8, adapt_tmp, adaptation);

  int    nIG = adaptation[0]+adaptation[1]+adaptation[2]+adaptation[3];
  int    nAG = adaptation[4]+adaptation[5]+adaptation[6]+adaptation[7];
  int    Naggr = (int) NAGGR;
  int    Npatho = pow(2, nIG) * pow(Naggr, nAG);

  /* Landscape */
  //char *nomfile_area     = "AREA";
  //char *nomfile_habitat1 = "HABITAT1";
  //char *nomfile_habitat2 = "HABITAT2";
  //double *y_areaTot;
  //double *y_habitat1Tot;
  //double *y_habitat2Tot;
  //int *areaTot;
  //int *area;  
  //int **habitatTot;
  //int **habitat;
  //y_areaTot = init_d1(NPoly);
  //areaTot = init_i1(NPoly);
  //area = init_i1(NPoly);
  //y_habitat1Tot = init_d1(NPoly);
  //y_habitat2Tot = init_d1(NPoly);
  //habitatTot = init_i2(2,NPoly);
  //habitat = init_i2(2,NPoly);
  //
  //
  //lire(120, NPoly, (char*)Rcpp::as<std::string>(nomfile_area).c_str(), pdelim, y_areaTot);
  //lire(120, NPoly, (char*)Rcpp::as<std::string>(nomfile_habitat1).c_str(), pdelim, y_habitat1Tot);
  //lire(120, NPoly, (char*)Rcpp::as<std::string>(nomfile_habitat2).c_str(), pdelim, y_habitat2Tot);
  //for (poly = 0 ; poly < NPoly ; poly++) {
  //  areaTot[poly] = (int) y_areaTot[poly];
  //  habitatTot[0][poly] = (int) y_habitat1Tot[poly];
  //  habitatTot[1][poly] = (int) y_habitat2Tot[poly];
  //}
  //for (poly=0;poly<NPoly;poly++) {
  //  area[poly] = areaTot[poly];
  //  habitat[0][poly] = habitatTot[0][poly];
  //  habitat[1][poly] = habitatTot[1][poly];
  //  if (area[poly]==0)
  //  	printf("CAREFULL, one of the areas is 0: check if NPoly, NPoly and idLAN are correct\n");
  //}
 
  /* Crop rotation */
  //char *nomfile_rotation = "ROTATION";
  //double *y_rotation;
  //int *rotation;
  //y_rotation = init_d1(nYears+1);
  //rotation = init_i1(nYears+1);
  //lire(120, nYears+1, (char*)Rcpp::as<std::string>(nomfile_rotation).c_str(), pdelim, y_rotation);
  //for (year = 0 ; year < (nYears+1) ; year++) 
  //     rotation[year] = (int) y_rotation[year];

  /*-----------*/  
  /* Dispersal */
  /*-----------*/
  //int N_elements_dispTot  = NPoly*NPoly;
  //char *nomfile_disppatho = "DISPPATHO";
  //char *nomfile_disphote  = "DISPHOTE";
  //double *y_disppathoTot;
  double **disppathoTot;
  double **disppatho;  
  //double *y_disphoteTot;
  double **disphoteTot;
  double **disphote;
  //y_disppathoTot = init_d1(N_elements_dispTot);
  disppathoTot = init_d2(NPoly,NPoly);
  disppatho = init_d2(NPoly,NPoly);
  //y_disphoteTot = init_d1(N_elements_dispTot);
  disphoteTot = init_d2(NPoly,NPoly);
  disphote = init_d2(NPoly,NPoly);

  //lire(120, N_elements_dispTot, (char*)Rcpp::as<std::string>(nomfile_disppatho).c_str(), pdelim, y_disppathoTot);
  //lire(120, N_elements_dispTot, (char*)Rcpp::as<std::string>(nomfile_disphote).c_str(), pdelim, y_disphoteTot);
  for (i=0 ; i<NPoly ; i++){
    for (j=0 ; j<NPoly ; j++) {
      //disppathoTot[j][i] = y_disppathoTot[j+i*NPoly];
      //disphoteTot[j][i] = y_disphoteTot[j+i*NPoly];
      disppathoTot[j][i] = dispP_tmp[j+i*NPoly];
      disphoteTot[j][i] = dispP_tmp[j+i*NPoly];
    }
  }

  /*--------------------------------------*/  
  /* Write and Print the model parameters */
  /*--------------------------------------*/   
//  FILE *fP;
//  fP = fopen("parameters.txt","w");
  //print_param(fP, seed, (char*)Rcpp::as<std::string>(nomfile_disppatho).c_str(), (char*)Rcpp::as<std::string>(nomfile_disphote).c_str(), "nomfile_habitat1", "nomfile_habitat2", nYears, nTSpY, pSurv, eff, Tlat, repro, Tspo, NPoly, NPoly, Nhote, Npatho, Naggr, C0, pI0, Cmax, delta, reproH, area, habitat, rotation, mort, strat, propRR, resistance, adaptation, MGeff, QReff, taumut, costInfect, costAggr, beta);
    //print_param(stdout, seed, (char*)Rcpp::as<std::string>(nomfile_disppatho).c_str(), (char*)Rcpp::as<std::string>(nomfile_disphote).c_str(), "nomfile_habitat1", "nomfile_habitat2", nYears, nTSpY, pSurv, eff, Tlat, repro, Tspo, NPoly, NPoly, Nhote, Npatho, Naggr, C0, pI0, Cmax, delta, reproH, area, habitat, rotation, mort, strat, propRR, resistance, adaptation, MGeff, QReff, taumut, costInfect, costAggr, beta);
//  fclose(fP);
  
  Rprintf("Pathogen dispersal matrix [poly][poly]\n");
  for (i=0 ; i<NPoly ; i++){
    for (j=0 ; j<NPoly ; j++)
      disppatho[i][j] = disppathoTot[i][j];
    /*if (NPoly<11)
        print_d1(stdout, NPoly, disppatho[i], "");*/
  } 
  Rprintf("Host dispersal matrix [poly][poly]\n");
  for (i=0 ; i<NPoly ; i++){
    for (j=0 ; j<NPoly ; j++)
      disphote[i][j] = disphoteTot[i][j];
    /*if (NPoly<11)
         print_d1(stdout, NPoly, disphote[i],""); */
  }

  /*-----------------------------------------*/  
  /* Mutation, pathotypes and aggressiveness */
  /*-----------------------------------------*/
  int ********aggrToPatho;
  int **pathoToAggr; 
  double **mutkernelMG;
  double **mutkernelQR;
  double **infect;
  double **aggr;
  
  aggrToPatho = init_i8(2,2,2,2,Naggr,Naggr,Naggr,Naggr);
  pathoToAggr = init_i2(Npatho,8);
  mutkernelMG = init_d2(2,2);
  mutkernelQR = init_d2(Naggr,Naggr);
  infect      = init_d2(2,2);
  aggr        = init_d2(Naggr,2);
  
  init_aggrFormula(Npatho, Naggr, adaptation, aggrToPatho, pathoToAggr);
  init_mutkernel(mutkernelMG, mutkernelQR, taumut, Naggr);
  init_infectAggr(Naggr, MGeff, QReff, costInfect, costAggr, beta, infect, aggr);
  
  /* -------------- */
  /* Epidemic model */
  /* -------------- */
  Rprintf("\n*** SPATIOTEMPORAL MODEL SIMULATING THE SPREAD AND EVOLUTION OF A PATHOGEN IN A LANDSCAPE ***\n\n");
  poLayer = poDS->GetLayerByName( Rcpp::as<std::string>(layername_res).c_str());
  dynepi(gen, nYears, nTSpY, NPoly, Nhote, Npatho, area, habitat, rotation, C0, pI0, pSurv, Cmax, mort, delta, reproH, disphote, disppatho, eff, Tlat, repro, Tspo, mutkernelMG,  mutkernelQR, pathoToAggr, aggrToPatho, Naggr, infect, aggr, strat, resistance, adaptation, printOn, poLayer, khote, sighote, shote, kpatho, sigpatho, spatho);//timesStep
       
  /* ----------- */
  /* Free memory */
  /* ----------- */

#ifdef GDALV2
  GDALClose( poDS );
#else
  OGRDataSource::DestroyDataSource( poDS );
#endif

  
 // free(areaTot);
  free(area); 

  free(Cmax);
  free(delta);
  free(reproH);

 // free(y_rotation);
  free(rotation);
 // free(y_disphoteTot);
 // free(y_disppathoTot);
 // free(y_areaTot);
 // free(y_habitat1Tot); 
 // free(y_habitat2Tot); 
 // free_i2(habitatTot,2);
  free_i2(habitat,2);
  free_d2(disphote, NPoly);
  free_d2(disphoteTot, NPoly);  
  free_d2(disppatho, NPoly);
  free_d2(disppathoTot, NPoly);
  free_i8(aggrToPatho, 2,2,2,2,Naggr,Naggr,Naggr);
  free_i2(pathoToAggr, Npatho);
  free_d2(mutkernelMG, 2);
  free_d2(infect, 2);
  free_d2(aggr, Naggr);    
  free_d2(mutkernelQR, Naggr);  
  free_d2(mort,NPoly);
  free_i2(resistance,Nhote);
  
}




