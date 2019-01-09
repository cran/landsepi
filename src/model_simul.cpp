/* 
 * Part of the landsepi R package.
 * Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inra.fr>
 *                    Julien Papaix <julien.papaix@inra.fr>
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
/* The following functions are nested in a { for poly } */


/* -------------------------------------- */
/*        MODE OF REPRODUCTION            */
/* -------------------------------------- */
/* Split reproducing infections (I) between clonal and sexual modes */
void split_IclonalIsex(const gsl_rng *gen, int Npatho, int Nhost, int **I, int **Iclonal_poly, int **Isex_poly, double probSex){
    int patho, host;
    /* Random multinomial draw WITHOUT replacement */
    for (patho = 0 ; patho < Npatho ; patho++){
        for (host = 0 ; host < Nhost ; host++){
            Isex_poly[patho][host]    = gsl_ran_binomial(gen, probSex, I[patho][host]); 
            Iclonal_poly[patho][host] = I[patho][host] - Isex_poly[patho][host];
        } /* for host */
    } /* for patho */
    //printf("$$$$  function split_IclonalIsex  $$$$\n");
    //print_i2(stdout, Npatho, Nhost, I, "I_poly [patho][host]");
    //print_i2(stdout, Npatho, Nhost, Iclonal_poly, "Iclonal_poly [patho][host]");
    //print_i2(stdout, Npatho, Nhost, Isex_poly, "Isex_poly [patho][host]");
}


/* -------------------------------------- */
/*        CLONAL REPRODUCTION             */
/* -------------------------------------- */
/* Production of propagules through clonal reproduction */
/*  (i.e. same genotype as parental individual)         */
/*       Update S in a given poly                       */
void reproClonal(const gsl_rng *gen, int t, int Nhost, int Npatho, int *S, int **I, double repro, int **resistance, int **pathoToAggr, double **aggr, int activeQR){
    /* S and I are the numbers of individuals in a given poly */ 
    int patho, host;
    int ag_R, qr_R;
    double Sprod_tmp, Sprod_exp;
    int Sprod;
    
    //printf("$$$$  function reproClonal  $$$$\n");
    for (patho = 0 ; patho < Npatho ; patho++){
        /* Poisson draw of the number of produced propagules, depending on the effect of aggressiveness on repro */	
        Sprod_exp = 0.0;
        for (host = 0 ; host < Nhost ; host++) {
            Sprod_tmp = repro * I[patho][host];
            ag_R = pathoToAggr[patho][ID_R];        /* indicate if the pathogen has a aggressiveness gene on reproduction rate */
        
        /* indicate if the cultivar has an active QR on reproduction rate */
        if (t < activeQR){
            qr_R = resistance[0][ID_R];          /* if QR is inactive (APR), the cultivar is similar to cultivar 0 */
        }else{
            qr_R = resistance[host][ID_R];
        }
        
        Sprod_tmp *= aggr[ag_R][qr_R];
        Sprod_exp += Sprod_tmp;           /* expectation of the number of produced propagules */
        } /* for host */	
        
        Sprod = gsl_ran_poisson(gen, Sprod_exp);	
        S[patho] += Sprod;
        //printf("# Patho=%d : Sprod=%d\n", patho, Sprod);
    } /* for patho */
        //print_i1(stdout, Npatho, S, "S_poly [patho]");
}          


/* -------------------------------------- */
/*        SEXUAL REPRODUCTION             */
/* -------------------------------------- */
/* Production of propagules through sexual reproduction (or recombination) */
/* Sex/recombination occurs between parents located in the same field and on the same host */
/* Update S in a given poly */
void reproSex(const gsl_rng *gen, int t, int Nhost, int Npatho, int *S, int **I, double repro, int **resistance, int **pathoToAggr, int ********aggrToPatho, double **aggr, int activeQR){
    /* S and I are the numbers of individuals in a given poly */ 
    int i, j, s, patho, host, ag_R1, ag_R2, qr_R, locus;
    int Itot, coupleTot, pathoS, randNum, randBin;
    int *pathoItot_tmp, *pathoItot, pathoParents[2];
    double Sprod_exp;
    int Sprod;
    int binFormula[NLOCI], aggrS[NLOCI];
    
    //printf("$$$$  function reproSex  $$$$\n");
    for (host = 0 ; host < Nhost ; host++){
        /* indicate if the cultivar has an active QR on reproduction rate */
        if (t < activeQR){
            qr_R = resistance[0][ID_R];          /* if QR is inactive (APR), the cultivar is similar to cultivar 0 */
        }else{
            qr_R = resistance[host][ID_R];
        }
        
        /* Total number of parents and couples */
        Itot = 0;
        for (patho = 0 ; patho < Npatho ; patho++)
            Itot += I[patho][host];
        coupleTot = Itot / 2;
        
        /* Sexual reproduction needs at least two parents */
        if (coupleTot > 0){
            /* Vectors containing parent genotypes */
            pathoItot_tmp = init_i1(Itot);
            pathoItot = init_i1(Itot);
            j=0;
            for (patho = 0 ; patho < Npatho ; patho++){
                for (i = 0 ; i < I[patho][host] ; i++)
                    pathoItot_tmp[j++] = patho;
            } /* for patho */
            sample(gen, pathoItot_tmp, Itot, pathoItot, Itot);
            
            for (j = 0; j < coupleTot; j++){   /* i.e. for each couple */
                pathoParents[0] = pathoItot[j];
                pathoParents[1] = pathoItot[j+coupleTot];
                
                /* Propagule production */
                /* Reproduction rate per parent is the mean of both parents        */
                /* i.e. reproduction rate of the couple is the sum of both parents */
                ag_R1 = pathoToAggr[pathoParents[0]][ID_R];        /* indicate if the pathogen has a aggressiveness gene on reproduction rate */        
                ag_R2 = pathoToAggr[pathoParents[1]][ID_R];        /* indicate if the pathogen has a aggressiveness gene on reproduction rate */
                Sprod_exp = repro * (aggr[ag_R1][qr_R] + aggr[ag_R2][qr_R]);
                Sprod = gsl_ran_poisson(gen, Sprod_exp);
                //printf("# Couple j=%d :  Sprod=%d  :  ", j, Sprod);
                if (pathoParents[0] == pathoParents[1]){
                    S[pathoParents[0]] += Sprod;
                    //printf("same parent genotype [%d]\n", pathoParents[0]);
                }else{  
                    //printf("parent1=[%d]  parent2=[%d]\n", pathoParents[0], pathoParents[1]);
                    /* Random loci segregation for each propagule */
                    for (s = 0 ; s < Sprod ; s++){
                        randNum = pow(2,NLOCI) * gsl_rng_uniform(gen);
                        randBin = as_binary((int) randNum);							
                        for (locus=0; locus<NLOCI; locus++){	/* i.e. for each locus */
                            binFormula[locus] = randBin % 10;	/* Algorithm to decompose randBin in a vector of size NLOCI */
                            randBin /= 10;
                            aggrS[locus] = pathoToAggr[ pathoParents[binFormula[locus]] ][locus];
                        } /* for locus */
                        /* Increment S with the corresponding pathogen genotype */
                        pathoS = aggrToPatho[aggrS[0]][aggrS[1]][aggrS[2]][aggrS[3]][aggrS[4]][aggrS[5]][aggrS[6]][aggrS[7]];
                        S[pathoS]++;
                        //printf("patho[%d]=%d\n", s, pathoS);   
                    } /* for s */
                } /* else pathoParents[0] != pathoParents[1] */
            } /* for j */
            //printf("Itot=%d   coupleTot=%d\n", Itot, coupleTot);
            //print_i1(stdout, Itot, pathoItot_tmp, "pathoItot_tmp");
            //print_i1(stdout, Itot, pathoItot, "pathoItot");
                    
            free(pathoItot_tmp);
            free(pathoItot);
        } /* if coupleTot>1 */
    } /* for host */
    //print_i1(stdout, Npatho, S, "S_poly [patho]");
}          


/* -------------------------------------------------- */
/*      MUTATION OF THE PROPAGULES AT A GIVEN LOCUS   */
/* -------------------------------------------------- */  
/* Update SpathoMut with mutation on "trait_mut" through multinomial draw */
void mutation_locus(const gsl_rng *gen, int patho, int Npatho, int Naggr, int trait_mut, double **mutkernel, int **pathoToAggr, int ********aggrToPatho, int **SpathoMut) {
    
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

/* Old code for mutation_locus */
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
/*               MUTATION                 */
/* -------------------------------------- */
/* Update S with mutation of propagules in a given poly */
void mutation(const gsl_rng *gen, int Npatho, int Naggr, int *S, int *adaptation, double **mutkernelMG,  double **mutkernelQR, int **pathoToAggr, int ********aggrToPatho){
    /* S is the numbers of individuals in a given poly */ 
    int patho, patho_mut, id_IG, id_AG;
    int **SpathoMut;
    SpathoMut = init_i2(Npatho,Npatho);  
    
    /* Initialisation SpathoMut at 0 */		
    init_SpathoMut(Npatho, SpathoMut);
    
    for (patho = 0 ; patho < Npatho ; patho++){
        SpathoMut[patho][patho] = S[patho];
        
        /* Mutation of the different traits listed in the adaptation formula */
        for (id_IG = 0; id_IG < 4; id_IG++) {
            if (adaptation[id_IG])
                mutation_locus(gen, patho, Npatho, 2, id_IG, mutkernelMG, pathoToAggr, aggrToPatho, SpathoMut);
        }
        for (id_AG = 4; id_AG < 8; id_AG++) {
            if (adaptation[id_AG])
                mutation_locus(gen, patho, Npatho, Naggr, id_AG, mutkernelQR, pathoToAggr, aggrToPatho, SpathoMut); 
        }
    } /* for patho */
        
    /* Update S after mutation */    
    for (patho_mut = 0; patho_mut < Npatho; patho_mut++) {
        S[patho_mut] = 0;
        for (patho = 0; patho<Npatho; patho++)
            S[patho_mut] += SpathoMut[patho][patho_mut];
    } /* for patho_mut */ 
    //printf("$$$$  function mutation  $$$$\n");
    //print_i1(stdout, Npatho, S, "S_poly [patho]");
    free_i2(SpathoMut,Npatho);  
}     


/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*    functions affecting all poly    */


/* -------------------------------------- */
/*               DISPERSAL                */
/* -------------------------------------- */
/* Dispersal of new host and new pathogen propagules in the landscape */
/* Update Hjuv and S */
void dispersal(const gsl_rng *gen, int Npoly, int Nhost, int Npatho, int **H, int **Hjuv, int **S, double *reproH, double **disphost, double **disppatho){
    /* H, Hjuv, S are the numbers of individuals in a given poly */ 
    int poly, polyE, patho, host;
    int ***Sdisp;
    int ***Hjuvtmp;
    Sdisp = init_i3(Npatho,Npoly,Npoly);      
    Hjuvtmp = init_i3(Nhost,Npoly,Npoly); 
    
    /* Production and dispersal of the host & dispersal of pathogen propagules */
    for (poly = 0 ; poly < Npoly ; poly++){
        /* Pathogen dispersal */
        for (patho = 0 ; patho < Npatho ; patho++)    
            gsl_ran_multinomial (gen, Npoly, S[poly][patho], disppatho[poly], (unsigned int*) Sdisp[patho][poly]);
        /* host reproduction: production and dispersal of Hjuv */
        for (host = 0 ; host < Nhost ; host++)
            gsl_ran_multinomial (gen, Npoly, (int) (reproH[host] * H[poly][host]), disphost[poly], (unsigned int*) Hjuvtmp[host][poly]);
    } /* for poly */
        
    for (poly = 0 ; poly < Npoly ; poly++){
        /* Update the number of spores (S) and Hjuv landing in each field */
        for (patho = 0 ; patho < Npatho ; patho++){
            S[poly][patho] = 0;
            for (polyE = 0 ; polyE < Npoly ; polyE++)
                S[poly][patho] += Sdisp[patho][polyE][poly];	
        } /* for patho */
            
        for (host = 0 ; host < Nhost ; host++){
            Hjuv[poly][host] = 0;
            for (polyE = 0 ; polyE < Npoly ; polyE++)
                Hjuv[poly][host] += Hjuvtmp[host][polyE][poly];
        } /* for host */
    } /* for poly */
        
    //  printf("# Dispersal\n");
    //	print_i3(stdout, Npatho,Npoly,Npoly,Sdisp,"Sdisp");
    //  print_i2(stdout, Npoly,Npatho,S,"S (after dispersal)");
            
    free_i3(Hjuvtmp,Nhost,Npoly);  
    free_i3(Sdisp,Npatho,Npoly);
}  


/* --------------------------------------------- */
/*       BOTTLENECK AT THE END OF THE SEASON     */
/* --------------------------------------------- */
void bottleneck(const gsl_rng *gen, int Npoly, int t, int Nhost, int Npatho, double pSurv, double *Tspo, int ***L, int ***I, double **aggr, int **resistance, int **pathoToAggr, int ***eqIsurv, int *activeQR){
    int poly, patho, host;
    int Tspo_exp;
    int ag_S, qr_S;
    
    for (patho = 0 ; patho < Npatho ; patho++){
        for (host = 0 ; host < Nhost ; host++){
            /* calculate the mean sporulation period */
            ag_S = pathoToAggr[patho][ID_S];        /* indicate if the pathogen has a aggressiveness gene on sporulation duration */
            for (poly = 0; poly < Npoly; poly++){
                /* Reduce the number of infected hosts (bottleneck) */
                eqIsurv[poly][patho][host] = gsl_ran_binomial(gen, pSurv, L[poly][patho][host] + I[poly][patho][host]); 
                
                /* indicate if the cultivar has an active QR on sporulation duration */
                if (t < activeQR[poly]){
                    qr_S = resistance[0][ID_S];          /* if QR is inactive (APR), the cultivar is similar to cultivar 0 */
                }else{
                    qr_S = resistance[host][ID_S];
                }
                Tspo_exp = (int) (Tspo[0] * aggr[ag_S][qr_S]); 
                
                /* Calculate the equivalent number of infectious hosts */
                eqIsurv[poly][patho][host] *= Tspo_exp;
            } /* for poly */
        } /* for host */
    } /* for patho */
}


/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/* The following functions are nested in a { for poly } */

/* -------------------------- */
/*         HOST DYNAMIC       */
/* -------------------------- */
/* Compute host reproduction, death and growth and updtate the number of H in the concerned poly */
void host_dynamic(const gsl_rng *gen, int Nhost, int Npatho, int area, double *Cmax, int *H, int *Hjuv, int **L, int **I, int **R, int *N, double *deathH, double *growthH, double khost, double sighost, double shost){ 
  /* H, Hjuv, L, I and R are the number of host in a given poly */    
  int patho, host;
  double f1host;
  int availSites, siteaccess;
  int K[Nhost];
  int L_host[Nhost];
  int I_host[Nhost];
  int R_host[Nhost];
  int H2M, H2Mjuv, Hnewind, H2H;

  for (host = 0 ; host < Nhost ; host++){
       
       /* Calculation of totals for L, I, R and K */
       K[host] = (int) (area * Cmax[host]);    /* carrying capacity of the cultivar in the concerned paddock */
       L_host[host] = 0;
       I_host[host] = 0;
       R_host[host] = 0;
       for (patho = 0 ; patho < Npatho ; patho++){
            L_host[host] += L[patho][host];
            I_host[host] += I[patho][host];
            R_host[host] += R[patho][host];
       }
       N[host] = H[host] + L_host[host] + I_host[host] + R_host[host];        /* Number of occupied sites */
       
       /* HOST MORTALITY: H2M */
       /* ------------------- */
       H2M = gsl_ran_binomial(gen, (double) deathH[host], H[host]);
       H2Mjuv = gsl_ran_binomial(gen, (double) deathH[host], Hjuv[host]);    
       Hjuv[host] -= H2Mjuv;   /* update Hjuv */
       
       /* HOST REPRODUCTION: Hnewind */
       /* -------------------------- */
       /* Hjuv settlement in the field */
       availSites = K[host] - N[host] - H2M;                      /* Number of available sites */
       if (availSites < 0)  /* security */
            availSites = 0;
       sigmoid((double) shost, (double) khost, (double) sighost, availSites / (double) K[host], &f1host);
       siteaccess = gsl_ran_binomial(gen, f1host, availSites);
       if (siteaccess < Hjuv[host])
            Hnewind = siteaccess;
       else
            Hnewind = Hjuv[host];
       
       /* HOST GROWTH: H2H */
       /* ---------------- */
       H2H = growthH[host] * (H[host] - H2M) * (1 - ((N[host] - H2M + Hnewind)/(double) K[host]));
       if (H2H < 0) {   /* security */
            Rprintf("CAREFUL !!!!!  H2H < 0   one of the areas may be 0: check if Npoly, NpolyTot and idLAN are correct\n");
            H2H = 0;
       } else if ((N[host] - H2M + Hnewind + H2H) > K[host]) {
            Rprintf("CAREFUL !!!!!  H2H too big\n");
            H2H = K[host] - (N[host] - H2M + Hnewind);
       }
       
       /* UPDATE NUMBER OF HOSTS */
       /* ---------------------- */
       H[host] = H[host] - H2M + Hnewind + H2H;
       N[host] = N[host] - H2M + Hnewind + H2H;
       
       /*  printf("# Host dynamic: death, reproduction and growth of host %d\n", host);
        printf("death : H2M = %d   H2Mjuv = %d\n", H2M, H2Mjuv);
        printf("repro : Hnewind = %d\n", Hnewind);
        printf("growth: H2H = %d\n", H2H);
        printf("final H: H = %d\n", H[host]);
        */       
  } /* for host */
}


/* -------------------------------------------------------------- */
/*         CONTAMINATION : spore deposition on healthy sites      */
/* -------------------------------------------------------------- */
/* Calculation of the number of contaminated sites */
void contamination(const gsl_rng *gen, int Nhost, int Npatho, int *H, int *S, int *N, int **Hcontaminated, double kpatho, double sigpatho, double spatho){ 
     /* H, S are the number of individuals in a given poly */ 
     int patho, host;
     double f1patho;
     int Stot=0;
     int Htot=0;
     int S_host[Nhost+1];
     int *S_host_patho;
     double probaH[Nhost+1];
     double *probaP;
     double probaHtot=0;
     double probaPtot=0;
     int Hcontaminable;        /* Contaminable site, where a spore may deposit */
     int *Hcontaminated_tmp;
     probaP = init_d1(Npatho+1);
     S_host_patho = init_i1(Npatho+1);
     Hcontaminated_tmp = init_i1(Npatho+1);

 /* Calculation of total for H and S */
 for (host = 0 ; host < Nhost ; host++)
      Htot += H[host];
 for (patho = 0 ; patho < Npatho ; patho++)
      Stot += S[patho];

  /* Probability for H to belong to each host */
  if (Htot == 0){
       for (host=0; host<Nhost; host++)
            probaH[host] = 0;
  }
  else {
       for (host=0; host<Nhost; host++)
            probaH[host] = (double) H[host] /(double) Htot;
  }
  for (host=0 ; host<Nhost; host++)
       probaHtot += probaH[host];
  probaH[Nhost] = 1 - probaHtot;   
  
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
  gsl_ran_multinomial (gen, Nhost+1, Stot, probaH, (unsigned int*)S_host); 
  
  for (host=0; host<Nhost; host++){
     /* distribution of the spores among the different pathotypes */ 
     gsl_ran_multinomial (gen, Npatho+1, S_host[host], probaP, (unsigned int*)S_host_patho); 
       
     /* Calculation of the number of contaminable sites */
     if (N[host] > 0)
          sigmoid((double) spatho, (double) kpatho, (double) sigpatho, (H[host]/(double) N[host]), &f1patho);
     else
          f1patho = 0.0;
     Hcontaminable = gsl_ran_binomial (gen, f1patho, H[host]);
     
     /* distribution of the contaminable sites among the different pathotypes */ 
     gsl_ran_multinomial (gen, Npatho+1, Hcontaminable, probaP, (unsigned int*)Hcontaminated_tmp);
     
     /* The true number of contaminated sites is the minimum between sites and spores */
     for (patho=0; patho<Npatho; patho++){
          if (Hcontaminated_tmp[patho] < S_host_patho[patho])
               Hcontaminated[patho][host] = Hcontaminated_tmp[patho];
          else
               Hcontaminated[patho][host] = S_host_patho[patho];
     } /* for patho */     
  } /* for host */
     
  //printf("# Contamination and infection of healthy hosts\n");
  //print_i2(stdout, Npatho, Nhost, Hcontaminated,"Hcontaminated [host][patho]");

 free(probaP);
 free(S_host_patho);
 free(Hcontaminated_tmp);

}  

       
/* ----------------------------------------------------------- */
/*         INFECTIOUS CYCLE : transitions H -> L -> I -> R     */
/* ----------------------------------------------------------- */    
/* Calculate the number of contaminated sites that become infected, infectious or removed and update H, L, I, R */
void infection(const gsl_rng *gen, int t, int nTSpY, int Nhost, int Npatho, int *H, int **Hcontaminated, int **L, int **I, int **R, int ***L2I, int ***I2R, double eff, double *Tlat, double *Tspo, int **resistance, int **pathoToAggr, double **infect, double **aggr, int activeQR){ 
  /* H, Hcontaminated, L, I, R, L2I and I2R are the number of individuals in a given poly */   
  int patho, host, h2l, id_MG;
  int ig, mg, ag_E, ag_L, ag_S, qr_E, qr_L, qr_S;
  double eff_exp, Tlat_exp, Tspo_exp;
  int H2L;
  double Tlat_alpha[2];
  double Tspo_alpha[2];
  int lag1, lag2;
  
  for (patho = 0 ; patho < Npatho ; patho++){
    for (host = 0 ; host < Nhost ; host++){
         
      /* infection of healthy sites: H2L */
      /* ------------------------------- */
      eff_exp = eff;
      /* Interaction major genes - Infectivity genes */
      for (id_MG=0; id_MG<4; id_MG++) {
        ig = pathoToAggr[patho][id_MG];        /* indicate if the pathogen has an infectivity gene */
        mg = resistance[host][id_MG];          /* indicate if the cultivar has an major gene */
        eff_exp *= infect[ig][mg];
      } /* for id_MG */
      /* Interaction Quantitative resistance - Aggressiveness gene */
      ag_E = pathoToAggr[patho][ID_E];        /* indicate if the pathogen has a aggressiveness gene on infection efficiency */
      if (t < activeQR){
          qr_E = resistance[0][ID_E];          /* if QR is inactive (APR), the cultivar is similar to cultivar 0 */
      }else{
          qr_E = resistance[host][ID_E];       /* indicate if the cultivar has an active QR on infection efficiency */
      }
      eff_exp *= aggr[ag_E][qr_E];
      H2L = gsl_ran_binomial(gen, eff_exp, Hcontaminated[patho][host]);
     
      /* Latent and infectious periods */
      /* ----------------------------- */
      /* find parameters of gamma distributions from mean and variance */     	    	
      Tlat_exp = Tlat[0];
      ag_L = pathoToAggr[patho][ID_L];        /* indicate if the pathogen has a aggressiveness gene on latent period */
      if (t < activeQR){
          qr_L = resistance[0][ID_L];          /* if QR is inactive (APR), the cultivar is similar to cultivar 0 */
      }else{
          qr_L = resistance[host][ID_L];       /* indicate if the cultivar has an active QR on latent period */
      }
      Tlat_exp /= (aggr[ag_L][qr_L] + 0.001*(aggr[ag_L][qr_L] == 0));  /* security to avoid division by 0 */
      Tlat_exp += 0.001*(Tlat_exp == 0);       /* Security to avoid problem in alpha calculation */
      find_paramGamma(Tlat_exp, Tlat[1], Tlat_alpha);
      
      Tspo_exp = Tspo[0];
      ag_S = pathoToAggr[patho][ID_S];        /* indicate if the pathogen has a aggressiveness gene on sporulation duration */
      if (t < activeQR){
          qr_S = resistance[0][ID_S];          /* if QR is inactive (APR), the cultivar is similar to cultivar 0 */
      }else{
          qr_S = resistance[host][ID_S];       /* indicate if the cultivar has an active QR on sporulation duration */
      }
      Tspo_exp *= aggr[ag_S][qr_S];
      Tspo_exp += 0.001*(Tspo_exp == 0);       /* Security to avoid problem in alpha calculation */				
      find_paramGamma(Tspo_exp, Tspo[1], Tspo_alpha);
      
      for (h2l=0; h2l<H2L; h2l++) { /* recently infected hosts */
          /* Latent period */		
          lag1 = (int) gsl_ran_gamma(gen, Tlat_alpha[0], Tlat_alpha[1]);
           //printf("lat = %d (mean=%.2f)  ", lag1, Tlat[0] * (2 - aggr[ag_L][qr_L]));
           if ((t+lag1) < nTSpY)     	
                L2I[patho][host][t+lag1]++;
           
           /* Sporulation period */		       		
           lag2 = lag1 + (int) gsl_ran_gamma(gen, Tspo_alpha[0], Tspo_alpha[1]);	
           //printf("spo = %d (mean=%.2f)  ", lag2, Tspo[0] * aggr[ag_S][qr_S]);
           if ((t+lag2) < nTSpY)
                I2R[patho][host][t+lag2]++;
      } /* for h2l */
        
      /* Update H, L, I, R */
      /* ----------------- */ 
      H[host]        -= H2L;
      L[patho][host] += H2L;
      L[patho][host] -= L2I[patho][host][t];
      I[patho][host] += L2I[patho][host][t];
      I[patho][host] -= I2R[patho][host][t];
      R[patho][host] += I2R[patho][host][t];  

    } /* for patho */
  } /* for host */

//print_i3(stdout, Npatho, Nhost, nTSpY, L2I, "L2I(poly) [Npatho][Nhost][t]");
//print_i3(stdout, Npatho, Nhost, nTSpY, I2R, "I2R(poly) [Npatho][Nhost][t]"); 
}


 
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
 
/* --------------------------- */
/*   DYNAMIC OF THE EPIDEMIC   */
/* --------------------------- */

void dynepi(const gsl_rng *gen, int nYears, int nTSpY, int Npoly, int Nhost, int Npatho, int *area, int **habitat, int *rotation, double *C0, double pI0, double pSurv, double probSex, double *Cmax, double **deathH, double *growthH, double *reproH, double **disphost, double **disppatho, double eff, double *Tlat, double repro, double *Tspo, double **mutkernelMG,  double **mutkernelQR, int **pathoToAggr, int ********aggrToPatho, int Naggr, int **activeQR, double **infect, double **aggr, char *strat, int **resistance, int *adaptation/*, int printOn*/, OGRLayer *poLayer, double khost, double sighost, double shost, double kpatho, double sigpatho, double spatho){//, NumericVector timesStep,

  int year, t, poly;
  char name_fH[15], name_fHjuv[15], name_fS[15], name_fL[15], name_fI[15], name_fR[15];
  int **S;
  int **H;
  int **Hjuv;
  int **Hcontaminated;      /* Contaminated sites (where a spore is deposited) */
  int ***L;
  int ***I;
  int ***R;
  int N[Nhost];
  int ****L2I;
  int ****I2R;
  int ***eqIsurv;    /* equivalent number of I that survive and produce spores for the next season */ 
  int **Iclonal_poly; /* Number of I exhibiting clonal reproduction in a given poly */
  int **Isex_poly;    /* Number of I exhibiting sexual reproduction in a given poly */

  S = init_i2(Npoly,Npatho); 
  H = init_i2(Npoly,Nhost);
  Hjuv = init_i2(Npoly,Nhost);
  Hcontaminated = init_i2(Npatho,Nhost);
  L = init_i3(Npoly,Npatho,Nhost);
  I = init_i3(Npoly,Npatho,Nhost);
  R = init_i3(Npoly,Npatho,Nhost);
  L2I = init_i4(Npoly,Npatho,Nhost,nTSpY);
  I2R = init_i4(Npoly,Npatho,Nhost,nTSpY); 
  eqIsurv = init_i3(Npoly, Npatho, Nhost);
  Iclonal_poly = init_i2(Npatho, Nhost);
  Isex_poly    = init_i2(Npatho, Nhost);
  
  /* Initialisation (t=0) */ 
  init_HHjuvSLIR(Npoly, Nhost, Npatho, H, Hjuv, S, L, I, R);
  init_L2I2R(Npoly, Npatho, Nhost, nTSpY, L2I, I2R);
  intro_H(Npoly, Nhost, H, area, habitat, C0, strat, rotation[0]); 
  intro_I(gen, nTSpY, Npoly, Nhost, Npatho, H, I, I2R, pI0, Tspo, resistance, infect, aggr, ID_E, ID_S);
       
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
	for (t = 1 ; t < nTSpY ; t++){
	     
	     /* Writing model output for timestep t */
	     //if(t==timesStep[INDtime]){
	       //sortie_layer(poLayer,H,Npoly,Nhost,t,year);
	       write_HHjuvSLIR(Npoly, Npatho, Nhost, t, H, Hjuv, S, L, I, R, fH, fHjuv, fS, fL, fI, fR/*, printOn*/);
	       //INDtime=INDtime+1;
	       //}	     
	     init_S(Npoly, Npatho, S);                                              /* Re-initialisation at 0 */
	     for (poly = 0 ; poly < Npoly ; poly++){ 
	         //printf("___________________  poly %d  ___________________\n", poly);
	         split_IclonalIsex(gen, Npatho, Nhost, I[poly], Iclonal_poly, Isex_poly, probSex);
	         reproClonal(gen, t, Nhost, Npatho, S[poly], Iclonal_poly, repro, resistance, pathoToAggr, aggr, activeQR[year-1][poly]);
	         reproSex(gen, t, Nhost, Npatho, S[poly], Isex_poly, repro, resistance, pathoToAggr, aggrToPatho, aggr, activeQR[year-1][poly]);
	         mutation(gen, Npatho, Naggr, S[poly], adaptation, mutkernelMG, mutkernelQR, pathoToAggr, aggrToPatho);
	     } /* for poly */
	     dispersal(gen, Npoly, Nhost, Npatho, H, Hjuv, S, reproH, disphost, disppatho);
	     
	     for (poly = 0 ; poly < Npoly ; poly++){     
	          //      printf("poly  %d\n", poly);
	       host_dynamic(gen, Nhost, Npatho, area[poly], Cmax, H[poly], Hjuv[poly], L[poly], I[poly], R[poly], N, deathH[poly], growthH, khost, sighost, shost);
	       contamination(gen, Nhost, Npatho, H[poly], S[poly], N, Hcontaminated, kpatho, sigpatho, spatho);
	       infection(gen, t, nTSpY, Nhost, Npatho, H[poly], Hcontaminated, L[poly], I[poly], R[poly], L2I[poly], I2R[poly], eff, Tlat, Tspo, resistance, pathoToAggr, infect, aggr, activeQR[year-1][poly]);
	     } /* for poly */
	}  /* for t */

     /* last time-step of the season: bottleneck before starting a new season */
     //sortie_layer(poLayer,H,Npoly,Nhost,nTSpY,year);	
     write_HHjuvSLIR(Npoly, Npatho, Nhost, nTSpY, H, Hjuv, S, L, I, R, fH, fHjuv, fS, fL, fI, fR/*, printOn*/);  /* Writing model output for last timestep */
     bottleneck(gen, Npoly, nTSpY, Nhost, Npatho, pSurv, Tspo, L, I, aggr, resistance, pathoToAggr, eqIsurv, activeQR[year-1]);     /* Calculation of the equivalent number of I */
         
     init_HHjuvSLIR(Npoly, Nhost, Npatho, H, Hjuv, S, L, I, R);                                              /* Re-initialisation at 0 */
     init_L2I2R(Npoly, Npatho, Nhost, nTSpY, L2I, I2R);                                                      /* Re-initialisation at 0 */
    
    /* Generate S issued from eqIsurv = (remaining L+I) * Tspo */
    for (poly = 0 ; poly < Npoly ; poly++){ 
        split_IclonalIsex(gen, Npatho, Nhost, eqIsurv[poly], Iclonal_poly, Isex_poly, probSex);
        reproClonal(gen, nTSpY, Nhost, Npatho, S[poly], Iclonal_poly, repro, resistance, pathoToAggr, aggr, activeQR[year-1][poly]);
        reproSex(gen, nTSpY, Nhost, Npatho, S[poly], Isex_poly, repro, resistance, pathoToAggr, aggrToPatho, aggr, activeQR[year-1][poly]);
        mutation(gen, Npatho, Naggr, S[poly], adaptation, mutkernelMG, mutkernelQR, pathoToAggr, aggrToPatho);
    } /* for poly */
    dispersal(gen, Npoly, Nhost, Npatho, H, Hjuv, S, reproH, disphost, disppatho);
     
     /* re-plantation --> regenerate H */
     intro_H(Npoly, Nhost, H, area, habitat, C0, strat, rotation[year]);
     /* infection of newly planted hosts to generate the inoculum of the next season */
     for (poly = 0 ; poly < Npoly ; poly++){ 
          contamination(gen, Nhost, Npatho, H[poly], S[poly], H[poly], Hcontaminated, kpatho, sigpatho, spatho);  /* N = H[poly] in beginning of next season */
          infection(gen, 0, nTSpY, Nhost, Npatho, H[poly], Hcontaminated, L[poly], I[poly], R[poly], L2I[poly], I2R[poly], eff, Tlat, Tspo, resistance, pathoToAggr, infect, aggr, activeQR[year][poly]);
     }
      
	fclose(fH);
	fclose(fHjuv);
	fclose(fL);
	fclose(fI);
	fclose(fS);
	fclose(fR);	    
  } /* for year */ 

  free_i2(Hcontaminated,Npatho);  
  free_i4(L2I, Npoly, Npatho, Nhost);
  free_i4(I2R, Npoly, Npatho, Nhost);	    	    
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

/*void sortie_layer(OGRLayer *poLayer,int **H,int Npoly, int Nhost, int timeStep, int year){

  char nomH[256];
  double *resH;
  resH=init_d1(Npoly);
     
  for (int i=0; i<Nhost; i++){
    for (int j=0; j<Npoly; j++){
      resH[j] = H[j][i];
    }
    sprintf(nomH, "H%i_%i-%i", i, year, timeStep);
    addField(poLayer, nomH, resH);
  }
 
 free(resH);
 
}*/



/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/

void model_landsepi(Rcpp::List timeP, Rcpp::List landscape, Rcpp::List dispersal, Rcpp::List inits ,int val_seed, Rcpp::List hostP, Rcpp::List pathoP, Rcpp::List evolP) {
    
    int    i,j,poly,host,year;
    char   *pdelim = ",";
    //  int    printOn = 0;
    
    
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
    
    
    /*-------------------*/ 
    /*  Time parameters  */ 
    /*-------------------*/
    int nYears = Rcpp::as<int>(timeP["nYears"]);
    int nTSpY = Rcpp::as<int>(timeP["nTSpY"]);
    //NumericVector timesStep=Rcpp::as<NumericVector>(timeP["timesStep"]);
    
    
    /*-------------------------------------*/
    /*  Landscape & deployment parameters  */
    /*-------------------------------------*/
    CharacterVector shapefilename = Rcpp::as<CharacterVector>(landscape["shapefilename"]);
    CharacterVector layername_hab = Rcpp::as<CharacterVector>(landscape["layername_hab"]);
    CharacterVector layername_res = Rcpp::as<CharacterVector>(landscape["layername_res"]);
 //   char * strat = (char*) Rcpp::as<std::string>(Rcpp::as<CharacterVector>(landscape["strat"])).c_str();
 
    /* avoid a bug when importing a string from R */
    int strat_tmp = Rcpp::as<int>(timeP["strat_tmp"]);
    char strat[3];
    if (strat_tmp==1)
        strcpy(strat, "MO");
    if (strat_tmp==2)
        strcpy(strat, "MI");
    if (strat_tmp==3)
        strcpy(strat, "RO");
    if (strat_tmp==4)
        strcpy(strat, "PY");
 
    
    IntegerVector rotation_tmp = Rcpp::as<IntegerVector>(landscape["rotation"]);
    double propRR = Rcpp::as<double>(landscape["propRR"]);
    //fprintf(stderr,"\ntype rotation: %s\n",typeid(rotation[0]).name());

    int * rotation=init_i1(nYears+1);
    IntegerVector::iterator it;
    int roti;
    for (roti=0, it = rotation_tmp.begin() ; it != rotation_tmp.end(); ++it, roti++){
        rotation[roti]=*it;
    }
    
    int Npoly = 0;
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
    
    //GIntBig Npoly = poDS->GetFeatureCount(); 
    poLayer->ResetReading();
    OGRFeature *poFeature;
    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        OGRGeometry *poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL && (wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon || wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon)) Npoly++;
        else Rcerr<<"WARNING unknow geometry type : "<<poGeometry->getGeometryType()<<endl;
    }
    
    int ipoly=0;
    int * area = init_i1(Npoly);
    int **habitat;
    habitat = init_i2(2,Npoly);
    poLayer->ResetReading();
    while( (poFeature = poLayer->GetNextFeature()) != NULL ){
        habitat[0][ipoly] = poFeature->GetFieldAsInteger(0);
        habitat[1][ipoly] = poFeature->GetFieldAsInteger(1);
        area[ipoly] = poFeature->GetFieldAsInteger(2);
        ipoly = ipoly+1;
    }
    
    
    /*------------------------*/
    /*  Dispersal parameters  */
    /*------------------------*/
    NumericVector dispP_tmp=Rcpp::as<NumericVector>(dispersal["dispP"]);
    NumericVector dispH_tmp=Rcpp::as<NumericVector>(dispersal["dispH"]);
    double **disppathoTot;
    double **disppatho;  
    double **disphostTot;
    double **disphost;
    disppathoTot = init_d2(Npoly,Npoly);
    disppatho = init_d2(Npoly,Npoly);
    disphostTot = init_d2(Npoly,Npoly);
    disphost = init_d2(Npoly,Npoly);
    
    for (i=0 ; i<Npoly ; i++){
        for (j=0 ; j<Npoly ; j++) {
            disppathoTot[j][i] = dispP_tmp[j+i*Npoly];
            disphostTot[j][i] = dispP_tmp[j+i*Npoly];
        }
    }
    //Rprintf("Pathogen dispersal matrix [poly][poly]\n");
    for (i=0 ; i<Npoly ; i++){
        for (j=0 ; j<Npoly ; j++)
            disppatho[i][j] = disppathoTot[i][j];
        /*if (Npoly<11)
         print_d1(stdout, Npoly, disppatho[i], "");*/
    } 
    //Rprintf("Host dispersal matrix [poly][poly]\n");
    for (i=0 ; i<Npoly ; i++){
        for (j=0 ; j<Npoly ; j++)
            disphost[i][j] = disphostTot[i][j];
        /*if (Npoly<11)
         print_d1(stdout, Npoly, disphost[i],""); */
    }
    
    
    /*----------------------*/
    /*  Initial conditions  */
    /*----------------------*/
    double pI0 = Rcpp::as<double>(inits["pI0"]);
    
    
    /*-------------------*/
    /*  Host parameters  */
    /*-------------------*/
    int Nhost = Rcpp::as<int>(hostP["Nhost"]);
    NumericVector C0_tmp=Rcpp::as<NumericVector>(hostP["C0"]);
    NumericVector Cmax_tmp=Rcpp::as<NumericVector>(hostP["Cmax"]);
    NumericVector growthH_tmp=Rcpp::as<NumericVector>(hostP["growthH"]);
    NumericVector reproH_tmp=Rcpp::as<NumericVector>(hostP["reproH"]);
    NumericVector deathH_tmp=Rcpp::as<NumericVector>(hostP["deathH"]);
    double *C0 = init_d1(Nhost);
    double *Cmax = init_d1(Nhost);
    double *growthH = init_d1(Nhost);
    double *reproH = init_d1(Nhost);
    double **deathH; 
    for (host=0; host<Nhost; host++) {
        C0[host]     = (double) C0_tmp[host];
        Cmax[host]   = (double) Cmax_tmp[host];
        growthH[host]  = (double) growthH_tmp[host];
        reproH[host] = (double) reproH_tmp[host];	
    }
    //  if (strcmp(strat,"MI") == 0) {    /* Adjustment of C0 and Cmax for mixtures */
    //       C0[0]   = C0[0] * (1-propRR);
    //       C0[1]   = C0[1] * propRR;
    //       Cmax[1] = CMAX1 * (1-propRR);
    //       Cmax[2] = CMAX1 * propRR;
    //  }
    deathH = init_d2(Npoly,Nhost);
    for (poly = 0 ; poly < Npoly ; poly++){
        for (host = 0 ; host < Nhost ; host++)
            deathH[poly][host] = deathH_tmp[host];
    }
    
    double khost = Rcpp::as<double>(hostP["khost"]);
    double sighost = Rcpp::as<double>(hostP["sighost"]);
    double shost = Rcpp::as<double>(hostP["shost"]);
    IntegerVector resistance_tmp = Rcpp::as<IntegerVector>(hostP["resistance"]);
    int **resistance;  
    resistance = init_i2(Nhost,NLOCI);
    for (int ihost=0; ihost<Nhost; ihost++){
        for (int itrait=0; itrait<NLOCI; itrait++){
            resistance[ihost][itrait] = resistance_tmp[itrait+ihost*NLOCI];
        }
    }  
    
    
    /*---------------------*/
    /* Pathogen parameters */
    /*---------------------*/
    double pSurv = Rcpp::as<double>(pathoP["pSurv"]);
    double probSex = Rcpp::as<double>(pathoP["probSex"]);
    double eff = Rcpp::as<double>(pathoP["eff"]);
    double repro = Rcpp::as<double>(pathoP["repro"]);
    double Tlat[2], Tspo[2]; 
    Tlat[0] = Rcpp::as<double>(pathoP["TlatEXP"]);
    Tlat[1] = Rcpp::as<double>(pathoP["TlatVAR"]);
    Tspo[0] = Rcpp::as<double>(pathoP["TspoEXP"]);
    Tspo[1] = Rcpp::as<double>(pathoP["TspoVAR"]);
    double kpatho = Rcpp::as<double>(pathoP["kpatho"]);
    double sigpatho = Rcpp::as<double>(pathoP["sigpatho"]);
    double spatho = Rcpp::as<double>(pathoP["spatho"]);
    
    
    /*------------------------*/
    /*  Evolution parameters  */
    /*------------------------*/
    int Naggr = Rcpp::as<int>(evolP["Naggr"]);
    IntegerVector adaptation_tmp = Rcpp::as<IntegerVector>(evolP["adaptation"]);
    int adaptation[NLOCI];  
    for (int itrait=0; itrait<NLOCI; itrait++){
        adaptation[itrait] = adaptation_tmp[itrait];
    }
    int    nIG = adaptation[0]+adaptation[1]+adaptation[2]+adaptation[3];
    int    nAG = adaptation[4]+adaptation[5]+adaptation[6]+adaptation[7];
    int    Npatho = pow(2, nIG) * pow(Naggr, nAG);
    
    double costInfect = Rcpp::as<double>(evolP["costinfect"]);
    double costAggr = Rcpp::as<double>(evolP["costaggr"]);
    double taumut = Rcpp::as<double>(evolP["taumut"]);
    double MGeff = Rcpp::as<double>(evolP["MGeff"]);
    double QReff = Rcpp::as<double>(evolP["QReff"]);
    double beta = Rcpp::as<double>(evolP["beta"]);
    /* Time to quantitative resistance expression */
    double timeToQR[2];
    timeToQR[0] = Rcpp::as<double>(evolP["timeToQR_exp"]);
    timeToQR[1] = Rcpp::as<double>(evolP["timeToQR_var"]);
    int **activeQR;
    double timeToQR_alpha[2];
    activeQR = init_i2(nYears+1, Npoly);
    init_activeQR(gen, nYears, Npoly, timeToQR, timeToQR_alpha, activeQR);
    
    
    /*--------------------------------------*/  
    /* Write and Print the model parameters */
    /*--------------------------------------*/   
    FILE *fP;
    fP = fopen("parameters.txt","w");
    print_param(fP, seed, nYears, nTSpY, Npoly, Nhost, strat, C0, Cmax, growthH, reproH, pI0, pSurv, kpatho, spatho, sigpatho, eff, repro, Tlat, Tspo, Npatho, Naggr, resistance, adaptation, MGeff, QReff, timeToQR, taumut, probSex, costInfect, costAggr, beta);
    fclose(fP);
    
    
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
    pathoToAggr = init_i2(Npatho,NLOCI);
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
    dynepi(gen, nYears, nTSpY, Npoly, Nhost, Npatho, area, habitat, rotation, C0, pI0, pSurv, probSex, Cmax, deathH, growthH, reproH, disphost, disppatho, eff, Tlat, repro, Tspo, mutkernelMG,  mutkernelQR, pathoToAggr, aggrToPatho, Naggr, activeQR, infect, aggr, strat, resistance, adaptation/*, printOn*/, poLayer, khost, sighost, shost, kpatho, sigpatho, spatho);//timesStep
    
    /* ----------- */
    /* Free memory */
    /* ----------- */
    
#ifdef GDALV2
    GDALClose( poDS );
#else
    OGRDataSource::DestroyDataSource( poDS );
#endif
    
    
    free(area); 
    free(C0);
    free(Cmax);
    free(growthH);
    free(reproH);
    free(rotation);
    free_i2(habitat,2);
    free_d2(disphost, Npoly);
    free_d2(disphostTot, Npoly);  
    free_d2(disppatho, Npoly);
    free_d2(disppathoTot, Npoly);
    free_i8(aggrToPatho, 2,2,2,2,Naggr,Naggr,Naggr);
    free_i2(pathoToAggr, Npatho);
    free_d2(mutkernelMG, 2);
    free_d2(infect, 2);
    free_d2(aggr, Naggr);    
    free_d2(mutkernelQR, Naggr); 
    free_i2(activeQR, nYears);
    free_d2(deathH,Npoly);
    free_i2(resistance,Nhost);
  
}




