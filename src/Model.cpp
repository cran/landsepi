/*
 * Part of the landsepi R package.
 * Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
 *                    Julien Papaix <julien.papaix@inrae.fr>
 *                    Jean-François Rey <jean-francois.rey@inrae.fr>
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
/*  (int) a/b is the euclidian quotient  */
/*  (int) a%b is the euclidian remainder */
/*  ++variable == variable++ == variable=variable+1 */
/*  ++*pointer == (*pointer)++ == *pointer=*pointer+1    CAREFULL   *pointer++   DOES NOT WORK */

/* CULTIVAR GENOTYPE : mg are major resistant genes, qr are traits of quantitative resistance  */
/* PATHOGEN GENOTYPE : ig are infectivity genes,     ag are aggressiveness genes               */
/* Cultivar 0 is susceptible (no mg, no qr), Cultivars > 0 have one or more resistance sources */
/* Pathogen 0 is avirulent   (no ig, no ag), Pathogen > 0 have some ig or ag                   */

#include "Model.hpp"

// Model constructor
Model::Model(const int& Nyears, const int& time_steps_per_year, const int& Npoly, const int& Nhost, const int& Npatho,
             const int& Ngene, const std::vector<double>& area, const Vector2D<int>& rotation,
             const gsl_rng* random_generator, const std::vector<Cultivar>& cultivars, const std::vector<Gene>& genes,
             const Basic_patho& basic_patho, const Treatment& treatment, const std::vector<Croptype>& croptypes, const double& sigmoid_kappa_host,
             const double& sigmoid_sigma_host, const double& sigmoid_plateau_host, const double& pI0,
             const Vector2D<double>& disp_patho, const Vector2D<double>& disp_patho_sex, const Vector2D<double>& disp_host, const int& seed)
  : Nyears(Nyears),
    time_steps_per_year(time_steps_per_year),
    Npoly(Npoly),
    Nhost(Nhost),
    Npatho(Npatho),
    Ngene(Ngene),
    area(area),
    rotation(rotation),
    random_generator(random_generator),
    cultivars(cultivars),
    genes(genes),
    basic_patho(basic_patho),
    treatment(treatment),
    croptypes(croptypes),
    sigmoid_kappa_host(sigmoid_kappa_host),
    sigmoid_sigma_host(sigmoid_sigma_host),
    sigmoid_plateau_host(sigmoid_plateau_host),
    pI0(pI0),
    disp_patho(disp_patho),
    disp_patho_sex(disp_patho_sex),
    disp_host(disp_host) {
  gsl_rng_set(random_generator, seed); // This function initializes (or `seeds') the random number generator.
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/* The following functions are nested in a { for poly } */

/* -------------------------------------- */
/*        MODE OF REPRODUCTION            */
/* -------------------------------------- */

/* Split reproducing infections (I) between clonal and sexual modes using prob at step t */
std::array<Vector2D<int>, 2> Model::split_IclonalIsex(const int& t, const Vector2D<int>& I) {
  /* Random multinomial draw WITHOUT replacement */
  
  // Number of I exhibiting clonal reproduction in a given poly
  Vector2D<int> Iclonal_poly(this->Npatho, std::vector<int>(this->Nhost));
  // Number of I exhibiting sexual reproduction in a given poly
  Vector2D<int> Isex_poly(this->Npatho, std::vector<int>(this->Nhost));
  //Rcpp::Rcout << "Repro sex probability is" << " " << this->basic_patho.repro_sex_prob[t]<< "\n"; 
  for(int patho = 0; patho < this->Npatho; patho++) {
    for(int host = 0; host < this->Nhost; host++) {
      Isex_poly[patho][host] = ran_binomial(this->basic_patho.repro_sex_prob[t], I[patho][host]);
      Iclonal_poly[patho][host] = I[patho][host] - Isex_poly[patho][host];
    }
  }
  return {Iclonal_poly, Isex_poly};
}

/* -------------------------------------- */
/*        CLONAL REPRODUCTION             */
/* -------------------------------------- */
/* Production of propagules through clonal reproduction */
/*  (i.e. same genotype as parental individual)         */
/*       Update P in a given poly                       */
void Model::reproClonal(const int& t, std::vector<int>& P, const Vector2D<int>& I, const std::vector<int>& activeR) {
  // Rcpp::Rcout  << "\n" << "Repro clonal" << " ";
  /* P and I are the numbers of individuals in a given poly */
  for(int patho = 0; patho < this->Npatho; patho++) {
    // Indicate if the pathogen has a aggressiveness gene on propagule production rate
    const std::vector<int> aggr = this->switch_patho_to_aggr(patho);
    // Poisson draw of the number of produced propagules, depending on the effect of aggressiveness on
    // propagule_prod_rate
    double Pprod_exp = 0.0;
    for(int host = 0; host < this->Nhost; host++) {
      double Pprod_tmp = this->basic_patho.propagule_prod_rate * I[patho][host];
      /* Interaction with resistance genes */
      for(int g = 0; g < this->Ngene; g++) {
        if(this->genes[g].target_trait == "PR") {
          // Indicate if the cultivar has an active resistance gene
          const bool Rgene = this->get_resistance(g, host, t, activeR[g]);
          Pprod_tmp *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
        }
      }
      Pprod_exp += Pprod_tmp; // Expectation of the number of produced propagules
    }
    const int Pprod = ran_poisson(Pprod_exp);
    P[patho] += Pprod;
  }
}

/* -------------------------------------- */
/*        SEXUAL REPRODUCTION             */
/* -------------------------------------- */
/* Production of propagules through sexual reproduction (or recombination) */
/* (propagules' genotypes is sampled from a normal distribution with mean equal to the mean of the parents' traits 
 and fixed variance) */
/* Sex/recombination occurs between parents located in the same field and on the same host */
/* Update P in a given poly */

void Model::reproSex(const int& t, std::vector<int>& P, const Vector2D<int>& I, const std::vector<int>& activeR, const std::vector<int>& Nlevels_aggressiveness, const int& Nquali_gene){
  int Nquanti = this->Ngene - Nquali_gene; // number of quantitative genes
  // P and I are the numbers of individuals in a given poly 
  for(int host = 0; host < this->Nhost; host++) {
    int Itot = 0; // Total number of parents in a given host
    for(int patho = 0; patho < this->Npatho; patho++) {
      Itot += I[patho][host];
    }
    const int coupleTot = floor(Itot / 2); // Total number of couples in a given host
    
    // resistance define whether a locus (of the genotype of the host) carries the resistance gene
    // Note that, since the sampling for defining the propagules genotype is done from a multivariate
    // gaussian distribution centered in the mean of the parents trait ON THE S HOST,
    // resistance is always set to 0 (no matter the real host considered)
    std::vector<bool> resistance;
    resistance =  std::vector<bool>(this->Ngene,0);
  
    // Sexual reproduction needs at least two parents 
    if(coupleTot > 0) {
      int counter = 0; // counter counts for the number of individuals with the same genotype
      std::vector<int> pathoItot_tmp(Itot); // Vectors containing  the genotype of each parent (ordered)
      for(int patho = 0; patho < this->Npatho; patho++) {
        for(int i = 0; i < I[patho][host]; i++) {
          pathoItot_tmp[i+counter] = patho;
        }
        counter += I[patho][host];
      }
      
      // pathoItot is a vector of length Itot cointaining the genotype of each parent (randomized)
      const std::vector<int> pathoItot(sample(this->random_generator, pathoItot_tmp));
      
      // couple_set is a square matrix of size Npatho, it contains the number of couples
      // composed by a parent with genotype = r (row of couple_set) and the other with genotype = c (column of couple_set)  
      Vector2D<int> couple_set(this->Npatho, std::vector<int>(this->Npatho,0));
      for(int j = 0; j < coupleTot; j++) {
        couple_set[pathoItot[j]][pathoItot[j + coupleTot]]+=1;
      }
      // Here the sexual reproduction starts: first the number of propagules for the couples (r,c)
      // is defined, then the genotype for each propagule is defined
      
      for(int c = 0; c< this->Npatho; c++) {
        for(int r = c; r< this->Npatho; r++) {
          if(r!=c){
            couple_set[r][c]+=couple_set[c][r];
          }
          // I check if the couple (r,c) exists
          if(couple_set[r][c]!=0){
            // pathoParents is a 2 elements vector containing the genotypes of the couple of parents
            std::vector<int> pathoParents = {r, c};
            
            // Propagule production 
            // Reproduction rate per parent is the mean of both parents        
            // i.e. reproduction rate of the couple is the sum of both parents
            
            double Pprod_exp = 2* this->basic_patho.propagule_prod_rate;
            
            Vector2D<int> aggr_parents(2, std::vector<int>(this->Ngene)); 
            Vector2D<double> trait_parents(2, std::vector<double>(this->Ngene)); 
            for(int k = 0; k < 2; k++) {
              // parents' indices of aggressiveness 
              aggr_parents[k] = this->switch_patho_to_aggr(pathoParents[k]); 
              // trait values for resistance (in the host) associated with parents' aggressiveness 
              // NOTE the trait value is always computed on the S host
              trait_parents[k] = this->switch_aggr_to_trait(aggr_parents[k],resistance); 
            }
            // Interaction with resistance genes 
            for(int g = 0; g < this->Ngene; g++) {
              if(this->genes[g].target_trait == "PR") {
                // Indicate if the cultivar has an active resistance gene
                const bool Rgene = this->get_resistance(g, host, t, activeR[g]);
                Pprod_exp *= 0.5*(this->genes[g].aggressiveness_matrix[aggr_parents[0][g]][Rgene] +
                  this->genes[g].aggressiveness_matrix[aggr_parents[1][g]][Rgene]);
              }
            }
            
            double Pprod = 0.0;
            // Pprod is the total number of propagules produced by all the couples (r,c)
            Pprod = ran_poisson(Pprod_exp*couple_set[r][c]); 
            
            // If the two parents has the same genotype (r == c), all the produced propagule will have the same genotype too
            
            if( r == c ){ // if we are on the diagonal of the matrix couple_set
              P[r]+= Pprod ;
            } else {
            // trait_prop contains the trait values (relative to the quantitatite resistance genes) 
            // for each propagules produced by couples (r,c)
            std::vector<std::vector<double>> trait_prop(Pprod, std::vector<double>(Nquanti)); 
            
            // Print resistance and parents genotype
            
#ifdef DEBUG
              Rcpp::Rcout  << "\n" << "Resistance" << " ";
              for(int i=0; i < this->Ngene; i++){
                Rcpp::Rcout << resistance[i]<< " ";
              }
              
              Rcpp::Rcout  << "\n" << "Parent 1" << " ";
              for(int locus = 0; locus < this->Ngene; locus++) {
                Rcpp::Rcout << aggr_parents[0][locus] << " " << trait_parents[0][locus] << " " ;
              }
              Rcpp::Rcout << "\n" << "Parent 2" << " ";
              for(int locus = 0; locus < this->Ngene; locus++) {
                Rcpp::Rcout << aggr_parents[1][locus] << " "<< trait_parents[1][locus] << " " ;
              }
#endif
        
            // compute the trait values (relative to the quantitative resistance genes) all the propagules
            if(Nquanti>0){
              std::vector<double> mu_vec(Nquanti);
              for(int locus = 0; locus < Nquanti; locus++) {
                // mu_vec is a vector with the average value of each trait of the two parents
                // NOTE: the traits are always evaluated with respect to the S host
                mu_vec[locus] = (trait_parents[0][locus+Nquali_gene] + trait_parents[1][locus+Nquali_gene])/2.0;
              } 
              
              // Variance covariance matrix contains the variance and covariance of each qualitative trait
              // (it is fixed for all the pathogen population).
              
              double wild_type_value = 0.0;
              Vector2D<double> cov;
              cov = Vector2D<double>(Nquanti, std::vector<double>(Nquanti,0.0));
  
              for(int g = 0; g < Nquanti; g++) {
                   if(this->genes[g+Nquali_gene].target_trait == "IR") {
                     wild_type_value = this->basic_patho.infection_rate;
                   } else if (this->genes[g+Nquali_gene].target_trait == "LAT"){
                     wild_type_value = this->basic_patho.latent_period_mean;
                 } else if (this->genes[g+Nquali_gene].target_trait == "IP"){
                     wild_type_value = this->basic_patho.infectious_period_mean;
                   } else if (this->genes[g+Nquali_gene].target_trait == "PR"){
                     wild_type_value = this->basic_patho.propagule_prod_rate;
                   } 
                 cov[g][g] = pow((this->genes[g+Nquali_gene].recombination_sd * wild_type_value),2);
               }

              // the traits values for quantitative resistance associated to each propagule Pprod is sampled
              // from a multivariate gaussian centered in the mean parents' trait values'
              trait_prop = ran_multisample_multivariate_gaussian(Pprod, mu_vec, cov);
              
            }
            // aggr_prop is a vector containing the index of aggressiveness associated to all 
            //  the resistance genes (MG + QTL) for a single propagule
            std::vector<int> aggr_prop(this->Ngene);
            
            // Random loci segregation for each propagule (only for qualitative resistance gene) 
            for(int s = 0; s < Pprod; s++) {
#ifdef DEBUG
                Rcpp::Rcout << "\n" << "Propagule " << s << "\n";
#endif
              const int randNum = static_cast<int>(pow(2, static_cast<double>(this->Ngene)) * rng_uniform());
              int randBin = as_binary(randNum);
              
              // QUALITATIVE RESISTANCE
              // If the parents have the same index of aggressiveness associated to a major gene,
              // the propagule "s" will have the same index of aggressiveness, otherwise, the propagule 
              // will have the index of aggressiveness either of one or of the other parent
              
              for(int locus = 0; locus < Nquali_gene; locus++) { // i.e. for each locus of qualitative genes
                // Algorithm to decompose randBin in a vector of size Ngene
                if(aggr_parents[0][locus] == aggr_parents[1][locus]) {
                  aggr_prop[locus] = aggr_parents[0][locus];
#ifdef DEBUG
                    Rcpp::Rcout << " " << aggr_prop[locus] ;
#endif
                } else {
                  const int binFormula = randBin % 10;
                  randBin /= 10;
                  aggr_prop[locus] = aggr_parents[binFormula][locus];
#ifdef DEBUG
                    Rcpp::Rcout << " " <<  aggr_prop[locus] ;
#endif
                }
              }
              
              // QUANTITATIVE RESISTANCE
              // The value of the traits for quantitative resistance have been already determined (trait_new_sample),
              // here we convert them in aggressiveness index and we combine the aggressiveness indices 
              // related to qualitative and quantitative resistances (in aggr_prop) 
  
              if(Nquanti > 0){
                std::vector<double> trait_prop_temp(Nquali_gene);
                for (int g = 0; g < Nquali_gene; g++){
                  // trait_prop_temp is a temporary vector containing the traits values for all resistances (MG + QTL)
                  // (the trait values related to MGs are "fake" but we need it for the swithc function)
                  trait_prop_temp[g] = trait_parents[0][g];
                }
                
                trait_prop_temp.insert(trait_prop_temp.end(), trait_prop[s].begin(),trait_prop[s].end());     
                std::vector<int> aggr_prop_temp(this->Ngene);
                // aggr_prop_temp is a temporary vector containing the aggressiveness indices
                // (the aggressiveness indices associated to MGs are "fake" and will be replaced in aggr_prop)
                aggr_prop_temp = this-> switch_trait_to_aggr(trait_prop_temp, resistance);
                aggr_prop.insert(aggr_prop.begin() + Nquali_gene, aggr_prop_temp.begin() + Nquali_gene, aggr_prop_temp.end());
                
#ifdef DEBUG
                  Rcpp::Rcout << "\n" << " traits values propagule" ;
                  for(int i=0; i < Nquanti; i++){
                    Rcpp::Rcout << " " << trait_prop[s][i] << " ";
                  }
                  
                  Rcpp::Rcout << "\n" << " aggressiveness indices (QTL)";
                  for(int j=0 + Nquali_gene; j < this-> Ngene; j++){
                    Rcpp::Rcout << " "<< aggr_prop_temp[j] << " ";
                  }
                  
                  Rcpp::Rcout << "\n" << " aggressiveness indices (MG + QTL)";
                  for(int j=0; j < this->Ngene; j++){
                    Rcpp::Rcout << " "<< aggr_prop[j] << "\n ";
                  }
#endif
               
              }
              
              // Increment P with the corresponding pathogen genotype 
              const int patho_new = this->switch_aggr_to_patho(aggr_prop);
#ifdef DEBUG
                Rcpp::Rcout << " " << "pathotype" << " " <<  patho_new << "\n" ;

#endif
              P[patho_new]++ ;
              }
            }
          }
        }
      }
    }
  }
}

/* ------------------------------------------------------- */
/*  DISTRIBUTION OF PRIMARY INOCULUM BETWEEN THE SEASONS   */
/* ------------------------------------------------------- */
/* Exponential distribution of sexual spores (P_sex_temp) between the growing seasons */
/* (sexual spore may stay dormant multiple seasons before germination)*/
/* Update P_stock in a given poly */

void Model::between_season_pr_inoc(std::vector<int>& P_sex_tmp, Vector2D<int>& P_stock, int& year){
  // sex_propagule_release_mean is the average number of cropping seasons after which a sexual propagule is released.
  // sex_propagule_viability_limit is the maximum number of cropping seasons up to which a sexual propagule is viable
  for(int patho = 0; patho<this->Npatho; patho++){
    for(int counter = 0; counter<P_sex_tmp[patho];counter++){
      int lag_germination = static_cast<int>(ran_exponential_trunc(this->basic_patho.sex_propagule_release_mean,
                                                                   this->basic_patho.sex_propagule_viability_limit));
      if(lag_germination >= 0 && lag_germination < this->basic_patho.sex_propagule_viability_limit) {
        P_stock[patho][(year - 1 + lag_germination)%this->basic_patho.sex_propagule_viability_limit]+=1;
      }
    }
  }
}

/* ----------------------------------------------------- */
/*  DISTRIBUTION OF PRIMARY INOCULUM WITHIN THE SEASON   */
/* ----------------------------------------------------- */
/* If distr=1: Uniform distribution of sexual/asexual spores that will germinate */
/*  in the considered season (P_stock_germ) within the growing season */
/* if distr=0: all the spore germinate the first day of the following season*/
/* Update P (which can be P_sex or P_asex) in a given poly */

void Model::in_season_pr_inoc(std::vector<int>& P_stock_germ, Vector2D<int>& P, const bool& distr){
  for(int patho = 0; patho<this->Npatho; patho++){
    for(int counter = 0; counter<P_stock_germ[patho];counter++){
      // the timestep t when the spore will germinate is sampled from a uniform distribution between [0,120)
      // int lag_sporulation_seas = static_cast<int>((this->time_steps_per_year - 0)*rng_uniform()+0);
      int lag_sporulation_seas = 0; // all the spores germinate the first day of the following year
      if(distr > 0){
        lag_sporulation_seas = static_cast<int>((this->time_steps_per_year - 0)*rng_uniform()+0);
      }
      if(lag_sporulation_seas >= 0 && lag_sporulation_seas < this->time_steps_per_year){
        P[patho][lag_sporulation_seas]+=1;
      }
    }
  }
}

/* -------------------------------------------------- */
/*      MUTATION OF THE PROPAGULES AT A GIVEN LOCUS   */
/* -------------------------------------------------- */

/* Update PpathoMut with mutation on "trait_mut" through multinomial draw */
void Model::mutation_locus(const int& patho, const int& trait_mut, Vector2D<int>& PpathoMut) {
  const int Nlevels = this->genes[trait_mut].Nlevels_aggressiveness;
  Vector2D<int> PaggrMut(this->Npatho, std::vector<int>(Nlevels));
  
  /* Mutation of trait_mut */
  for(int patho_old = 0; patho_old < this->Npatho; patho_old++) {
    const std::vector<int> aggr_old = this->switch_patho_to_aggr(patho_old);
    const int aggr_to_mutate = aggr_old[trait_mut]; // Aggressiveness index before mutation
    PaggrMut[patho_old] = 
      ran_multinomial(PpathoMut[patho][patho_old], this->genes[trait_mut].mutkernel[aggr_to_mutate]);
    PpathoMut[patho][patho_old] = 0; // Re-initialisation of PpathoMut
  }
  
  /* Update PpathoMut with mutants */
  for(int patho_old = 0; patho_old < this->Npatho; patho_old++) {
    // Aggressiveness index relative to the different traits
    const std::vector<int> aggr_old = this->switch_patho_to_aggr(patho_old);
    for(int aggr_mut = 0; aggr_mut < Nlevels; aggr_mut++) {
      std::vector<int> aggr_new = aggr_old;
      aggr_new[trait_mut] = aggr_mut;
      const int id_patho_mut = this->switch_aggr_to_patho(aggr_new); // Pathotype index after mutation
      PpathoMut[patho][id_patho_mut] += PaggrMut[patho_old][aggr_mut]; // Add in index of new pathotype
    }
  }
}

/* -------------------------------------- */
/*               MUTATION                 */
/* -------------------------------------- */

/* Update P with mutation of propagules in a given poly */
void Model::mutation(std::vector<int>& P) {
  /* P is the numbers of individual propagules in a given poly */
  //Rcpp::Rcout  << "\n" << "Mutation" << " ";
  Vector2D<int> PpathoMut(this->Npatho, std::vector<int>(this->Npatho, 0));
  for(int patho = 0; patho < this->Npatho; patho++) {
    PpathoMut[patho][patho] = P[patho];
    /* Mutation of the different genes */
    for(int g = 0; g < this->Ngene; g++) {
      this->mutation_locus(patho, g, PpathoMut);
    }
  }
  
  /* Update P after mutation */
  for(int patho_mut = 0; patho_mut < this->Npatho; patho_mut++) {
    P[patho_mut] = 0;
    for(int patho = 0; patho < this->Npatho; patho++) {
      P[patho_mut] += PpathoMut[patho][patho_mut];
    }
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*    functions affecting all poly    */

/* --------------------------------------- */
/*               DISPERSAL                 */
/* --------------------------------------- */
/* Dispersal of new host and new pathogen propagules in the landscape */

/* Update Hjuv and P */
void Model::dispersal(const Vector2D<int>& H, Vector2D<int>& Hjuv, Vector2D<int>& P, const Vector2D<double>& disp_matrix) {
  /* H, Hjuv, P are the numbers of individuals in a given poly */
  Vector3D<int> Pdisp(this->Npatho, Vector2D<int>(this->Npoly, std::vector<int>(this->Npoly)));
  Vector3D<int> Hjuvtmp(this->Nhost, Vector2D<int>(this->Npoly, std::vector<int>(this->Npoly)));
  
  /* Production and dispersal of the host & dispersal of pathogen propagules */
  for(int poly = 0; poly < this->Npoly; poly++) {
    /* Pathogen dispersal */
    for(int patho = 0; patho < this->Npatho; patho++) {
      Pdisp[patho][poly] = ran_multinomial(P[poly][patho], disp_matrix[poly]);
    }
    /* Host reproduction: production and dispersal of Hjuv */
    for(int host = 0; host < this->Nhost; host++) {
      Hjuvtmp[host][poly] = ran_multinomial(
        static_cast<int>(this->cultivars[host].reproduction_rate * H[poly][host]), this->disp_host[poly]);
    }
  }
  
  for(int poly = 0; poly < this->Npoly; poly++) {
    /* Update the number of propagules (P) and Hjuv landing in each field */
    for(int patho = 0; patho < this->Npatho; patho++) {
      P[poly][patho] = 0;
      for(int polyE = 0; polyE < this->Npoly; polyE++) {
        P[poly][patho] += Pdisp[patho][polyE][poly];
      }
    }
    
    for(int host = 0; host < this->Nhost; host++) {
      Hjuv[poly][host] = 0;
      for(int polyE = 0; polyE < this->Npoly; polyE++) {
        Hjuv[poly][host] += Hjuvtmp[host][polyE][poly];
      }
    }
  }
}

/* --------------------------------------------- */
/*       BOTTLENECK AT THE END OF THE SEASON     */
/* --------------------------------------------- */
Vector3D<int> Model::bottleneck(const int& t, const Vector3D<int>& L, const Vector3D<int>& I,
                                const Vector2D<int>& activeR) {
  Vector3D<int> eqIsurv(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost)));
  
  for(int patho = 0; patho < this->Npatho; patho++) {
    const std::vector<int> aggr = this->switch_patho_to_aggr(patho);
    for(int host = 0; host < this->Nhost; host++) {
      for(int poly = 0; poly < this->Npoly; poly++) {
        /* Reduce the number of infected hosts (bottleneck) */
        eqIsurv[poly][patho][host] =
        ran_binomial(this->basic_patho.survival_prob, L[poly][patho][host] + I[poly][patho][host]);
        /* Calculate the mean infectious period */
        double infectious_period_mean_new = this->basic_patho.infectious_period_mean;
        for(int g = 0; g < this->Ngene; g++) {
          if(this->genes[g].target_trait == "IP") {
            // Indicate if the cultivar has an active resistance gene
            const bool Rgene = this->get_resistance(g, host, t, activeR[patho][g]);
            infectious_period_mean_new *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
          }
        }
        /* Calculate the equivalent number of infectious hosts */
        eqIsurv[poly][patho][host] *= static_cast<int>(infectious_period_mean_new);
      }
    }
  }
  return eqIsurv;
}


/* -------------------------- */
/*         HOST DYNAMIC       */
/* -------------------------- */

/* Compute host reproduction, death and growth and updtate the number of H in the concerned poly */
void Model::host_dynamic(const int& poly, const int& year, const int& t, std::vector<int>& H, std::vector<int>& Hjuv,
                         const Vector2D<int>& L, const Vector2D<int>& I, const Vector2D<int>& R, std::vector<int>& N, std::vector<int>& Nspray) {
  /* H, Hjuv, L, I and R are the number of host in a given poly */
  
  // If there is no rotation (same croptype each year), only take the first year rotation
  const int id_croptype = (this->rotation[poly].size() == 1) ? this->rotation[poly][0] : this->rotation[poly][year];
  
  for(std::pair<int, double> cultivar_prop : this->croptypes[id_croptype].cultivar_proportion) {
    const int id_host = cultivar_prop.first;
    const double prop = cultivar_prop.second;
    
    /* Calculation of totals for L, I, R and K */
    // Carrying capacity of the cultivar in the concerned paddock
    const int K = static_cast<int>(this->area[poly] * this->cultivars[id_host].max_density * prop); 
    int L_host = 0, I_host = 0, R_host = 0;
    for(int patho = 0; patho < this->Npatho; patho++) {
      L_host += L[patho][id_host];
      I_host += I[patho][id_host];
      R_host += R[patho][id_host];
    }
    N[id_host] = H[id_host] + L_host + I_host + R_host; // Number of occupied sites
    // N[id_host] = H[id_host] + L_host + I_host ; // Number of occupied sites IN THE CASE OF INFINITIVE HOST GROWTH
    
    /* HOST MORTALITY: H2M */
    /* ------------------- */
    const int H2M = ran_binomial(this->cultivars[id_host].death_rate, H[id_host]);
    const int H2Mjuv = ran_binomial(this->cultivars[id_host].death_rate, Hjuv[id_host]);
    Hjuv[id_host] -= H2Mjuv; // Update Hjuv
    
    /* HOST REPRODUCTION: Hnewind */
    /* -------------------------- */
    /* Hjuv settlement in the field */
    int availSites = K - N[id_host] - H2M; // Number of available sites
    if(availSites < 0) { /* Security */
    availSites = 0;
    }
    const double f1host = sigmoid(this->sigmoid_plateau_host, this->sigmoid_kappa_host, this->sigmoid_sigma_host,
                                  availSites / static_cast<double>(K));
    int siteaccess = ran_binomial(f1host, availSites);
    const int Hnewind = (siteaccess < Hjuv[id_host]) ? siteaccess : Hjuv[id_host];
    
    /* HOST GROWTH: H2H */
    /* ---------------- */
    int H2H = static_cast<int>(this->cultivars[id_host].growth_rate * (H[id_host] - H2M) *
    (1 - ((N[id_host] - H2M + Hnewind) / static_cast<double>(K))));
    if(H2H < 0) { /* Security */
    Rcpp::Rcerr << "hostID" << id_host << " growrate " << this->cultivars[id_host].growth_rate << " H " << H[id_host] 
                << " H2M " << H2M << " N " << N[id_host] << " HnewId " << Hnewind << " K " << static_cast<double>(K) 
                << std::endl;
    Rprintf("CAREFUL ! H2H < 0, one of the areas may be 0: check if Npoly, NpolyTot and idLAN are correct\n");
    H2H = 0;
    } else if((N[id_host] - H2M + Hnewind + H2H) > K) {
      Rprintf("CAREFUL ! H2H too big\n");
      H2H = K - (N[id_host] - H2M + Hnewind);
    }
    
    /* UPDATE NUMBER OF HOSTS */
    /* ---------------------- */
    H[id_host] = H[id_host] - H2M + Hnewind + H2H;
    N[id_host] = N[id_host] - H2M + Hnewind + H2H;
    
    /* SAVE Nspray (i.e. the number of host at treatment dates) */
    /* This information is needed to compute the effect of pesticide treatments */
    /* ------------------------------------------------------------------------- */ 
    if (std::find(this->treatment.treatment_timesteps.begin(), this->treatment.treatment_timesteps.end(), t) != this->treatment.treatment_timesteps.end()) {
      Nspray[id_host] = N[id_host];
    }
  }
}
/* ------------------------------------------------------------------ */
/*         CONTAMINATION : propagule deposition on healthy sites      */
/* ------------------------------------------------------------------ */

/* Calculation of the number of contaminated sites */
Vector2D<int> Model::contamination(const std::vector<int>& H, const std::vector<int>& P, const std::vector<int>& N) {
  /* H, P are the number of individuals in a given poly */
  Vector2D<int> Hcontaminated(this->Npatho, std::vector<int>(this->Nhost));
  
  /* Calculation of total for H and P */
  const int Htot = std::accumulate(H.begin(), H.end(), 0);
  const int Ptot = std::accumulate(P.begin(), P.end(), 0);
  
  std::vector<double> probaH(this->Nhost + 1); // Probability for H to belong to each host
  double probaHtot = 0;
  for(int host = 0; host < this->Nhost; host++) {
    probaH[host] = (Htot == 0) ? 0 : static_cast<double>(H[host]) / static_cast<double>(Htot);
    // (cond) ? res1 : res2;
    probaHtot += probaH[host];
  }
  probaH[this->Nhost] = 1 - probaHtot;
  
  std::vector<double> probaP(this->Npatho + 1); // Probability to belong to the different pathotypes
  double probaPtot = 0;
  for(int patho = 0; patho < this->Npatho; patho++) {
    // If there is no propagule, set everything to 0, else put each proportion
    probaP[patho] = (Ptot == 0) ? 0 : static_cast<double>(P[patho]) / static_cast<double>(Ptot);
    probaPtot += probaP[patho];
  }
  probaP[this->Npatho] = 1 - probaPtot;
  
  // Distribution of the propagules among the different cultivars
  const std::vector<int> P_host = ran_multinomial(Ptot, probaH);
  
  for(int host = 0; host < this->Nhost; host++) {
    // Distribution of the propagules among the different pathotypes
    const std::vector<int> P_host_patho = ran_multinomial(P_host[host], probaP);
    /* Calculation of the number of contaminable sites */
    const double f1patho = (N[host] > 0)
      ? sigmoid(this->basic_patho.sigmoid_plateau, this->basic_patho.sigmoid_kappa,
                this->basic_patho.sigmoid_sigma, (H[host] / static_cast<double>(N[host])))
        : 0.0;
    
    // Contaminable site, where a propagule may deposit
    const int Hcontaminable = ran_binomial(f1patho, H[host]);
    // Distribution of the contaminable sites among the different pathotypes
    const std::vector<int> Hcontaminated_tmp = ran_multinomial(Hcontaminable, probaP);
    /* The true number of contaminated sites is the minimum between sites and propagules */
    for(int patho = 0; patho < this->Npatho; patho++) {
      Hcontaminated[patho][host] = std::min(Hcontaminated_tmp[patho], P_host_patho[patho]);
    }
  }
  return Hcontaminated;
}

/* ----------------------------------------------------------- */
/*         INFECTIOUS CYCLE : transitions H -> L -> I -> R     */
/* ----------------------------------------------------------- */

/* Calculate the number of contaminated sites that become infected, infectious or removed and update H, L, I, R */
void Model::infection(const int& t, std::vector<int>& H, const Vector2D<int>& Hcontaminated, Vector2D<int>& L,
                      Vector2D<int>& I, Vector2D<int>& R, Vector3D<int>& L2I, Vector3D<int>& I2R,
                      const std::vector<int>& activeR, const std::vector<int>& N,  const std::vector<int>& Nspray) {
  /* H, Hcontaminated, L, I, R, L2I and I2R are the number of individuals in a given poly */
  for(int patho = 0; patho < this->Npatho; patho++) {
    const std::vector<int> aggr = this->switch_patho_to_aggr(patho);
    for(int host = 0; host < this->Nhost; host++) {
      /* Infection of healthy sites: H2L */
      /* ------------------------------- */
      double infection_rate_exp = this->basic_patho.infection_rate;
      /* Interaction with resistance genes */
      for(int g = 0; g < this->Ngene; g++) {
        if(this->genes[g].target_trait == "IR") {
          // Indicate if the cultivar has an active resistance gene
          bool Rgene = this->get_resistance(g, host, t, activeR[g]);
          infection_rate_exp *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
        }
      }
      
      /* Interaction with pesticide treatment*/
      
      if (std::find(this->treatment.treatment_cultivars.begin(), this->treatment.treatment_cultivars.end(), host) != this->treatment.treatment_cultivars.end()) {
        infection_rate_exp *= this->get_treat_effect(N[host], Nspray[host], t);
      }
      const int H2L = ran_binomial(infection_rate_exp, Hcontaminated[patho][host]);
      
      /* Latent and infectious periods */
      /* ----------------------------- */
      /* Find parameters of gamma distributions from mean and variance */
      double latent_period_mean_new = this->basic_patho.latent_period_mean;
      double infectious_period_mean_new = this->basic_patho.infectious_period_mean;
      std::array<double, 2> latent_period_alpha;
      std::array<double, 2> infectious_period_alpha;
      for(int g = 0; g < this->Ngene; g++) {
        /* Latent period */
        if(this->genes[g].target_trait == "LAT") {
          // Indicate if the cultivar has an active resistance gene
          bool Rgene = this->get_resistance(g, host, t, activeR[g]);
          latent_period_mean_new /= (this->genes[g].aggressiveness_matrix[aggr[g]][Rgene] +
            0.001 * (this->genes[g].aggressiveness_matrix[aggr[g]][Rgene] == 0));
        }
        /* Security to avoid problem in alpha calculation */
        latent_period_mean_new += 0.001 * (latent_period_mean_new == 0);
        
        /* Infectious period */
        if(this->genes[g].target_trait == "IP") {
          // Indicate if the cultivar has an active resistance gene
          bool Rgene = this->get_resistance(g, host, t, activeR[g]);
          infectious_period_mean_new *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
        }
        // Security to avoid problem in alpha calculation
        infectious_period_mean_new += 0.001 * (infectious_period_mean_new == 0);
      }
      latent_period_alpha = find_paramGamma(latent_period_mean_new, this->basic_patho.latent_period_var);
      infectious_period_alpha = find_paramGamma(infectious_period_mean_new, this->basic_patho.infectious_period_var);
      
      /* Recently infected hosts */
      for(int h2l = 0; h2l < H2L; h2l++) {
        /* Latent period */
        const int lag1 = static_cast<int>(ran_gamma(latent_period_alpha[0], latent_period_alpha[1]));
        if((t + lag1) < this->time_steps_per_year) {
          L2I[patho][host][t + lag1]++;
        }
        /* Infectious period */
        const int lag2 = lag1 + static_cast<int>(ran_gamma(infectious_period_alpha[0], infectious_period_alpha[1]));
        if((t + lag2) < this->time_steps_per_year) {
          I2R[patho][host][t + lag2]++;
        }
      }
      
      /* Update H, L, I, R */
      /* ----------------- */
      H[host] -= H2L;
      L[patho][host] += H2L;
      L[patho][host] -= L2I[patho][host][t];
      I[patho][host] += L2I[patho][host][t];
      I[patho][host] -= I2R[patho][host][t];
      R[patho][host] += I2R[patho][host][t];
    }
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

/* --------------------------- */
/*   DYNAMIC OF THE EPIDEMIC   */
/* --------------------------- */

void Model:: dynepi() {
  using namespace std::chrono;
  auto start_tot = high_resolution_clock::now();  
  char name_fH[20], name_fHjuv[20], name_fP[20], name_fL[20], name_fI[20], name_fR[20], name_fPbefinter[20];
  Vector2D<int> Hcontaminated; // Contaminated sites (where a propagule is deposited)
  Vector2D<int> Hjuv, P;
  Vector3D<int> L(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost)));
  Vector3D<int> I(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost)));
  Vector3D<int> R(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost)));
  std::vector<int> N(this->Nhost);
  Vector2D<int> Nspray(this->Npoly, std::vector<int>(this->Nhost));
  std::vector<int> Nlevels_aggressiveness(this->Ngene);
  Vector4D<int> L2I;
  Vector4D<int> I2R;
  Vector2D<int> P_sex_tmp; // P_sex_tmp contains the propagules produced sexually
  Vector2D<int> P_asex_tmp; // P_asex_tmp contains the propagules produced asexually
  Vector2D<int> P_bef_interseas; 
  // P_bef_interseas contains the propagules produced sexually AND asexually before the interseason
  Vector3D<int> P_stock(this->Npoly,Vector2D<int>(this->Npatho, std::vector<int>(this->basic_patho.sex_propagule_viability_limit)));
  // P_stock contains the propagules that will germinates in each year within the "sex_propagule_viability_limit"
  Vector3D<int> P_sex(this->Npoly,Vector2D<int>(this->Npatho, std::vector<int>(this->time_steps_per_year)));
  // P_sex cointains the propagules issued from sexual reproduction after the bottleneck
  // that will germinate at each time of the time_steps_per_year
  Vector3D<int> P_asex(this->Npoly,Vector2D<int>(this->Npatho, std::vector<int>(this->time_steps_per_year)));
  // P_asex cointains the propagules issued from asexual reproduction after the bottleneck 
  // that will germinate at each time of the time_steps_per_year
  Vector2D<int> P_sex_daily (this->Npoly, std::vector<int>(this->Npatho));
  // P_sex_daily contains the propagule issued from sexual reproduction after the bottleneck
  // that germinate in a specific time of the time_steps_per_year
  Vector2D<int> P_asex_daily (this->Npoly, std::vector<int>(this->Npatho));
  // P_asex_daily contains the propagule issued from asexual reproduction after the bottleneck
  // that germinate in a specific time of the time_steps_per_year
  
  /* Initialisation (t=0) */
  this->init_HjuvLIR(Hjuv, L, I, R);
  this->init_P(P, P_sex_tmp, P_asex_tmp, P_sex, P_asex, P_sex_daily, P_asex_daily, P_stock); // initialised a 0
  this->init_Nspray(Nspray); // initialised a 0
  this->init_Nlevels_aggressiveness(Nlevels_aggressiveness);
  this->init_L2I2R(L2I, I2R);
  Vector2D<int> H = this->intro_H(0);
  Vector2D<int> activeR = this->init_activeR();
  intro_I(H, I, I2R, activeR);
  
  /* Compute the number of Qualitative gene in Gene and check if they are all at the beginning of the sequence */
  /* (otherwise, the reproSex function would not work properly) */
  
  for(int locus = 0; locus < this->Ngene; locus++) {
    Nlevels_aggressiveness[locus] = genes[locus].Nlevels_aggressiveness; 
  }
  
  int Nquali_gene = std::count(Nlevels_aggressiveness.begin(),Nlevels_aggressiveness.end(),2); // number of qualitative genes in Ngene
  
  for(int i = 0; i < Nquali_gene; i++){
    if(Nlevels_aggressiveness[i] != 2){
      Rcpp::Rcout << "Attention: qualitative genes are not at the beginning of the gene sequences. \n"; 
    }  
  }

#ifdef DEBUG
  const int row = 1; // to print just a specific poly
#endif
  
  for(int year = 1; year < (this->Nyears + 1); year++) {
    auto start = high_resolution_clock::now(); 
    Rprintf("----------------------------- YEAR %d -----------------------------\n", year);
    /* Create the files to write the output */
    sprintf(name_fH, "H-%02d.bin", year);
    sprintf(name_fHjuv, "Hjuv-%02d.bin", year);
    sprintf(name_fP, "P-%02d.bin", year);
    sprintf(name_fL, "L-%02d.bin", year);
    sprintf(name_fI, "I-%02d.bin", year);
    sprintf(name_fR, "R-%02d.bin", year);
    sprintf(name_fPbefinter, "Pbefinter-%02d.bin", year);
    
    FILE* fH = fopen(name_fH, "wb");
    FILE* fHjuv = fopen(name_fHjuv, "wb");
    FILE* fL = fopen(name_fL, "wb");
    FILE* fI = fopen(name_fI, "wb");
    FILE* fP = fopen(name_fP, "wb");
    FILE* fR = fopen(name_fR, "wb");
    FILE* fPbefinter = fopen(name_fPbefinter, "wb");
    
    /* Loop for all the timesteps of the cropping season */
    for(int t = 1; t < this->time_steps_per_year; t++) {

#ifdef DEBUG
      Rcpp::Rcout <<"----------------------------- T "<< t << " -----------------------------" << std::endl;
#endif
      
      // Writing model output for timestep t
      this->write_HHjuvPLIR(H, Hjuv, P, L, I, R, fH, fHjuv, fP, fL, fI, fR);
      P = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0)); // Re-initialisation at 0
      P_sex_tmp = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0)); // Re-initialisation at 0
      P_asex_tmp = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0)); // Re-initialisation at 0
      for(int poly = 0; poly < this->Npoly; poly++) {
        if(this->basic_patho.repro_sex_prob[t] > 0 && this->basic_patho.repro_sex_prob[t] < 1){
          // the split between individuals doing sexual and clonal reproduction only takes place if
          // 0 < repro_sex_prob < 1
          const std::array<Vector2D<int>, 2> splited_I(this->split_IclonalIsex(t, I[poly]));
          this->reproClonal(t, P_asex_tmp[poly], splited_I[0], activeR[poly]); 
          this->mutation(P_asex_tmp[poly]); // assumption: mutation only takes place after clonal reproduction
          this->reproSex(t, P_sex_tmp[poly],  splited_I[1], activeR[poly], Nlevels_aggressiveness, Nquali_gene);
          // note: we assume that if repro sex takes place in season the propagule are treated like clonal propagule (except for the mutation)
          // they are added to the P matrix and they disperse like clonal propagule
        } else if(this->basic_patho.repro_sex_prob[t] == 0){
          // if repro_sex_prob = 0 all the individuals do clonal reproduction
          this->reproClonal(t, P_asex_tmp[poly], I[poly], activeR[poly]); // assumption: mutation only takes place after clonal reproduction
          this->mutation(P_asex_tmp[poly]);
        } else {
          // if repro_sex_prob = 1 all the individuals do sexual reproduction
          this->reproSex(t, P_sex_tmp[poly], I[poly], activeR[poly], Nlevels_aggressiveness, Nquali_gene);
          // note: we assume that if repro sex takes place in season the propagule are treated like clonal propagule (except for the mutation)
          // they are added to the P matrix and they disperse like clonal propagule
        }
      }
      
      // Get the asexual spores produced after the bottleneck that will germinate at time step t  
      this->get_P_daily(P_asex_daily, P_asex, t); 
      
#ifdef DEBUG
        Rcpp::Rcout << "P asex germinated at time" << t << "\n ";
        for(unsigned int i = 0; i < P_asex_daily[row].size(); i++){
          Rcpp::Rcout <<P_asex_daily[row][i]<< " ";
        }
        Rcpp::Rcout << "\n";
        
        Rcpp::Rcout << "P asex produced at time" << t << "\n ";
        for(unsigned int i = 0; i < P_asex_tmp[row].size(); i++){
          Rcpp::Rcout <<P_asex_tmp[row][i]<< " ";
        }
        Rcpp::Rcout << "\n";
#endif      
      // sum the asexual spores germinated on day t (primary inoculum) with the sexual 
      // spore produced on day t (secondary inouculum)
      
      P_asex_daily = this->get_sum_Vector2D(P_asex_daily, P_asex_tmp);
#ifdef DEBUG     
      Rcpp::Rcout << "P asex TOT" << t << "\n ";
      for(unsigned int i = 0; i < P_asex_daily[row].size(); i++){
        Rcpp::Rcout <<P_asex_daily[row][i]<< " ";
      }
      Rcpp::Rcout << "\n";
#endif     
      this->dispersal(H, Hjuv, P_asex_daily, this->disp_patho); // dispersion of asexual spore 
      
      // Get the sexual spores that will germinate at time step t  
      this->get_P_daily(P_sex_daily, P_sex, t); 
      // sum the sexual spores germinated on day t (primary inoculum)
      // with the sexual spore produced on day t (secondary inoculum)
      
      P_sex_daily = this->get_sum_Vector2D(P_sex_daily, P_sex_tmp);
      
      this->dispersal(H, Hjuv, P_sex_daily,this->disp_patho_sex); // dispersion of sexual spore
      
      // Update the number of spores (sex + asex) after dispersal
      P = this->get_sum_Vector2D(P_sex_daily, P_asex_daily);
      
      for(int poly = 0; poly < this->Npoly; poly++) {
        //Rprintf("----------------------------- Poly %d -----------------------------\n",poly);
        this->host_dynamic(poly, year - 1, t, H[poly], Hjuv[poly], L[poly], I[poly], R[poly], N, Nspray[poly]);
        Hcontaminated = this->contamination(H[poly], P[poly], N);
        this->infection(t, H[poly], Hcontaminated, L[poly], I[poly], R[poly], L2I[poly], I2R[poly],
                        activeR[poly], N,  Nspray[poly]);
      }
    }
    
    /* Last time-step of the season: bottleneck before starting a new season */
    
    // Writing model output for last timestep
    this->write_HHjuvPLIR(H, Hjuv, P, L, I, R, fH, fHjuv, fP, fL, fI, fR);
    // Calculation of the equivalent number of I that survive and produce propagules for the next season
    const Vector3D<int> eqIsurv = this->bottleneck(this->time_steps_per_year, L, I, activeR);
    
    /* Re-initialisation at 0 */
    this->init_HjuvLIR(Hjuv, L, I, R);
    this->init_L2I2R(L2I, I2R);
    P = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
    P_sex_tmp = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
    P_asex_tmp = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
    P_bef_interseas = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
    P_sex = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->time_steps_per_year, 0)));
    P_asex = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->time_steps_per_year, 0)));
    

    /* Generate P issued from eqIsurv = (remaining L+I) * Tspo */
    
    for(int poly = 0; poly < this->Npoly; poly++) {
      //Rprintf("----------------------------- Poly %d -----------------------------\n",poly);
      if(this->basic_patho.repro_sex_prob[time_steps_per_year] > 0 && this->basic_patho.repro_sex_prob[time_steps_per_year] < 1){
        // the split between individuals doing sexual and clonal reproduction only takes place if
        // 0 < repro_sex_prob < 1
        const std::array<Vector2D<int>, 2> splited_I(this->split_IclonalIsex(this->time_steps_per_year, eqIsurv[poly]));
        this->reproClonal(this->time_steps_per_year, P_asex_tmp[poly], splited_I[0], activeR[poly]);
        this->mutation(P_asex_tmp[poly]); // assumption: mutation only takes place after clonal reproduction
        this->reproSex(this->time_steps_per_year, P_sex_tmp[poly],splited_I[1], activeR[poly], Nlevels_aggressiveness, Nquali_gene);

#ifdef DEBUG
        if (poly==row){
          Rcpp::Rcout << "P sex produced" << "\n ";
          for(unsigned int i = 0; i < P_sex_tmp[poly].size(); i++){
            Rcpp::Rcout <<P_sex_tmp[poly][i]<< " ";
          }
          Rcpp::Rcout << "\n";
        }
        if (poly==row){
          for(unsigned int z = 0; z < P_stock[0][0].size(); z++){
            for(unsigned int c = 0; c < P_stock[0].size(); c++){
              Rcpp::Rcout <<P_stock[poly][c][z] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
        // Distribution of sexual spores between-seasons 
        this->between_season_pr_inoc(P_sex_tmp[poly], P_stock[poly], year);

#ifdef DEBUG
        if (poly==row){
          for(unsigned int z = 0; z < P_stock[0][0].size(); z++){
            Rcpp::Rcout << "YEAR of GERMINATION" << " " <<  z << "\n ";
            for(unsigned int c = 0; c < P_stock[0].size(); c++){
              Rcpp::Rcout <<P_stock[poly][c][z] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
        
        // Get the sexual spores that germinate in the following season
        std::vector<int> P_stock_germ = this->get_P_stock_germ(P_stock[poly], year);

#ifdef DEBUG
        if (poly==row){
          Rcpp::Rcout << "Germination in the following year" << "\n ";
          for(unsigned int count = 0; count < P_stock_germ.size(); count++){
            Rcpp::Rcout << P_stock_germ[count] << " ";
          }
          Rcpp::Rcout << "\n";
        }
        if (poly==row){
          for(unsigned int z = 0; z < P_stock[0][0].size(); z++){
            Rcpp::Rcout << "YEAR of GERMINATION" << " " <<  z << "\n ";
            for(unsigned int c = 0; c < P_stock[0].size(); c++){
              Rcpp::Rcout <<P_stock[poly][c][z] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
        // Distribution of sexual spores that germinate in 
        // the following season within the season time-steps
        this->in_season_pr_inoc(P_stock_germ, P_sex[poly], 1);
        
        // Distribution of asexual spores that germinate in 
        // the following season within the season time-steps
        this->in_season_pr_inoc(P_asex_tmp[poly], P_asex[poly], this->basic_patho.clonal_propagule_gradual_release);
#ifdef DEBUG
        if (poly==row){
          for(unsigned int r = 0; r < P_sex[0].size(); r++){
            for(unsigned int c = 0; c < P_sex[0][0].size(); c++){
              Rcpp::Rcout << P_sex[poly][r][c] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
      } else if(this->basic_patho.repro_sex_prob[time_steps_per_year] == 0){
        // if repro_sex_prob = 0, only clonal reproduction takes place
        this->reproClonal(this->time_steps_per_year, P_asex_tmp[poly], eqIsurv[poly], activeR[poly]);
        this->mutation(P_asex_tmp[poly]); // assumption: mutation only takes place after clonal reproduction
        
        // Distribution of asexual spores that germinate in 
        // the following season within the season time-steps
        this->in_season_pr_inoc(P_asex_tmp[poly], P_asex[poly],this->basic_patho.clonal_propagule_gradual_release);
        
#ifdef DEBUG
  if (poly==row){
    Rcpp::Rcout << "P asex produced" << "\n ";
    for(unsigned int i = 0; i < P_asex_tmp[poly].size(); i++){
      Rcpp::Rcout <<P_asex_tmp[poly][i]<< " ";
    }
    Rcpp::Rcout << "\n";
  }

        if (poly==row){
          Rcpp::Rcout << "Repro sex \n";
          for(unsigned int r = 0; r < P_sex[0].size(); r++){
            for(unsigned int c = 0; c < P_sex[0][0].size(); c++){
              Rcpp::Rcout << P_sex[poly][r][c] << " ";
            }
            Rcpp::Rcout << "\n";
          }
          Rcpp::Rcout << "Repro asex \n";
          for(unsigned int r = 0; r < P_asex[0].size(); r++){
            for(unsigned int c = 0; c < P_asex[0][0].size(); c++){
              Rcpp::Rcout << P_asex[poly][r][c] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
        
      } else {
        // if repro_sex_prob = 1, only sexual reproduction takes place
        this->reproSex(this->time_steps_per_year, P_sex_tmp[poly],eqIsurv[poly], activeR[poly], Nlevels_aggressiveness, Nquali_gene);
        
#ifdef DEBUG   
          if (poly==row){
          Rcpp::Rcout << "P sex produced" << "\n ";
            for(unsigned int i = 0; i < P_sex_tmp[poly].size(); i++){
              Rcpp::Rcout <<P_sex_tmp[poly][i]<< " ";
            }
           Rcpp::Rcout << "\n";
          }
        
 
        if (poly==row){
          for(unsigned int z = 0; z < P_stock[0][0].size(); z++){
            for(unsigned int c = 0; c < P_stock[0].size(); c++){
              Rcpp::Rcout <<P_stock[poly][c][z] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
        // Distribution of sexual spores between-seasons 
        this->between_season_pr_inoc(P_sex_tmp[poly], P_stock[poly], year);
#ifdef DEBUG
        if (poly==row){
          for(unsigned int z = 0; z < P_stock[0][0].size(); z++){
            Rcpp::Rcout << "YEAR of GERMINATION" << " " <<  z << "\n ";
            for(unsigned int c = 0; c < P_stock[0].size(); c++){
              Rcpp::Rcout <<P_stock[poly][c][z] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
        // Get the sexual spores that germinate in the following season
        std::vector<int> P_stock_germ = this->get_P_stock_germ(P_stock[poly], year);
#ifdef DEBUG
        if (poly==row){
          Rcpp::Rcout << "Germination in the following year" << "\n ";
          for(unsigned int count = 0; count < P_stock_germ.size(); count++){
            Rcpp::Rcout << P_stock_germ[count] << " ";
          }
          Rcpp::Rcout << "\n";
        }
        if (poly==row){
          for(unsigned int z = 0; z < P_stock[0][0].size(); z++){
            Rcpp::Rcout << "YEAR of GERMINATION" << " " <<  z << "\n ";
            for(unsigned int c = 0; c < P_stock[0].size(); c++){
              Rcpp::Rcout <<P_stock[poly][c][z] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
       
        // Distribution of sexual spores that germinate in 
        // the following season within the season time-steps
        this->in_season_pr_inoc(P_stock_germ, P_sex[poly],1);
#ifdef DEBUG
        if (poly==row){
          Rcpp::Rcout << "Repro sex \n";
          for(unsigned int r = 0; r < P_sex[0].size(); r++){
            for(unsigned int c = 0; c < P_sex[0][0].size(); c++){
              Rcpp::Rcout << P_sex[poly][r][c] << " ";
            }
            Rcpp::Rcout << "\n";
          }
          Rcpp::Rcout << "Repro asex \n";
          for(unsigned int r = 0; r < P_asex[0].size(); r++){
            for(unsigned int c = 0; c < P_asex[0][0].size(); c++){
              Rcpp::Rcout << P_asex[poly][r][c] << " ";
            }
            Rcpp::Rcout << "\n";
          }
        }
#endif
      }
      
    }
    
    // Writing P_before_interseason (sum of propagules produced by sexual AND 
    // asexual reproduction before the interseason)
    P_bef_interseas = this->get_sum_Vector2D(P_sex_tmp, P_asex_tmp);
    this->write_Pbefinter(P_bef_interseas, fPbefinter);
    
    //Spore dispersal
    // Get the asexual spores that will germinate in the first day (t=0) of the following season 
    this->get_P_daily(P_asex_daily, P_asex, 0); 
    
    this->dispersal(H, Hjuv, P_asex_daily, this->disp_patho); // dispersion of asexual spore produced in the fall
    // and that will germinate the first day of the following year
    
    // Get the sexual spores that will germinate in the first day (t=0) of the following season 
    this->get_P_daily(P_sex_daily, P_sex, 0); 
    
    this->dispersal(H, Hjuv, P_sex_daily,this->disp_patho_sex); // dispersion of sexual spore produced in the fall 
    // and that will germinate the first day of the following year

    
    // Update the number of spore (sex + asex) after dispersal
    P = this->get_sum_Vector2D(P_sex_daily, P_asex_daily);
    
    /* Re-plantation --> regenerate H */
    H = this->intro_H(year);
    activeR = this->init_activeR();
    // Nspray is re-initializated at 0, i.e. there is no fungicide left on plant tissue
    this->init_Nspray(Nspray);
    
    /* Infection of newly planted hosts to generate the inoculum of the next season */
    for(int poly = 0; poly < this->Npoly; poly++) {
      /* N = H[poly] in beginning of next season */ 
      Hcontaminated = this->contamination(H[poly], P[poly], H[poly]);
      this->infection(0, H[poly], Hcontaminated, L[poly], I[poly], R[poly], L2I[poly], I2R[poly], activeR[poly],
                      N,  Nspray[poly]);
    }
    
    fclose(fH);
    fclose(fHjuv);
    fclose(fL);
    fclose(fI);
    fclose(fP);
    fclose(fR);
    fclose(fPbefinter);
    
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
#ifdef DEBUG
    Rcpp::Rcout << "computational time" << " " << duration.count()<< " " << "seconds. \n"; 
#endif
  }
  auto stop_tot = high_resolution_clock::now(); 
  auto duration_tot = duration_cast<seconds>(stop_tot - start_tot); 
  Rcpp::Rcout << "total computational time" << " " << duration_tot.count()<< " " << "seconds. \n"; 
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void model_landsepi(Rcpp::List time_param, Rcpp::NumericVector area_vector, Rcpp::NumericMatrix rotation_matrix, Rcpp::NumericMatrix croptypes_cultivars_prop, Rcpp::List dispersal, Rcpp::List inits, int seed,
                    Rcpp::List cultivars_param, Rcpp::List basic_patho_param, Rcpp::List genes_param, 
                    Rcpp::List treatment_param) {
  
// To active Debug mode
// CLI or in projet build configuration
// add in 'R CMD INSTALL' this parameters : -d --configure-args="--enable-debug"
// or in Makevars.win PKG_CPPFLAGS=-UNDEBUG -DDEBUG=1
#ifdef DEBUG
  Rcpp::Rcerr << "### C++ Debug mode actived" <<std::endl;
#endif

  /*------------------*/
  /* Seed generation  */
  /*------------------*/
  /* Makoto Matsumoto and Takuji Nishimura, “Mersenne Twister: A 623-dimensionally equidistributed uniform
   * pseudorandom number generator”. ACM Transactions on Modeling and Computer Simulation, Vol. 8, No. 1 (Jan. 1998),
   * Pages 3–30 */
  const gsl_rng_type* gen_type = gsl_rng_mt19937;
  // This function returns a pointer to a newly-created instance of a random number generator of type gen_type
  const gsl_rng* gen = gsl_rng_alloc(gen_type);
  
  /*-------------------*/
  /*  Time parameters  */
  /*-------------------*/
  const int Nyears = Rcpp::as<int>(time_param["Nyears"]);
  const int nTSpY = Rcpp::as<int>(time_param["nTSpY"]);
  
  /*-------------------------------------*/
  /*  Landscape & deployment parameters  */
  /*-------------------------------------*/
  
  std::vector<double> area = Rcpp::as<std::vector<double>>(area_vector);
  Vector2D<int> rotation(0);
  for( int i = 0; i < rotation_matrix.nrow() ; i++) {
    Rcpp::NumericVector poly_rotation = rotation_matrix.row(i);
    //poly_rotation.push_back(poFeature_rotation->GetFieldAsInteger(year));
    rotation.push_back(Rcpp::as<std::vector<int>>(poly_rotation));
  }
  
  const int Npoly = area.size();
  
  std::vector<Croptype> croptypes(0);
  
  for(int i=0 ; i < croptypes_cultivars_prop.nrow() ; i++) {
    if( croptypes_cultivars_prop(i,0) >= croptypes.size()) croptypes.push_back(Croptype());
    croptypes[croptypes_cultivars_prop(i,0)].cultivar_proportion.push_back(
        std::pair<int, double>(
            croptypes_cultivars_prop(i,1),
            croptypes_cultivars_prop(i,2) ));                        
  }
  
  /*------------------------*/
  /*  Dispersal parameters  */
  /*------------------------*/
  const std::vector<double> disp_patho_tmp = Rcpp::as<std::vector<double>>(dispersal["disp_patho"]);
  const std::vector<double> disp_patho_sex_tmp = Rcpp::as<std::vector<double>>(dispersal["disp_patho_sex"]);
  const std::vector<double> disp_host_tmp = Rcpp::as<std::vector<double>>(dispersal["disp_host"]);
  Vector2D<double> disp_patho(Npoly, std::vector<double>(Npoly));
  Vector2D<double> disp_patho_sex(Npoly, std::vector<double>(Npoly));
  Vector2D<double> disp_host(Npoly, std::vector<double>(Npoly));
  
  // Matrix as Vector was created by columns
  for(int i = 0; i < Npoly; i++) {
    for(int j = 0; j < Npoly; j++) {
      disp_patho[j][i] = disp_patho_tmp[j + i * Npoly];
      disp_patho_sex[j][i] = disp_patho_sex_tmp[j + i * Npoly];
      disp_host[i][j] = disp_host_tmp[j + i * Npoly];
    }
  }
  
  /*----------------------*/
  /*  Initial conditions  */
  /*----------------------*/
  const double pI0 = Rcpp::as<double>(inits["pI0"]);
  
  /*-------------------*/
  /*  Host parameters  */
  /*-------------------*/
  //const int Nhost = Rcpp::as<int>(cultivars_param["Nhost"]);
  const std::vector<double> initial_density = Rcpp::as<std::vector<double>>(cultivars_param["initial_density"]);
  const std::vector<double> max_density = Rcpp::as<std::vector<double>>(cultivars_param["max_density"]);
  const std::vector<double> growth_rate = Rcpp::as<std::vector<double>>(cultivars_param["growth_rate"]);
  const std::vector<double> reproduction_rate = Rcpp::as<std::vector<double>>(cultivars_param["reproduction_rate"]);
  const std::vector<double> death_rate = Rcpp::as<std::vector<double>>(cultivars_param["death_rate"]);
  const int Nhost = initial_density.size();
  Rcpp::List cultivars_genes_list = Rcpp::as<Rcpp::List>(cultivars_param["cultivars_genes_list"]);
  std::vector<Cultivar> cultivars;
  std::vector<int> total_genes_id(0); // Contains all the genes_id used by the cultivars
  for(int i = 0; i < Nhost; i++) {
    const std::vector<int> genes_id = Rcpp::as<std::vector<int>>(cultivars_genes_list[i]);
    //Rcpp::Rcerr << "host " << i << std::endl;
    //for(int t=0; t<genes_id.size();t++) Rcpp::Rcerr << "\tGenes " << genes_id[t] << std::endl;
    cultivars.push_back(Cultivar(initial_density[i], max_density[i], growth_rate[i], reproduction_rate[i],
                                 death_rate[i], genes_id));
    total_genes_id.insert(total_genes_id.end(), genes_id.begin(), genes_id.end()); // Appends a vector to another
  }
  
  // Remove duplicate value in "total_genes_id"
  std::sort(total_genes_id.begin(), total_genes_id.end());
  auto last = std::unique(total_genes_id.begin(), total_genes_id.end());
  total_genes_id.erase(last, total_genes_id.end());
  
  double sigmoid_kappa_host = Rcpp::as<double>(cultivars_param["sigmoid_kappa_host"]);
  // security to avoid kappa = 0
  sigmoid_kappa_host += 1E-6 * (sigmoid_kappa_host == 0);
  const double sigmoid_sigma_host = Rcpp::as<double>(cultivars_param["sigmoid_sigma_host"]);
  const double sigmoid_plateau_host = Rcpp::as<double>(cultivars_param["sigmoid_plateau_host"]);
  
  /*---------------------*/
  /* Pathogen parameters */
  /*---------------------*/
  const double infection_rate = Rcpp::as<double>(basic_patho_param["infection_rate"]);
  const double propagule_prod_rate = Rcpp::as<double>(basic_patho_param["propagule_prod_rate"]);
  const double latent_period_mean = Rcpp::as<double>(basic_patho_param["latent_period_mean"]);
  double latent_period_var = Rcpp::as<double>(basic_patho_param["latent_period_var"]);
  // security to avoid variance = 0
  latent_period_var += 1E-6 * (latent_period_var == 0);
  const double infectious_period_mean = Rcpp::as<double>(basic_patho_param["infectious_period_mean"]);
  double infectious_period_var = Rcpp::as<double>(basic_patho_param["infectious_period_var"]);
  // security to avoid variance = 0
  infectious_period_var += 1E-6 * (infectious_period_var == 0);
  const double survival_prob = Rcpp::as<double>(basic_patho_param["survival_prob"]);
  const std::vector<double> repro_sex_prob = Rcpp::as< std::vector<double> >(basic_patho_param["repro_sex_prob"]);
  double sigmoid_kappa = Rcpp::as<double>(basic_patho_param["sigmoid_kappa"]);
  // security to avoid kappa = 0
  sigmoid_kappa += 1E-6 * (sigmoid_kappa == 0);
  const double sigmoid_sigma = Rcpp::as<double>(basic_patho_param["sigmoid_sigma"]);
  const double sigmoid_plateau = Rcpp::as<double>(basic_patho_param["sigmoid_plateau"]);
  const int sex_propagule_viability_limit = Rcpp::as<int>(basic_patho_param["sex_propagule_viability_limit"]);
  const double sex_propagule_release_mean = Rcpp::as<double>(basic_patho_param["sex_propagule_release_mean"]);
  const bool clonal_propagule_gradual_release = Rcpp::as<bool>(basic_patho_param["clonal_propagule_gradual_release"]);
  
  Basic_patho basic_patho(infection_rate, propagule_prod_rate, latent_period_mean, latent_period_var,
                          infectious_period_mean, infectious_period_var, survival_prob, repro_sex_prob, sigmoid_kappa,
                          sigmoid_sigma, sigmoid_plateau, sex_propagule_viability_limit, sex_propagule_release_mean, clonal_propagule_gradual_release);
#ifdef DEBUG
  Rcpp::Rcerr << basic_patho.to_string() << std::endl;
#endif
  
  Treatment treatment(Rcpp::as<double>(treatment_param["treatment_reduction_rate"]),
                      Rcpp::as<double>(treatment_param["treatment_efficiency"]),
                      Rcpp::as< std::vector<int> >(treatment_param["treatment_timesteps"]),
                      Rcpp::as< std::vector<int> >(treatment_param["treatment_cultivars"]),
                      Rcpp::as<double>(treatment_param["treatment_cost"])
                      );
#ifdef DEBUG
  Rcpp::Rcerr << treatment.to_string() << std::endl;
#endif
  
  /*------------------------*/
  /*  Evolution parameters  */
  /*------------------------*/
  const std::vector<double> time_to_activ_mean = Rcpp::as<std::vector<double>>(genes_param["time_to_activ_mean"]);
  const std::vector<double> time_to_activ_var = Rcpp::as<std::vector<double>>(genes_param["time_to_activ_var"]);
  const std::vector<int> Nlevels_aggressiveness = Rcpp::as<std::vector<int>>(genes_param["Nlevels_aggressiveness"]);
  const std::vector<std::string> target_trait = Rcpp::as<std::vector<std::string>>(genes_param["target_trait"]);
  const std::vector<double> mutation_prob = Rcpp::as<std::vector<double>>(genes_param["mutation_prob"]);
  const std::vector<double> efficiency = Rcpp::as<std::vector<double>>(genes_param["efficiency"]);
  const std::vector<double> fitness_cost = Rcpp::as<std::vector<double>>(genes_param["fitness_cost"]);
  const std::vector<double> tradeoff_strength = Rcpp::as<std::vector<double>>(genes_param["tradeoff_strength"]);
  const std::vector<double> recombination_sd = Rcpp::as<std::vector<double>>(genes_param["recombination_sd"]);

  
  // ! BE CAREFUL ! 
  // total_genes_id -> GENE ID in DB
  // genes_param -> liste sorted by GENE ID order ?
  std::vector<Gene> genes(0);
  //for( int g=0; g < total_genes_id.size(); g++) {
  for(int g : total_genes_id) { // Only keep the gene used by the cultivars (Reduce Npatho)
    genes.push_back(Gene(time_to_activ_mean[g], time_to_activ_var[g], Nlevels_aggressiveness[g], target_trait[g],
                         mutation_prob[g], efficiency[g], fitness_cost[g], tradeoff_strength[g], recombination_sd[g]));
  }
  const int Ngene = genes.size();
  
  int Npatho = 1;
  for(int g = 0; g < Ngene; g++) {
    Npatho *= genes[g].Nlevels_aggressiveness;
  }
  
  // Create the model
  Model model(Nyears, nTSpY, Npoly, Nhost, Npatho, Ngene, area, rotation, gen, cultivars, genes, basic_patho, treatment, 
              croptypes, sigmoid_kappa_host, sigmoid_sigma_host, sigmoid_plateau_host, pI0, disp_patho, disp_patho_sex, disp_host,
              seed);
  
  /*--------------------------------------*/
  /* Write and Print the model parameters */
  /*--------------------------------------*/
    model.print_param(seed, mutation_prob, efficiency, fitness_cost, tradeoff_strength);

  /* -------------- */
  /* Epidemic model */
  /* -------------- */
  Rprintf("\n*** SPATIOTEMPORAL MODEL SIMULATING THE SPREAD AND EVOLUTION OF A PATHOGEN IN A LANDSCAPE ***\n\n");
  model.dynepi();
}
