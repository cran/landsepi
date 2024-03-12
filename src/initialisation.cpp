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

/****************************************************************/
/*                    initialisation.c                          */
/****************************************************************/
#include "Model.hpp"

/* Switch from indices of aggressiveness to index pathotype */
int Model::switch_aggr_to_patho(const std::vector<int>& aggr) {
    int index_patho = 0;

    for(int i = 0; i < this->Ngene; i++) {
        int prod = 1;
        for(int j = i + 1; j < this->Ngene; j++) {
            prod *= this->genes[j].Nlevels_aggressiveness;
        }
        index_patho += aggr[i] * prod;
    }
    return index_patho;
}

/* Switch from index pathotype to indices of aggressiveness (vector of dimension Ngene) */
std::vector<int> Model::switch_patho_to_aggr(const int& index_patho) {
    std::vector<int> aggr(this->Ngene);

    int remainder = index_patho;
    for(int i = 0; i < this->Ngene; i++) {
        int prod = 1;
        for(int j = i + 1; j < this->Ngene; j++) {
            prod *= this->genes[j].Nlevels_aggressiveness;
        }
        aggr[i] = remainder / prod; /* Quotient */
        remainder = remainder % prod;
    }
    return aggr;
}


/* Switch from indices of aggressiveness (in the patho) to value of the traits for resistance (in the host)
 * activeR is a vector of the same dimension of aggressiveness which define, for a given cultivar, which resistance genes are active at a
 * given timestep t */
  
std::vector<double> Model::switch_aggr_to_trait(const std::vector<int>& aggr, const std::vector<bool>& activeR) {
    std::vector<double> trait(aggr.size());
    double basic_value = 0.0;
    for(int g = 0; g < this->Ngene; g++) {
        if(this->genes[g].target_trait == "IR") {
            basic_value = this->basic_patho.infection_rate;
            trait[g] = basic_value * this->genes[g].aggressiveness_matrix[aggr[g]][activeR[g]];  
        } else if (this->genes[g].target_trait == "LAT"){
            basic_value = this->basic_patho.latent_period_mean;
            trait[g] = basic_value / (this->genes[g].aggressiveness_matrix[aggr[g]][activeR[g]]+
              0.001 * (this->genes[g].aggressiveness_matrix[aggr[g]][activeR[g]] == 0));  
        } else if (this->genes[g].target_trait == "IP"){
            basic_value = this->basic_patho.infectious_period_mean;
            trait[g] = basic_value * this->genes[g].aggressiveness_matrix[aggr[g]][activeR[g]];  
        } else if (this->genes[g].target_trait == "PR"){
            basic_value = this->basic_patho.propagule_prod_rate;
            trait[g] = basic_value * this->genes[g].aggressiveness_matrix[aggr[g]][activeR[g]];  
        }
     
        }
    return trait;
}

/* Switch from value of the traits for resistance (in the host) to indices of aggressiveness (in the patho)
 * activeR is a vector of th same dimension on trait which define, for a given cultivar, which resistance genes are active at a
 * given timestep t */

std::vector<int> Model::switch_trait_to_aggr(const std::vector<double>& trait, const std::vector<bool>& activeR) {
    std::vector<int> aggr(trait.size());
    double aggr_coeff = 0.0;
    for(int g = 0; g < this->Ngene; g++) {
        if(this->genes[g].target_trait == "IR") {
            aggr_coeff = trait[g]/this->basic_patho.infection_rate;
        } else if (this->genes[g].target_trait == "LAT"){
          aggr_coeff = this->basic_patho.latent_period_mean/trait[g];
        } else if (this->genes[g].target_trait == "IP"){
            aggr_coeff = trait[g]/this->basic_patho.infectious_period_mean;
        } else if (this->genes[g].target_trait == "PR"){
            aggr_coeff = trait[g]/this->basic_patho.propagule_prod_rate;
        }
        for(int i = 0; i < genes[g].Nlevels_aggressiveness - 1; i++){
            if (abs(abs(aggr_coeff) - this->genes[g].aggressiveness_matrix[i][activeR[g]]) < abs(abs(aggr_coeff) - this->genes[g].aggressiveness_matrix[i+1][activeR[g]])){
                aggr[g] = i;
                break;
            } else {
                aggr[g] = i+1;
            }
            }
                
        }
    return aggr;
}

/* Initialisation of activeR: matrix indicating for each poly and each gene, at which time step resistance starts to be
 * active */
Vector2D<int> Model::init_activeR() {
    Vector2D<int> activeR(this->Npoly, std::vector<int>(this->Ngene));
    // note: if Ngene==0, activeR is of size 0 and will never be used
    for(int poly = 0; poly < this->Npoly; poly++) {
        for(int gene = 0; gene < this->Ngene; gene++) {
            const double exp = this->genes[gene].age_of_activ_mean;
            const double var = this->genes[gene].age_of_activ_var;
            if(exp > 0 && var > 0) {
                const std::array<double, 2> alpha = find_paramGamma(exp, var);
                activeR[poly][gene] = static_cast<int>(ran_gamma(alpha[0], alpha[1]));
            } else if (var == 0){
                activeR[poly][gene] = (int) exp;  /* Deterministic delay if variance==0 */
            } else {
                activeR[poly][gene] = 0;
            }
        }
    }
    return activeR;
}

/* Indicates if a cultivar carries a given gene and this one is active at a given timestep t */
bool Model::get_resistance(const int& index_gene, const int& host, const int& t, const int& activeR) {
    std::vector<int> genes_id = this->cultivars[host].genes_id;
    if(std::find(genes_id.begin(), genes_id.end(), index_gene) != genes_id.end()) {
        if(t >= activeR) {
            return 1;
        }
    }
    return 0;
}

/* Indicates the efficacy of the pesticide treatment at time t         */
/* Following the model of Arneson et al. (2002), the temporal dynamics */
/*  of fungicide residue after application is influenced by fungicide  */ 
/*  characteristics and plant growth.                                  */

double Model::get_treat_effect(const int& Nt, const int& Nspray, const int& t, const int& t_lastspray) {
  
  /* Evaluating fungicide CONCENTRATION at time t after fungicide application Ct */
  /* C1 is the degradation of fungicide concentration due to time */
  /* C2 is the reduction of fungicide concentration due to plant growth, */
  /* since new plant tissue is not covered by fungicide */
  /* Fungicide concentration at the time of application (C0) is considered equal to 1 */
 
  int t_after_spray = 0; // time after fungicide application
  
  double treat_effect = 1; 
  
  /* Evaluating the variable "time after spray" */
  if(t_lastspray == 0){
    treat_effect = 1; // No treatment effect before the first spray application 
  }else{
    t_after_spray = t - t_lastspray;
    double C1 = 1*exp(- this->treatment.treatment_degradation_rate*t_after_spray);
    double C2 = 1;
    if (Nt > 0){
      C2 = std::min(1.0,static_cast<double>(Nspray)/static_cast<double>(Nt));
    }
    double C = C1*C2;
    
    
    /* Evaluating the fungicide EFFICIENCY, which is predicted by a logistic function */
    /* of the fungicide concentration C. logistic_a and logistic_b are shape parameters, the range */
    /* of their values is between parentheses */ 
    
    double logistic_a = 4.0; // 3.5; [3.5 - 4.5]
    double logistic_b = 8.5; // 9.0; [8.0 - 9.0]
    treat_effect = 1 - this->treatment.treatment_efficiency/(1+exp(logistic_a-logistic_b*C));
  }
  return treat_effect;
}

/* Extract the number of primary sexual propagules (P_stock_poly) that will */
/* be released the following season for a given poly and for each pathogen */

std::vector<int> Model::get_P_stock_release(Vector2D<int>& P_stock_poly, const int& year){
  // year_idx is an index of the position of the current year into the P_stock 3D matrix
  int year_idx = (year-1) % this->basic_patho.sex_propagule_viability_limit; 
  std::vector<int> P_stock_release(P_stock_poly.size(),0);
  for(int patho = 0; patho<this->Npatho; patho++){
      P_stock_release[patho] = P_stock_poly[patho][year_idx];
      P_stock_poly[patho][year_idx] = 0;
    }
  return P_stock_release;
}
/* Extract the number of propagules among the pool of primary inoculum (P_primary) that will be released at time-step t */
/*   for each poly and for each pathogen */
void Model::get_P_daily(Vector2D<int>& P_daily, Vector3D<int>& P_primary, const int& t) {   
  for(unsigned int r = 0; r < P_primary.size(); r++){
    for(unsigned int c = 0; c < P_primary[0].size(); c++){
      P_daily[r][c] = P_primary[r][c][t];
    }
  }
}

/* Sum of two Vector2D */
Vector2D<int> Model::get_sum_Vector2D(Vector2D<int>& M1, Vector2D<int>& M2) {
  Vector2D<int> M_sum(M1.size(), std::vector<int>(M1[0].size()));
  for(unsigned int r = 0; r < M1.size(); r++){
    for(unsigned int c = 0; c < M1[0].size(); c++){
      M_sum[r][c] = M1[r][c] + M2[r][c];
    }
  }
  return M_sum;
}

/* Initialisation of H, L, I, and R at 0 */
void Model::init_HjuvLIR(Vector2D<int>& Hjuv, Vector3D<int>& L, Vector3D<int>& I, Vector3D<int>& R) {
    Hjuv = Vector2D<int>(this->Npoly, std::vector<int>(this->Nhost, 0));
    L = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost, 0)));
    I = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost, 0)));
    R = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost, 0)));
}

/* Initialisation of P, P_sex_secondary, P_clonal_secondary, P_sex_primary, P_sex_primary_tmp, 
  P_clonal_primary, P_clonal_primary_tmp, P_sex_daily, P_clonal_daily, P_stock  at 0 */
void Model::init_P( Vector2D<int>& P, Vector2D<int>& P_sex_secondary, Vector2D<int>& P_clonal_secondary, 
                    Vector3D<int>& P_sex_primary, Vector2D<int>& P_sex_primary_tmp, Vector3D<int>& P_clonal_primary, 
                    Vector2D<int>& P_clonal_primary_tmp, Vector2D<int>& P_sex_daily, Vector2D<int>& P_clonal_daily,
                    Vector3D<int>& P_stock) {
  P = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
  P_sex_secondary = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
  P_clonal_secondary = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
  P_sex_primary = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->time_steps_per_year, 0)));
  P_sex_primary_tmp = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
  P_clonal_primary = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->time_steps_per_year, 0)));
  P_clonal_primary_tmp = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
  P_sex_daily = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
  P_clonal_daily = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
  P_stock = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->basic_patho.sex_propagule_viability_limit, 0)));
}

/* Initialisation of Nspray and t_lastspray at 0 */
void Model::init_Nspray_t_lastspray(Vector2D<int>& Nspray, Vector2D<int>& t_lastspray) {
  Nspray = Vector2D<int>(this->Npoly, std::vector<int>(this->Nhost, 0));
  t_lastspray = Vector2D<int>(this->Npoly, std::vector<int>(this->Nhost, 0));
}

/* Initialisation of TFI at 0 */
void Model::init_TFI(Vector3D<int>& TFI) {
  TFI = Vector3D<int>(this->Nyears, Vector2D<int>(this->Npoly, std::vector<int>(this->Nhost, 0)));
}

/* Initialisation of Nlevels_aggressiveness at 0 */
void Model::init_Nlevels_aggressiveness(std::vector<int>& Nlevels_aggressiveness) {
    Nlevels_aggressiveness = std::vector<int>(this->Ngene, 0);
}

/* Initialise L2I and I2R with 0 */
void Model::init_L2I2R(Vector4D<int>& L2I, Vector4D<int>& I2R) {
    L2I = Vector4D<int>(
        this->Npoly,
        Vector3D<int>(this->Npatho, Vector2D<int>(this->Nhost, std::vector<int>(this->time_steps_per_year, 0))));
    I2R = Vector4D<int>(
        this->Npoly,
        Vector3D<int>(this->Npatho, Vector2D<int>(this->Nhost, std::vector<int>(this->time_steps_per_year, 0))));
}

/* Plantation of H in the beginning of a season */
Vector2D<int> Model::intro_H(const int& year) {
    Vector2D<int> H(this->Npoly, std::vector<int>(this->Nhost));
    for(int poly = 0; poly < this->Npoly; poly++) {
        // If there is no rotation (same croptype each year)
        int id_croptype = (this->rotation[poly].size() == 1) ? this->rotation[poly][0] : this->rotation[poly][year];

        for(std::pair<int, double> cultivar_prop : this->croptypes[id_croptype].cultivar_proportion) {
            int id_host = cultivar_prop.first;
            double prop = cultivar_prop.second;
            H[poly][id_host] = static_cast<int>(this->area[poly] * this->cultivars[id_host].initial_density * prop);
        }
    }
    return H;
}

/* Pathogen introduction : infectious sites from pathotype 0 in cultivar 0 */
void Model::intro_I(Vector2D<int>& H, Vector3D<int>& I, Vector4D<int>& I2R, const Vector2D<int>& activeR) {
    const int t = 0; // Introduction at t=0
    
    for(unsigned int poly=0; poly<this->pI0.size(); poly++)  {
      for(unsigned int patho=0; patho<this->pI0[poly].size(); patho++)    {
        for(unsigned int host=0; host<this->pI0[poly][patho].size(); host++)      {
          // Infection of hosts (no consideration of resistance genes)
          I[poly][patho][host] = ran_binomial(this->pI0[poly][patho][host], H[poly][host]);
          H[poly][host] -= I[poly][patho][host];
          
          
          
          // Calculation of infectious period for infected sites at t=0 
          const std::vector<int> aggr = switch_patho_to_aggr(patho);
          for(int i = 0; i < I[poly][patho][host]; i++) {
            double infectious_period_mean_new = this->basic_patho.infectious_period_mean;
            for(int g = 0; g < this->Ngene; g++) {
              if(this->genes[g].target_trait == "IP") {
                // Indicate if the cultivar has an active resistance gene
                bool Rgene = get_resistance(g, host, t, activeR[poly][g]);
                infectious_period_mean_new *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
              }
            }
            // Security to avoid problem in alpha calculation
            infectious_period_mean_new += 0.001 * (infectious_period_mean_new == 0);
            std::array<double, 2> infectious_period_alpha =
              find_paramGamma(infectious_period_mean_new, this->basic_patho.infectious_period_var);
            int lag = static_cast<int>(ran_gamma(infectious_period_alpha[0], infectious_period_alpha[1]));
            lag += 1 * (lag == 0); // Security to avoid infectious duration of 0 day
            if(lag < this->time_steps_per_year) {
              I2R[poly][patho][host][lag]++;
            }
          }
        
        }
      }
    }
    /*
     With a multinomial draw instead 
     for(int poly=0;poly<pI0_mat.size();poly++)  {
        for(int patho=0;patho<pI0_mat[poly].size();patho++)    {
           pI0_mat[poly][patho]=tmp;
           pos++;
           I[poly][patho] = ran_multinomial(pI0_mat[poly][patho], H[poly]);
           H[poly] -= I[poly][patho];
        }
     }
     */
}
