// ///////////////////////////////////////////////////////////////////////////// //
//
// TMB OPTIMIZATION MODEL FOR ITALY MORTALITY MODEL
// 
// Author: Nat Henry
// Created: 21 August 2020
// Purpose: TMB objective function for Italy all-cause mortality model with
//   space, time (year/week), and age random effects
// 
// /////////////////////////////////////////////////////////////////////////////


#include <TMB.hpp>
#include <math.h>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix; 


// HELPER FUNCTIONS ----------------------------------------------------------->

// Function for preparing a LCAR RE structure matrix given a neighborhood graph 
//  and rho parameter
// Adapted from Laura Dwyer-Lindgren
template<class Type> 
SparseMatrix<Type> lcar_strmat(SparseMatrix<Type> graph, Type rho) {
  SparseMatrix<Type> K = rho * graph; 
  for (int i = 0; i < K.rows(); i++)
    K.coeffRef(i,i) += (1 - rho);
  return K; 
}

// Robust Inverse Logit that sets min and max values to avoid numerical instability
template<class Type>
Type invlogit_robust(Type x){
  if (x < -20.723){
    x = -20.723; // corresponds to p=1e-9
  } else if ( x > 20.723 ){
    x = 20.723;  // corresponds to p=1-1e-9
  }
  return 1 / (1 + exp( -1.0 * x ));
}

// Transformation from (-Inf, Inf) to (-1, 1) for all rho parameters
template<class Type>
Type rho_transform(Type rho){
  return (exp(rho) - 1) / (exp(rho) + 1);
}


// OBJECTIVE FUNCTION --------------------------------------------------------->

template<class Type>
Type objective_function<Type>::operator() () {

  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // INPUT DATA --------------------------------------------------------------->

    DATA_INTEGER(flag);

    // OPTION: Holdout number
    // Any observation where `idx_holdout` is equal to `holdout` will be 
    //   excluded from this model fit
    // Holdouts are 1-indexed, so if holdout==0 all data will be used
    DATA_INTEGER(holdout);

    // Input mortality data
    DATA_VECTOR(y_i); // Number of deaths in a given location/year/week/age group
    DATA_VECTOR(n_i); // Population size for the corresponding group
    DATA_VECTOR(days_exp_i); // Days of exposure in a given week of observations
    
    // Fixed effects
    // Matrix of (n observations) * (m covariates including intercept)
    DATA_MATRIX(X_ij);

    // Indices
    DATA_IVECTOR(idx_loc);     // Index for the location of the observation
    DATA_IVECTOR(idx_year);    // Index for the time period (year)
    DATA_IVECTOR(idx_week);    // Index for the time period (week) 
    DATA_IVECTOR(idx_age);     // Index for the age group
    DATA_IVECTOR(idx_holdout); // Holdout index for each BH observation

    // Adjacency matrix for locations
    DATA_SPARSE_MATRIX(loc_adj_mat);


  // INPUT PARAMETERS --------------------------------------------------------->

    // Fixed effects
    PARAMETER_VECTOR(beta_covs); // Vector of fixed effects on covariates
    PARAMETER_VECTOR(beta_ages); // Vector of fixed effects on age group

    // Correlated random effect surface
    // 3-dimensional array of size: (# locations) x (# years) x (# ages)
    PARAMETER_ARRAY(Z_sta);
    // 3-dimensional array of size: (# years) x (# weeks) x (# ages)
    // PARAMETER_ARRAY(Z_twa);

    // Nugget
    // Vector of random effects, same length as number of observations
    PARAMETER_VECTOR(nugget);

    // Random effect autocorrelation parameters, transformed scale
    PARAMETER(rho_loc_trans_sta);  // By location
    PARAMETER(rho_year_trans_sta); // By year
    PARAMETER(rho_age_trans_sta);  // By age group
    // PARAMETER(rho_year_trans_twa); // By year
    // PARAMETER(rho_week_trans_twa); // By week
    // PARAMETER(rho_age_trans_twa);  // By age group

    // Variance of space-time-age-year random effect
    PARAMETER(log_sigma_sta);
    // PARAMETER(log_sigma_loc_sta);
    // PARAMETER(log_sigma_year_sta);
    // PARAMETER(log_sigma_age_sta);
    // PARAMETER(log_sigma_year_twa);
    // PARAMETER(log_sigma_week_twa);
    // PARAMETER(log_sigma_age_twa);
    PARAMETER(log_sigma_nugget);


  // TRANSFORM DATA AND PARAMETER OBJECTS ------------------------------------->

    // Basic indices
    int num_obs = y_i.size();
    int num_covs = beta_covs.size();

    // Transform some of our parameters
    // - Convert rho from (-Inf, Inf) to (-1, 1)
    Type rho_loc_sta = rho_transform(rho_loc_trans_sta);
    Type rho_year_sta = rho_transform(rho_year_trans_sta);
    Type rho_age_sta = rho_transform(rho_age_trans_sta);
    // Type rho_year_twa = rho_transform(rho_year_trans_twa);
    // Type rho_week_twa = rho_transform(rho_week_trans_twa);
    // Type rho_age_twa = rho_transform(rho_age_trans_twa);

    // Convert from log-sigma (-Inf, Inf) to sigmas (must be positive)
    Type sigma_sta = exp(log_sigma_sta);
    // Type sigma_loc_sta = exp(log_sigma_loc_sta);
    // Type sigma_year_sta = exp(log_sigma_year_sta);
    // Type sigma_age_sta = exp(log_sigma_age_sta);
    // Type sigma_year_twa = exp(log_sigma_year_twa);
    // Type sigma_week_twa = exp(log_sigma_week_twa);
    // Type sigma_age_twa = exp(log_sigma_age_twa);
    Type sigma_nugget = exp(log_sigma_nugget);

    // Create the LCAR covariance matrix 
    SparseMatrix<Type> loc_structure = lcar_strmat(loc_adj_mat, rho_loc_sta);

    // Vectors of fixed and structured random effects for all data points
    vector<Type> fes_i(num_obs);
    vector<Type> struct_res_i(num_obs);

    // Create a vector to hold individual data estimates
    vector<Type> logit_prob_i(num_obs);
    logit_prob_i.setZero();


  // Instantiate joint negative log-likelihood -------------------------------->

    Type jnll = 0.0;


  // JNLL CONTRIBUTION FROM PRIORS -------------------------------------------->

    // N(0, 3) prior for fixed effects
    for(int j = 0; j < num_covs; j++){
      PARALLEL_REGION jnll -= dnorm(beta_covs(j), Type(0.0), Type(3.0), true);
    }

    // N(0, 3) prior for sigmas
    PARALLEL_REGION jnll -= dnorm(sigma_sta, Type(0.0), Type(3.0), true);
    // PARALLEL_REGION jnll -= dnorm(sigma_loc_sta, Type(0.0), Type(3.0), true);
    // PARALLEL_REGION jnll -= dnorm(sigma_year_sta, Type(0.0), Type(3.0), true);
    // PARALLEL_REGION jnll -= dnorm(sigma_age_sta, Type(0.0), Type(3.0), true);
    // PARALLEL_REGION jnll -= dnorm(sigma_year_twa, Type(0.0), Type(3.0), true);
    // PARALLEL_REGION jnll -= dnorm(sigma_week_twa, Type(0.0), Type(3.0), true);
    // PARALLEL_REGION jnll -= dnorm(sigma_age_twa, Type(0.0), Type(3.0), true);
    PARALLEL_REGION jnll -= dnorm(sigma_nugget, Type(0.0), Type(3.0), true);

    // Evaluation of separable space-year-age random effect surface
    // PARALLEL_REGION jnll += SEPARABLE( \
    //     SCALE(AR1(rho_age_sta), sigma_age_sta), SEPARABLE( \ 
    //     SCALE(AR1(rho_year_sta), sigma_year_sta), \
    //     SCALE(GMRF(loc_structure, false), sigma_loc_sta) \
    // ))(Z_sta);
    PARALLEL_REGION jnll += SCALE(SEPARABLE(\
        AR1(rho_age_sta), SEPARABLE(AR1(rho_year_sta), GMRF(loc_structure, false))\
    ), sigma_sta)(Z_sta);

    // Evaluation of separable year-week-age random effect surface
    // PARALLEL_REGION jnll += SEPARABLE( \ 
    //     SCALE(AR1(rho_age_twa), sigma_age_twa), SEPARABLE( \
    //     SCALE(AR1(rho_week_twa), sigma_week_twa), \
    //     SCALE(AR1(rho_year_twa), rho_year_twa) \
    // ))(Z_twa);

    // Evaluation of nugget
    for(int i = 0; i < num_obs; i++){
      PARALLEL_REGION jnll -= dnorm(nugget[i], Type(0.0), sigma_nugget, true);
    }

    if(flag == 0) return jnll;


  // JNLL CONTRIBUTION FROM DATA ---------------------------------------------->

    // Determine fixed effect component for all observations
    fes_i = X_ij * beta_covs.matrix();

    for(int i=0; i < num_obs; i++){
      if(idx_holdout[i] != holdout){
        // Determine structured random effect component for this observation
        // struct_res_i[i] = Z_sta[idx_loc[i], idx_year[i], idx_age[i]] + \
        //     Z_twa[idx_year[i], idx_week[i], idx_age[i]];
        struct_res_i[i] = Z_sta[idx_loc[i], idx_year[i], idx_age[i]];
        
        // Determine logit probability for this observation
        logit_prob_i[i] = fes_i[i] + struct_res_i[i] + nugget[i];

        PARALLEL_REGION jnll -= dbinom_robust(\
            y_i[i]*7/days_exp_i[i], n_i[i], logit_prob_i[i], true);
      }
    }    


  // RETURN JNLL -------------------------------------------------------------->

    return jnll;

} // END objective function
