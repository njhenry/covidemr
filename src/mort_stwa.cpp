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
#include <stdio.h>
using namespace density;
using Eigen::SparseMatrix;


// HELPER FUNCTIONS ----------------------------------------------------------->

// Function for preparing a LCAR precision matrix based on an adjacency matrix
//  and a distribution
// Taken from MacNab 2011, adapted from the `ar.matrix` package:
// https://rdrr.io/cran/ar.matrix/src/R/Q.lCAR.R
template<class Type>
SparseMatrix<Type> lcar_q_from_graph(SparseMatrix<Type> graph, Type sigma, Type rho){
  SparseMatrix<Type> D(graph.rows(), graph.cols());
  SparseMatrix<Type> I(graph.rows(), graph.cols());
  for(int ii=0; ii < D.rows(); ii ++){
    D.insert(ii, ii) = graph.col(ii).sum();
    I.insert(ii, ii) = 1.0;
  }
  SparseMatrix<Type> Q = 1/sigma * (rho * (D - graph) + (1 - rho) * I);
  return Q;
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
    DATA_IVECTOR(idx_fourier); // Index for the Fourier transform group
    DATA_IVECTOR(idx_holdout); // Holdout index for each BH observation

    // Adjacency matrix for locations
    DATA_SPARSE_MATRIX(loc_adj_mat);

    // Which of the correlated random effects structures should be used?
    DATA_INTEGER(use_Z_stwa);
    DATA_INTEGER(use_Z_sta);
    DATA_INTEGER(use_Z_fourier);
    DATA_INTEGER(use_nugget);
    DATA_INTEGER(harmonics_level);


  // INPUT PARAMETERS --------------------------------------------------------->

    // Fixed effects
    PARAMETER_VECTOR(beta_covs); // Vector of fixed effects on covariates
    PARAMETER_VECTOR(beta_ages); // Vector of fixed effects on age group

    // Random effect autocorrelation parameters, transformed scale
    PARAMETER(rho_loc_trans);  // By location
    PARAMETER(rho_year_trans); // By year
    PARAMETER(rho_week_trans); // By week
    PARAMETER(rho_age_trans);  // By age group

    // Variance of space-time-age-year random effect
    PARAMETER(log_sigma_loc);
    PARAMETER(log_sigma_year);
    PARAMETER(log_sigma_week);
    PARAMETER(log_sigma_age);
    PARAMETER(log_sigma_nugget);

    // Correlated random effect surfaces
    // Multiple possible surfaces can be added together
    // -> 4-dimensional array of size: (# locations) by (# years) by (# weeks) by (# ages)
    PARAMETER_ARRAY(Z_stwa);
    // -> 3-dimensional array of size: (# locations) by (# years) by (# ages)
    PARAMETER_ARRAY(Z_sta);

    // Harmonics matrix
    // Dimensions: (# separately-fit harmonics) by (2 x harmonics series level)
    PARAMETER_ARRAY(Z_fourier);

    // Nugget
    // Vector of random effects, same length as number of observations
    PARAMETER_VECTOR(nugget);


  // DATA CHECKS -------------------------------------------------------------->

    // Warn if both Z_stwa and Z_sta are being used
    if(use_Z_sta == 1 && use_Z_stwa == 1){
      printf("Warning: both STA and STWA random effects are in use.");
    }
    if(use_Z_fourier == 1 && Z_fourier.cols() / harmonics_level != 2){
      printf("Warning: incorrect number of columns for Z harmonics.");
    }

  // TRANSFORM DATA AND PARAMETER OBJECTS ------------------------------------->

    // Basic indices
    int num_obs = y_i.size();
    int num_covs = beta_covs.size();

    // Transform some of our parameters
    // - Convert rho from (-Inf, Inf) to (-1, 1)
    Type rho_loc = rho_transform(rho_loc_trans);
    Type rho_year = rho_transform(rho_year_trans);
    Type rho_week = rho_transform(rho_week_trans);
    Type rho_age = rho_transform(rho_age_trans);

    // Convert from log-sigma (-Inf, Inf) to sigmas (must be positive)
    Type sigma_loc = exp(log_sigma_loc);
    Type sigma_year = exp(log_sigma_year);
    Type sigma_week = exp(log_sigma_week);
    Type sigma_age = exp(log_sigma_age);
    Type sigma_nugget = exp(log_sigma_nugget);

    // Create the LCAR covariance matrix
    SparseMatrix<Type> loc_Q = lcar_q_from_graph(\
        loc_adj_mat, sigma_loc, rho_loc\
    );

    // Vectors of fixed and structured random effects for all data points
    vector<Type> fes_i(num_obs);
    vector<Type> struct_res_i(num_obs);

    // Create a vector to hold individual data estimates
    vector<Type> logit_prob_i(num_obs);

    // Harmonic frequency of 1 year (52 weeks)
    Type year_freq = 2.0 * 3.1415926535 / 52.0;


  // Instantiate joint negative log-likelihood -------------------------------->

    Type jnll = 0.0;


  // JNLL CONTRIBUTION FROM PRIORS -------------------------------------------->

    // N(0, 3) prior for fixed effects
    // Skip the intercept (index 0)
    for(int j = 1; j < num_covs; j++){
      PARALLEL_REGION jnll -= dnorm(beta_covs(j), Type(0.0), Type(3.0), true);
    }

    // N(0, 3) prior for age effects
    // Skip the first age group (index 0)
    for(int j = 1; j < beta_ages.size(); j++){
      PARALLEL_REGION jnll -= dnorm(beta_ages(j), Type(0.0), Type(3.0), true);
    }

    // N(0, 3) prior for sigmas
    // TODO: Try messing with variance
    PARALLEL_REGION jnll -= dnorm(sigma_loc, Type(0.0), Type(3.0), true);
    PARALLEL_REGION jnll -= dnorm(sigma_year, Type(0.0), Type(3.0), true);
    PARALLEL_REGION jnll -= dnorm(sigma_week, Type(0.0), Type(3.0), true);
    PARALLEL_REGION jnll -= dnorm(sigma_age, Type(0.0), Type(3.0), true);
    PARALLEL_REGION jnll -= dnorm(sigma_nugget, Type(0.0), Type(0.1), true);

    if(use_Z_fourier == 1){
      // N(0, 1) prior for harmonics
      for(int i = 0; i < Z_fourier.rows(); i++){
        for(int j = 0; j < Z_fourier.cols(); j++){
          PARALLEL_REGION jnll -= dnorm(Z_fourier(i, j), Type(0.0), Type(1.0), true);
        }
      }
    }

    if(use_Z_stwa == 1){
      // Evaluation of separable space-year-week-age random effect surface
      // Rescale AR1 in age and time; spatial RE has already been scaled
      PARALLEL_REGION jnll += SEPARABLE(\
          SCALE(AR1(rho_age), sigma_age),\
          SEPARABLE(\
              SCALE(AR1(rho_week), sigma_week),\
              SEPARABLE(\
                  SCALE(AR1(rho_year), sigma_year),\
                  GMRF(loc_Q, false)\
              )\
          )\
      )(Z_stwa);
    }

    if(use_Z_sta == 1){
      // Evaluation of separable space-year-age random effect surface
      // Rescale AR1 in age and time; spatial RE has already been scaled
      PARALLEL_REGION jnll += SEPARABLE(\
          SCALE(AR1(rho_age), sigma_age),\
          SEPARABLE(\
              SCALE(AR1(rho_year), sigma_year),\
              GMRF(loc_Q, false)\
          )\
      )(Z_sta);
    }

    // Evaluation of nugget
    if(use_nugget == 1){
      for(int i = 0; i < num_obs; i++){
        PARALLEL_REGION jnll -= dnorm(nugget(i), Type(0.0), sigma_nugget, true);
      }
    }

    if(flag == 0) return jnll;


  // JNLL CONTRIBUTION FROM DATA ---------------------------------------------->

    // Determine fixed effect component for all observations
    fes_i = X_ij * beta_covs.matrix();

    for(int i=0; i < num_obs; i++){
      if(idx_holdout(i) != holdout){
        // Determine structured random effect component for this observation
        struct_res_i(i) = 0;
        if(use_Z_stwa == 1){
          struct_res_i(i) += Z_stwa(idx_loc(i), idx_year(i), idx_week(i), idx_age(i));
        }
        if(use_Z_sta == 1){
          struct_res_i(i) += Z_sta(idx_loc(i), idx_year(i), idx_age(i));
        }
        if(use_Z_fourier == 1){
          for(int lev=1; lev <= harmonics_level; lev++){
            struct_res_i(i) += Z_fourier(idx_fourier(i), 2*lev-2) * sin(lev * (idx_week(i) + 1.0) * year_freq ) + \
              Z_fourier(idx_fourier(i), 2*lev-1) * cos(lev * (idx_week(i) + 1.0) * year_freq);
          }
        }
        if(use_nugget == 1){
          struct_res_i(i) += nugget(i);
        }

        // Determine logit probability for this observation
        logit_prob_i(i) = fes_i(i) + beta_ages(idx_age(i)) + struct_res_i(i);

        // Use the dbinom_robust PDF function, parameterized using logit(prob)
        PARALLEL_REGION jnll -= dbinom_robust(\
            y_i(i), n_i(i) * days_exp_i(i) / 7.0, logit_prob_i(i), true\
        );
      }
    }


  // RETURN JNLL -------------------------------------------------------------->

    return jnll;

} // END objective function
