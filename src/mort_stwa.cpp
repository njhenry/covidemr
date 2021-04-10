// ///////////////////////////////////////////////////////////////////////////////////////
//
// TMB OPTIMIZATION MODEL FOR ITALY MORTALITY MODEL
//
// Author: Nat Henry
// Created: 21 August 2020
// Purpose: TMB objective function for Italy all-cause mortality model with
//   space, time (year/week), and age random effects
//
// MODEL FORMULATION
//
// D_i \sim Poisson(N_i * p_i)
//
// logit(p_i) = \sum_{k=1}^{4}\alpha_i * \mathbb{I}[age_i = k] + \vec{\beta} X +
//   Z_{prov_i, year_i, age_i} + f_{prov_i, age_i}(week_i) + \epsilon_i
//
// ///////////////////////////////////////////////////////////////////////////////////////


#include <TMB.hpp>
using Eigen::SparseMatrix;


// HELPER FUNCTIONS --------------------------------------------------------------------->

// Transformation from (-Inf, Inf) to (-1, 1) for all rho parameters
template<class Type>
Type rho_transform(Type rho){
  return (exp(rho) - 1) / (exp(rho) + 1);
}

// Function for preparing a CAR precision matrix (Q) based on an adjacency matrix
//  and a neighborhood autocorrelation parameter.
//
//  param W: Sparse adjacency matrix W = {w_ij}, where w_ij = w_ji = 1 if areal units i
//    and j are neighbors and 0 otherwise.
//  param rho: Neighborhood autocorrelation parameter in [-1, 1]
//
//  For more information, see:
//  Banerjee, Carlin, and Gelfand (2015). Hierarchical Modeling and Analysis for Spatial
//    Data, 2nd Edition. Section 4.3: Conditionally autoregressive models
//
template<class Type>
SparseMatrix<Type> car_precision_from_graph(SparseMatrix<Type> W, Type rho){
  SparseMatrix<Type> D(W.rows(), W.cols());
  for(int ii=0; ii < D.rows(); ii++){
    D.insert(ii, ii) = W.row(ii).sum();
  }
  SparseMatrix<Type> Q = (D - rho * W);
  return Q;
}


// OBJECTIVE FUNCTION ------------------------------------------------------------------->

template<class Type>
Type objective_function<Type>::operator() () {

  // INPUT DATA ------------------------------------------------------------------------->

    // OPTION: Holdout number
    // Any observation where `idx_holdout` is equal to `holdout` will be
    //   excluded from this model fit
    // Holdouts are 1-indexed, so if holdout==0 all data will be used
    DATA_INTEGER(holdout);

    // Input mortality data
    DATA_VECTOR(y_i); // Number of deaths in a given location/year/week/age group
    DATA_VECTOR(n_i); // Population size for the corresponding group
    DATA_VECTOR(days_exp_i); // Days of exposure in a given week of observations

    // Covariate matrix: dimensions (n observations) by (m covariates including intercept)
    DATA_MATRIX(X_ij);

    // Indices
    DATA_IVECTOR(idx_loc);     // Index for the location of the observation
    DATA_IVECTOR(idx_year);    // Index for the time period (year)
    DATA_IVECTOR(idx_week);    // Index for the time period (week)
    DATA_IVECTOR(idx_age);     // Index for the age group
    DATA_IVECTOR(idx_fourier); // Index for the Fourier transform group
    DATA_IVECTOR(idx_holdout); // Holdout index for each observation

    // Adjacency matrix capturing province neighborhood structure
    DATA_SPARSE_MATRIX(adjacency_matrix);

    // Which of the correlated random effects structures should be used?
    DATA_INTEGER(use_Z_sta);
    DATA_INTEGER(use_Z_fourier);
    DATA_INTEGER(use_nugget);

    // Number of harmonic terms used to fit seasonality
    DATA_INTEGER(harmonics_level);

    // Is random effects automatic process normalization being used? If so, should only
    //  priors be returned?
    DATA_INTEGER(auto_normalize);
    DATA_INTEGER(early_return);


  // INPUT PARAMETERS ------------------------------------------------------------------->

    // Fixed effects
    PARAMETER_VECTOR(beta_covs); // Vector of fixed effects on covariates
    PARAMETER_VECTOR(beta_ages); // Vector of fixed effects on age group

    // Random effect autocorrelation parameters, transformed scale
    PARAMETER(rho_year_trans); // By year
    PARAMETER(rho_age_trans);  // By age group
    PARAMETER(rho_loc_trans);  // Neighborhood correlation by province

    // Log precision of space-time-age-year random effect
    PARAMETER(log_tau_sta);
    // Log precision of the nugget
    PARAMETER(log_tau_nugget);

    // Correlated random effect surfaces
    // -> 3-dimensional array of size: (# locations) by (# years) by (# ages)
    PARAMETER_ARRAY(Z_sta);

    // Harmonics matrix
    // Dimensions: (# separately-fit harmonics) by (2 x harmonics series level)
    PARAMETER_ARRAY(Z_fourier);

    // Nugget
    // Vector of random effects, same length as number of observations
    PARAMETER_VECTOR(nugget);


  // TRANSFORM DATA AND PARAMETER OBJECTS ----------------------------------------------->

    // Basic indices
    int num_obs = y_i.size();
    int num_covs = beta_covs.size();
    int num_fourier_groups = Z_fourier.rows();

    // Transform some of our parameters
    // - Convert rho from (-Inf, Inf) to (-1, 1)
    Type rho_year = rho_transform(rho_year_trans);
    Type rho_age = rho_transform(rho_age_trans);
    Type rho_loc = rho_transform(rho_loc_trans);

    // Convert from log-tau (-Inf, Inf) to tau (must be positive)
    Type tau_sta = exp(log_tau_sta);
    Type sd_sta = exp(log_tau_sta * Type(-0.5));
    Type tau_nugget = exp(log_tau_nugget);
    Type sd_nugget = exp(log_tau_nugget * Type(-0.5));

    // Construct CAR precision matrix
    SparseMatrix<Type> Q_car = car_precision_from_graph(adjacency_matrix, rho_loc);

    // Vectors of fixed and structured random effects for all data points
    vector<Type> fix_effs(num_obs);
    vector<Type> ran_effs(num_obs);
    ran_effs.setZero();

    // Create a vector to hold data-specific estimates of the mortality rate
    //   per person-week
    vector<Type> weekly_mort_rate_i(num_obs);

    // Harmonic frequency of 1 year (52 weeks)
    Type year_freq = 2.0 * 3.1415926535 / 52.0;


  // Instantiate joint negative log-likelihood (JNLL) ----------------------------------->

    parallel_accumulator<Type> jnll(this);


  // JNLL CONTRIBUTION FROM PRIORS ------------------------------------------------------>

    // N(mean=0, sd=3) prior for fixed effects
    // Skip the intercept (index 0)
    for(int j = 1; j < num_covs; j++){
      jnll -= dnorm(beta_covs(j), Type(0.0), Type(3.0), true);
    }

    // N(mean=0, sd=3) prior for age effects
    // Skip the first age group (index 0)
    for(int j = 1; j < beta_ages.size(); j++){
      jnll -= dnorm(beta_ages(j), Type(0.0), Type(3.0), true);
    }

    if(use_Z_sta){
      // Wide gamma priors for tau precision parameters
      jnll -= dlgamma(tau_sta, Type(1.0), Type(1000.0), true);
      // Evaluate separable prior against the space-time-age random effects:
      // Spatial effect = CAR model using province neighborhood structure
      // Time effect = AR1 by year
      // Age effect = AR1 by age group
      jnll += density::SCALE(
        density::SEPARABLE(
          density::AR1(rho_age),
          density::SEPARABLE(
            density::AR1(rho_year),
            density::GMRF(Q_car, !auto_normalize)
          )
        ), sd_sta
      )(Z_sta);
    }

    if(use_nugget){
      // Wide gamma priors for tau precision hyperparameters
      jnll -= dlgamma(tau_nugget, Type(1.0), Type(1000.0), true);
      // Evaluate prior on each nugget
      jnll -= dnorm(nugget, Type(0.0), sd_nugget, true).sum();
    }

    if(use_Z_fourier){
      // N(mean=0, sd=3) prior for harmonic terms
      for(int i = 0; i < num_fourier_groups; i++){
        for(int j = 0; j < Z_fourier.cols(); j++){
          jnll -= dnorm(Z_fourier(i,j), Type(0.0), Type(3.0), true);
        }
      }
    }

    // Option to return the JNLL early for automatic process normalization
    if(auto_normalize && early_return) return jnll;


  // JNLL CONTRIBUTION FROM DATA -------------------------------------------------------->

    // Determine fixed effect component for all observations
    fix_effs = X_ij * beta_covs.matrix();

    for(int i=0; i < num_obs; i++){
      if(idx_holdout(i) != holdout){
        // Random effects and seasonality terms
        if(use_Z_sta){
          ran_effs(i) += Z_sta(idx_loc(i), idx_year(i), idx_age(i));
        }
        if(use_nugget){
          ran_effs(i) += nugget(i);
        }
        if(use_Z_fourier){
          for(int lev=1; lev <= harmonics_level; lev++){
            ran_effs(i) += (
              Z_fourier(idx_fourier(i), 2*lev-2) * sin(lev * (idx_week(i) + 1.0) * year_freq) +
              Z_fourier(idx_fourier(i), 2*lev-1) * cos(lev * (idx_week(i) + 1.0) * year_freq)
            );
          }
        }
        // Estimate weekly mortality rate based on log-linear mixed effects model
        weekly_mort_rate_i(i) = exp(beta_ages(idx_age(i)) + fix_effs(i) + ran_effs(i));
        // Use the dpois PDF function centered around:
        //  lambda = population * weekly mort rate * (observed days / 7)
        jnll -= dpois(y_i(i), n_i(i) * weekly_mort_rate_i(i) * days_exp_i(i) / 7.0, true);
      }
    }


  // RETURN JNLL ------------------------------------------------------------------------>

    return jnll;

} // END objective function
