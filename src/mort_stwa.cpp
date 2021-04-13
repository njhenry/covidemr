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
using namespace density;
using Eigen::SparseMatrix;

// HELPER FUNCTIONS --------------------------------------------------------------------->

// Transformation from (-Inf, Inf) to (-1, 1) for all rho parameters
template<class Type>
Type rho_transform(Type rho){
  return (exp(rho) - 1) / (exp(rho) + 1);
}

// Function for preparing a precision matrix corresponding to a Leroux CAR spatial model,
//   based on a spatial adjacency matrix. Adapted from the `ar.matrix` package:
//   https://rdrr.io/cran/ar.matrix/src/R/Q.lCAR.R
//
// Parameter W: A sparse adjacency matrix indicating spatial neighborhood structure.
//   Defined as W = {w_ij}, where w_ij = w_ji = 1 if areal units i and j are neighbors,
//   and 0 otherwise.
// Parameter phi: A mixing parameter indicating the relative contribution of spatial and
//   IID variation, strictly between 0 and 1.
//
template<class Type>
SparseMatrix<Type> lcar_precision_from_adjacency(SparseMatrix<Type> W, Type sigma, Type phi){
  SparseMatrix<Type> D(W.rows(), W.cols());
  SparseMatrix<Type> I(W.rows(), W.cols());
  for(int ii=0; ii < W.rows(); ii ++){
    D.insert(ii, ii) = W.row(ii).sum();
    I.insert(ii, ii) = 1.0;
  }
  SparseMatrix<Type> Q = 1.0 / sigma * (phi * (D - W) + (1 - phi) * I);
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

    // Adjacency matrix capturing spatial neighborhood structure
    DATA_SPARSE_MATRIX(adjacency_matrix);
    // Rank deficiency of the ICAR graph (1 if the adjacency graph is fully connected)
    DATA_SCALAR(icar_rank_deficiency);

    // Which of the correlated random effects structures should be used?
    DATA_INTEGER(use_Z_sta);
    DATA_INTEGER(use_Z_fourier);
    DATA_INTEGER(use_nugget);

    // Number of harmonic terms used to fit seasonality
    DATA_INTEGER(harmonics_level);

    DATA_INTEGER(auto_normalize);
    DATA_INTEGER(early_return);


  // INPUT PARAMETERS ------------------------------------------------------------------->

    // Fixed effects
    PARAMETER_VECTOR(beta_covs); // Vector of fixed effects on covariates
    PARAMETER_VECTOR(beta_ages); // Vector of fixed effects on age group

    // Random effect autocorrelation parameters, transformed scale
    PARAMETER(rho_year_trans); // By year
    PARAMETER(rho_age_trans);  // By age group

    // Log precision of space-time-age-year random effect
    PARAMETER(log_tau_loc);
    // Log precision of the nugget
    PARAMETER(log_tau_nugget);

    // Mixing parameter controlling spatial vs. nonspatial correlation by province
    PARAMETER(logit_phi_loc);

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
    int num_locs = Z_sta.dim(0);
    int num_years = Z_sta.dim(1);
    int num_ages = Z_sta.dim(2);
    int num_fourier_groups = Z_fourier.rows();

    // Transform some of our parameters
    // - Convert rho from (-Inf, Inf) to (-1, 1)
    Type rho_year = rho_transform(rho_year_trans);
    Type rho_age = rho_transform(rho_age_trans);

    // Convert from log-tau (-Inf, Inf) to tau (must be positive)
    Type tau_loc = exp(log_tau_loc);
    Type sigma_loc = exp(log_tau_loc * Type(-0.5));
    Type tau_nugget = exp(log_tau_nugget);
    Type sigma_nugget = exp(log_tau_nugget * Type(-0.5));

    // Convert mixing parameter to the space (0, 1)
    Type phi_loc = exp(logit_phi_loc)/(Type(1.0) + exp(logit_phi_loc));

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

    Type jnll = 0.0;


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
      jnll -= dlgamma(tau_loc, Type(1.0), Type(10.0), true);
      // Evaluate separable prior against the space-time-age random effects:
      // Spatial effect = CAR model using province neighborhood structure
      // Time effect = AR1 by year
      // Age effect = AR1 by age group
      SparseMatrix<Type> Q_loc = lcar_precision_from_adjacency(adjacency_matrix, sigma_loc, phi_loc);
      jnll += SEPARABLE(AR1(rho_age), SEPARABLE(AR1(rho_year), GMRF(Q_loc, bool(1-auto_normalize))))(Z_sta);
      // SEPARABLE is calculating the density of Q_sta if Q_space was full rank. We need
      //   to subtract the difference in density caused by the rank deficiency of the
      //   ICAR precision matrix.
      // jnll -= 0.5 * icar_rank_deficiency * (
      //   (num_years - 1) * log(1 - rho_year * rho_year) - log(2 * PI) +
      //   (num_ages - 1) * log(1 - rho_age * rho_age) - log(2 * PI)
      // );
    }

    if(use_nugget){
      // Wide gamma priors for tau precision hyperparameters
      jnll -= dlgamma(tau_nugget, Type(1.0), Type(1000.0), true);
      // Evaluate prior on each nugget
      jnll -= dnorm(nugget, Type(0.0), sigma_nugget, true).sum();
    }

    if(use_Z_fourier){
      // N(mean=0, sd=3) prior for harmonic terms
      for(int i = 0; i < num_fourier_groups; i++){
        for(int j = 0; j < Z_fourier.cols(); j++){
          jnll -= dnorm(Z_fourier(i,j), Type(0.0), Type(3.0), true);
        }
      }
    }

    if(early_return) return jnll;


  // JNLL CONTRIBUTION FROM DATA -------------------------------------------------------->

    // // Soft sum-to-zero constraint on spatial REs for identifiability
    // jnll -= dnorm(Z_sta.sum(), Type(0.0), Type(0.001) * Z_sta.size(), true);

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
