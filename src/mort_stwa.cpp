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

// Function for rescaling a precision matrix to have standard deviation sigma
//
// Parameter Q: Unscaled precision matrix
// Parameter sigma: Standard deviation to scale to
//
template<class Type>
SparseMatrix<Type> scale_precision(SparseMatrix<Type> Q, Type sigma){
  SparseMatrix<Type> Q_scaled = Q / (sigma * sigma);
  return Q_scaled;
}

// Function to create an IID precision matrix (AKA a scaled identity matrix)
//
// Parameter dim: Number of rows (and columns) for the precision matrix
// Parameter sigma: Standard deviation of the iid process
//
template<class Type>
SparseMatrix<Type> iid_precision(int dim, Type sigma = 1.0){
  SparseMatrix<Type> I(dim, dim);
  for(int ii=0; ii < dim; ii++){
    I.insert(ii, ii) = 1.0;
  }
  SparseMatrix<Type> I_scaled = scale_precision(I, sigma);
  return I_scaled;
}

// Function to create a precision matrix corresponding with an autoregressive process of
//   order 1 (AR1 process). Adapted from the `ar.matrix` package:
//   https://rdrr.io/cran/ar.matrix/src/R/Q.AR1.R
//
// Corresponds to the process x_t = rho * x_(t-1) + epsilon
//
// Parameter steps: Number of steps (in time, age, etc) for the autoregressive process
// Paramter rho: Correlation between step t and step t+1
// Parameter sigma: Standard deviation of the random noise term epsilon
//
template<class Type>
SparseMatrix<Type> ar1_precision(int steps, Type rho, Type sigma = 1.0){
  SparseMatrix<Type> Q(steps, steps);
  for(int ii=1; ii < steps; ii++){
    Q.insert(ii, ii) = 1.0 + rho * rho;
    Q.insert(ii - 1, ii) = rho * -1.0;
    Q.insert(ii, ii - 1) = rho * -1.0;
  }
  Q.insert(0, 0) = 1.0;
  Q.insert(steps, steps) = 1.0;
  SparseMatrix<Type> Q_scaled = scale_precision(Q, sigma);
  return Q_scaled;
}

// Function for preparing a precision matrix corresponding to a BYM2 spatial model,
//   based on a scaled ICAR precision matrix. For more details, see:
//   Riebler et al. (2016). An intuitive Bayesian sptial model for disease mapping that
//     accounts for scaling. Statistical methods in medical research, 25(4):1145-65.
//
// Parameter Q_icar: A precision matrix corresponding to an intrinsic correlated
//   autoregressive (ICAR) model in space, scaled to have generalized variance 1
// Parameter phi: A mixing parameter indicating the relative contribution of spatial and
//   IID variation, strictly between 0 and 1.
// Parameter sigma: Standard deviation of the LCAR process
//
template<class Type>
SparseMatrix<Type> bym2_precision(SparseMatrix<Type> Q_icar, Type phi, Type sigma = 1.0){
  SparseMatrix<Type> I = iid_precision(Q_icar.rows(), Type(1.0));
  SparseMatrix<Type> Q_bym2 = phi * Q_icar + (1 - phi) * I;
  SparseMatrix<Type> Q_bym2_scaled = scale_precision(Q_bym2, sigma);
  return Q_bym2_scaled;
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

    // Precision matrix for an ICAR spatial model
    DATA_SPARSE_MATRIX(Q_icar);
    // Rank deficiency of the ICAR graph (1 if the adjacency graph is fully connected)
    DATA_SCALAR(icar_rank_deficiency);

    // Which of the correlated random effects structures should be used?
    DATA_INTEGER(use_Z_sta);
    DATA_INTEGER(use_Z_fourier);
    DATA_INTEGER(use_nugget);

    // Information about the seasonality harmonics fit
    // Are the parameters for the harmonics stationary?
    DATA_INTEGER(fourier_stationary);
    // Number of groupings for harmonics fit
    DATA_INTEGER(num_fourier_groups);
    // Number of harmonic terms used to fit seasonality
    DATA_INTEGER(harmonics_level);


  // INPUT PARAMETERS ------------------------------------------------------------------->

    // Fixed effects
    PARAMETER_VECTOR(beta_covs); // Vector of fixed effects on covariates
    PARAMETER_VECTOR(beta_ages); // Vector of fixed effects on age group

    // Random effect autocorrelation parameters, transformed scale
    PARAMETER(rho_year_trans); // By year
    PARAMETER(rho_age_trans);  // By age group

    // Log precision of space-time-age-year random effect
    PARAMETER(log_tau_sta);
    // Log precision of the nugget
    PARAMETER(log_tau_nugget);

    // Mixing parameter controlling spatial vs. nonspatial correlation by province
    PARAMETER(logit_phi_loc);

    // Hyperparameters used for nonstationary Fourier fits - used across all AR1 groups
    // phi (autocorrelation hyperparameter) for AR1 process
    PARAMETER(logit_phi_fourier);
    // log sigma (logged precision of innovation term in AR1 process)
    PARAMETER(log_tau_fourier);

    // Correlated random effect surfaces
    // -> 3-dimensional array of size: (# locations) by (# years) by (# ages)
    PARAMETER_ARRAY(Z_sta);

    // Vector of harmonics coefficients
    // IF stationary: Length = (2x harmonics series level) repeated (# unique fits) times
    // Eg: <coef1 group 1, coef2 group 1, coef1 group 2, coef2 group2, ...>
    //
    // IF nonstationary: Length = (2x harmonics series level), repeated (# unique fits)
    //   times, repeated (# years) times
    // Eg: <coef1 group1 yr1, coef2 group1 yr1, coef1 group2 yr1, coef2 group2 yr1,
    //      coef1 group1 yr2, coef2 group1 yr2, coef1 group2 yr2, coef2 group2 yr2, ...>
    PARAMETER_VECTOR(Z_fourier);

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
    int num_f_params = 2 * num_fourier_groups * harmonics_level;

    // Transform some of our parameters
    // - Convert rho from (-Inf, Inf) to (-1, 1)
    Type rho_year = rho_transform(rho_year_trans);
    Type rho_age = rho_transform(rho_age_trans);

    // Convert from log-tau (-Inf, Inf) to tau (must be positive)
    Type tau_sta = exp(log_tau_sta);
    Type sigma_sta = exp(log_tau_sta * Type(-0.5));
    Type tau_nugget = exp(log_tau_nugget);
    Type sigma_nugget = exp(log_tau_nugget * Type(-0.5));

    // Convert mixing parameter to the space (0, 1)
    Type phi_loc = invlogit(logit_phi_loc);

    // Transform harmonics hyperparameters (used when fourier terms are nonstationary)
    Type phi_fourier = invlogit(logit_phi_fourier);
    Type tau_fourier = exp(log_tau_fourier);
    Type sigma_fourier = exp(log_tau_fourier * Type(-0.5));

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

    if(use_Z_sta){
      // Gamma(1, 1E3) priors for tau precision parameters
      jnll -= dlgamma(tau_sta, Type(1.0), Type(1000.0), true);
      // Spatial effect = BYM2 (scaled CAR) model using province neighborhood structure
      SparseMatrix<Type> Q_loc = bym2_precision(Q_icar, phi_loc);
      // Time effect = AR1 by year
      SparseMatrix<Type> Q_year = ar1_precision(num_years, rho_year);
      // Age effect = AR1 by age group
      SparseMatrix<Type> Q_age = ar1_precision(num_ages, rho_age);

      // Evaluate separable prior against the space-time-age random effects
      jnll += SEPARABLE(GMRF(Q_age), SEPARABLE(GMRF(Q_year), GMRF(Q_loc)))(Z_sta);
      // SEPARABLE is calculating the density of Q_sta if Q_space was full rank. We need
      //   to subtract the difference in density caused by the rank deficiency of the
      //   ICAR precision matrix.
      jnll -= 0.5 * icar_rank_deficiency * (
        (num_years - 1) * log(1 - rho_year * rho_year) - log(2 * PI) +
        (num_ages - 1) * log(1 - rho_age * rho_age) - log(2 * PI)
      );
      // Sum-to-zero constraint on each layer of spatial REs for identifiability
      vector<Type> sum_res(num_ages * num_years);
      sum_res.setZero();
      for(int age_i = 0; age_i < num_ages; age_i++){
        for(int year_i = 0; year_i < num_years; year_i++){
          for(int loc_i = 0; loc_i < num_locs; loc_i++){
            sum_res(age_i * year_i + age_i) += Z_sta(loc_i, year_i, age_i);
          }
        }
      }
      jnll -= dnorm(sum_res, Type(0.0), Type(0.001) * num_locs, true).sum();
    }

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

    if(use_nugget){
      // Gamma(1, 1E3) priors for tau precision hyperparameters
      jnll -= dlgamma(tau_nugget, Type(1.0), Type(1000.0), true);
      // Evaluate prior on each nugget
      jnll -= dnorm(nugget, Type(0.0), sigma_nugget, true).sum();
    }

    if(use_Z_fourier){
      if(fourier_stationary){
        // CASE: Fourier terms are fixed across years
        // N(mean=0, sd=3) prior for harmonic terms
        jnll -= dnorm(Z_fourier, Type(0.0), Type(3.0), true).sum();
      } else {
        // CASE: Fourier terms are nonstationary (autoregressive) across years
        // Gamma(1, 1E3) priors for tau precision hyperparameters
        jnll -= dlgamma(tau_fourier, Type(1.0), Type(1000.0), true);
        // Evaluate all time series of Fourier terms against the same AR1 precision matrix
        vector<Type> coefficients_by_year(num_years);
        vector<Type> fourier_intercepts(num_f_params);
        coefficients_by_year.setZero();
        SparseMatrix<Type> Q_f_ar1 = ar1_precision(num_years, phi_fourier, sigma_fourier);
        for(int hp_i = 0; hp_i < num_f_params; hp_i++){
          // Populate coefficients across all years
          for(int yr_i = 0; yr_i < num_years; yr_i++){
            coefficients_by_year(yr_i) = Z_fourier(yr_i * num_f_params + hp_i);
          }
          // Evaluate density against the mean-zero time series
          fourier_intercepts(hp_i) = coefficients_by_year.mean();
          jnll += GMRF(Q_f_ar1)(coefficients_by_year - fourier_intercepts(hp_i));
        }
        // N(mean=0, sd=3) for the mean across all years
        jnll -= dnorm(fourier_intercepts, Type(0.0), Type(3.0), true).sum();
      }
    }


  // JNLL CONTRIBUTION FROM DATA -------------------------------------------------------->

    // Determine fixed effect component for all observations
    fix_effs = X_ij * beta_covs.matrix();

    // Convenience variable for indexing Fourier harmonic coefficients
    int f_term_i = 0;

    for(int i=0; i < num_obs; i++){
      if(idx_holdout(i) != holdout){
        // Random effects and seasonality terms
        if(use_Z_sta){
          ran_effs(i) += Z_sta(idx_loc(i), idx_year(i), idx_age(i)) * sigma_sta;
        }
        if(use_nugget){
          ran_effs(i) += nugget(i);
        }
        if(use_Z_fourier){
          for(int lev=1; lev <= harmonics_level; lev++){
            // Stationary case: Fourier terms indexed across groups > harmonic levels
            f_term_i = 2 * harmonics_level * idx_fourier(i) + 2 * lev;
            if(!fourier_stationary){
              // Nonstationary case: Also indexed across years
              f_term_i += num_f_params * idx_year(i);
            }
            ran_effs(i) += (
              Z_fourier(f_term_i) * sin(lev * (idx_week(i) + 1.0) * year_freq) +
              Z_fourier(f_term_i + 1) * cos(lev * (idx_week(i) + 1.0) * year_freq)
            );
          }
        }
        // Estimate weekly mortality rate based on log-linear mixed effects model
        weekly_mort_rate_i(i) = exp(beta_ages(idx_age(i)) + fix_effs(i) + ran_effs(i));
        // Use the dpois PDF function centered around:
        //  lambda = population * weekly mort rate * (observed days / 7)
        if(weekly_mort_rate_i(i) > 0.0){
          jnll -= dpois(y_i(i), n_i(i) * weekly_mort_rate_i(i) * days_exp_i(i) / 7.0, true);
        }
      }
    }


  // RETURN JNLL ------------------------------------------------------------------------>

    return jnll;

} // END objective function
