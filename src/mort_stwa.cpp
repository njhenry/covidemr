// ///////////////////////////////////////////////////////////////////////////// //
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
// /////////////////////////////////////////////////////////////////////////////


#include <TMB.hpp>
using std::isfinite;
using std::printf;

// HELPER FUNCTIONS ----------------------------------------------------------->

// Transformation from (-Inf, Inf) to (-1, 1) for all rho parameters
template<class Type>
Type rho_transform(Type rho){
  return (exp(rho) - 1) / (exp(rho) + 1);
}

template<class Type>
Type check_finite(Type x){
  if(isfinite(x)){
    return x;
  }
  printf("NAN EVAL ");
  return 0;
}

// OBJECTIVE FUNCTION --------------------------------------------------------->

template<class Type>
Type objective_function<Type>::operator() () {

  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // INPUT DATA --------------------------------------------------------------->

    // OPTION: Holdout number
    // Any observation where `idx_holdout` is equal to `holdout` will be
    //   excluded from this model fit
    // Holdouts are 1-indexed, so if holdout==0 all data will be used
    DATA_INTEGER(holdout);

    // Input mortality data
    DATA_VECTOR(y_i); // Number of deaths in a given location/year/week/age group
    DATA_VECTOR(n_i); // Population size for the corresponding group
    DATA_VECTOR(days_exp_i); // Days of exposure in a given week of observations

    // Covariate matrix: dimensions (n observations) * (m covariates including intercept)
    DATA_MATRIX(X_ij);

    // Indices
    DATA_IVECTOR(idx_loc);     // Index for the location of the observation
    DATA_IVECTOR(idx_year);    // Index for the time period (year)
    DATA_IVECTOR(idx_week);    // Index for the time period (week)
    DATA_IVECTOR(idx_age);     // Index for the age group
    DATA_IVECTOR(idx_fourier); // Index for the Fourier transform group
    DATA_IVECTOR(idx_holdout); // Holdout index for each observation

    // Precision matrix for ICAR priors
    DATA_SPARSE_MATRIX(Q_icar);
    // Rank deficiency of precision matrix
    DATA_SCALAR(Q_rank_deficiency);

    // Which of the correlated random effects structures should be used?
    DATA_INTEGER(use_Z_sta);
    DATA_INTEGER(use_Z_fourier);
    DATA_INTEGER(use_nugget);

    DATA_INTEGER(harmonics_level);


  // INPUT PARAMETERS --------------------------------------------------------->

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

    // Correlated random effect surfaces
    // -> 3-dimensional array of size: (# locations) by (# years) by (# ages)
    PARAMETER_ARRAY(Z_sta);

    // Harmonics matrix
    // Dimensions: (# separately-fit harmonics) by (2 x harmonics series level)
    PARAMETER_ARRAY(Z_fourier);

    // Nugget
    // Vector of random effects, same length as number of observations
    PARAMETER_VECTOR(nugget);


  // TRANSFORM DATA AND PARAMETER OBJECTS ------------------------------------->

    // Basic indices
    int num_obs = y_i.size();
    int num_covs = beta_covs.size();
    int num_locs = Z_sta.dim(0);
    int num_years = Z_sta.dim(1);
    int num_ages = Z_sta.dim(2);

    // Transform some of our parameters
    // - Convert rho from (-Inf, Inf) to (-1, 1)
    Type rho_year = rho_transform(rho_year_trans);
    Type rho_age = rho_transform(rho_age_trans);

    // Convert from log-tau (-Inf, Inf) to tau (must be positive)
    Type tau_sta = exp(log_tau_sta);
    Type sd_sta = exp(log_tau_sta * Type(-0.5));
    Type tau_nugget = exp(log_tau_nugget);
    Type sd_nugget = exp(log_tau_nugget * Type(-0.5));

    // Vectors of fixed and structured random effects for all data points
    vector<Type> fix_effs(num_obs);
    vector<Type> ran_effs(num_obs);
    ran_effs.setZero();

    // Create a vector to hold data-specific estimates of the mortality rate
    //   per person-week
    vector<Type> weekly_mort_rate_i(num_obs);

    // Harmonic frequency of 1 year (52 weeks)
    Type year_freq = 2.0 * 3.1415926535 / 52.0;


  // Instantiate joint negative log-likelihood -------------------------------->

    Type jnll = 0.0;


  // JNLL CONTRIBUTION FROM PRIORS -------------------------------------------->

    // N(mean=0, sd=3) prior for fixed effects
    // Skip the intercept (index 0)
    for(int j = 1; j < num_covs; j++){
      jnll -= check_finite(dnorm(beta_covs(j), Type(0.0), Type(3.0), true));
    }

    // N(mean=0, sd=3) prior for age effects
    // Skip the first age group (index 0)
    for(int j = 1; j < beta_ages.size(); j++){
      jnll -= check_finite(dnorm(beta_ages(j), Type(0.0), Type(3.0), true));
    }

    // Wide gamma priors for tau precision parameters
    if(use_nugget == 1){
      jnll -= check_finite(dlgamma(tau_nugget, Type(1.0), Type(10.0), true));

      // Soft sum-to-zero constraint on nugget effects for each year, location, and age group
      for(int age_i = 0; age_i < num_ages; age_i++){
        for(int year_i = 0; year_i < num_years; year_i++){
          for(int loc_i = 0; loc_i < num_locs; loc_i++){
            Type sum_nuggets = 0.0;
            for(int i=0; i < num_obs; i++){
              if((idx_loc(i)==loc_i) && (idx_age(i)==age_i) && (idx_year(i)==year_i)){
                sum_nuggets += nugget(i);
              }
            }
            jnll -= check_finite(dnorm(sum_nuggets, Type(0.0), Type(0.1), true));
          }
        }
      }
    }

    if(use_Z_sta == 1){
      // Wide gamma priors for tau precision parameters
      jnll -= check_finite(dlgamma(tau_sta, Type(1.0), Type(10.0), true));

      // Evaluate separable prior against the space-time-age random effects:
      // Spatial effect = ICAR by province
      // Time effect = AR1 by year
      // Age effect = AR1 by age group
      jnll += check_finite(SCALE(
        SEPARABLE(AR1(rho_age), SEPARABLE(AR1(rho_year), GMRF(Q_icar))),
        sd_sta
      )(Z_sta));

      // Soft sum-to-zero constraint on each layer of spatial random effects
      for(int age_i = 0; age_i < num_ages; age_i++){
        for(int year_i = 0; year_i < num_years; year_i++){
          Type sum_res = 0.0;
          for(int loc_i = 0; loc_i < num_locs; loc_i++){
            sum_res += Z_sta(loc_i, year_i, age_i);
          }
          jnll -= check_finite(dnorm(sum_res, Type(0.0), Type(0.001) * num_locs, true));
        }
      }

      // Adjust normalizing constant to account for rank deficiency of the ICAR precision
      //  matrix:
      // 1) Calculate log(generalized variance) of outer product matrix
      Type kronecker_log_genvar = (
        log(1 - rho_age * rho_age) * num_locs * num_years +
        log(1 - rho_year * rho_year) * num_locs * num_ages
      );
      // 2) Adjust normalizing constant
      jnll -= check_finite(Q_rank_deficiency * 0.5 * (kronecker_log_genvar - log(2 * PI)));
    }

    if(use_Z_fourier == 1){
      for(int i = 0; i < Z_fourier.rows(); i++){
        for(int j = 0; j < Z_fourier.cols(); j++){
          // N(mean=0, sd=3) prior for harmonic terms
          jnll -= check_finite(dnorm(Z_fourier(i,j), Type(0.0), Type(3.0), true));
        }
      }
    }

    // Evaluation of prior on nugget
    if(use_nugget == 1){
      jnll -= check_finite(dnorm(nugget, Type(0.0), sd_nugget, true).sum());
    }


  // JNLL CONTRIBUTION FROM DATA ---------------------------------------------->

    // Determine fixed effect component for all observations
    fix_effs = X_ij * beta_covs.matrix();

    for(int i=0; i < num_obs; i++){
      if(idx_holdout(i) != holdout){

        // Random effects and seasonality terms
        if(use_Z_sta == 1){
          ran_effs(i) += Z_sta(idx_loc(i), idx_year(i), idx_age(i));
        }
        if(use_nugget == 1){
          ran_effs(i) += nugget(i);
        }
        if(use_Z_fourier == 1){
          for(int lev=1; lev <= harmonics_level; lev++){
            ran_effs(i) += (
              Z_fourier(idx_fourier(i), 2*lev-2) * sin(lev * (idx_week(i) + 1.0) * year_freq) +
              Z_fourier(idx_fourier(i), 2*lev-1) * cos(lev * (idx_week(i) + 1.0) * year_freq)
            );
          }
        }

        // Estimate weekly mortality rate based on log-linear mixed effects model
        weekly_mort_rate_i(i) = exp(beta_ages(idx_age(i)) + fix_effs(i) + ran_effs(i));
        if(weekly_mort_rate_i(i) < 0.000001 ){
          weekly_mort_rate_i(i) = 0.000001;
        }

        // Use the dpois PDF function centered around:
        //  lambda = population * weekly mort rate * (observed days / 7)
        jnll -= check_finite(dpois(
            y_i(i), n_i(i) * weekly_mort_rate_i(i) * days_exp_i(i) / 7.0, true
        ));
      }
    }


  // RETURN JNLL -------------------------------------------------------------->

    return jnll;

} // END objective function
