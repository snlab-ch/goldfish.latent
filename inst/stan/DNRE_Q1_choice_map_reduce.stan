// Copyright (C) 2025, Alvaro Uzaheta - SNlab-ETH Zurich
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT
// License for more details.
//
// You should have received a copy of the MIT License along with this
// program. If not, see <https://opensource.org/licenses/MIT>.
//
// DNRE_Q1_choice_map_reduce.stan
// Multinomial choice model with P fixed effect covariates and one random effect
// using map_reduce for within-chain parallelization.

functions {
  real partial_sum_lpmf(
    array[] int start_group,
    int start,
    int end,
    array[] int end_group,
    matrix X_choice,
    vector Z_choice,
    array[] int chose_choice,
    array[] int start_choice,
    array[] int end_choice,
    vector beta_choice,
    vector gamma
 ) {
    real log_lik = 0.0;
    int g_index;

    for (g in 1:(end - start + 1)) {
      g_index = g + start - 1;
      
      for (t in start_group[g]:end_group[g_index]) {
        int size_slice = end_choice[t] - start_choice[t] + 1;
        array[size_slice] int event_slice =
          linspaced_int_array(size_slice, start_choice[t], end_choice[t]);
        
        vector[size_slice] xb_choice = X_choice[event_slice] * beta_choice +
          Z_choice[event_slice] .* gamma[g_index];
        
        log_lik += xb_choice[chose_choice[t]] - log_sum_exp(xb_choice);
      }
    }
    return log_lik;
  }
}

data {
  int<lower=1> N_choice;       // Total number of choices in all events
  int<lower=1> T_choice;       // Number of events
  int<lower=0> P_choice;       // Number of fixed-effect covariates

  matrix[N_choice, P_choice] X_choice;

  // the starting, ending and chosen index observation for each event
  array[T_choice] int<lower=1, upper=N_choice> start_choice;
  array[T_choice] int<lower=1, upper=N_choice> end_choice;
  array[T_choice] int<lower=1, upper=N_choice> chose_choice;

  // random effects vars
  int<lower=1> A;              // Number of actors/groups
  int<lower=1> Q_choice;       // Number of random effects (must be 1)
  array[Q_choice] int<lower=1> V1_choice;
  array[A] int<lower=1, upper=T_choice> start_group;

  int<lower=1> grain_size;     // Grain size for map_reduce
}

transformed data {
  array[A] int<lower=1, upper=T_choice> end_group;
  for (g in 1:(A - 1)) {
    end_group[g] = start_group[g + 1] - 1;
  }
  end_group[A] = T_choice;

  vector[N_choice] Z_choice = to_vector(X_choice[, V1_choice]);
}

parameters {
  vector[P_choice] beta_choice; // Fixed effects
  real<lower=0> sigma;          // Variance of the random effect
  vector[A] gamma_raw;          // Uncentered random effects
}

transformed parameters {
  vector[A] gamma = sigma * gamma_raw; // Centered random effects
}

model {
  // Priors
  target += std_normal_lpdf(beta_choice);
  target += exponential_lpdf(sigma | 1);
  target += std_normal_lpdf(gamma_raw);

  // Likelihood
  target += reduce_sum(partial_sum_lpmf,
                       start_group,
                       grain_size,
                       end_group,
                       X_choice,
                       Z_choice,
                       chose_choice,
                       start_choice,
                       end_choice,
                       beta_choice,
                       gamma);
}
