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
// DNRE_Q0_rate_interaction_map_reduce.stan
// Rate model with P fixed effect covariates, random intercepts
// and interactions for a categorical variable coded from 1 to C
// using map_reduce for within-chain parallelization.

functions {
  real partial_sum_lpmf(
    array[] int start_rate,
    int start,
    int end,
    array[] int end_rate,
    array[] int chose_rate,
    vector timespan,
    array[] int is_dependent,
    real log_crude_rate,
    matrix X_rate,
    array[] int sender,
    array[] int interaction,
    array[] vector beta_rate,
    vector gamma
 ) {
    real log_lik = 0.0;
    int t_index;
    int size_slice;
    int start_event;
    int end_event;
    int interaction_event;
    int chose_event;

    for (t in 1:(end - start + 1)) {
      t_index = t + start - 1;
      start_event = start_rate[t];
      end_event = end_rate[t_index];
      interaction_event = interaction[t_index];
      chose_event = chose_rate[t_index] - start_event + 1;
      size_slice = end_event - start_event + 1;
      array[size_slice] int event_slice =
        linspaced_int_array(size_slice, start_event, end_event);
      
      vector[size_slice] xb_rate =
        X_rate[event_slice] * beta_rate[interaction_event] +
        gamma[sender[event_slice]] + log_crude_rate;
      
      if (timespan[t_index] > 0)
        log_lik += (is_dependent[t_index] ? xb_rate[chose_event] : 0) -
          timespan[t_index] * exp(log_sum_exp(xb_rate));
    }

    return log_lik;
  }
}

data {
  int<lower=1> N_rate;       // Total number of choices in all events
  int<lower=1> T_rate;       // Number of events
  int<lower=0> P_rate;       // Number of fixed-effect covariates

  matrix[N_rate, P_rate] X_rate;

  // the starting, ending and chosen index observation for each event
  array[T_rate] int<lower=1, upper=N_rate> start_rate;
  array[T_rate] int<lower=1, upper=N_rate> end_rate;
  array[T_rate] int<lower=0, upper=N_rate> chose_rate;

  // random effects vars
  int<lower=1> A;              // Number of actors/groups
  array[N_rate] int<lower=1, upper=A> sender;
  // array[A] int<lower=1, upper=T_choice> start_group;

  // rate information: time span and right-censored events
  vector<lower=0>[T_rate] timespan;
  array[T_rate] int<lower=0, upper=1> is_dependent;

  // interaction var for event
  int<lower=2> C; // number of categories
  array[T_rate] int<lower=1, upper=C> interaction;

  int<lower=1> grain_size;     // Grain size for map_reduce
}

transformed data {
  // array[A] int<lower=1, upper=T_choice> end_group;
  // for (g in 1:(A - 1)) {
  //   end_group[g] = start_group[g + 1] - 1;
  // }
  // end_group[A] = T_choice;

  real log_crude_rate = log(T_rate / (N_rate * mean(timespan)));
}

parameters {
  array[C] vector[P_rate] beta_rate; // Fixed effects
  real<lower=0> sigma;          // Variance of the random effect
  vector[A] gamma_raw;          // Uncentered random effects
}

transformed parameters {
  vector[A] gamma = sigma * gamma_raw; // Centered random effects
}

model {
  // Priors
  for (c in 1:C) {
    // target += normal_lpdf(beta_rate[c, 1] | 0, 4);
    target += std_normal_lpdf(beta_rate[c]);
  }
  target += exponential_lpdf(sigma | 1);
  target += std_normal_lpdf(gamma_raw);

  // Likelihood
  target += reduce_sum(partial_sum_lpmf,
                       start_rate,
                       grain_size,
                       end_rate,
                       chose_rate,
                       timespan,
                       is_dependent,
                       log_crude_rate,
                       X_rate,
                       sender,
                       interaction,
                       beta_rate,
                       gamma);
}
