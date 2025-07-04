// MCM.stan
// Multinomial choice model with P covariates
// Update to avoid warnings from deprecated syntax
// trying changing log_sum_exp in transformed quantities for speed up
//   from mc-stan discourse: speeding up a hierarchical negative binomial
//     regression with a huge number of groups
// Copyright (C) 2024, SNlab- ETH Zurich
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// version 2, as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
data {
  int N_choice; // number of events * choice set (actors - 1)
  int T_choice; // number of events
  int P_choice; // number of covariates

  matrix[N_choice, P_choice] X_choice;
  // position receiver chose
  array[T_choice] int<lower = 1, upper = N_choice> chose_choice; 

  // start and end position for each event
  array[T_choice] int<lower = 1, upper = N_choice> start_choice;
  array[T_choice] int<lower = 1, upper = N_choice> end_choice;
}
parameters {
  vector[P_choice] beta_choice;
}
model {
  // priors
  target += normal_lpdf(beta_choice | 0, 4);

  // helper for likelihood
  vector[N_choice] xb_choice;
  xb_choice = X_choice * beta_choice;

  for(t in 1:T_choice)
    target += xb_choice[chose_choice[t]] -
      log_sum_exp(xb_choice[start_choice[t]:end_choice[t]]);
}
