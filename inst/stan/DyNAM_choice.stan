// MCM.stan
// Multinomial choice model with P covariates
// Update to avoid warnings from deprecated syntax
// trying changing log_sum_exp in transformed quantities for speed up
//   from mc-stan discourse: speeding up a hierarchical negative binomial
//     regression with a huge number of groups
# Copyright (C) 2024, SNlab- ETH Zurich
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# version 2, as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
data {
  int N; // number of events * choice set (actors - 1)
  int T; // number of events
  int P; // number of covariates

  array[T] int<lower = 1, upper = N> chose; // position receiver chose
  matrix[N, P] X;

  array[T] int<lower = 1, upper = N> start; // the starting observation for each event
  array[T] int<lower = 1, upper = N> end; // the ending observation for each event
}
parameters {
  vector[P] beta;
}
model {
  // priors
  target += normal_lpdf(beta | 0, 4);

  // helper for likelihood
  vector[N] xb;
  xb = X * beta;

  for(t in 1:T) target += xb[chose[t]] - log_sum_exp(xb[start[t]:end[t]]);
}
