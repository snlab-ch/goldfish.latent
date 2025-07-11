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
// DNCHMMRE_choice_mapReduce.stan
// The choice model with P covariates
// Continuous-time HMM with K states
// Random effects in the state-dependent parameters
// Map-reduce using group level
functions {
   real partial_sum(
    array[] int startGroup,
    int start, int end,
    array[] int endGroup,
    matrix Xchoice,
    matrix Zchoice,
    array[] int choseChoice,
    array[] int startChoice, array[] int endChoice,
    array[] real timesSender,
    int K,
    int Qchoice,
    matrix ta,
    vector pi1,
    array[] vector betaChoice,
    array[] matrix gammaChoice) {

    array[K] vector[K] log_theta_tr;
    matrix[K, K] theta_exp;
    vector[K] lp;
    vector[K] lp_p1;
    real logLik = 0;
    
    array[K] vector[Qchoice] gammaChoice_group;
    vector[K] log_omega;

    for (g in 1:(end - start + 1)) {  
      for (k in 1:K)
        for (q in 1:Qchoice)
          gammaChoice_group[k, q] = gammaChoice[k, g + start - 1, q];

      for (t in startGroup[g]:endGroup[g + start - 1]) {
        int size_slice = endChoice[t] - startChoice[t] + 1;
        array[size_slice] int event_slice =
          linspaced_int_array(size_slice, startChoice[t], endChoice[t]);
        
        vector[size_slice] xbChoice;

        for (k in 1:K) {
          xbChoice = Xchoice[event_slice] * betaChoice[k] + Zchoice[event_slice] * gammaChoice_group[k];
          log_omega[k] = xbChoice[choseChoice[t] - startChoice[t] + 1] - log_sum_exp(xbChoice);
        }
        
        // forward algorithm
        if (t == startGroup[g]) {
          // first observation
          lp = log(pi1) + log_omega;
        } else {
          // compute log probabilities of transition given time in state
          theta_exp = matrix_exp(timesSender[t] * ta);
          for (n in 1:K) {
            for (n_from in 1:K) {
              log_theta_tr[n, n_from] = log(theta_exp[n_from, n]);
            }
          }

          for (n in 1:K) {
            lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp) + log_omega[n];
          }
          lp = lp_p1;
        }
      }
      logLik += log_sum_exp(lp);
    }
    return logLik;
  }
}
data {
  int<lower = 2> K;  // number of states
  // real offsetInt; // kappa for Blackwell Integrated CHMM
  //
  int Nchoice; // number of events * actors in the choice set
  int Tchoice; // number of events
  int Pchoice; // number of covariates
  int Qchoice; // number of random effects
  int G; // number of levels/random effects

  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position receiver
  matrix[Nchoice, Pchoice] Xchoice;

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;

  // time span between sender events
  array[Tchoice] real<lower = 0> timesSender;

  // random effects vars
  array[Qchoice] int<lower = 1> V1choice; // index random effects
  array[G] int<lower = 1, upper = Tchoice> startGroup; // start event group
  array[Nchoice] int<lower = 1, upper = G> groupChoice; // group random effect

  //
  int grain_size;
}
transformed data {
  //row_vector[K] v_ones = rep_row_vector(1.0, K);
  vector[K - 1] v_ones_Km1 = rep_vector(1.0, K - 1);
  //matrix[K, K] m_ones = rep_matrix(1.0, K, K);
  //real log_crude_rate = log(Tchoice / ((Nchoice + Tchoice) * mean(timespan))); // Trate / sum(timespan) / mean(actors) 

  array[G] int<lower = 1, upper = Tchoice> endGroup; 
  for (g in 1:(G - 1)) {
    endGroup[g] = startGroup[g + 1] - 1;
  }
  endGroup[G] = Tchoice;

  matrix[Nchoice, Qchoice] Zchoice = Xchoice[, V1choice];
}
parameters {
  matrix<lower = 0>[K, K - 1] theta; // rates of transition between states
  array[K] vector[Pchoice] betaChoice; // parms for each state
  simplex[K] pi1;
  
  // parms random effects
  array[K] cholesky_factor_corr[Qchoice] L_OmegaChoice; // correlation RE
  array[K] vector<lower = 0>[Qchoice] L_sigmaChoice; // scale RE  real<lower=0>

  array[K] matrix[Qchoice, G] gamma_rawChoice;  // group Random Effect
}
transformed parameters {
  matrix[K, K] ta;

  {
    // \lambda_i = \sum_{j \neq i} \lambda_{ij}
    vector[K] theta_row_sum = theta * v_ones_Km1;
    for (j in 1:K) {
      for (i in 1:K) {
        if (i == j) {
          // row sums set to zero
          ta[i, j] = - theta_row_sum[i];
        } else {
          ta[i, j] = theta[i, (j > i? j - 1 : j)];
        }
      }
    }
  }

  array[K] matrix[G, Qchoice] gammaChoice;
  for (k in 1:K)
    gammaChoice[k] = (diag_pre_multiply(L_sigmaChoice[k], L_OmegaChoice[k]) * gamma_rawChoice[k])';
}
model {
  // target += lognormal_lpdf(to_vector(theta) | 0.28, 0.5);  
  target += gamma_lpdf(to_vector(theta) | 2, 0.1); // Blackwell 2016  

  target += dirichlet_lpdf(pi1 | rep_vector(1.0, K)); // prior initial state

  for (k in 1:K) {
    target += std_normal_lpdf(betaChoice[k]);
    target += lkj_corr_cholesky_lpdf(L_OmegaChoice[k] | 2); // bigger than one
    target += gamma_lpdf(L_sigmaChoice[k] | 2, 2.0); // prior over sd, instead of cauchy
    target += std_normal_lpdf(to_vector(gamma_rawChoice[k]));
  }

  target += reduce_sum(partial_sum, startGroup, grain_size, 
    endGroup, Xchoice, Zchoice, choseChoice, startChoice, endChoice,
    timesSender, K, Qchoice, ta, pi1, betaChoice, gammaChoice);
}
generated quantities {
  array[K] matrix[Qchoice, Qchoice] SigmaChoice;
  for (k in 1:K)
    SigmaChoice[k] = diag_pre_multiply(L_sigmaChoice[k], L_OmegaChoice[k]) *
      diag_pre_multiply(L_sigmaChoice[k], L_OmegaChoice[k])';
}
