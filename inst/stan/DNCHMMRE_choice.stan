// Copyright (C) 2024, Alvaro Uzaheta-SNlab- ETH Zurich
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
// DNCHMMRE_choice.stan
// The choice model with P covariates
// Continuous-time HMM with K states
// Random effects in the state-dependent parameters
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
}
transformed data {
  //row_vector[K] v_ones = rep_row_vector(1.0, K);
  vector[K - 1] v_ones_Km1 = rep_vector(1.0, K - 1);
  //matrix[K, K] m_ones = rep_matrix(1.0, K, K);
  //real log_crude_rate = log(Tchoice / ((Nchoice + Tchoice) * mean(timespan))); // Trate / sum(timespan) / mean(actors) 

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
    target += gamma_lpdf(L_sigmaChoice[k] | 2, 1.0 / 2.0); // prior over sd, instead of cauchy
    target += std_normal_lpdf(to_vector(gamma_rawChoice[k]));
  }

  // auxiliar variables
  vector[G] log_prob;
  array[K] vector[K] log_theta_tr;
  matrix[K, K] theta_exp;
  // real sojourn;
  // real kappa;
  vector[K] lp;
  vector[K] lp_p1;

  int n_event_group;
  int last_event_group;

  // compute log probabilities of observed data given state
  array[Tchoice] vector[K] log_omega;
  {
    vector[Nchoice] xbChoice;
    for (n in 1:K) {
      xbChoice = Xchoice * betaChoice[n] + rows_dot_product(Zchoice, gammaChoice[n, groupChoice, ]);
      for (t in 1:Tchoice) {
        log_omega[t, n] = xbChoice[choseChoice[t]] -
          log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]);
      }
    }
  }

  // loop over random effects levels
  for (g in 1:G) {
    last_event_group = (g < G ? startGroup[g + 1] : Tchoice + 1);    
    n_event_group = last_event_group - startGroup[g];
    
    // forward algorithm by group
    lp = log(pi1) + log_omega[startGroup[g]]; // first observation

    if (n_event_group > 1) {
    
    for (t in (startGroup[g] + 1):last_event_group) { // looping over observations
      // compute log probabilities of transition given time in state

      theta_exp = matrix_exp(timesSender[t] * ta);
      for (n in 1:K) {
        // -\lambda_i * \Delta_t
        // sojourn = timespan[t] * ta[n, n];
        for (n_from in 1:K) {
          // transpose the tpm and take natural log of entries
          log_theta_tr[n, n_from] = log(theta_exp[n_from, n]);
        }
      }

      for (n in 1:K) { // looping over states
  /*       lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp) +
          xbChoice[n][choseChoice[t]] -
          log_sum_exp(xbChoice[n][startChoice[t]:endChoice[t]]);
  */
        lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp) + log_omega[t, n];
      }
      lp = lp_p1;
    }
    }
    log_prob[g] = log_sum_exp(lp);
  }

  target += log_prob;
}
generated quantities {
  array[K] matrix[Qchoice, Qchoice] SigmaChoice;
  for (k in 1:K)
    SigmaChoice[k] = diag_pre_multiply(L_sigmaChoice[k], L_OmegaChoice[k]) *
      diag_pre_multiply(L_sigmaChoice[k], L_OmegaChoice[k])';
}
