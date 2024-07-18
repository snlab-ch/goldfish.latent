// DyNAMSR_choice.stan
// The choice model with P covariates
// Switching regime dynamics with kR states
data {
  int<lower = 2> kR;  // number of states
  real offsetInt; // prior for rate of transition between states
  //
  int Nchoice; // number of events * actors in the choice set
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position receiver
  matrix[Nchoice, Pchoice] Xchoice;

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;

  // time span between events
  array[Tchoice] real<lower = 0> timespan;
}
transformed data {
  row_vector[kR] v_ones = rep_row_vector(1.0, kR);
  vector[kR - 1] v_ones_kRm1 = rep_vector(1.0, kR - 1);
  matrix[kR, kR] m_ones = rep_matrix(1.0, kR, kR);
  real alpha = offsetInt / kR;
}
parameters {
  matrix<lower = 0>[kR, kR - 1] theta; // rates of transition between states
  array[kR] vector[Pchoice] betaChoice; // parms for each state
}
transformed parameters {
  simplex[kR] pi1;
  matrix[kR, kR] ta;

  {
    // \lambda_i = \sum_{j \neq i} \lambda_{ij}
    vector[kR] theta_row_sum = theta * v_ones_kRm1;
    for (j in 1:kR) {
      for (i in 1:kR) {
        if (i == j) {
          // row sums set to zero
          ta[i, j] = - theta_row_sum[i];
        } else {
          ta[i, j] = theta[i, (j > i? j - 1 : j)];
        }
      }
    }
  }
  // compute stationary distribution from transition prob matrix
  pi1 = to_vector(v_ones / (m_ones - ta));
}
model {
  target += lognormal_lpdf(to_vector(theta) | 0.28, 0.5);  

  for (n in 1:kR) {
    target += std_normal_lpdf(betaChoice[n]);
  }

  array[kR] vector[kR] log_theta_tr;
  real sojourn;
  vector[kR] lp;
  vector[kR] lp_p1;

  array[kR] vector[Nchoice] xbChoice;
  for (n in 1:kR)
    xbChoice[n] = Xchoice * betaChoice[n];

  // forward algorithm implementation

  for(n in 1:kR) // first observation
    lp[n] = log(pi1[n]) + xbChoice[n][choseChoice[1]] -
        log_sum_exp(xbChoice[n][startChoice[1]:endChoice[1]]);

  for (t in 2:Tchoice) { // looping over observations
    // compute log probabilities of transition given time in state
    for (n in 1:kR) {
      // -\lambda_i * \Delta_t
      sojourn = timespan[t] * ta[n, n];
      for (n_from in 1:kR) {
        // transpose the tpm and take natural log of entries
        if (n == n_from) {
          // the probability of staying in the same state: 
          // 1 - \lambda_i * exp(-\lambda_i * \Delta_t)
          log_theta_tr[n, n_from] = log(1 + ta[n, n] * exp(sojourn));
        } else {
          log_theta_tr[n, n_from] = log(ta[n_from, n]) + sojourn;
        }
      }
    }

    for (n in 1:kR) { // looping over states
      lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp) +
        xbChoice[n][choseChoice[t]] -
        log_sum_exp(xbChoice[n][startChoice[t]:endChoice[t]]);
    }
    lp = lp_p1;
  }

  target += log_sum_exp(lp);
}
// generated quantities {
//   array[Tchoice] int<lower=1, upper=kR> zstar;
//   real logp_zstar;
//
//   array[kR] vector[Nchoice] xbChoice;
//   for (n in 1:kR)
//     xbChoice[n] = Xchoice * betaChoice[n];
//
//   { // Viterbi algorithm
//     array[Tchoice, kR] int bpointer; // backpointer to the most likely previous state on the most probable path
//     array[Tchoice, kR] real delta; // max prob for the sequence up to t
//     // that ends with an emission from state k
//     for(k in 1:kR) // first observation
//       delta[1, k] = log(pi1[k]) + xbChoice[k][choseChoice[1]] -
//           log_sum_exp(xbChoice[k][startChoice[1]:endChoice[1]]);
//
//     for (t in 2:Tchoice) {
//       for (j in 1:kR) { // i = current (t)
//         delta[t, j] = negative_infinity();
//         for (i in 1:kR) { // i = previous (t-1)
//           real logp;
//           logp = delta[t-1, i] + log(theta[i, j]) +
//             xbChoice[j][choseChoice[t]] -
//             log_sum_exp(xbChoice[j][startChoice[t]:endChoice[t]]);
//             if (logp > delta[t, j]) {
//               bpointer[t, j] = i;
//               delta[t, j] = logp;
//             }
//         }
//       }
//     }
//     logp_zstar = max(delta[Tchoice]);
//     for (j in 1:kR)
//       if (delta[Tchoice, j] == logp_zstar)
//         zstar[Tchoice] = j;
//     for (t in 1:(Tchoice - 1)) {
//       zstar[Tchoice - t] = bpointer[Tchoice - t + 1, zstar[Tchoice - t + 1]];
//     }
//   }
// }
