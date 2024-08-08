// DyNAMSR_choice.stan
// The choice model with P covariates
// Switching regime dynamics with kR states
data {
  int<lower = 2> K;  // number of states
  // real offsetInt; // kappa for Blackwell Integrated CHMM
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
  row_vector[K] v_ones = rep_row_vector(1.0, K);
  vector[K - 1] v_ones_Km1 = rep_vector(1.0, K - 1);
  matrix[K, K] m_ones = rep_matrix(1.0, K, K);
  real log_crude_rate = log(Tchoice / ((Nchoice + Tchoice) * mean(timespan))); // Trate / sum(timespan) / mean(actors) 
}
parameters {
  matrix<lower = 0>[K, K - 1] theta; // rates of transition between states
  array[K] vector[Pchoice] betaChoice; // parms for each state
}
transformed parameters {
  simplex[K] pi1;
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
  // compute stationary distribution from transition prob matrix
  pi1 = to_vector(v_ones / (m_ones - ta));
}
model {
  // target += lognormal_lpdf(to_vector(theta) | 0.28, 0.5);  
  target += gamma_lpdf(to_vector(theta) | 2, 0.1); // Blackwell 2016  

  for (n in 1:K) {
    target += std_normal_lpdf(betaChoice[n]);
  }

  array[K] vector[K] log_theta_tr;
  matrix[K, K] theta_exp;
  // real sojourn;
  // real kappa;
  vector[K] lp;
  vector[K] lp_p1;

  // compute log probabilities of observed data given state
  array[Tchoice] vector[K] log_omega;
  {
    vector[Nchoice] xbChoice;
    for (n in 1:K) {
      xbChoice = Xchoice * betaChoice[n];
      for (t in 1:Tchoice) {
        log_omega[t, n] = xbChoice[choseChoice[t]] -
          log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]);
      }
    }
  }

  /* array[kR] vector[Nchoice] xbChoice;
  for (n in 1:kR)
    xbChoice[n] = Xchoice * betaChoice[n]; */

  // forward algorithm implementation

/*   for(n in 1:kR) // first observation
    lp[n] = log(pi1[n]) + xbChoice[n][choseChoice[1]] -
        log_sum_exp(xbChoice[n][startChoice[1]:endChoice[1]]);
 */
  lp = log(pi1) + log_omega[1]; // first observation

  for (t in 2:Tchoice) { // looping over observations
    // compute log probabilities of transition given time in state

    theta_exp = matrix_exp(timespan[t] * ta);
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
