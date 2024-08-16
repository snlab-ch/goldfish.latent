// DNCHMM_both.stan
// The DyNAM (rate and choice) model with P covariates
// Switching regime dynamics follow a continuous time HMM with kR states
data {
  int<lower = 2> K;  // number of states
  //real offsetInt; // crude event rate: 1 / mean(timespan) = Trate / sum(timespan) 
  // data rate sub-model
  int Nrate; // number of events * present actors
  int Trate; // number of events
  int Prate; // number of covariates

  array[Trate] int<lower = 0, upper = Nrate> choseRate; // position sender
  matrix[Nrate, Prate] Xrate; // effects statistics rate model

  // the starting and ending index observation for each event
  array[Trate] int<lower = 1, upper = Nrate> startRate;
  array[Trate] int<lower = 1, upper = Nrate> endRate;

  // inter-event times and indicator right-censoring
  array[Trate] real<lower = 0> timespan;
  array[Trate] int<lower = 0, upper = 1> isDependent;
  
  // data choice sub-model
  int Nchoice; // number of events * actors in the choice set
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position receiver
  matrix[Nchoice, Pchoice] Xchoice; // effects statistics choice model

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;
}
transformed data {
  row_vector[K] v_ones = rep_row_vector(1.0, K);
  vector[K - 1] v_ones_Km1 = rep_vector(1.0, K - 1);
  matrix[K, K] m_ones = rep_matrix(1.0, K, K);
  //real delta = mean(timespan);
  real log_crude_rate = log(Trate / (Nrate * mean(timespan))); // Trate / sum(timespan) / mean(actors) 
  //real mu_rate = log(2 / (delta * K * sqrt(4 + (K - 1)^2)));
  //real sigma_rate = sqrt(log(1 + (K - 1)^2 / 4));
}
parameters {
  matrix<lower = 0>[K, K - 1] theta; // rates of transition between states
  array[K] vector[Prate] betaRate; // parms for each state rate model
  array[K] vector[Pchoice] betaChoice; // parms for each state choice model
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
  //target += lognormal_lpdf(to_vector(theta) | mu_rate, sigma_rate);  
  // target += gamma_lpdf(to_vector(theta) | 2, 0.1); // Blackwell 2016  
  target += gamma_lpdf(to_vector(theta) | 2, 0.5); // distribution concentrated below 1
  // target += exponential_lpdf(to_vector(theta) | 1); //


  for (n in 1:K) {
    // Bayesian Survival Analysis Using rstanarm N(0, 20)
    target += normal_lpdf(betaRate[n][1] | 0, 10);
    target += std_normal_lpdf(betaRate[n][2 : ]);
    target += std_normal_lpdf(betaChoice[n]);
  }

  array[K] vector[K] log_theta_tr;
  matrix[K, K] theta_exp;
  // real sojourn;
  // real kappa;
  vector[K] lp;
  vector[K] lp_p1;

  // compute log probabilities of observed data given state
  array[Trate] vector[K] log_omega;
  {
    vector[Nrate] xbRate;
    vector[Nchoice] xbChoice;
    real acum;
    int ttChoice;
    for (n in 1:K) {
      xbRate = Xrate * betaRate[n] + log_crude_rate;
      xbChoice = Xchoice * betaChoice[n];
      ttChoice = 1;
      for (t in 1:Trate) {
        acum = -timespan[t] * exp(log_sum_exp(xbRate[startRate[t]:endRate[t]]));
        if (isDependent[t]) {
          log_omega[t, n] = xbRate[choseRate[t]] + acum + // rate log-likelihood
            xbChoice[choseChoice[ttChoice]] - // choice log-likelihood
              log_sum_exp(xbChoice[startChoice[ttChoice]:endChoice[ttChoice]]); 
          ttChoice += 1;
        } else {
          log_omega[t, n] = acum; // rate log-likelihood: right-censored event
        }
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
