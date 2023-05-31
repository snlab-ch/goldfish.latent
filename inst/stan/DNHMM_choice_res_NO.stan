// DyNAMSR_choice.stan
// The rate model with intercept model with P covariates
// Switching regime dynamics with kR states
functions {
  // Normalize
  vector normalize(array[] real x) {
    return to_vector(x) / sum(x);
  }
}
data {
  vector[2] alpha;
  //
  int Nchoice; // number of events * present actors
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position receiver
  matrix[Nchoice, Pchoice] Xchoice;

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;

  // resolution HMM
  int Nres;
  array[Nres + 1] int<lower = 1, upper = Tchoice + 1> resA;
}
transformed data {
  int<lower = 2> K;
  K = %%kR%%;
  int T = Nres;
}
parameters {
  array[K] simplex[K] theta; // transition probabilities
  array[K] vector[Pchoice] betaChoice;
}
transformed parameters {
  simplex[K] pi1;
  matrix[K, K] ta;

  for (j in 1:K)
    for (i in 1:K)
      ta[i, j] = theta[i, j];
  // compute stationary distribution from transition prob matrix
  pi1 = to_vector((to_row_vector(rep_vector(1.0, K)) /
    (diag_matrix(rep_vector(1.0, K)) - ta + rep_matrix(1, K, K))));
}
model {
  vector[K] alphas;
  for (n in 1:K) {
    alphas = rep_vector(alpha[2], K);
    alphas[n] = alpha[1];
    target += dirichlet_lpdf(theta[n] | alphas); // prior for trans probs
    target += std_normal_lpdf(betaChoice[n]);
  }

  array[K] vector[K] log_theta_tr;
  vector[K] lp;
  vector[K] lp_p1;

  array[T] vector[K] log_omega;
  {
    vector[Nchoice] xbChoice;
    real acum;
    for (n in 1:K) {
      xbChoice = Xchoice * betaChoice[n];
      for (t in 1:Nres) {
        acum = 0;
        for (event in resA[t]:(resA[t + 1] - 1))
          acum += xbChoice[choseChoice[event]] -
        log_sum_exp(xbChoice[startChoice[event]:endChoice[event]]);

        log_omega[t, n] = acum;
      }
    }
  }

  // transpose the tpm and take natural log of entries
  for (n_from in 1:K)
    for (n in 1:K)
      log_theta_tr[n, n_from] = log(theta[n_from, n]);

  // forward algorithm implementation

  lp = log(pi1) + log_omega[1]; // first observation

  for (t in 2:Nres) { // looping over observations
    for (n in 1:K) // looping over states
      lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp) + log_omega[t, n];

    lp = lp_p1;
  }

  target += log_sum_exp(lp);
}
generated quantities {
  array[T] int zstar;
  array[K] vector[T] gamma;

  // viterbi and forward-backward alg adapted from BayesHMM package
  {
  // compute log likelihood
  array[K] vector[T] log_omega;
  vector[Nchoice] xbChoice;
  real acum;
  for (n in 1:K) {
    xbChoice = Xchoice * betaChoice[n];
    for (t in 1:Nres) {
      acum = 0;
      for (event in resA[t]:(resA[t + 1] - 1))
        acum += xbChoice[choseChoice[event]] -
          log_sum_exp(xbChoice[startChoice[event]:endChoice[event]]);

      log_omega[n, t] = acum;
    }
  }
  // viterbi
  real logp_zstar;
  array[T, K] int bpointer;
  array[T, K] real delta;

  for (j in 1:K)
      delta[1, j] = log_omega[j, 1] + log(pi1[j]);

  for (t in 2:T) {
      for (j in 1:K) { // j = current (t)
        delta[t, j] = negative_infinity();
        for (i in 1:K) { // i = previous (t-1)
          real logp;
          logp = delta[t-1, i] + log(ta[i, j]) + log_omega[j, t];
          if (logp > delta[t, j]) {
            bpointer[t, j] = i;
            delta[t, j] = logp;
          }
        }
      }
    }

    logp_zstar = max(delta[T]);

    for (j in 1:K) {
      if (delta[T, j] == logp_zstar) {
        zstar[T] = j;
      }
    }

    for (t in 1:(T - 1)) {
      zstar[T - t] = bpointer[T - t + 1, zstar[T - t + 1]];
    }

    // forward-backward alg
    array[K] vector[T] logalpha;
    array[K] vector[T] p_alpha;

    array[K] vector[T] logbeta;
    array[K] vector[T] loggamma;
    array[K] vector[T] p_beta;
    vector[K] accumulator;

    // forward
    for(j in 1:K) {
      logalpha[j, 1] = log(pi1[j]) + log_omega[j, 1];
    }

    for (t in 2:T) {
      for (j in 1:K) { // j = current (t)
        for (i in 1:K) { // i = previous (t-1)
          accumulator[i] = logalpha[i, t-1] + log(ta[i, j]) + log_omega[j, t];
        }
        logalpha[j, t] = log_sum_exp(accumulator);
      }
    }

    for (t in 1:T) {
      p_alpha[, t] = to_array_1d(softmax(to_vector(logalpha[, t])));
    }

    // backward
    for (j in 1:K) {
      logbeta[j, T] = 1;
    }

    for (tforward in 0:(T-2)) {
      int t;
      t = T - tforward;

      for (j in 1:K) { // j = previous (t-1)
        for (i in 1:K) { // i = next (t)
          accumulator[i] = logbeta[i, t] + log(ta[j, i]) + log_omega[i, t];
        }
        logbeta[j, t-1] = log_sum_exp(accumulator);
      }
    }

    for (t in 1:T) {
      p_beta[, t] = to_array_1d(softmax(to_vector(logbeta[, t])));
    }

    for(t in 1:T) {
      loggamma[, t] = to_array_1d(
        to_vector(p_alpha[, t]) .* to_vector(p_beta[, t])
      );
    }

    for(t in 1:T) {
      gamma[, t] = to_array_1d(normalize(loggamma[, t]));
    }
  }
}
