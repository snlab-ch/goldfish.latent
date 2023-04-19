// DyNAMSR_rate.stan
// The rate model with intercept model with P covariates
// Switching regime dynamics with kR states
data {
  int<lower = 2> kR;
  vector[2] alpha;
  //
  int Nrate; // number of events * present actors
  int Trate; // number of events
  int Prate; // number of covariates

  array[Trate] int<lower = 0, upper = Nrate> choseRate; // position sender
  matrix[Nrate, Prate] Xrate;

  // the starting and ending index observation for each event
  array[Trate] int<lower = 1, upper = Nrate> startRate;
  array[Trate] int<lower = 1, upper = Nrate> endRate;

  //
  array[Trate] real<lower = 0> timespan;
  array[Trate] int<lower = 0, upper = 1> isDependent;

  real offsetInt;
}
parameters {
  array[kR] simplex[kR] theta; // transition probabilities

  ordered[kR] intercept;
  array[kR] vector[Prate - 1] betaWOInt;

}
transformed parameters {
  //array[kR] vector[Prate] beta;
  //for (p in 1:Prate)
  //  for (n in 1:kR)
  //    beta[n, p] = betaOrd[p, n];

  array[kR] vector[Prate] beta;
  for (n in 1:kR) {
    beta[n][1] = intercept[n];
    beta[n][2:Prate] = betaWOInt[n];
  }

  simplex[kR] pi1;
  matrix[kR, kR] ta;
  //array[kR] vector[Trate] logAlpha;

  for (j in 1:kR)
    for (i in 1:kR)
      ta[i, j] = theta[i, j];

  // compute stationary distribution from transition prob matrix
  pi1 = to_vector((to_row_vector(rep_vector(1.0, kR)) /
    (diag_matrix(rep_vector(1.0, kR)) - ta + rep_matrix(1, kR, kR))));
}
model {
  vector[kR] alphas;
  for (n in 1:kR) {
    alphas = rep_vector(alpha[2], kR);
    alphas[n] = alpha[1];
    target += dirichlet_lpdf(theta[n] | alphas); // prior for trans probs
    target += normal_lpdf(beta[n] | 0, 4);
  }

  array[kR] vector[kR] log_theta_tr;
  vector[kR] lp;
  vector[kR] lp_p1;

  array[kR] vector[Nrate] xb;
  for (n in 1:kR)
    xb[n] = Xrate * beta[n] + offsetInt;

  // transpose the tpm and take natural log of entries
  for (n_from in 1:kR)
    for (n in 1:kR)
      log_theta_tr[n, n_from] = log(theta[n_from, n]);

  // forward algorithm implementation

  for(n in 1:kR) // first observation
    lp[n] = log(pi1[n]) + isDependent[1] * xb[n][choseRate[1]] -
        timespan[1] * exp(log_sum_exp(xb[n][startRate[1]:endRate[1]]));

  for (t in 2:Trate) { // looping over observations
    for (n in 1:kR) // looping over states
      lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp) +
      isDependent[t] * xb[n][choseRate[t]] -
        timespan[t] * exp(log_sum_exp(xb[n][startRate[t]:endRate[t]]));

    lp = lp_p1;
  }

  target += log_sum_exp(lp);
}
generated quantities {
  array[Trate] int<lower=1, upper=kR> zstar;
  real logp_zstar;

  array[kR] vector[Nrate] xb;
  for (n in 1:kR)
    xb[n] = Xrate * beta[n] + offsetInt;

  { // Viterbi algorithm
    array[Trate, kR] int bpointer; // backpointer to the most likely previous state on the most probable path
    array[Trate, kR] real delta; // max prob for the sequence up to t
    // that ends with an emission from state k
    for(k in 1:kR) // first observation
      delta[1, k] = log(pi1[k]) + isDependent[1] * xb[n][choseRate[1]] -
        timespan[1] * exp(log_sum_exp(xb[n][startRate[1]:endRate[1]]));

    for (t in 2:Trate) {
      for (j in 1:kR) { // i = current (t)
        delta[t, j] = negative_infinity();
        for (i in 1:kR) { // i = previous (t-1)
          real logp;
          logp = delta[t-1, i] + log(theta[i, j]) +
            isDependent[t] * xb[n][choseRate[t]] -
            timespan[t] * exp(log_sum_exp(xb[n][startRate[t]:endRate[t]]));
            if (logp > delta[t, j]) {
              bpointer[t, j] = i;
              delta[t, j] = logp;
            }
        }
      }
    }
    logp_zstar = max(delta[Trate]);
    for (j in 1:kR)
      if (delta[Trate, j] == logp_zstar)
        zstar[Trate] = j;
    for (t in 1:(Trate - 1)) {
      zstar[Trate - t] = bpointer[Trate - t + 1, zstar[Trate - t + 1]];
    }
  }
}
