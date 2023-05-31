// DNHMM_rate_res_HMM.stan
// The rate model with intercept with P covariates
// offset by the log crude rate
// HMM dynamics with kR states
// the resolution of the HMM is different that the events
data {
  vector[2] alpha;
  //
  int Nrate; // number of events * present actors
  int Trate; // number of events
  int Prate; // number of covariates

  matrix[Nrate, Prate] Xrate;

  // waiting times between events and right-censoring
  array[Trate] real<lower = 0> timespan;
  array[Trate] int<lower = 0, upper = 1> isDependent;

  real offsetInt; // log crude rate: log(Trate / sum(timespan) / avgActors)

  array[Trate] int<lower = 0, upper = Nrate> choseRate; // position sender
  // the starting and ending index observation for each event
  array[Trate] int<lower = 1, upper = Nrate> startRate;
  array[Trate] int<lower = 1, upper = Nrate> endRate;

  // resolution HMM
  int Nres;
  array[Nres + 1] int<lower = 1, upper = Trate + 1> resA;
}
transformed data {
  int<lower = 2> kR;
  kR = %%kR%%;
}
parameters {
  array[kR] simplex[kR] theta; // transition probabilities

  ordered[kR] intercept;
  array[kR] vector[Prate - 1] betaWOInt;

}
transformed parameters {
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
    target += normal_lpdf(beta[n][1] | 0, 10);
    target += std_normal_lpdf(beta[n][2:]);
  }

  array[kR] vector[Nrate] xb;
  for (n in 1:kR)
    xb[n] = Xrate * beta[n] + offsetInt;

  // Using Stan implementation of HMM
  matrix[kR, Nres] log_omega;
  real acum;
  for (t in 1:Nres) {
    for (n in 1:kR) {
      acum = 0;
      for (event in resA[t]:(resA[t + 1] - 1))
        acum += isDependent[event] * xb[n][choseRate[event]] - timespan[event] *
          exp(log_sum_exp(xb[n][startRate[event]:endRate[event]]));

      log_omega[n, t] = acum;
    }
  }

  target += hmm_marginal(log_omega, ta, pi1);
}
generated quantities {
  matrix[kR, Nres] hidden_probs;
  array[Nres] int y_sim;
  {
    array[kR] vector[Nrate] xb;
    for (n in 1:kR)
      xb[n] = Xrate * beta[n] + offsetInt;

    matrix[kR, Nres] log_omega;
    real acum;
    for (t in 1:Nres) {
      for (n in 1:kR) {
        acum = 0;
        for (event in resA[t]:(resA[t + 1] - 1))
          acum += isDependent[event] * xb[n][choseRate[event]] - timespan[event] *
            exp(log_sum_exp(xb[n][startRate[event]:endRate[event]]));

        log_omega[n, t] = acum;
      }
    }

    hidden_probs = hmm_hidden_state_prob(log_omega, ta, pi1);
    y_sim = hmm_latent_rng(log_omega, ta, pi1);
  }
}
