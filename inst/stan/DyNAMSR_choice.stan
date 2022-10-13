// DyNAMSR_choice.stan
// The rate model with intercept model with P covariates
// Switching regime dynamics with kR states
data {
  int<lower = 2> kR;
  //
  int Nchoice; // number of events * present actors
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position receiver
  matrix[Nchoice, Pchoice] Xchoice;

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;
}
parameters {
  array[kR] simplex[kR] theta; // transition probabilities

  ordered[kR] intercept;
  array[kR] vector[Pchoice - 1] betaWOInt;
}
transformed parameters {
  array[kR] vector[Pchoice] betaChoice;
  for (n in 1:kR) {
    betaChoice[n][1] = intercept[n];
    betaChoice[n][2:Pchoice] = betaWOInt[n];
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
  // priors for beta
  for (n in 1:kR)
    target += normal_lpdf(betaChoice[n] | 0, 4);

  array[kR] vector[kR] log_theta_tr;
  vector[kR] lp;
  vector[kR] lp_p1;

  array[kR] vector[Nchoice] xbChoice;
  for (n in 1:kR)
    xbChoice[n] = Xchoice * betaChoice[n];

  // transpose the tpm and take natural log of entries
  for (n_from in 1:kR)
    for (n in 1:kR)
      log_theta_tr[n, n_from] = log(theta[n_from, n]);

  // forward algorithm implementation

  for(n in 1:kR) // first observation
    lp[n] = log(pi1[n]) + xbChoice[n][choseChoice[1]] -
        log_sum_exp(xbChoice[n][startChoice[1]:endChoice[1]]);

  for (t in 2:Tchoice) { // looping over observations
    for (n in 1:kR) // looping over states
      lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp) +
        xbChoice[n][choseChoice[t]] -
        log_sum_exp(xbChoice[n][startChoice[t]:endChoice[t]]);

    lp = lp_p1;
  }

  target += log_sum_exp(lp);
}
