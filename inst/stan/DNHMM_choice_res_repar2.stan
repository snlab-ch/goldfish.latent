// DyNAMSR_choice.stan
// The rate model with intercept model with P covariates
// Switching regime dynamics with kR states
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
  int<lower = 2> kR;
  kR = %%kR%%;
}
parameters {
  array[kR] simplex[kR] theta; // transition probabilities

  real int0;
  array[kR - 1] real<lower = 0.5> intk;

  array[kR] vector[Pchoice - 1] betaWOInt;
}
transformed parameters {
  // array[kR] vector[Pchoice] betaChoice;
  // for (p in 1:Pchoice)
  //   for (n in 1:kR)
  //     betaChoice[n, p] = betaOrd[p, n];

  array[kR] vector[Pchoice] betaChoice;
  betaChoice[1][1] = int0;
  for (n in 1:kR) {
    if (n > 1)  betaChoice[n][1] = betaChoice[n-1][1] + intk[n-1];
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
  vector[kR] alphas;
  target += std_normal_lpdf(int0);
  target += normal_lpdf(intk | 1, 1);
  for (n in 1:kR) {
    alphas = rep_vector(alpha[2], kR);
    alphas[n] = alpha[1];
    target += dirichlet_lpdf(theta[n] | alphas); // prior for trans probs
    target += std_normal_lpdf(betaWOInt[n]);
  }

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

  for(n in 1:kR) { // first observation
    lp[n] = log(pi1[n]);
    for (event in resA[1]:(resA[2] - 1))
      lp[n] += xbChoice[n][choseChoice[event]] -
        log_sum_exp(xbChoice[n][startChoice[event]:endChoice[event]]);
  }

  for (t in 2:Nres) { // looping over observations
    for (n in 1:kR) { // looping over states
      lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp);
      for (event in resA[t]:(resA[t + 1] - 1))
        lp_p1[n] += xbChoice[n][choseChoice[event]] -
          log_sum_exp(xbChoice[n][startChoice[event]:endChoice[event]]);
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
