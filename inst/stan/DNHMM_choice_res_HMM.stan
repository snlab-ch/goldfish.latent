// DyNAMSR_choice.stan
// The rate model with intercept model with P covariates
// HMM dynamics with kR states
// the resolution of the HMM is different that the events
data {
  vector[2] alpha;
  //
  int Nchoice; // number of events * present actors
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  matrix[Nchoice, Pchoice] Xchoice;

  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position receiver
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
  vector[kR] alphas;
  for (n in 1:kR) {
    alphas = rep_vector(alpha[2], kR);
    alphas[n] = alpha[1];
    target += dirichlet_lpdf(theta[n] | alphas); // prior for trans probs
    target += std_normal_lpdf(betaChoice[n]);
  }

  matrix[kR, Nres] log_omega;


    array[kR] vector[Nchoice] xbChoice;
    for (n in 1:kR)
      xbChoice[n] = Xchoice * betaChoice[n];

    // Using Stan implementation of HMM

    real acum;
    for (t in 1:Nres) {
      for (n in 1:kR) {
        acum = 0;
        for (event in resA[t]:(resA[t + 1] - 1))
          acum += xbChoice[n][choseChoice[event]] -
            log_sum_exp(xbChoice[n][startChoice[event]:endChoice[event]]);

        log_omega[n, t] = acum;
      }
    }


  target += hmm_marginal(log_omega, ta, pi1);
}
// generated quantities {
//   matrix[kR, Nres] hidden_probs;
//   array[Nres] int y_sim;
//   {
//     array[kR] vector[Nchoice] xbChoice;
//     for (n in 1:kR)
//       xbChoice[n] = Xchoice * betaChoice[n];
//
//     matrix[kR, Nres] log_omega;
//     real acum;
//     for (t in 1:Nres) {
//       for (n in 1:kR) {
//         acum = 0;
//         for (event in resA[t]:(resA[t + 1] - 1))
//           acum += xbChoice[n][choseChoice[event]] -
//             log_sum_exp(xbChoice[n][startChoice[event]:endChoice[event]]);
//
//         log_omega[n, t] = acum;
//       }
//     }
//
//     hidden_probs = hmm_hidden_state_prob(log_omega, ta, pi1);
//     y_sim = hmm_latent_rng(log_omega, ta, pi1);
//   }
// }
