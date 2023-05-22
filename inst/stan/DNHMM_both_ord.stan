// DyNAMSR_both_ord.stan
// The rate model with intercept model with P covariates
//
data {
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
  int Nchoice; // number of events * present actors
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position sender
  matrix[Nchoice, Pchoice] Xchoice;

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;
}
parameters {
  vector[Prate] betaRate;
  vector[Pchoice] betaChoice;
}
model {
  // priors
  target += normal_lpdf(betaRate | 0, 4);
  target += normal_lpdf(betaChoice | 0, 4);

  // helper for likelihood
  {
    vector[Nrate] xbRate;
    xbRate = Xrate * betaRate;

    vector[Nchoice] xbChoice;
    xbChoice = Xchoice * betaChoice;

    for(t in 1:Tchoice)
        target += xbChoice[choseChoice[t]] -
          log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]) +
          xbRate[choseRate[t]] -
          log_sum_exp(xbRate[startRate[t]:endRate[t]]);
  }
}
generated quantities {
  int<lower=1,upper=N> viterbi[T];
  real stateProbs[T,N];
  vector[N] lp;
  vector[N] lp_p1;

  // Viterbi algorithm (most likely state sequence)
  {
    real max_logp;
    int back_ptr[T, N];
    real best_logp[T, N];
    for (t in 1:T) {
      if(t==1 || ID[t]!=ID[t-1]) {
        for(n in 1:N)
        best_logp[t, n] = gamma_lpdf(steps[t] | shape[n], rate[n]);
      } else {
        for (n in 1:N) {
          best_logp[t, n] = negative_infinity();
          for (j in 1:N) {
            real logp;
            logp = best_logp[t-1, j] + log_theta[t,j,n];
            if(steps[t]>0)
            logp = logp + gamma_lpdf(steps[t] | shape[n], rate[n]);
            if(angles[t]>(-pi()))
            logp = logp + von_mises_lpdf(angles[t] | loc[n], kappa[n]);
            if (logp > best_logp[t, n]) {
              back_ptr[t, n] = j;
              best_logp[t, n] = logp;
            }
          }
        }
      }
    }
    for(t0 in 1:T) {
      int t = T - t0 + 1;
      if(t==T || ID[t+1]!=ID[t]) {
        max_logp = max(best_logp[t]);
        for (n in 1:N)
        if (best_logp[t, n] == max_logp)
        viterbi[t] = n;
      } else {
        viterbi[t] = back_ptr[t+1, viterbi[t+1]];
      }
    }
  }
  // forward-backward algorithm (state probabilities)
  {
    real logalpha[T,N];
    real logbeta[T,N];
    real llk;
    // log alpha probabilities
    for(t in 1:T) {
      if(t==1 || ID[t]!=ID[t-1]) {
        for(n in 1:N)
        lp[n] = -log(N);
      }
      for (n in 1:N) {
        lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t,n]) + lp);
        if(steps[t]>=0)
        lp_p1[n] = lp_p1[n] + gamma_lpdf(steps[t] | shape[n], rate[n]);
        if(angles[t]>=(-pi())) {
          lp_p1[n] = lp_p1[n] + von_mises_lpdf(angles[t] | loc[n], kappa[n]);
        }
        logalpha[t,n] = lp_p1[n];
      }
      lp = lp_p1;
    }
    // log beta probabilities
    for(t0 in 1:T) {
      int t = T - t0 + 1;
      if(t==T || ID[t+1]!=ID[t]) {
        for(n in 1:N)
        lp_p1[n] = 0;
      } else {
        for(n in 1:N) {
          lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t+1,n]) + lp);
          if(steps[t+1]>=0)
          lp_p1[n] = lp_p1[n] + gamma_lpdf(steps[t+1] | shape[n], rate[n]);
          if(angles[t+1]>=(-pi()))
          lp_p1[n] = lp_p1[n] + von_mises_lpdf(angles[t+1] | loc[n], kappa[n]);
        }
      }
      lp = lp_p1;
      for(n in 1:N)
      logbeta[t,n] = lp[n];
    }
    // state probabilities
    for(t0 in 1:T) {
      int t = T - t0 + 1;
      if(t==T || ID[t+1]!=ID[t])
      llk = log_sum_exp(logalpha[t]);
      for(n in 1:N)
      stateProbs[t,n] = exp(logalpha[t,n] + logbeta[t,n] - llk);
    }
  }
}

