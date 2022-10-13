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
